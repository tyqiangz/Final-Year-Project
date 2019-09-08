import csv
import math
from sympy import *
import numpy as np
from scipy import special
import networkx as nx
from networkx.algorithms import bipartite
import matplotlib.pyplot as plt

# input:
#   r: length of first row
#   wi: Hamming weight of first row
# output:
#   numpy array of the first row of a circulant block
def genFirstRow(r, wi):
    firstRow = np.zeros(r, dtype=np.int32)
    # pick wi numbers from range(wi)
    randomArray = np.random.choice(r, wi, replace=False)

    for i in randomArray:
        firstRow[i] = 1

    return firstRow

# input:  the first row of the binary circulant matrix
# output: the binary circulant matrix
def genCirculant(firstRow):
    rows = firstRow.size
    firstRow = np.array(firstRow, dtype = np.int32)
    M = np.zeros((rows, rows), dtype = np.int32)
    M[0,:] = firstRow

    for i in range(1,rows):
        for j in range(rows):
            M[i,j] = M[i-1, (j-1)%rows]
    return M

# input: the first row of a circulant matrix
# output: the first row of the transpose of the circulant matrix
def genTransposePoly(firstRow):
    transposePoly = np.zeros(len(firstRow), dtype = np.int32)
    transposePoly[0] = firstRow[0]
    transposePoly[1:] = np.flip(firstRow[1:])
    return transposePoly

# input: the first row of 2 circulant matrices
# output: the first row of the sum of the circulant matrices
def genSumPoly(firstRowA, firstRowB):
    return firstRowA + firstRowB

# input: the first row of 2 circulant matrices
# output: the first row of the product of the circulant matrices
def genProdPoly(firstRowA, firstRowB):
    r = len(firstRowA)
    prodPoly = np.zeros(r, dtype = np.int32)
    convolvePoly = np.convolve(firstRowA, firstRowB)

    print("ConvolvePoly:\n", convolvePoly)

    if len(convolvePoly) > r:
        prodPoly[0:r] = convolvePoly[0:r]
        prodPoly[0: len(convolvePoly) - r] += convolvePoly[r:]
    else:
        prodPoly[0: len(convolvePoly)] = convolvePoly[0: len(convolvePoly)]
        
    return prodPoly
# input: the first row of a binary circulant matrix
# output: the first row of the inverse of the binary circulant matrix
def genInvPoly(firstRow):
    r = len(firstRow)

    # convert a numpy array to a polynomial
    invPoly = convertNumpyToSympy(firstRow)

    #defining a polynomial ring, F_2 / (x**r - 1)
    FField = "x**" + str(r) + "-" + str(1)
    FField = poly(FField, domain=FF(2))

    #finding the inverse of the polynomial in the specified polynomial ring
    invPoly = invert(poly(invPoly, domain=FF(2)), poly(FField, domain=FF(2)))
    
    return convertSympyToNumpy(invPoly)

# input: a numpy array
# output: a Sympy polynomial
# Assumptions: The coefficients are non-negative
def convertNumpyToSympy(f):
    polynomial = ""
    polynomial = polynomial + str(int(f[0]))
    
    for i in range(1, f.size):
        if f[i] != 0:
            polynomial += " + " + str(int(f[i])) + "*x**" + str(i)
    return polynomial

# input: a Sympy polynomial
# output: a numpy array
def convertSympyToNumpy(f):
    return np.array(f.all_coeffs())

# input:
#   firstRow: the first row of the QC-MDPC matrix
#   n: length of the QC-MDPC code/ no. of bit nodes in the QC-MDPC code
#   r: no. of rows/cols of each circulant block of the QC-MDPC code
#   w: sum of the weights (wi) of the first row of all circulant blocks
#
# output: A QC-MDPC matrix constructed uing the first row
#   note that the weight distribution is randomised.

# Note that firstRow is of length n = n0 * r
# Define f_i(x) to be the Hall's polynomial of H_i
def genParityCheck(n, r, w):
    n0 = int(n / r)
    wi = int(w / r)
    
    H = np.zeros((r, r * n0), dtype=np.int32)
    H_i = np.zeros((r, r), dtype=np.int32)

    for i in range(n0):
        firstRow = genFirstRow(r, wi)
        H_i = genCirculant(firstRow)

        if i == 0:
            H = H_i
        else :
            H = np.concatenate((H, H_i), axis=1)

    print("Parity-Check matrix H:\n", H, "\n")

    filename = str(n) + "_" + str(r) + "_" + str(w) + "_" + "parityCheckMatrix.csv"

    np.savetxt(filename, H, delimiter = ",", fmt = "%d")

    return H

# input: an (n,r,w)-QC-MDPC matrix, H
# output: the number of 4 cycles in the Tanner graph of H
def count4Cycles(H, n, r, w):
    num4Cycles = 0
    
    HHT = np.matmul(H, (np.transpose(H)))

    for i in range(r-1):
        for j in range(i+1,r):
            if HHT[i, j] >= 2:
                num4Cycles += special.binom(int(HHT[i, j]), 2)
    print("Number of 4 cycles in parity-check matrix H:", int(num4Cycles))
    
    return num4Cycles

# Input: A parity-check matrix H (not necessarily QC-MDPC)
# Output: A plot of the Tanner graph of H
def drawTanner(H):
    edges = []
    pos = {}
    
    # Give names to the nodes in the two node sets
    C = [ "c{}".format(i) for i in range(1,H.shape[0]+1) ]
    B = [ "b{}".format(i) for i in range(1,H.shape[1]+1) ]

    # Create the graph and add each set of nodes
    G = nx.Graph()
    G.add_nodes_from(B, bipartite=1)
    G.add_nodes_from(C, bipartite=0)

    # Find the non-zero indices in the biadjacency matrix to connect 
    # those nodes
    for i in range(H.shape[1]):
        for j in range(H.shape[0]):
            if H[j,i] == 1:
                edges.append((B[i], C[j]))

    G.add_edges_from( edges )

    # Setting the cartesian coordinates of the nodes, the nodes spread across [-1,1] on the x-axis
    # Bit nodes at the top of the plot, hence y-coor is 1
    count = 0
    for i in B:
        pos[i] = [-1 + (2 * count) / (len(B) - 1), 1]
        count += 1

    # Check nodes below bit nodes, hence y-coor is -1
    count = 0
    for i in C:
        pos[i] = [-1 + (2 * count) / (len(C) - 1), -1]
        count += 1

    #plot graph G
    nx.draw_networkx(G, pos = pos, node_color = "white", with_labels = True)
    plt.show()
    

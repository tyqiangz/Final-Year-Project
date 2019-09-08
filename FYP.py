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
    invPoly = np.zeros(r, dtype = np.int32)

    # convert a numpy array to a polynomial
    inv = convertNumpyToSympy(firstRow)

    #defining a polynomial ring, F_2 / (x**r - 1)
    FField = "x**" + str(r) + "-" + str(1)
    FField = poly(FField, domain=FF(2))

    #finding the inverse of the polynomial in the specified polynomial ring
    inv = invert(poly(inv, domain=FF(2)), poly(FField, domain=FF(2)))
    temp = convertSympyToNumpy(inv)

    invPoly[0 : len(temp)] = temp
    
    return invPoly

# input: numpy array with integer values
# output: numpy array modulo 2
def convertBinary(v):
    for i in range(len(v)):
        v[i] = v[i] % 2

    return v

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
    v = np.array(f.all_coeffs())

    #Sympy stores the coefficients of polynomial f(x) in decreasing powers of x
    #thus the order is reversed to obtain numpy array of increasing powers of x
    v = np.flip(v)
    
    return v

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
def genQCMDPC(n, r, w):
    n0 = int(n / r)
    wi = int(w / r)
    
    H = np.zeros((r, r * n0), dtype = np.int32)
    H_i = np.zeros((r, r), dtype = np.int32)

    for i in range(n0):
        firstRow = genFirstRow(r, wi)
        H_i = genCirculant(firstRow)

        if i == 0:
            H = H_i
        else :
            H = np.concatenate((H, H_i), axis=1)

    filename = str(n) + "_" + str(r) + "_" + str(w) + "_" + "parityCheckMatrix.csv"

    np.savetxt(filename, H, delimiter = ",", fmt = "%d")

    return H

def genGenQCMDPC(H):
    r, n = H.shape
    n0 = int(n / r)
    H_i = np.zeros((r, r), dtype = np.int32)
    G = np.eye(n - r, dtype = np.int32)
    block0 = np.zeros((r, r), dtype = np.int32)

    #extract the generator polynomials for h_{n_0-1} from H
    lastPoly = H[0, n - r : n]

    #compute generator polynomial of inverse of h_{n_0-1}
    invLastPoly = genInvPoly(lastPoly)
    print(invLastPoly)

    #compute the first circulant block for G
    temp = genProdPoly(invLastPoly, H[0, 0:r])
    temp1 = convertBinary(temp)
    temp2 = genTransposePoly(temp1)
    block0 = genCirculant(temp2)

    #compute the subsequent circulant block for G and concatenating them
    for i in range(1, n0 - 1):
        temp = genProdPoly(invLastPoly, H[0, i*r : (i+1)*r])
        temp1 = convertBinary(temp)
        temp2 = genTransposePoly(temp1)
        block = genCirculant(temp2)
        block0 = np.concatenate((block0, block), axis = 0)

    #concatenate identity matrix of order (n - r) and the stack of circulant matrices to form G
    G = np.concatenate((G, block0), axis = 1)
    
    return G

# input: an (n,r,w)-QC-MDPC matrix, H
# output: the number of 4 cycles in the Tanner graph of H
def count4Cycles(H, n, r, w):
    num4Cycles = 0
    firstRow = np.zeros(r, dtype = np.int32)
    temp1 = np.zeros(r, dtype = np.int32)
    temp2 = np.zeros(r, dtype = np.int32)
    temp3 = np.zeros(r, dtype = np.int32)
    
    #compute the Hall's polynomial of HH^T
    for i in range( int(n / r) ):
        temp1 = H[0, i*r : (i+1)*r]
        temp2 = genTransposePoly(temp1)
        temp3 = genProdPoly(temp1, temp2)
        firstRow = genSumPoly(temp3, firstRow)

    HHT = genCirculant(firstRow)
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

# input:
#   H: Parity-check matrix
#   c: word to be decoded
#   N: cutoff for number of bit-flipping iterations
# output:
#   if decoding is successful, return the decoded word
#   else return 0
def bitFlipping(H, c, N):
    print("\nStarting Bit-Flipping Algorithm...")
    rows, cols = H.shape
    
    #The Tanner graph of the parity-check matrix is represented using adjancency lists
    bitsAdjList = [[] for i in range(cols)]
    checksAdjList = [[] for i in range(rows)]
    
    #bit nodes and check nodes, and initialise other parameters
    bits = np.copy(c)
    checks = np.zeros(rows)
    t = 0           #no. of rounds
    numOnes = 0     #count no. of ones in the check nodes, if numOnes = 0, the codeword is decoded
    flipOrNot = 0   #will use majority voting to decide if a bit is to be flipped
    checkSum = 0    #just an intermediate variable for computation
    
    for i in range(rows):
        for j in range(cols):
            if H[i, j] == 1:
                bitsAdjList[j].append(i)
                checksAdjList[i].append(j)
    
    while t < N:
        #print("Bit-Flipping Decoding Round", t, ":")
        for i in range(rows):
            for j in checksAdjList[i]:
                checkSum += bits[j]
            checks[i] = checkSum % 2
            checkSum = 0
            if checks[i] == 1:
                numOnes += 1

        if numOnes == 0:
            break
        else:
            numOnes = 0
        
        for i in range(cols):
            for j in bitsAdjList[i]:
                if checks[j] == 1:
                    flipOrNot += 1
            if 2*flipOrNot > len(bitsAdjList[i]):
                bits[i] = (bits[i] + 1) % 2
            flipOrNot = 0

        t += 1
                
    if t < N:
        print("Decoded in", t, "step(s).")
        
        return bits
    else:
        print("Cannot decode")
    
    return 0

# input: 
#   H: Parity-check matrix
#   y: word to be decoded
#   N: cutoff for the number of sum-product iterations
#   p: probability of a bit being digit 0
# output:
#   if decoding is successful, return the decoded word
#   else return 0

def sumProduct(H, y, N, p):
    print("\nStarting Sum-Product Algorithm...")
    rows, cols = H.shape
    
    #initialise variables
    t = 0
    checkSum = 0
    numOnes = 0
    prodTerm = 1

    #initialise arrays and matrices
    bitsAdjList = [[] for i in range(cols)]
    checksAdjList = [[] for i in range(rows)]
    estimate = np.zeros(2)
    bits = np.zeros(rows, dtype = np.int32)
    checks = np.zeros(rows, dtype = np.int32)
    M = np.zeros((rows, cols))
    E = np.zeros((rows, cols))
    r = np.zeros(cols)
    L = np.zeros(cols)
    x = np.zeros(cols, dtype = np.int32)
    
    estimate = [math.log((1-p)/p), math.log(p/(1-p))]   #L(C_j|Y_j=y_j)
    
    #setup the Tanner-graph with adjancency lists
    for i in range(rows):
        for j in range(cols):
            if H[i, j] == 1:
                bitsAdjList[j].append(i)
                checksAdjList[i].append(j)

    #Step 1: check if y is a codeword with 1 round of bit-flipping
    for i in range(rows):
        for j in checksAdjList[i]:
            checkSum += y[j]
        if (checkSum % 2) == 1:
            numOnes += 1
        checkSum = 0
        
    if numOnes == 0:
        return y
        
    checkSum = 0
    numOnes = 0

    #Step 2
    for j in range(cols):
        r[j] = estimate[int(y[j])]
        for i in bitsAdjList[j]:
            M[i,j] = r[j]
    while (t < N):
        t += 1
        
        #Step 3: compute E[i,j]

        #print("M:", "\n", M, "\n")
        for i in range(rows):
            prodTerm = 1
            for j in checksAdjList[i]:
                prodTerm *= np.tanh(0.5*M[i,j])
            for j in checksAdjList[i]:
                divisor = np.tanh(0.5*M[i,j])
                try:
                    E[i,j] = math.log((1+(prodTerm/divisor))/(1-(prodTerm/divisor)))
                except ValueError:
                    print("Cannot compute E[i,j]")
                    return -1
        #print("E:", "\n", E, "\n")
        
        #Step 4:
        L = np.copy(r)
        #L = L * (N-t) / N
            
        for j in range(cols):
            for s in bitsAdjList[j]:
                L[j] += E[s,j]
            if L[j] > 0:
                x[j] = 0
            elif L[j] < 0:
                x[j] = 1
        
        #Step 4: check if x is a codeword with 1 round of bit-flipping    
        checkSum = 0
        numOnes = 0
        for i in range(rows):
            for j in checksAdjList[i]:
                checkSum += x[j]
            if (checkSum % 2) == 1:
                numOnes += 1
            checkSum = 0
        if numOnes == 0:
            print("Decoded in", t, "step(s).")
            return x

        #Step 5:
        for j in range(cols):
            for i in bitsAdjList[j]:
                M[i,j] = L[j] - E[i,j]

    print("Cannot decode")
    return 0


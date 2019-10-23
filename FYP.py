import csv
import math
from sympy import *
import numpy as np
from scipy import special
from scipy.stats import binom
import networkx as nx
from networkx.algorithms import bipartite
import matplotlib.pyplot as plt
import random

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

# input: integer array v 
# output: integer array of v modulo 2
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
    wi = int(w / n0)
    
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

    if (n0 == 1):
        return G
    else :
        #concatenate identity matrix of order (n - r) and the stack of circulant matrices to form G
        G = np.concatenate((G, block0), axis = 1)
    
    return G

# input: an (n,r,w)-QC-MDPC matrix, H
# output: the number of 4 cycles in the Tanner graph of H
def count4Cycles(H):
    r, n = H.shape
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
    
    print("First row of HHT:\n", firstRow)
    
    for i in range(1,r):
        if (firstRow[i] >= 2):
            num4Cycles += (r-i) * special.binom(int(firstRow[i]), 2)
                
    print("Number of 4 cycles in parity-check matrix H:", int(num4Cycles))
    
    return int(num4Cycles)


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
#   H: Parity-check matrix (not necessarily QC-MDPC)
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
#   H: Parity-check matrix (not necessarily QC-MDPC)
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
        print("L:\n", L)
        
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

# input:
#   k: dimension of vector
#   t: Hamming weight of vector
# output: random vector of Hamming weight t and dimension k
def genRandomVector(k, t):
    randomVector = np.zeros(k, dtype = np.int32)

    # pick t random positions out of k positions
    randomPositions = np.random.choice(k, t, replace=False)

    # assign the random positions the value 1
    for j in randomPositions:
        randomVector[j] = 1
    
    return randomVector

# input:
#   G: Generator Matrix of a QC-MDPC matrix
#   m: Plaintext
#   e: error vector
# output:
#   ciphertext: encrypted message
def encryptMcEliece(G, m, e):
    rows, cols = G.shape
    n = cols
    r = cols - rows
    
    #encryption follows research paper, ciphertext = m*G + e
    ciphertext = np.copy(np.add(np.matmul(m, G), e))
    
    ciphertext = convertBinary(ciphertext)

    #Uncomment this to check if the encryption is correct
    #print("Plaintext : ", convertBinary(np.matmul(m, G)))
    #print("Error     : ", randomError)
    #print("Ciphertext: ", ciphertext)
    
    return ciphertext

# input: 
#   H: QC-MDPC matrix
#   y: ciphertext
#   method: either 'BF' or 'SP', representing Bit-Flipping and Sum-Product resp.
#   N: cutoff for no. of decoding iterations
#   p: probability of error (only for method = 'SP'. If method = 'BF',
#       it doesn't matter what value p is
# output:
#   decryptedText: decrypted text
#       (decrypted text can just be an integer if decryption fails)
#       bitFlipping return 0 if decoding fails due to exceeding max iteration
#       SumProduct returns 0 if decoding fails due to exceeding max iteration
#       SumProduct returns -1 if decoding fails due to E[i,j] computation error
def decryptMcEliece(H, y, method, N, p):
    r, n = H.shape

    if (method == 'BF'):
        #Decryption
        decryptedText = bitFlipping(H, y, N)
        
        if type(decryptedText) == int:
            print("Cannot decode by Bit-Flipping algorithm")
        else :
            decryptedText = decryptedText[0: n - r]
            
    elif (method == 'SP'):
        #Decryption
        decryptedText = sumProduct(H, y, N, p)
        
        if type(decryptedText) == int:
            print("Cannot decode by Sum-Product algorithm")
        else :
            decryptedText = decryptedText[0: n - r]
            
    return decryptedText

# input: 
#   plaintext: the original message
#   decryptedText: decrypted text
# ouput: returns true if (plaintext == decryptedText) element-wise, returns false otherwise
def decryptSuccess(plaintext, decryptedText):
    status = np.array_equal(plaintext, decryptedText)
    if (status == True):
        print("Decryption success!")
    else:
        print("Decryption failure!")
        
    return status

# input: vector v
# output: distance spectrum of v and distance spectrum with multiplicity of v
def genDistSpec(v):
    maxDist = len(v) // 2
    distSpec = []
    
    vT = genTransposePoly(v)
    temp = genProdPoly(v, vT)

    #first element of temp is always w(v) ( Hamming weight of v) which is not relevant
    #distance greater than maxDist have already been accounted for, cyclically
    distSpecMult = temp[1: maxDist+1]

    #if v is of even length, then there is double-counting for distance of length [len(v) / 2]
    if (len(v) % 2 == 0):
        distSpecMult[len(distSpecMult) - 1] /= 2

    for i in range( len(distSpecMult) ):
        if (distSpecMult[i] != 0):
            distSpec.append(i+1)
    
    return distSpec, distSpecMult

# input: code parameters n, r
# output: the first row of an (n,r)-QCMDPC code with no 4 cycles
# assumptions: n / r is a positive integer
def genNo4Cycles(n, r):
    n0 = int(n / r)
    pos = []
    w = []
    error = False
    # Hall's polynomial of resultant (n,r,w)-QC code
    h = np.zeros(n, dtype = np.int32)
    
    for i in range(n0):
        # list to keep track of the available positions for bits to placed
        pos.append([i for i in range(r)])
        # list to keep track of the weight of each circulant block
        w.append(0)

        # pick a position in each block
        randPos = random.choice(pos[i])
        h[randPos + i * r] = 1
        pos[i].remove(randPos)
        w[i] += 1

    print("w:", w)
        
    while(True):
        h_prev = np.copy(h)
        for i in range(n0):
            # will add 2 bits at once to ensure in each block has odd weight
            for j in range(2):
                randPos = random.choice(pos[i])
                h[randPos + i * r] = 1
                pos[i].remove(randPos)
                w[i] += 1
            print("w:", w)
            # check if the QC-code has 4-cycles
            H = genCirculant(h)
            if (count4Cycles(H) > 0):
                w[i] -=2
                return h, sum(w)

# input: code parameters n, r, w
# output: the first row of an (n,r,w)-QCMDPC code with no 4 cycles
#           if such a code is attainable in 100 random trials, else return False
# assumptions: w must be odd
def genQCNo4(n,r,w):
    for i in range(100):
        h, weight = genNo4Cycles(n,r)
        if (weight >= w):
            break

    if (weight < w):
        print("(n,r,w)-QC-code without 4 cycles is not obtainable in 100 trials.")
        return 0
    
    return h
# input:
#   H: (n,r,w)-QCMDPC matrix
#   G: generator matrix of H
#   method: 'BF' or 'SP'
#   N: number of decoding iterations
#   p: probability of success
#output: the number of errors the decoding method can correct using H
def decodeMax(H, G, method, N, p):
    r, n = H.shape
    d = r // 2
    t = 0
    decryptedText = []

    while (True):
        # generate a random message m of weight d
        m = genRandomVector(r, d)

        # generate a random error vector e of weight t
        e = genRandomVector(n, t)

        # encrypt the message m
        y = encryptMcEliece(G, m, e)

        # decrypt the ciphertext
        decryptedText = decryptMcEliece(H, y, method, N, p)
        
        if (type(decryptedText) == int):
            print(method, "can correct", t-1, "errors")
            return t
        else:
            print(method, "can correct", t, "errors")
        t += 1
        
def XBar(S,t,n,w):
    numer = 0
    denom = 0

    for i in range(1, t+1, 2):
        rho = rhoL(n, w, t, i)
        numer += (i - 1) * rho
        denom += rho
 
    return S * numer / denom

def rhoL(n, w, t, ell):
    return special.binom(w, ell) * special.binom(n - w, t - ell) / special.binom(n, t)
    
def counterDist(S, XBar, n, w, d, t):
    pi1prime = (S + XBar) / (d * t)
    pi0prime = ((w - 1) * S - XBar) / (d * (n - t))

    return pi1prime, pi0prime
    
def threshold(d, pi1prime, pi0prime, n, t):
    '''
    Input: d, pi1prime, pi0prime, n, t
    Output: Threhold required for SBSBF
    '''
    T = 0
    while(t * binom.pmf(k=T, n=d, p=pi1prime) <
          (n - t) * binom.pmf(k=T, n=d, p=pi0prime)):
        T += 1
        
    return T

def sampling(H, y):
    s = convertBinary(np.matmul(H, y))
    print("s:", s)
    # extract nonzero entries of s
    unsatEqn = [idx for idx, s_j in enumerate(s) if s_j == 1]
    # pick a random unsatisfied eqn
    i = random.choice(unsatEqn)

    print("i:",i)
    
    # extract nonzero entries of ith row of H
    ones = [bit for bit, h_ij in enumerate(H[i,:]) if h_ij == 1]
    # pick a random index of nonzero entry of ith row of H
    j = random.choice(ones)

    return j

def SBSBF(H, y, w, t):
    '''
    Step-by-Step Bit Flipping algorithm
    Input: Parity-check matrix H, ciphertext y
    output: decrypted text y'
    '''
    r, n = H.shape
    d = w / 2

    print(np.transpose(H).shape)
    print((y[np.newaxis, :]).shape)
    s = np.matmul(y, np.transpose(H))
    print("s:", s)
    s = convertBinary(s)
    print("s:", s)
    
    while (np.count_nonzero(s) > 0):
        # syndrome weight
        S = sum(s==1)
        print("S:", S)
            
        X_bar = XBar(S,t,n,w)
        pi1prime, pi0prime = counterDist(S, X_bar, n, w, d, t)
        tau = threshold(d, pi1prime, pi0prime, n, t)
        tau = (d + 1) / (2 * d)
        j = sampling(H, y)

        print("j:", j)
        print("unsatified pc eqn involving j bit:", np.matmul(s, H[:,j]))
        if (np.matmul(s, H[:,j]) > tau * d):
            y[j] = y[j] ^ 1

        # syndrome
        s = np.matmul(y, np.transpose(H))
        s = convertBinary(s)

        print("syndrome:", s)

    print("Decrypted text:", y)        
    return y

def transProb(n, d, t, w, S, pi1prime, pi0prime, sigma, T):
    i = 0
    p_sigma_neg = t * sigma * binom.pmf(sigma, d, pi1prime) / (w * S)
    p_sigma_pos = (n - t) * sigma * binom.pmf(sigma, d, pi0prime) / (w * S)
    
    while(i < T):
        p += t * i * binom.pmf(i, d, pi1prime) / (w * S)
        p += (n - t) * i * binom.pmf(i, d, pi0prime) / (w * S)

    return p_sigma_neg, p_sigma_pos, p

def newTransProb(d, pi1prime, pi0prime, p_sigma_neg, p_sigma_pos, p, T):
    pL = 0
    prob1 = 0
    prob2 = 0

    for i in range(T):
        prob1 += binom.pmf(i, d, pi1prime)
        prob2 += binom.pmf(i, d, pi0prime)

    pL = (prob1 ** t) * (prob2 ** (n - t))
    return p_sigma_neg * (1 - pL) / (1 - p), p_sigma_pos * (1 - pL) / (1 - p), pL

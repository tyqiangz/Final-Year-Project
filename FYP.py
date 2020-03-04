import csv
import math
from sympy import *
import numpy as np
from scipy import special
from scipy.stats import binom, hypergeom
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

    w = sum(H[0, :])
    
    filename = str(n) + "_" + str(r) + "_" + str(w) + "_" + "GeneratorMatrix.csv"

    np.savetxt(filename, G, delimiter = ",", fmt = "%d")
    
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

def BF(H, y, N):
    '''
    H: parity-check matrix
    y: word to be decoded
    N: max no. of decoding iterations
    '''
    r, n = H.shape
    d = sum(H[0,:] == 1) // 2
    i = 0
    s = convertBinary(np.matmul(y, np.transpose(H)))

    while ( (sum(s==1) != 0) and i < N ):
        s = convertBinary(np.matmul(y, np.transpose(H)))
        for j in range(n):
            sigma_j = np.matmul(s, H[:,j])
            if (sigma_j >= 0.5 * d):
                y[j] = (1 - y[j]) % 2
        i += 1

    s = convertBinary(np.matmul(y, np.transpose(H)))
    if (sum(s==1) == 0):
        return y
    else:
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
        #print("L:\n", L)
        
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
#   method: 'BF' or 'SP' or  'SBSBF'
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
########################## Step-By-Step Bit-Flipping Decoder ##########################

# input: parity-check matrix H, syndrome s, threshold T, method:
# method 1: random select an unsatisfied eqn, then select a nonzero bit of that eqn
# method 2: random selection amongst j where sigma_j >= T
# method 3: return argmax_j sigma_j
# output: a bit j such that sigma_j >= T, else return 'F'
def sampling(H, s, T, method):
    r, n = H.shape

    # method 1: random select an unsatisfied eqn, then select a nonzero bit of that eqn
    # loop until a nonzero bit j where sigma_j >= T is found
    # bound the loop by n iterations, returns F if such a bit cannot be found
    if (method == 1):
        for q in range(n):
            # extract nonzero entries of s
            unsatEqn = [idx for idx, s_j in enumerate(s) if s_j == 1]
            
            # pick a random unsatisfied eqn
            i = random.choice(unsatEqn)
            
            # extract nonzero entries of ith row of H
            ones = [bit for bit, h_ij in enumerate(H[i,:]) if h_ij == 1]
            
            # pick a random index of nonzero entry of ith row of H
            j = random.choice(ones)

            if ( sum((s + H[:, j]) == 2) >= T ):
                return j

        return 'F'
    
    # method 2: random selection amongst j where sigma_j >= T
    if (method == 2):
        # compute all the sigma_j for all j = 0, 1, ..., n-1
        sigmaJ = np.matmul(s, H)
        
        # if sigma_j >= threshold T, then add to a list of bits to be flipped
        toFlip = []
        for k in range(n):
            if (sigmaJ[k] >= T):
                toFlip += [k]

        # randomly pick a bit to flip if it exists
        if (len(toFlip) > 0):
            return random.choice(toFlip)
        else:
            return 'F'
        
    # method 3: return argmax_j sigma_j
    # in the case where there are many positions k = argmax_j sigma_j, randomly return any of them
    if (method == 3):
        # compute all the sigma_j for all j = 0, 1, ..., n-1
        sigmaJ = np.matmul(s, H)
        j = np.argmax(sigmaJ)

        if (sigmaJ[j] < T):
            return 'F'
        else:
            maxIndices = []
            for i in range(len(sigmaJ)):
                if (sigmaJ[i] == sigmaJ[j]):
                    maxIndices += [i]

            return random.choice(maxIndices)
        
def SBSBF(H, y, w, t, N, codeword, samp_method):
    '''
    Step-by-Step Bit Flipping algorithm
    Input: Parity-check matrix H, ciphertext y, weight of each parity-check eqn w,
    number of errors t, number of iterations N
    output: decrypted text y'
    '''
    print("Starting Step-by-Step Bit-Flipping Algorithm...")
    
    r, n = H.shape
    d = w // 2
    iteration = 1
    flipped = 1

    s = convertBinary(np.matmul(y, np.transpose(H)))
    
    while (np.count_nonzero(s) > 0 and iteration <= N and flipped == 1):
        flipped = 0
        iteration += 1
        # syndrome weight
        S = sum(s==1)

        t = sum( convertBinary(np.array(codeword) + np.array(y)) == 1)
        X_bar = XBar(S,t,n,w)
        pi1prime, pi0prime = counterDist(S, X_bar, n, w, d, t)
            
        T = threshold(d, pi1prime, pi0prime, n, t)
        T = max(floor(d / 2) + 1, T)

        j = sampling(H, s, T, samp_method)
        # if no such bit can be sampled, then exit algorithm, decoding failure
        if (j == 'F'):
            print("Cannot sample a bit")
            return 0
        else:
            y[j] = y[j] ^ 1
            flipped = 1
            
        # syndrome
        s = np.matmul(y, np.transpose(H))
        s = convertBinary(s)

    print("Decrypted text:\n", y)
    print("Codeword:\n", codeword)
    if (sum(s == 1) == 0):
        return y
    else:
        print("Cannot decode")
        return 0

#################### Markov Chain Monte Carlo Simulation Algorithms ####################

def XBar(S,t,n,w):
    numer = 0
    denom = 0

    for i in range(1, t+1, 2):
        rho = rhoL(n, w, t, i)
        if (not(np.isnan(rho)) and not(np.isinf(rho))):
            numer += (i - 1) * rho
            denom += rho
    
    return S * numer / denom

def rhoL(n, w, t, ell):
    return hypergeom.pmf(ell, n, w, t)

def rhoBar(n, w, t):
    temp = 0
    for ell in range(1,t+1,2):
        temp += rhoL(n,w,t,ell)
    return temp
        
def g_0(n, w, t, rhobar, k):
    if (k % 2 == 1):
        return 0
    else:
        return rhoL(n, w, t, k) / (1- rhobar)

def g_1(n, w, t, rhobar, k):
    if (k % 2 == 0):
        return 0
    else:
        return rhoL(n, w, t, k) / rhobar

# return h(l) for l = 0,1,...,r
def convolveH(n, w, t, rhobar):
    d = w // 2
    r = n // 2

    G0 = [g_0(n, w, t, rhobar, i) for i in range(w+1)]
    G1 = [g_1(n, w, t, rhobar, i) for i in range(w+1)]

    G0conv = [[1]]
    G1conv = [[1]]

    for i in range(r):
        G0conv += [np.convolve(G0conv[i], G0)]
        G1conv += [np.convolve(G1conv[i], G1)]

    # h(0) = g_0^{*r}(d * t)
    h = [G0conv[r][d * t]]

    # compute h(ell) for ell = 1, ..., r-1
    for i in range(1, r):
        temp = 0
        j = 0
        while (j >= 0 and j < len(G1conv[i]) and d * t - j >=0 and
               d * t - j < len(G0conv[r - i])):
            temp += G1conv[i][j] * G0conv[r - i][d * t - j]
            j += 1
        h += [temp]
        
    # h(r) = g_1^{*r}(d * t)
    h += [ G1conv[r][d * t] ]

    return h

# return a list of Pr(S = ell | H is regular) for ell = 0,1,...,r
def syndromeDist(n, w, t, rhobar):
    r = n // 2
    d = w // 2

    # if t = 1, Pr(S = d | H is regular) = 1 and 0 elsewhere 
    if (t == 1):
        return [0 for i in range(d)] + [1] + [0 for i in range(r - d)]
    
    h = convolveH(n, w, t, rhobar)
    sDist = []

    for ell in range(r+1):
        denom = 0
        for k in range(r):
            denom += binom.pmf(k, r, rhobar) * h[k]
        sDist += [binom.pmf(ell, r, rhobar) * h[ell] / denom]

    return sDist
    
def counterDist(S, XBar, n, w, d, t):
    pi1prime = (S + XBar) / (d * t)
    pi0prime = ((w - 1) * S - XBar) / (d * (n - t))

    return pi1prime, pi0prime
    
def threshold(d, pi1prime, pi0prime, n, t):
    '''
    Input: d, pi1prime, pi0prime, n, t
    Output: Threshold required for SBSBF
    '''
    if (pi1prime > 1):
        pi1prime = 1
        
    if (pi1prime < 0 or pi0prime > 1 or pi0prime < 0):
        print("Probabilities must be within 0 and 1")
        print("pi1prime: %f, pi0prime: %f" % (pi1prime, pi0prime))
        return 'F'
    
    if (pi1prime == 1):
        if ( t >= (n-t) * (pi0prime ** d) ):
            return d
        else:
            return 'F'
        
    numer = math.log((n - t) / t) + d * math.log((1 - pi0prime) / (1 - pi1prime))
    denom = math.log(pi1prime / pi0prime) + math.log((1 - pi0prime) / (1 - pi1prime))
    return math.ceil(numer / denom)

def p_sigmas(n, d, t, w, S, pi1prime, pi0prime, sigma):
        
    p_sigma_neg = t * sigma * binom.pmf(sigma, d, pi1prime) / (w * S)
    p_sigma_pos = (n - t) * sigma * binom.pmf(sigma, d, pi0prime) / (w * S)        
        
    return p_sigma_neg, p_sigma_pos

def q_sigmas(n, d, t, w, S, pi1prime, pi0prime, sigma):
        
    q_sigma_neg = t * binom.pmf(sigma, d, pi1prime) / n
    q_sigma_pos = (n - t) * binom.pmf(sigma, d, pi0prime) / n
        
    return q_sigma_neg, q_sigma_pos

def calcP(T, n, d, t, w, S, pi1prime, pi0prime):
    p = 0
    for i in range(T):
        sigma = i
        p_sigma_neg, p_sigma_pos = p_sigmas(n, d, t, w, S, pi1prime, pi0prime, sigma)
        p += p_sigma_neg + p_sigma_pos

    return p

def calcQ(T, n, d, t, w, S, pi1prime, pi0prime):
    q = 0
    for i in range(T):
        sigma = i
        q_sigma_neg, q_sigma_pos = q_sigmas(n, d, t, w, S, pi1prime, pi0prime, sigma)
        q += q_sigma_neg + q_sigma_pos

    return q

def pL(d, pi1prime, pi0prime, p, T, t, n):
    prob1 = 0
    prob2 = 0

    for i in range(T):
        prob1 += binom.pmf(i, d, pi1prime)
        prob2 += binom.pmf(i, d, pi0prime)

    pL = (prob1 ** t) * (prob2 ** (n - t))
    return pL

def p_sigmas_prime(p_sigma_neg, p_sigma_pos, pL, p):
    return p_sigma_neg * (1 - pL) / (1 - p), p_sigma_pos * (1 - pL) / (1 - p)

def q_sigmas_prime(q_sigma_neg, q_sigma_pos, pL, q):
    return q_sigma_neg * (1 - pL) / (1 - q), q_sigma_pos * (1 - pL) / (1 - q)

def q_maxs(n, d, t, pi1prime, pi0prime, sigma):
    temp1 = 0
    temp0 = 0

    for x in range(sigma):
        temp1 += binom.pmf(x, d, pi1prime)
        temp0 += binom.pmf(x, d, pi0prime)

    temp1_new = temp1 + binom.pmf(sigma, d, pi1prime)
    temp0_new = temp0 + binom.pmf(sigma, d, pi0prime)
    
    prob_max = (temp1_new ** t) * (temp0_new ** (n-t)) - (temp1 ** t) * (temp0 ** (n-t))

    denom = t * binom.pmf(sigma, d, pi1prime) + (n - t) * binom.pmf(sigma, d, pi0prime)

    q_max_minus = t * binom.pmf(sigma, d, pi1prime) / denom * prob_max
    q_max_plus = (n - t) * binom.pmf(sigma, d, pi0prime) / denom * prob_max
    
    return q_max_plus, q_max_minus

############################# DFR Algorithms #############################
def DFR_exp(H, G, w, t, N, trials, method, samp_method):
    '''
    H: parity-check matrix
    G: generator matrix
    w: row-weight of H
    t: number of errors
    N: max no. of iterations for decoder
    trials: no. of decoding trials
    method: decoding method 'BF' or 'SP' or 'SBSBF'
    samp_method: sampling method of a bit to flip
    '''
    r, n = H.shape
    d = w // 2

    DFR = 0
    
    for k in range(trials):
        print("Trial ", k)
        # generate a random message m of weight between 1 and r
        m = genRandomVector(r, random.randint(1, r))

        # generate a random error vector e of weight t
        e = genRandomVector(n, t)

        # encrypt the message m
        y = encryptMcEliece(G, m, e)

        codeword = convertBinary(np.array(y) + np.array(e))

        if (method == 'SBSBF'):
            # decrypt the ciphertext
            decryptedText = SBSBF(H, y, w, t, N, codeword, samp_method)
        elif (method == 'BF'):
            # decrypt the ciphertext
            decryptedText = BF(H, y, N=20)
                    
        # check if decryption is correct
        if (type(decryptedText) == int):
            status = decryptSuccess(m, decryptedText)
        else:
            status = decryptSuccess(m, decryptedText[0: n-r])

        if (status == False):
            DFR += 1

    print("Failed:", DFR, "Trials", trials, "DFR:", DFR / trials)

    return(DFR / trials)

def DFR_model(t_pass, t_fail, n, w, t_init, prob, samp_method):
    '''
    Input:
    t_pass: integer t such that SBSBF decoder always decodes for given code parameters
    t_fail: integer t such that SBSBF decoder always fail to decode for given code parameters
    n: codelength
    w: row-weight
    t_init: initial no. of errors
    prob: distribution of initial syndrome weight
    samp_method: sampling method of a bit to flip
    '''
    d = w // 2
    r = n // 2
    maxS = r
    DFR = {}

    while (prob[maxS] < 10 ** (-30)):
        maxS -= 1

    # syndrome weight and t_init must have same parity (odd/even)
    if ( (maxS + t_init) % 2 == 1 ):
        maxS -= 1
    
    print("maxS:", maxS)

    minS = int(np.floor(d / 2)) + 1
    print("minS:", minS)
    
    # initialisation of DFR
    for t in range(0, t_fail + 1):
        DFR[(0,t)] = 0
        
    for S in range(1, minS):
        DFR[(S,0)] = 0
        for t in range(1, t_fail):
            DFR[(S,t)] = 1
            
    for S in range(1, r + 1):
        DFR[(S, t_fail)] = 1
        
    for S in range(minS, r + 1):
        DFR[(S, t_pass)] = 0
    
    #print("DFR:\n", DFR)
    
    # computation in ascending manner
    for S in range(minS, maxS + 1):
        for t in range(t_pass + 1, t_fail):
            X_bar = XBar(S,t,n,w)
            pi1prime, pi0prime = counterDist(S, X_bar, n, w, d, t)

            T = threshold(d, pi1prime, pi0prime, n, t)

            # if a threshold can't be found, set as ceil((d + 1) / 2)
            if (T == 'F'):
                T = int(minS)

            T = max(minS, T)
            print("T:", T)
            
            p = calcP(T, n, d, t, w, S, pi1prime, pi0prime)
            PL = pL(d, pi1prime, pi0prime, p, T, t, n)
            q = calcQ(T, n, d, t, w, S, pi1prime, pi0prime)
            if (samp_method == 1):
                DFR[(S,t)] = PL
            elif (samp_method == 2):
                DFR[(S,t)] = PL
            elif (samp_method == 3):
                DFR[(S,t)] = pL(d, pi1prime, pi0prime, p, minS, t, n)
            
            if (samp_method == 1):
                for sigma in range(T, min(d + 1, S + 1)):
                    print("(S,t,sigma) = (%d,%d,%d)" % (S,t,sigma))
                    p_sigma_neg, p_sigma_pos = p_sigmas(n, d, t, w, S, pi1prime, pi0prime, sigma)
                    p_sigma_neg_prime, p_sigma_pos_prime = p_sigmas_prime(p_sigma_neg, p_sigma_pos, PL, p)
            
                    DFR[(S,t)] += p_sigma_neg_prime * DFR[S + d - 2 * sigma, t - 1] + p_sigma_pos_prime * DFR[S + d - 2 * sigma, t + 1]
            
            elif (samp_method == 2):
                for sigma in range(T, min(d + 1, S + 1)):
                    print("(S,t,sigma) = (%d,%d,%d)" % (S,t,sigma))
                    q_sigma_neg, q_sigma_pos = q_sigmas(n, d, t, w, S, pi1prime, pi0prime, sigma)
                    q_sigma_neg_prime, q_sigma_pos_prime = q_sigmas_prime(q_sigma_neg, q_sigma_pos, PL, q)

                    DFR[(S,t)] += q_sigma_neg_prime * DFR[S + d - 2 * sigma, t - 1] + q_sigma_pos_prime * DFR[S + d - 2 * sigma, t + 1]
                    
            elif (samp_method == 3):
                for sigma in range(minS, min(d + 1, S + 1)):
                    print("(S,t,sigma) = (%d,%d,%d)" % (S,t,sigma))
                    q_max_plus, q_max_minus = q_maxs(n, d, t, pi1prime, pi0prime, sigma)
                    
                    DFR[(S,t)] += q_max_minus * DFR[S + d - 2 * sigma, t - 1] + q_max_plus * DFR[S + d - 2 * sigma, t + 1]
    
    fail = 0

    if (t_init % 2 == 0):
        startS = 2
    else:
        startS = 1

    for S in range(startS, maxS + 1, 2):
        #print("DFR[(%d,%d)]: %f" % (S, t_init, DFR[(S,t_init)]) )
        if (math.isnan(DFR[(S,t_init)]) or math.isnan(prob[S])):
            continue
        fail += prob[S] * DFR[(S,t_init)]

    return(fail)

############################# BIKE Algorithms #############################
def BIKE1(H, s, u):
    '''
    Input: Parity-check matrix H, syndrome s, upper bound u
    Output: Error vector e
    '''
    S = sum(s==1)
    w = sum(H[0,:] == 1)
    r, n = H.shape
    d = w // 2
    delta = 5
    
    X_bar = XBar(S,t,n,w)
    pi1prime, pi0prime = counterDist(S, X_bar, n, w, d, t)
    T = threshold(d, pi1prime, pi0prime, n, t)

    J = [[]]
    for i in range(T):
        J += [[]]
    
    for j in range(n):
        ell = min(sum(H[:,j] == s), T)
        J[ell].append(j)

    e = np.zeros(n, dtype = np.int32)
    for bits in J[T]:
        e[bits] = 1

    s_prime = np.array(s) - convertBinary(np.matmul(e, np.transpose(H)))

############################### FYP Demo ################################

def demo(H, y):
    '''
    Simple demo for Bit-Flipping algorithm
    Input: Parity-check matrix H, ciphertext y
    output: decrypted text y'

    Assumptions: only 1 bit of error was introduced, thus set threshold = 1
    and decoding will finish in 1 step.
    '''

    print("H:\n", H)
    r, n = H.shape
    w = sum(H[0,:] == 1)
    d = w // 2
    iteration = 1
    flipped = 1

    s = convertBinary(np.matmul(y, np.transpose(H)))
    print("\ny=mH^T + e:", y, "#ciphertext")
    print("Threshold T:", d)
    print("\n######### Starting the Bit-Flipping Algorithm... #########\n")

    print("s = yH^T:", s)
    print("sigmaJ:", np.matmul(s, np.transpose(H)))
    while (np.count_nonzero(s) > 0 and flipped == 1):
        flipped = 0
        # syndrome weight
        T = 1

        for j in range(n):
            if (sum((s + H[:,j]) == 2)) >= T * d:
                print("FLIPPED position %d" % j)
                y[j] = y[j] ^ 1
                print("y:", y)
                s = convertBinary(np.matmul(y, np.transpose(H)))
                print("s = yH^T:", s)
                flipped = 1
                
        
        iteration += 1
        
        # syndrome
        s = np.matmul(y, np.transpose(H))
        s = convertBinary(s)

    print("Decrypted text:\n", y)
    if (sum(s == 1) == 0):
        return y[0: n-r]
    else:
        print("Cannot decode")
        return 0

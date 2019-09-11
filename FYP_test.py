from FYP import *

##################################### Setting parameters ###################################

# code parameters
n0 = 2
r = 1019
wi = 21

#decryption parameters
method = 'SP'
N = 10
p = 0.01
t = 1

d = r // 2

################################## Processing parameters ###################################

n = n0 * r
w = wi * n0

##################################### Testing functions ####################################

# generate a random (n,r,w)-QC-MDPC matrix H
H = genQCMDPC(n, r, w)

# Generate the corresponding generator matrix G
G = genGenQCMDPC(H)

# count the number of 4 cycles in the Tanner graph of H
count4Cycles(H, n, r, w)

# generate the distance spectrum of h (first row of H)
a, b = genDistSpec(H[0, :])
print("Distance spectrum of h:\n", a)
print("Distance spectrum with multiplicity of h:\n", b)

# generate a random message m of weight d
m = genRandomVector(r, d)

# generate a random error vector e of weight t
e = genRandomVector(n, t)

# encrypt the message m
y = encryptMcEliece(G, m, e)

# decrypt the ciphertext
decryptedText = decryptMcEliece(H, y, method, N, p)

# check if decryption is correct
decryptSuccess(m, decryptedText)

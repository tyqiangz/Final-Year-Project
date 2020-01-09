from FYP import *

##################################### Setting parameters ###################################

# code parameters
n0 = 2
r = 5
wi = 3

N = 20
k = r // 2
#decryption parameters
t = 1

################################## Processing parameters ###################################

n = n0 * r
w = wi * n0

##################################### Testing functions ####################################

# generate a random (n,r,w)-QC-MDPC matrix H
H = genQCMDPC(n, r, w)

# Generate the corresponding generator matrix G
G = genGenQCMDPC(H)
print("G:\n", G)
# generate a random message m of weight k
m = genRandomVector(r, k)

print("\nPlaintext m:", m)

# generate a random error vector e of weight t
e = genRandomVector(n, t)

# encrypt the message m
y = encryptMcEliece(G, m, e)

# decrypt the ciphertext
decryptedText = demo(H, y)

# check if decryption is correct
decryptSuccess(m, decryptedText)




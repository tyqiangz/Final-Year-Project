from FYP import *

##################################### Setting parameters ###################################

# code parameters
n0 = 2
r = 1289
wi = 39

#decryption parameters
t = 27

d = r // 2

################################## Processing parameters ###################################

n = n0 * r
w = wi * n0

##################################### Testing functions ####################################

print('########## Generate keys ##########')
print("(r,d,t):", (r,wi,t))
# generate a random (n,r,w)-QC-MDPC matrix H
H = genQCMDPC(n, r, w)

# Generate the corresponding generator matrix G
G = genGenQCMDPC(H)

print('########## DFR ##########')
DFR_exp(H, G, w, t, N=100, trials=10000, method='SBSBF', samp_method=2)
print("(r,d,t):", (r,wi,t))

from FYP import *

##################################### Setting parameters ###################################

# code parameters
n0 = 2
r = 2003
wi = 39

#decryption parameters
method = 'BF'
N = 10
p = 0.01
t = 21

d = r // 2

################################## Processing parameters ###################################

n = n0 * r
w = wi * n0

##################################### Testing functions ####################################

print('########## Generate keys ##########')
# generate a random (n,r,w)-QC-MDPC matrix H
H = genQCMDPC(n, r, w)

# Generate the corresponding generator matrix G
G = genGenQCMDPC(H)

print('########## DFR ##########')
DFR_Exp(H, G, w, t, N=100, trials=1000, method='SBSBF')
print("(n,w,t):", (n,w,t))

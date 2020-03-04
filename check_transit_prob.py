#!/usr/bin/env python
# coding: utf-8

# In[11]:


from FYP import *
import pandas as pd
import scipy as sp
from scipy.stats import hypergeom
import csv
from datetime import *


# In[2]:


def genDist(e,H):
    '''
    Input: 
    e: error vector
    H: parity-check matrix
    
    Output:
    sigmaJ: [sigma_0, sigma_1, ...], where sigma_j = |s \cap h_j|, where h_j is the jth col of H
    pmf0: distribution of sigma as given by paper for e=0
    pmf1: distribution of sigma as given by paper for e=1
    '''
    r, n = H.shape
    w = sum(H[0,:])
    d = int(w / 2)
    
    s = convertBinary(np.matmul(H, e))
    sigmaJ = np.matmul(s, H)
    
    S = sum(s == 1)
    t = sum(e == 1)
    X_bar = XBar(S,t,n,w)
    pi1prime, pi0prime = counterDist(S, X_bar, n, w, d, t)

    x = sp.linspace(0,d,d+1)
    pmf0 = binom.pmf(x,d,pi0prime)
    pmf1 = binom.pmf(x,d,pi1prime)
    
    return x, sigmaJ, pmf0, pmf1
    
def flip(e, H, samp_method, NUM_FLIPS):
    '''
    Input:
    e: error
    H: parity-check matrix
    samp_method: sampling method
    FLIPS: number of correct flips
    
    Output:
    error vector e after FLIPS number of flips
    '''
    r, n = H.shape
    w = sum(H[0,:])
    d = int(w / 2)
    
    t_init = sum(e == 1)
    t = t_init
    s = convertBinary(np.matmul(e, np.transpose(H)))
    S = sum(s == 1)
    
    while (t_init - t < NUM_FLIPS):
        s = convertBinary(np.matmul(e, np.transpose(H)))
        S = sum(s == 1)

        X_bar = XBar(S,t,n,w)
        pi1prime, pi0prime = counterDist(S, X_bar, n, w, d, t)
        
        if (samp_method == 1 or samp_method == 2):
            T = threshold(d, pi1prime, pi0prime, n, t)
        elif (samp_method == 3):
            T = int(np.floor(d / 2) + 1)

        j = sampling(H, s, T, samp_method)
        
        if (j == 'F'):
            #print('Cannot sample a bit to flip')
            return 'F'
        
        sigmaJ = np.matmul(s, H)
        #print('sigmaJ[j]', sigmaJ[j])
        
        e[j] = e[j] ^ 1
        t = sum(e == 1)
    
    #print('(S,t) = (%d,%d)' % (S,t))

    return e


# In[52]:


##################################### Setting parameters ###################################

# code parameters
r = 1019
d = 39

#decryption parameters
t_init = 20

################################## Processing parameters ###################################

n = 2 * r
w = 2 * d

##################################### Testing functions ####################################
print("(r,d,t):", (r,d,t_init))

H = genQCMDPC(n, r, w)


# In[57]:


NUM_FLIPS = 2
NUM_TRIALS = 100000
samp_method = 3

sampledBits0 = {}
sampledBits1 = {}
failed = {}

for i in range(NUM_TRIALS):
    print("Trial",i,end='\r')
    e = genFirstRow(n, t_init)
    e = flip(e, H, samp_method, NUM_FLIPS)
    if (type(e) == str):
        continue
        
    t = sum(e==1)
    s = convertBinary(np.matmul(H,e))
    S = sum(s==1)
    
    X_bar = XBar(S,t,n,w)
    pi1prime, pi0prime = counterDist(S, X_bar, n, w, d, t)
    
    if (samp_method == 1 or samp_method == 2):
        T = threshold(d, pi1prime, pi0prime, n, t)
    elif (samp_method == 3):
        T = int(np.floor(d/2) + 1)
    
    j = sampling(H, s, T, samp_method)
    
    if (str(j) == 'F'):
        if (S in failed):
            failed[S] += 1
        else:
            failed[S] = 1
    else:
        sigmaj = np.matmul(s, H[:,j])
        if (e[j] == 0):
            if (S in sampledBits0):
                sampledBits0[S] += [sigmaj]
            else:
                sampledBits0[S] = [sigmaj]
        elif (e[j] == 1):
            if (S in sampledBits1):
                sampledBits1[S] += [sigmaj]
            else:
                sampledBits1[S] = [sigmaj]

print("\nfailed\n", failed)


# In[58]:


with open('Ft_init_NUMFLIPS_samp_method_' + str(t_init) + '_' + str(NUM_FLIPS) + '_' + str(samp_method) + '.csv', 'w') as f:
    f.write("%s, %s\n" % ("Syndrome Weight (S)", "sigmas"))
    for key in failed.keys():
        f.write("%s,%s\n" % (key, failed[key]))
        
with open('0t_init_NUMFLIPS_samp_method_' + str(t_init) + '_' + str(NUM_FLIPS) + '_' + str(samp_method) + '.csv', 'w') as f:
    f.write("%s, %s\n" % ("Syndrome Weight (S)", "sigmas"))
    for key in sampledBits0.keys():
        f.write("%s," % key)
        for entry in sampledBits0[key]:
            f.write("%s," % entry)
        f.write("\n")
        
with open('1t_init_NUMFLIPS_samp_method_' + str(t_init) + '_' + str(NUM_FLIPS) + '_' + str(samp_method) + '.csv', 'w') as f:
    f.write("%s, %s\n" % ("Syndrome Weight (S)", "sigmas"))
    for key in sampledBits1.keys():
        f.write("%s," % key)
        for entry in sampledBits1[key]:
            f.write("%s," % entry)
        f.write("\n")


# In[59]:


########### generate the model transition probabilities ###########

t = t_init - NUM_FLIPS
print(t)
S = 372
X_bar = XBar(S,t,n,w)
pi1prime, pi0prime = counterDist(S, X_bar, n, w, d, t)

if (samp_method == 1 or samp_method == 2):
    T = threshold(d, pi1prime, pi0prime, n, t)
elif (samp_method == 3):
    T = int(np.floor(d/2) + 1)

p = calcP(T, n, d, t, w, S, pi1prime, pi0prime)
# failure probability
PL = pL(d, pi1prime, pi0prime, p, T, t, n)
q = calcQ(T, n, d, t, w, S, pi1prime, pi0prime)

transProb0, transProb1 = [[k for k in range(T, d+1)],[]], [[k for k in range(T, d+1)],[]]

if (samp_method == 1):
    for sigma in range(T, d+1):
        p_sigma_neg, p_sigma_pos = p_sigmas(n, d, t, w, S, pi1prime, pi0prime, sigma)
        p_sigma_neg_prime, p_sigma_pos_prime = p_sigmas_prime(p_sigma_neg, p_sigma_pos, PL, p)
        
        transProb0[1] += [p_sigma_pos_prime]
        transProb1[1] += [p_sigma_neg_prime]
elif (samp_method == 2):
    for sigma in range(T, d+1):
        q_sigma_neg, q_sigma_pos = q_sigmas(n, d, t, w, S, pi1prime, pi0prime, sigma)
        q_sigma_neg_prime, q_sigma_pos_prime = q_sigmas_prime(q_sigma_neg, q_sigma_pos, PL, q)
        
        transProb0[1] += [q_sigma_pos_prime]
        transProb1[1] += [q_sigma_neg_prime]
    
elif (samp_method == 3):
    for sigma in range(T, d+1):
        q_max_plus, q_max_minus = q_maxs(n, d, t, pi1prime, pi0prime, sigma)
        
        transProb0[1] += [q_max_plus]
        transProb1[1] += [q_max_minus]
    
######################## print results ##################################

print("Model failure prob: %f" % PL)
print("Simulation failure prob: %f\n" % failProb)

print("Model trans prob, e=0")
print(transProb0)

print("Model trans prob, e=1")
print(transProb1)

print("Total prob:", sum(transProb1[1])+sum(transProb0[1])+PL)


# In[ ]:





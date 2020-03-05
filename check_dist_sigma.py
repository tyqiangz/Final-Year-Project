# This program is used to check accuracy of the counter distributions

from FYP import *
import pandas as pd
import scipy as sp
from scipy.stats import hypergeom

def plot1Dist(ax, x, S, d, t, sigmaJ, pmf, title):
    '''
    plots a graph comparing the simulation distribution of sigma and that of the model
    
    input: 
    ax: axis to be plotted
    x: domain of the probability distribution of sigma
    S: syndrome weight
    d: block weight
    t: error weight
    sigmaJ: distribution of sigma
    pmf: probability distribution of sigma
    title: desired title for the graph to be plotted
    '''
    freqTable = np.array(np.unique(sigmaJ, return_counts=True)).T
    simu = ax.plot(freqTable[:,0], freqTable[:,1], label='Simulation')
    dist = pmf * sum(freqTable[:,1])
    model = ax.plot(x, dist, label='Model')
    
    max = d
    min = 0
    
    if (S < d):
        max = S
    if (S > d*(t-1)):
        min = (S + d*(1-t)) / 2
    
    print("Model\n", [x,dist])
    
    # add description to the plot
    ax.legend(loc="best")
    ax.set_xlabel('sigma')
    ax.set_ylabel('Frequency')
    ax.set_title(title)
    
    print("Frequency Table\n", freqTable)

def plotAllDist(x, sigmaJ, e, pmf0, pmf1, FLIPS, S, d, t):
    '''
    plots both graphs comparing the simulation distribution of sigma and that of the model
    
    input: 
    sigmaJ: distribution of sigma
    e: error vector
    pmf0: probability distribution of sigma given e=0
    pmf1: probability distribution of sigma given e=1
    FLIPS: number of flips to be made 
    S: syndrome weight
    d: block weight
    t: error weight
    title: desired title for the graph to be plotted
    '''
    ONES = [i for i in range(e.size) if e[i] == 1]
    ZEROS = [j for j in range(e.size) if e[j] == 0]

    if (FLIPS == 0):
        title = 'Distribution of sigma before bit-flipping'
    elif (FLIPS > 0):
        title = 'Distribution of sigma after ' + str(FLIPS) + ' flips'

    # number of columns = 2, so that the plots are side by side
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))
    plot1Dist(axes[0], x, S, d, t, sigmaJ[ZEROS], pmf0, title + ' (e=0)')
    plot1Dist(axes[1], x, S, d, t, sigmaJ[ONES], pmf1, title + ' (e=1)')
    
    plt.tight_layout(pad=1, w_pad=1)
    
def genDist(e,H):
    '''
    Input: 
    e: error vector
    H: parity-check matrix
    Output:
    x: domain of sigmaJ
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
    
def flip(e, H, samp_method, FLIPS):
    '''
    Input:
    e: error
    H: parity-check matrix
    samp_method: sampling method
    FLIPS: number of correct flips
    
    Output:
    e: error vector after FLIPS number of flips
    '''
    r, n = H.shape
    w = sum(H[0,:])
    d = int(w / 2)
    
    t_init = sum(e == 1)
    t = t_init
    s = convertBinary(np.matmul(e, np.transpose(H)))
    S = sum(s == 1)
    
    while (t_init - t < FLIPS):
        if (t - t_init > 10):
            return 0
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
            print('Cannot sample a bit to flip')
            return 'F'
        
        sigmaJ = np.matmul(s, H)
        #print('sigmaJ[j]', sigmaJ[j])
        
        e[j] = e[j] ^ 1
        t = sum(e == 1)
    
    #print('(S,t) = (%d,%d)' % (S,t))

    return e

##################################### Setting parameters ###################################

# code parameters
r = 9857
d = 71

#decryption parameters
t = 134

################################## Processing parameters ###################################

n = 2 * r
w = 2 * d

##################################### Testing functions ####################################
print("(r,d,t):", (r,d,t))

H = genQCMDPC(n, r, w)

samp_method = 1
flips = [0,5,10,20,40,70,100]

for FLIPS in flips:
    # generate random error e, length n, weight t
    e = genRandomVector(n, t)
    e = flip(e, H, samp_method, FLIPS)
    x, sigmaJ, pmf0, pmf1 = genDist(e,H)
    s = convertBinary(np.matmul(H,e))
    S = sum(s==1)
    t = sum(e==1)
    
    plotAllDist(x, sigmaJ, e, pmf0, pmf1, FLIPS, S, d, t)
    
    PIC_NAME = 'figures/r_d_t_S_samp_method_FLIPS_' + str(r) +'_' + str(d) + '_' + str(t-FLIPS) + '_' + str(S) + '_' + str(samp_method) + '_' + str(FLIPS) + '.png'
    plt.savefig(PIC_NAME)

samp_method = 2
flips = [0,5,10,20,40,70,100]

for FLIPS in flips:
    # generate random error e, length n, weight t
    e = genRandomVector(n, t)
    e = flip(e, H, samp_method, FLIPS)
    x, sigmaJ, pmf0, pmf1 = genDist(e,H)
    plotAllDist(x, sigmaJ, e, pmf0, pmf1, FLIPS)

    s = convertBinary(np.matmul(H,e))
    S = sum(s==1)
    
    PIC_NAME = 'figures/r_d_t_S_samp_method_FLIPS_' + str(r) +'_' + str(d) + '_' + str(t-FLIPS) + '_' + str(S) + '_' + str(samp_method) + '_' + str(FLIPS) + '.png'
    plt.savefig(PIC_NAME)

samp_method = 3
flips = [0,5,10,20,40,70,100]

for FLIPS in flips:
    # generate random error e, length n, weight t
    e = genRandomVector(n, t)
    e = flip(e, H, samp_method, FLIPS)
    x, sigmaJ, pmf0, pmf1 = genDist(e,H)
    plotAllDist(x, sigmaJ, e, pmf0, pmf1, FLIPS)

    s = convertBinary(np.matmul(H,e))
    S = sum(s==1)
    
    PIC_NAME = 'figures/r_d_t_S_samp_method_FLIPS_' + str(r) +'_' + str(d) + '_' + str(t-FLIPS) + '_' + str(S) + '_' + str(samp_method) + '_' + str(FLIPS) + '.png'
    plt.savefig(PIC_NAME)

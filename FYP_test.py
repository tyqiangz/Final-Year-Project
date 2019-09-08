from FYP import *
import networkx as nx, numpy as np
from networkx.algorithms import bipartite
import matplotlib.pyplot as plt

n0 = 1
r = 5
wi = 3

n = n0 * r
w = wi * r

H1 = genParityCheck(n,r,w)
H2 = genParityCheck(n,r,w)

print(genInvPoly(H1[0,:]))


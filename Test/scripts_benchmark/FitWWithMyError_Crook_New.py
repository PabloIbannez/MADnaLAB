#!/usr/bin/env python3

# inputs: 
# 1 -> file containing list of data

import sys
import math
import numpy as np
import glob
import os
from itertools import islice
import re
import subprocess
import pandas as pd

if (len(sys.argv) != 1+1):
    print("Usage: 1 input") 
    print("# 1 -> file containing list of data in the format f, beta, dbeta")
    print("Output: kbeta +- dkbeta, beta0 +- dbeta0, from formula cos(beta) = (cos(beta0)/kbeta)*f + cos(beta0), with beta0 = beta(f=0)")
    sys.exit(0)

def stampa(*argv):
    s = ""
    for arg in argv:
        s += str(arg)+" "
    print(s)

def acos(x):
    if x>1:
        x = 1
    if x<-1:
        x=-1
    return np.arccos(x)

lista = sys.argv[1]
#ntot = int(sys.argv[2])
ntot = 5000

dati0 = np.array(pd.read_csv(lista, sep=' ', header=None))

flist = dati0[:,0]
betalist = dati0[:,1]
betaerr = dati0[:,2]

xlist = flist
n = len(xlist)
ylist = np.zeros(n)
x = np.sum(xlist)
xx = np.sum(xlist*xlist)

kbeta = 0
dkbeta = 0
beta0 = 0
dbeta0 = 0
for i in range(0,ntot):
    for j in range(0,len(betalist)):
        ylist[j] = np.cos(betalist[j] + np.random.normal(0,betaerr[j]))
    xy = np.sum(xlist*ylist)
    y = np.sum(ylist)
    num = n*xy-x*y
    den = n*xx - x*x
    a = num/den
    b = (y-a*x)/n
    cos0loc = b
    beta0loc = acos(cos0loc)
    kbetaloc = cos0loc/a
    beta0 += beta0loc
    dbeta0 += beta0loc*beta0loc
    kbeta += kbetaloc
    dkbeta += kbetaloc*kbetaloc
beta0 /= ntot
dbeta0 /= ntot
kbeta /= ntot
dkbeta /= ntot
dbeta0 = np.sqrt(dbeta0-beta0*beta0)
dkbeta = np.sqrt(dkbeta - kbeta*kbeta)

stampa(kbeta,dkbeta,beta0,dbeta0)

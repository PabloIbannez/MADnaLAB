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
    print("# 1 -> file containing list of data in the format f, y, dy")
    print("Output: a +- da, y0 +- dy0, from formula y = a*f + y0")
    sys.exit(0)

def stampa(*argv):
    s = ""
    for arg in argv:
        s += str(arg)+" "
    print(s)

lista = sys.argv[1]
#ntot = int(sys.argv[2])
ntot = 5000

dati0 = np.array(pd.read_csv(lista, sep=' ', header=None))

xlist = dati0[:,0]
ylist_ave = dati0[:,1]
yerr = dati0[:,2]

n = len(xlist)
ylist = np.zeros(n)
x = np.sum(xlist)
xx = np.sum(xlist*xlist)

a = 0
da = 0
y0 = 0
dy0 = 0
for i in range(0,ntot):
    for j in range(0,n):
        ylist[j] = ylist_ave[j] + np.random.normal(0,yerr[j])
    xy = np.sum(xlist*ylist)
    y = np.sum(ylist)
    num = n*xy-x*y
    den = n*xx - x*x
    a_loc = num/den
    y0_loc = (y-a_loc*x)/n
    a += a_loc
    da += a_loc*a_loc
    y0 += y0_loc
    dy0 += y0_loc*y0_loc
a /= ntot
da /= ntot
da = np.sqrt(da - a*a)
y0 /= ntot
dy0 /= ntot
dy0 = np.sqrt(dy0 - y0*y0)

stampa(a,da,y0,dy0)

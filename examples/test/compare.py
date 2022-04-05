import sys
import numpy as np

import matplotlib.pyplot as plt

path1 = sys.argv[1]
path2 = sys.argv[2]

toCompare=["hrise","htwist","major_groove_depth","major_groove_width","minor_groove_depth","minor_groove_width","SBBS"]

for c in toCompare:
    d1 = np.genfromtxt(path1+"/"+c)[:,2::]
    d2 = np.genfromtxt(path2+"/"+c)[:,2::]

    plt.plot(d1[:,1],'.')
    plt.plot(d2[:,1],'.')
    plt.title(c)

    plt.show()

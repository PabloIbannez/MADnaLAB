import os
import sys

import libconf

from MADnaLAB import *

sequences = {"AA":["CGCGAAAAAAAAAACGCG",10.8748],
             "AC":["CGCGACACACACACCGCG",10.9725],
             "AG":["CGCGAGAGAGAGAGCGCG",10.7969],
             "AT":["CGCGATATATATATCGCG",10.8929],
             "CG":["CGCGCGCGCGCGCGCGCG",10.8929],
             "GG":["CGCGGGGGGGGGGGCGCG",10.8748]}

forces = [1.0,5.0,10.0,15.0,20.0]

nSimPerSeq = 100

with open("test_template.cfg","r") as inpt:
    conf = libconf.load(inpt)

#conf['SEQUENCES']['generate'] = {}
#conf['SEQUENCES']['generate']['list'] = {sequences["AA"][0]:nSimPerSeq}
#conf['SIMULATION']['options']['debyeLenght'] = sequences["AA"][1]
#conf['SIMULATION']['modifications']['pulling']['constantForce']['pullingForce'] = 10.0
#genMADna(conf)

os.system("rm -f *.coord")
os.system("rm -f *.top  ")
os.system("rm -f output*")
os.system("rm -f options.dat")

conf['SEQUENCES']['generate'] = {}
for seqName in sequences:
    os.system("rm -f -r "+seqName)
    os.system("mkdir "+seqName)
    os.chdir(seqName)

    for f in forces:
        os.system("mkdir results"+str(f))
        os.chdir("results"+str(f))
    
        os.system("rm -f *.coord")
        os.system("rm -f *.top  ")
        os.system("rm -f output*")

        conf['SEQUENCES']['generate']['list'] = {sequences[seqName][0]:nSimPerSeq}
        conf['SIMULATION']['options']['debyeLenght'] = sequences[seqName][1]
        conf['SIMULATION']['modifications']['pulling']['constantForce']['pullingForce'] = f
        
        genMADna(conf)
        
        os.system("MADnaLAB")
        
        os.system("rm -f *.dat")
        os.system("rm -f *.coord")
        os.system("rm -f *.top  ")
        os.system("rm -f *.list")
        
        os.chdir("..")
    
    os.chdir("..")

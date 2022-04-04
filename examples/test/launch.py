import os

with open("MADnaTestSimulationSets.dat") as f:
    for line in f:
        os.system("../../bin/MADnaLAB "+line)

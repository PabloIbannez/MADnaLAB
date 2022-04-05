import sys

from MADnaLAB.utils import *
from MADnaLAB.sequencesGeneration import *
from MADnaLAB.simulationPool import *
from MADnaLAB.simulationSets import *

simSetMaxSize = 50000

forces = [5.0,10.0,25.0,50.0,100.0]

K = 100.0

seqNum = 100
seqLen = 1000

sequenceList = []
for seq in range(seqNum):
    alias = "seq_"+str(seq)
    seq   = "A"*seqLen

    sequenceList.append({"alias":alias,"seq":seq})

simulationPool = []

simId=0
for seq in sequenceList:
    for f in forces:
        simulationPool.append({'alias':seq['alias'],
                               'seq':seq['seq'],
                               'simId':simId,
                               'externalForce': [{'force': [0.0,0.0,f], 'basePair': -2}],
                               'constraintsPositionOfBeads': [{'K':[K,K,K], 'basePair': 2}],
                               'filePath':"Simulations/f"+str(f)
                               })
        simId+=1

simPoolSet = splitSimulationPoolAccordingToMaxNumberOfBasePairs(simulationPool,simSetMaxSize)

cfg = loadCONF("wlc.cfg")
generateSimulationSetsFromListOfSimulationPool(cfg,simPoolSet)

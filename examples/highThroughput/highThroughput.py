import sys

import random

import MADnaLAB.utils
import MADnaLAB.sequencesGeneration
import MADnaLAB.simulationPool
import MADnaLAB.simulationSets

seqListFile = "seq.list"

forces = [5.0,10.0,25.0,50.0,100.0]
forcesTorques = [[2.0,1.0],[2.0,50.0],[2.0,100.0]]

K = 100.0

nCopies = 2

seqLength  = 40
nSequences = 10

seed = 123456789

with open(seqListFile,"w") as f:
    
    random.seed(seed)            
    
    B = ["A","T","C","G"]
    
    aliasIndex=0
    for n in range(nSequences):
        seq = ""
        for l in range(seqLength):
            seq = seq + random.choice(B)
        f.write("{:} {:}\n".format("rnd_"+str(aliasIndex),seq))
        aliasIndex+=1

sequenceList = []
with open(seqListFile,"r") as f:
    for line in f:
        alias = line.split()[0]
        seq   = line.split()[1]

        sequenceList.append({"alias":alias,"seq":seq})

simulationPool = []

simId=0
for seq in sequenceList:
    for f in forces:
        for c in range(nCopies):
            
            simulationPool.append({'alias':seq['alias'],
                                   'seq':seq['seq'],
                                   'simId':simId,
                                   'externalForceBtwCOM': [{'force': f, 'pairType': 'S', 'basePair1': 2, 'basePair2': -2}],
                                   'filePath':seq['alias']+"/"+str(f)
                                   })
            simId+=1
    
    for f,t in forcesTorques:
        for c in range(nCopies):
            
            simulationPool.append({'alias':seq['alias'],
                                   'seq':seq['seq'],
                                   'simId':simId,
                                   'externalForce': [{'force': [0.0,0.0,f], 'pairType': 'S', 'basePair': 2}],
                                   'externalTorque': [{'torque': [0.0,0.0,t], 'pairType': 'All', 'basePair': 2}],
                                   'constraintsPositionOfBeads': [{'K': [K, K, K], 'pairType': 'All', 'basePair': -2}],
                                   'filePath':seq['alias']+"/"+str(f)+"_"+str(t)
                                   })
            simId+=1


simPoolSet = MADnaLAB.simulationPool.splitSimulationPoolAccordingToMaxNumberOfBasePairs(simulationPool,int(sys.argv[1]))

cfg = MADnaLAB.utils.loadCONF("highThroughput.cfg")
MADnaLAB.simulationSets.generateSimulationSetsFromListOfSimulationPoolUsingMADnaModel(cfg,simPoolSet)

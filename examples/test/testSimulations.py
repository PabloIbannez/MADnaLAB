import os
import sys

import json

import libconf

from tqdm import tqdm

from MADnaLAB.utils import *
from MADnaLAB.sequencesGeneration import *
from MADnaLAB.simulationPool import *
from MADnaLAB.simulationSets import *

sequences = {"AA":["CGCGAAAAAAAAAACGCG",10.8748],
             "AC":["CGCGACACACACACCGCG",10.9725],
             "AG":["CGCGAGAGAGAGAGCGCG",10.7969],
             "AT":["CGCGATATATATATCGCG",10.8929],
             "CG":["CGCGCGCGCGCGCGCGCG",10.8929],
             "GG":["CGCGGGGGGGGGGGCGCG",10.8748]}

forces = [1.0,5.0,10.0,15.0,20.0]

nSimPerSeq = 100

#################################################

simulationPool = []

simId=0
for seq in sequences.keys():
    for f in forces:
        for c in range(nSimPerSeq):
            
            simulationPool.append({'alias':seq,
                                   'seq':sequences[seq][0],
                                   'simId':simId,
                                   'externalForceBtwCOM': [{'force': f, 'pairType': 'S', 'basePair1': 2, 'basePair2': -2}],
                                   'filePath':seq+"/"+str(f)
                                   })
            simId+=1

with open("simulationPool.json", "w") as f:
    json.dump(simulationPool,f)

simulationPoolList = splitSimulationPoolAccordingToEntry(simulationPool,"alias")

################################################################################################

conf = loadCONF("test.cfg")
    
mainName = conf['MAIN']['name']

simulationSetsName = mainName + 'SimulationSets'
os.makedirs(simulationSetsName, exist_ok=True)

print('[INFO] Generating simulation sets ...')
simSetPathsFile = open(mainName + 'SimulationSets.dat', 'w')

for simSet, simPool in enumerate(tqdm(simulationPoolList)):
    
    simSetName   = mainName + '_' + str(simSet)
    simSetFolder = simulationSetsName + '/' + simSetName + '/'
    
    #Create sim set folder
    os.makedirs(simSetFolder, exist_ok=True)
    
    ########################################################
    
    #!!!!!!
    conf['SIMULATION']['options']['debyeLength'] = sequences[simPool[0]["alias"]][1]
    #######
    
    OPT = generateOptionsFromCONF(conf)
    
    loadSimulationSetInfoIntoOptions(simSetName,simSetFolder,OPT)

    ########################################################

    seq_hsh = hashSequencesFromSimulationPool(simPool)
    generateMADnaTopologiesFromCONFAndHashedSequences(conf,seq_hsh)
    applyTransformationsFromCONFToTopologiesGeneratedFromHashedSequences(conf,seq_hsh)

    generateTopologyFromSimulationPoolMergingFromHashedSequences(simPool,seq_hsh,simSetFolder+simSetName)
    removeTopologiesFromHashedSequences(seq_hsh)
    
    ########################################################

    loadTopologyFileNamesIntoOptions(simSetFolder+simSetName,OPT)
    
    loadBoundariesFromCONFIntoOptions(conf,OPT)
    loadMeasuresFromCONFIntoOptions(conf,OPT)
    
    setUpModificationsFromSimulationPool(simPool,OPT)
    
    ########################################################
    
    generateOutputFoldersFromSimulationPoolAndOptions(simPool,OPT)
    
    ########################################################

    optionsFileName = simSetFolder + 'options_' + str(simSetName) + '.dat'
    simSetPathsFile.write(optionsFileName + '\n')
    
    writeOptionsToFile(OPT,optionsFileName)

simSetPathsFile.close()
print('[INFO] End')

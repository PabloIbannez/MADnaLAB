import sys
import os

from tqdm import tqdm

from MADnaLAB.optionsProcessing import *
from MADnaLAB.sequencesGeneration import *

#from options import *
#from common import *
#from sequencesGeneration import *
#from simulationPool import *

def generateSimulationSetsFromListOfSimulationPoolUsingMADnaModel(conf,simulationPoolList):
    
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
        
        OPT = generateOptionsFromCONF(conf)
        
        loadSimulationSetInfoIntoOptions(simSetName,simSetFolder,OPT)
    
        ########################################################

        seq_hsh = hashSequencesFromSimulationPool(simPool)
        generateMADnaTopologiesFromHashedSequences(seq_hsh)
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

import os

import numpy as np

from Topology import *

from MADnaLAB import *
from MADnaLAB.utils import *
from MADnaLAB.geometric import *

def generateOptionsFromCONF(conf):
    OPTIONS = DEFAULT_SIMULATION_OPTIONS.copy()

    OPTIONS['modelName'] = conf['MAIN']['modelName']

    if 'dt' in conf['SIMULATION']['options']:
        OPTIONS['dt'] = conf['SIMULATION']['options']['dt']

    for opt in conf['SIMULATION']['options']:
        if opt in TIME2STEPS_OPTIONS:
            if TIME2STEPS_OPTIONS[opt] in conf['SIMULATION']['options']:
                print('[ERROR] Both options', opt, 'and', TIME2STEPS_OPTIONS[opt], 'can not be set at the same time')
                sys.exit(0)
            else:
                dt = float(OPTIONS['dt'])
                OPTIONS[TIME2STEPS_OPTIONS[opt]] = int(float(conf['SIMULATION']['options'][opt]) / dt)
        else:
            OPTIONS[opt] = conf['SIMULATION']['options'][opt]
    
    OPTIONS['cutOffDstDH']       = OPTIONS['debyeLength']*DEBYE_FACTOR
    OPTIONS['VerletListDst']     = max(OPTIONS['cutOffDstDH'], OPTIONS['cutOffDstWCA']) + CUTOFF_VERLET_DIFF

    return OPTIONS

def loadMeasuresFromCONFIntoOptions(conf,OPTIONS):
    if 'measures' in conf['SIMULATION']:
        flag = 'measuresActive'
        OPTIONS[flag] = ''
        
        measuresTypes = ''

        dt = float(OPTIONS['dt'])
        
        OPTIONS['nStepsMeasure'] = int(float(conf['SIMULATION']['measures']['measureTime']) / dt)

        for m in conf['SIMULATION']['measures']['measuresList']:
            measuresTypes += m + ' '
        
        OPTIONS['measuresList'] = measuresTypes

def loadTopologyFileNamesIntoOptions(baseName,OPTIONS):
    OPTIONS['inputCoordPath']    = baseName + '.coord'
    OPTIONS['inputTopologyPath'] = baseName + '.top'

def loadBoundariesFromCONFIntoOptions(conf,OPTIONS):
    
    if 'boundaries' in conf['SIMULATION']:
        for boundaryType in conf['SIMULATION']['boundaries']:
            
            flag = 'boundary' + (boundaryType[0].capitalize() + boundaryType[1:]) + 'Active'
            OPTIONS[flag] = ''
            
            options = conf['SIMULATION']['boundaries'][boundaryType]

            if boundaryType == 'zPlates':
                if options['initialPlatesSeparation'] == 'auto':
                    
                    top = Topology(OPTIONS['inputCoordPath'],
                                   OPTIONS['inputTopologyPath'])

                    maxRadius = 0
                    for tinfo in top.propertiesLoaded['TYPES']:
                        t = tinfo[0].split()[0]
                        R = float(tinfo[0].split()[2])
                        maxRadius = max(maxRadius, R)
                        
                    maxHeight = 0
                    minHeight = 0
                    for c in top.coordLoaded:
                        pos = np.asarray([float(p) for p in c[1].split()])
                        z = pos[2]
                        maxHeight = max(maxHeight, z)
                        minHeight = min(minHeight, z)

                    maxSep = max(abs(maxHeight) + maxRadius * 1.5, abs(minHeight) + maxRadius * 1.5)
                    maxSep = maxSep * 2.0
                    options['initialPlatesSeparation'] = maxSep

                for opt in options:
                    OPTIONS[opt] = options[opt]

def loadSimulationSetInfoIntoOptions(simulationSetName,simulationSetFolder,OPTIONS):
    OPTIONS["simulationSetName"]   = simulationSetName
    OPTIONS["simulationSetFolder"] = simulationSetFolder

def generateOutputFoldersFromSimulationPoolAndOptions(simulationPool,OPTIONS):

    OPTIONS["outPutFolders"] = OPTIONS["simulationSetFolder"] + 'filePath_' +  OPTIONS["simulationSetName"] + '.dat'
        
    filePathFile = open(OPTIONS["outPutFolders"], 'w')
    for sim in simulationPool:
        filePathFile.write('{:} {:}\n'.format(sim['simId'], sim['filePath']))
        os.makedirs(sim['filePath'], exist_ok=True)

    filePathFile.close()

def writeOptionsToFile(OPTIONS,fl):
    dict2file(OPTIONS,fl)

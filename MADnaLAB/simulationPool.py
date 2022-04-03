import sys
import copy

from tqdm import tqdm

from MADnaLAB import *

def checkSimulationPool(simulationPool):
    
    if(type(simulationPool) != list):
        print("[ERROR] Simulation pool type has to be \"list\", but current type is:",type(simulationPool))
        sys.exit(0)
    else:
        for sim in simulationPool:
            if(type(sim) != dict):
                print("[ERROR] Elements of simulation pool have to be a \"dict\", but current type is:",type(sim))
                sys.exit(0)
            else:
                for info in SIMPOOL_COMPULSORY_INFO:
                    if(info not in sim.keys()):
                        print("[ERROR] Info:",info,"is not present for the simulation",sim)
                        sys.exit(0)

def generateSimulationPoolFromSequencesList(sequencesList):
    
    print('[INFO] Loading sequence list into simulation pool ...')
    simulationPool = []
    simId = 0

    for s in sequencesList:
        for c in range(s['copies']):
            simBuffer = {i:s[i] for i in s if i != 'copies'}
            simBuffer['simId'] = simId
            simulationPool.append(simBuffer)
            simId += 1
    
    for sim in simulationPool:
        filePath = sim['alias']
        sim['filePath'] = filePath

    checkSimulationPool(simulationPool)

    return simulationPool

def loadModificationsIntoSimulationPoolFromCONF(simulationPool,conf):
    
    if 'MODIFICATIONS' in conf:
        print('[INFO] Adding simulation modifications ...')
        for mod in conf['MODIFICATIONS']:
            for modificationType in conf['MODIFICATIONS'][mod]:
                print('[INFO] Adding', modificationType)
                name = mod + (modificationType[0].capitalize() + modificationType[1:])
                modificationInfo = conf['MODIFICATIONS'][mod][modificationType]
               
                if type(modificationInfo) != tuple:
                    modificationInfo = (modificationInfo,)

                for emod in modificationInfo:
                    emodDict = {name: [dict(emod)]}
                    if 'alias' in emod.keys():
                        for sim in simulationPool:
                            if sim['alias'] == emod['alias']:
                                if name not in sim:
                                    sim.update(copy.deepcopy(emodDict))
                                else:
                                    sim[name].append(copy.deepcopy(dict(emod)))
                    else:
                        for sim in simulationPool:
                            if name not in sim:
                                sim.update(copy.deepcopy(emodDict))
                            else:
                                sim[name].append(copy.deepcopy(dict(emod)))

def splitSimulationPoolAccordingToMaxNumberOfBasePairs(simulationPool,maxNumberOfBasePairs):
    
    if maxNumberOfBasePairs == None:
        simulationPoolList = [simulationPool]
    else:
        simulationPoolList = []
        
        listScore = len(simulationPool[0]['seq'])
        
        tmp = [simulationPool[0]]
        
        for sim in simulationPool[1:]:
            listScore += len(sim['seq'])
            if listScore > maxNumberOfBasePairs:
                simulationPoolList.append(copy.deepcopy(tmp))
                listScore = len(sim['seq'])
                tmp = []
                tmp.append(copy.deepcopy(sim))
            else:
                tmp.append(copy.deepcopy(sim))
        
        simulationPoolList.append(copy.deepcopy(tmp))

    totalSim = 0
    for l in simulationPoolList:
        for sim in l:
            totalSim += 1

    if totalSim != len(simulationPool):
        print('[ERROR] Error splitting simulation pool')
        sys.exit(0)

    print('[INFO] Simulation pool split into', len(simulationPoolList), 'sets')

    return simulationPoolList

def splitSimulationPoolAccordingToEntry(simulationPool,entry):

    #Check all simulations have entry
    for sim in simulationPool:
        if entry not in sim.keys():
            print('[ERROR] Entry:',entry,'is not present in the simulation:',sim)
            sys.exit(0)

    entry2simId = {}
    
    for sim in simulationPool:
        if sim[entry] not in entry2simId.keys():
            entry2simId[sim[entry]] = [sim['simId']]
        else:
            entry2simId[sim[entry]].append(sim['simId'])
        
    simulationPoolList = []

    for entry in entry2simId.keys():
        tmp = []
        for simId in entry2simId[entry]:
            for sim in simulationPool:
                if(sim["simId"] == simId):
                    tmp.append(copy.deepcopy(sim))
        
        simulationPoolList.append(copy.deepcopy(tmp))
    
    totalSim = 0
    for l in simulationPoolList:
        for sim in l:
            totalSim += 1

    if totalSim != len(simulationPool):
        print('[ERROR] Error splitting simulation pool')
        sys.exit(0)

    print('[INFO] Simulation pool split into', len(simulationPoolList), 'sets')

    return simulationPoolList


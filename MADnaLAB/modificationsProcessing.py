import numpy as np

from Topology import *

from MADnaLAB import *
from MADnaLAB.geometric import *
from MADnaLAB.basePairs import *

#def checkBasePairs(sim,modName,basePairName):
#    addedBasePairs = []
#    for info in sim[modName]:
#        for bp in info[basePairName]:
#            if bp not in addedBasePairs:
#                addedBasePairs.append(bp)
#            else:
#                print('[ERROR] Error processing \"{}\" in \"{}\". Base pair \"{}\" added before.'.format(basePairName,modName,bp))
#                sys.exit(1)

def basePairToList(sim,modName,basePairName):
    for info in sim[modName]:
        if(type(info[basePairName])==list):
                pass
        else:
            info[basePairName] = [info[basePairName]]

def constraintAuto(model,constraintType, sim, constraintInfo, top):
    if constraintType == 'constraintsDistanceBtwCOM':
        
        simId = sim['simId']
        
        bp1 = constraintInfo['basePair1']
        bp2 = constraintInfo['basePair2']
        
        types = getBaseType(constraintInfo)
        
        bp1_indices = pairBase2index(model, bp1, top, simId, types)
        bp2_indices = pairBase2index(model, bp2, top, simId, types)
        
        com1 = np.array([0.0, 0.0, 0.0])
        com1 += computeCenterOfMass(bp1_indices, top)
        
        com2 = np.array([0.0, 0.0, 0.0])
        com2 += computeCenterOfMass(bp2_indices, top)
        
        return np.linalg.norm(com2 - com1)

    elif constraintType == 'constraintsPositionOfCOM':
        
        simId = sim['simId']
        
        fbp = constraintInfo['basePair']
        
        types = getBaseType(constraintInfo)
        
        fbp_indices = pairBase2index(model, fbp, top, simId, types)
        
        com = np.array([0.0, 0.0, 0.0])
        com += computeCenterOfMass(fbp_indices, top)
        
        return com
    else:
        print('[ERROR] Auto for', constraintType, 'is not implemented')
        sys.exit(1)

def setUpModificationsFromSimulationPool(simulationPool,OPTIONS):

    model = OPTIONS['modelName']
                    
    top = Topology(OPTIONS['inputCoordPath'],
                   OPTIONS['inputTopologyPath'])
    
    addedModToPool = set()
    for sim in simulationPool:
        for mod in sim.keys():
            if mod.startswith('external'):
                addedModToPool.add(mod)
            if mod.startswith('constraints'):
                addedModToPool.add(mod)

    for mod in addedModToPool:
        if mod not in AVAILABLE_MODIFICATIONS:
            print("[ERROR] Modification:",mod,"is not avaible")
            sys.exit(1)

    simulationSetName   = OPTIONS["simulationSetName"]
    simulationSetFolder = OPTIONS["simulationSetFolder"]

    for mod in addedModToPool:
        flag = mod + 'Active'
        OPTIONS[flag] = ''

        filePath = simulationSetFolder + mod + '_' + simulationSetName + '.dat'

        modFile = open(filePath, 'w')
        OPTIONS[mod + 'FilePath'] = filePath

        for sim in simulationPool:
            if mod in sim.keys():
                simId = sim['simId']

                if mod == 'externalForceBtwCOM':
                    basePairToList(sim,mod,"basePair1")
                    #checkBasePairs(sim,mod,"basePair1")

                    basePairToList(sim,mod,"basePair2")
                    #checkBasePairs(sim,mod,"basePair2")

                    for info in sim[mod]:
                        baseInfo = str(len(info["basePair1"]))
                        for i1 in info["basePair1"]:
                            baseInfo+=" "+str(i1)
                        baseInfo += " "+str(len(info["basePair2"]))
                        for i2 in info["basePair2"]:
                            baseInfo+=" "+str(i2)
                        modFile.write('{:} {:} {} {:.4f}\n'.format(simId, 
                                                                   getBaseType(info), 
                                                                   baseInfo, 
                                                                   info['force']))
                
                if mod == 'externalTorqueBtwCOM':
                    basePairToList(sim,mod,"basePair1")
                    #checkBasePairs(sim,mod,"basePair1")

                    basePairToList(sim,mod,"basePair2")
                    #checkBasePairs(sim,mod,"basePair2")

                    for info in sim[mod]:
                        baseInfo = str(len(info["basePair1"]))
                        for i1 in info["basePair1"]:
                            baseInfo+=" "+str(i1)
                        baseInfo += " "+str(len(info["basePair2"]))
                        for i2 in info["basePair2"]:
                            baseInfo+=" "+str(i2)
                        modFile.write('{:} {:} {} {:.4f}\n'.format(simId, 
                                                                   getBaseType(info), 
                                                                   baseInfo, 
                                                                   info['torque']))

                elif mod == 'externalForce':
                    basePairToList(sim,mod,"basePair")
                    #checkBasePairs(sim,mod,"basePair")
                    for info in sim[mod]:
                        baseInfo = str(len(info["basePair"]))
                        for i in info["basePair"]:
                            baseInfo+=" "+str(i)
                        modFile.write('{:} {:} {} {:.4f} {:.4f} {:.4f}\n'.format(simId, 
                                                                                 getBaseType(info), 
                                                                                 baseInfo, 
                                                                                 info['force'][0], info['force'][1], info['force'][2]))

                elif mod == 'externalTorque':
                    basePairToList(sim,mod,"basePair")
                    #checkBasePairs(sim,mod,"basePair")
                    for info in sim[mod]:
                        baseInfo = str(len(info["basePair"]))
                        for i in info["basePair"]:
                            baseInfo+=" "+str(i)
                        modFile.write('{:} {:} {} {:.4f} {:.4f} {:.4f}\n'.format(simId, 
                                                                                 getBaseType(info), 
                                                                                 baseInfo, 
                                                                                 info['torque'][0], info['torque'][1], info['torque'][2]))

                elif mod == 'constraintsDistanceBtwCOM':
                    basePairToList(sim,mod,"basePair1")
                    #checkBasePairs(sim,mod,"basePair1")

                    basePairToList(sim,mod,"basePair2")
                    #checkBasePairs(sim,mod,"basePair2")

                    for info in sim[mod]:
                        baseInfo = str(len(info["basePair1"]))
                        for i1 in info["basePair1"]:
                            baseInfo+=" "+str(i1)
                        baseInfo += " "+str(len(info["basePair2"]))
                        for i2 in info["basePair2"]:
                            baseInfo+=" "+str(i2)
                        
                        if  info['distance'] == 'auto':
                            info['distance'] = constraintAuto(model, mod, sim, info, top)
                        
                        modFile.write('{:} {:} {} {:.4f} {:.4f}\n'.format(simId, 
                                                                          getBaseType(info), 
                                                                          baseInfo,
                                                                          info['distance'], 
                                                                          info['K']))

                elif mod == 'constraintsPositionOfCOM':
                    basePairToList(sim,mod,"basePair")
                    #checkBasePairs(sim,mod,"basePair")
                    for info in sim[mod]:
                        baseInfo = str(len(info["basePair"]))
                        for i in info["basePair"]:
                            baseInfo+=" "+str(i)

                        if  info['position'] == 'auto':
                            info['position'] = constraintAuto(model, mod, sim, info, top)
                        
                        modFile.write('{:} {:} {} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}\n'.format(simId, 
                                                                                                      getBaseType(info), 
                                                                                                      baseInfo,
                                                                                                      info['position'][0], info['position'][1], info['position'][2], 
                                                                                                      info['K'][0], info['K'][1], info['K'][2]))

                elif mod == 'constraintsPositionOfBeads':
                    basePairToList(sim,mod,"basePair")
                    #checkBasePairs(sim,mod,"basePair")
                    for info in sim[mod]:
                        baseInfo = str(len(info["basePair"]))
                        for i in info["basePair"]:
                            baseInfo+=" "+str(i)
                        modFile.write('{:} {:} {} {:.4f} {:.4f} {:.4f}\n'.format(simId, 
                                                                                 getBaseType(info), 
                                                                                 baseInfo, 
                                                                                 info['K'][0], info['K'][1], info['K'][2]))
        
        modFile.close()


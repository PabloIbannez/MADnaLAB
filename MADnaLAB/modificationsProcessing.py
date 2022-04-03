import numpy as np

from Topology import *

from MADnaLAB.geometric import *
from MADnaLAB.basePairs import *

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
        sys.exit(0)

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
                    for info in sim[mod]:
                        modFile.write('{:} {:} {:} {:} {:.4f}\n'.format(simId, 
                                                                        getBaseType(info), 
                                                                        info['basePair1'], info['basePair2'], 
                                                                        info['force']))

                elif mod == 'externalForce':
                    for info in sim[mod]:
                        modFile.write('{:} {:} {:} {:.4f} {:.4f} {:.4f}\n'.format(simId, 
                                                                                  getBaseType(info), 
                                                                                  info['basePair'], 
                                                                                  info['force'][0], info['force'][1], info['force'][2]))

                elif mod == 'externalTorque':
                    for info in sim[mod]:
                        modFile.write('{:} {:} {:} {:.4f} {:.4f} {:.4f}\n'.format(simId, 
                                                                                  getBaseType(info), 
                                                                                  info['basePair'], 
                                                                                  info['torque'][0], info['torque'][1], info['torque'][2]))

                elif mod == 'constraintsDistanceBtwCOM':
                    for info in sim[mod]:
                        if  info['distance'] == 'auto':
                            info['distance'] = constraintAuto(model, mod, sim, info, top)
                        
                        modFile.write('{:} {:} {:} {:} {:.4f} {:.4f}\n'.format(simId, 
                                                                               getBaseType(info), 
                                                                               info['basePair1'], info['basePair2'], 
                                                                               info['distance'], 
                                                                               info['K']))

                elif mod == 'constraintsPositionOfCOM':
                    for info in sim[mod]:
                        if  info['position'] == 'auto':
                            info['position'] = constraintAuto(model, mod, sim, info, top)
                        
                        modFile.write('{:} {:} {:} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}\n'.format(simId, 
                                                                                                       getBaseType(info), 
                                                                                                       info['basePair'], 
                                                                                                       info['position'][0], info['position'][1], info['position'][2], 
                                                                                                       info['K'][0], info['K'][1], info['K'][2]))

                elif mod == 'constraintsPositionOfBeads':
                    for info in sim[mod]:
                        modFile.write('{:} {:} {:} {:.4f} {:.4f} {:.4f}\n'.format(simId, 
                                                                                  getBaseType(info), 
                                                                                  info['basePair'], 
                                                                                  info['K'][0], info['K'][1], info['K'][2]))

                else:
                    print('[ERROR] Modification', mod, 'is not implemented')
                    sys.exit(0)
        
        modFile.close()


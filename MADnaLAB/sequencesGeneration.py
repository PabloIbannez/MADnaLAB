import sys
import os

import random

import hashlib

from tqdm import tqdm

from MADnaLAB import *
from MADnaLAB.basePairs import *
from MADnaLAB.geometric import *
from MADnaLAB.modificationsProcessing import *

from Topology import *

def generateSequencesListFromCONF(conf):
    
    print('[INFO] Generating sequences list ...')
    sequencesList = []
    for protocol in conf['SEQUENCES']['generate']:
        print('[INFO] Generating protocol:', protocol)
        if protocol == 'list':
            print('[INFO] Protocol "list"')
            dictTmp = [dict(d) for d in conf['SEQUENCES']['generate']['list']]
            for seqInfo in dictTmp:
                sequencesList.append(seqInfo)

        elif protocol == 'random':
            
            seqLength  = conf['SEQUENCES']['generate'][protocol]['sequencesLength']
            nsequences = conf['SEQUENCES']['generate'][protocol]['sequencesNumber']
            seqCopies  = conf['SEQUENCES']['generate'][protocol]['sequencesCopies']

            if 'seed' in conf['SEQUENCES']['generate'][protocol].keys():
                seed = conf['SEQUENCES']['generate'][protocol]['seed']
                random.seed(seed)
                print('[INFO] Protocol "random", with sequences length:', seqLength, ', number of sequences:', nsequences, ', sequences copies:', seqCopies, 'and seed:', seed)
            else:
                print('[INFO] Protocol "random", with sequences length:', seqLength, ', number of sequences:', nsequences, ' and sequences copies:', seqCopies)
            
            B = ['A', 'T', 'C', 'G']
            aliasIndex = 0
            
            for n in range(nsequences):
                seq = ''
                for l in range(seqLength):
                    seq = seq + random.choice(B)
                sequencesList.append({'alias':'rnd_' + str(aliasIndex),  'seq':seq,  'copies':seqCopies})
                aliasIndex += 1
        else:
            print('[ERROR] Generating protocol', protocol, ' not implemented')
            sys.exit(0)
    
    print('[INFO] sequences list ...')
    for s in sequencesList:
        print('      ', s['alias'], ':', s['seq'], ',', s['copies'])

    checkSequencesList(sequencesList)

    return sequencesList

def checkSequencesList(sequencesList):
    
    if(type(sequencesList) != list):
        print("[ERROR] Sequences list type has to be \"list\", but current type is:",type(sequencesList))
        sys.exit(0)
    else:
        for seq in sequencesList:
            if(type(seq) != dict):
                print("[ERROR] Elements of sequences list have to be a \"dict\", but current type is:",type(seq))
                sys.exit(0)
            else:
                for info in SEQ_COMPULSORY_INFO:
                    if(info not in seq.keys()):
                        print("[ERROR] Info:",info,"is not present for the sequence",seq)
                        sys.exit(0)

def hashSequencesFromSimulationPool(simulationPool):

    hashedSequences = {}
    
    for sim in simulationPool:
        seq = sim['seq']
        if seq not in hashedSequences:
            hsh = hashlib.md5(seq.encode('utf-8')).hexdigest()
            for _seq,_hsh in hashedSequences.items():
                if(hsh==_hsh):
                    print("[ERROR] WTF, hash collision !!!")
                    sys.exit(0)
            hashedSequences[seq] = hsh

    return hashedSequences

def removeTopologiesFromHashedSequences(hashedSequences):

    for hsh in hashedSequences.values():
        os.system('rm '+hsh+'.coord')
        os.system('rm '+hsh+'.top')

def generateMADnaTopologiesFromCONFAndHashedSequences(conf,hashedSequences):

    for seq,hsh in hashedSequences.items():
        command = 'python {:} {:} {:} {:}'.format(CREATE_MOLECULE_PATH, seq, hsh + '.coord', hsh + '.top')
        os.system(command)

def generateWLCTopologiesFromCONFAndHashedSequences(conf,hashedSequences):

    parameters = conf['MAIN']['model']['parameters']

    mass = parameters["M"]
    b  = parameters["b"]
    Kb = parameters["Kb"]
    Ka = parameters["Ka"]
    
    radius = b*0.5
    charge = 0.0

    for seq,hsh in hashedSequences.items():
        with open(hsh+".top","w") as f:
            f.write("[TYPES]\n")
            f.write("{:} {:} {:} {:}\n".format("A",mass,radius,charge))
            
            seqLen = len(seq)
            f.write("[STRUCTURE]\n")
            for i in range(seqLen):
                f.write("{:} {:} {:} {:} {:}\n".format(i,"A",i,0,0))
            
            f.write("[BONDS]\n")
            for i in range(seqLen-1):
                f.write("{:} {:} {:} {:}\n".format(i,i+1,b,Kb))
            
            f.write("[ANGLES]\n")
            for i in range(seqLen-2):
                f.write("{:} {:} {:} {:}\n".format(i,i+1,i+2,Ka))
        
        with open(hsh+".coord","w") as f:
            for i in range(seqLen):
                f.write("{:} {:} {:} {:}\n".format(i,0.0,0.0,i*b))
 

#def applyModeFromCONFToMADnaTopologiesGeneratedFromHashedSequences(conf,hashedSequences):
#
#    return
#    
#    #if 'mode' in conf['SEQUENCES']:
#    #    
#    #    mode = conf['SEQUENCES']['mode']
#
#    #    if mode.split()[0] == 'basic':
#    #        pass
#    #    if mode.kplit()[0] == 'fast': #(TODO)
#
#    #        for hsh in hashedSequences.values():
#
#    #            top = Topology(hsh + '.coord', hsh + '.top')
#
#    #            nBasis = set()
#    #            for s in top.propertiesLoaded['STRUCTURE']:
#    #                nBasis.add(s[2])
#
#    #            nBasisPairs = len(nBasis) // 2
#
#    #            phosphateIndexPair = []
#
#    #            for struct in top.propertiesLoaded['STRUCTURE']:
#    #                
#    #                index, tp, res, chain, mol, sim = struct
#    #                
#    #                basisPairIndex = res2basisPair(res, nBasisPairs)
#    #                
#    #                if tp == 'P':
#    #                    phosphateIndexPair.append([index, basisPairIndex])
#
#    #                sys.exit(0)

def applyTransformationsFromCONFToTopologiesGeneratedFromHashedSequences(conf,hashedSequences):

    model = conf['MAIN']['model']['name']

    if 'SEQUENCES' in conf:
        if 'transformations' in conf['SEQUENCES']:
            
            for trans in conf['SEQUENCES']['transformations']:
                if trans == 'alingAlongAxis':

                    axis = conf['SEQUENCES']['transformations']['alingAlongAxis']['axis']
                    bp1  = conf['SEQUENCES']['transformations']['alingAlongAxis']['basePair1']
                    bp2  = conf['SEQUENCES']['transformations']['alingAlongAxis']['basePair2']

                    for seq,hsh in hashedSequences.items():

                        top = Topology(hsh + '.coord', hsh + '.top')

                        sqLen = len(seq)
                        
                        res1 = pairBase2res(model, bp1, sqLen)
                        res2 = pairBase2res(model, bp2, sqLen)

                        allIndices = []
                        bp1Indices = []
                        bp2Indices = []

                        for struct in top.propertiesLoaded['STRUCTURE']:
                            
                            index, tp, res, chain, mol, sim = struct
                            allIndices.append(index)
                            if res in res1:
                                bp1Indices.append(index)
                            if res in res2:
                                bp2Indices.append(index)

                        comAll = computeCenterOfMass(allIndices, top)

                        for i in allIndices:
                            pos = np.asarray([float(p) for p in top.coordLoaded[i][1].split()])
                            pos = pos - comAll

                            pNew = [i, '{:.4f} {:.4f} {:.4f}'.format(pos[0], pos[1], pos[2])]

                            top.coordLoaded[i] = pNew.copy()
                        
                        com1 = computeCenterOfMass(bp1Indices, top)
                        com2 = computeCenterOfMass(bp2Indices, top)

                        r12 = com2 - com1

                        r12 = r12 / np.sqrt(np.sum(r12 ** 2))

                        axisN = axis / np.sqrt(np.sum(np.asarray(axis) ** 2))

                        angle = np.arccos(np.dot(r12, axisN))

                        rotAxis = np.cross(r12, axisN)
                        rotAxis = rotAxis / np.sqrt(np.sum(np.asarray(rotAxis) ** 2))

                        rotateStructure(rotAxis, angle, top)
                        comAll = computeCenterOfMass(allIndices, top)

                        top.write(hsh)

                else:
                    print('[ERROR] Transformation', trans, ' not implemented')
                    sys.exit(0)

def generateTopologyFromSimulationPoolMergingFromHashedSequences(simulationPool,hashedSequences,topologyPath):

    for sim in simulationPool:
        if sim["seq"] not in hashedSequences.keys():
            print("[ERROR] There is no equivalence between"
                  " simulation pool and hashed sequences."
                  " Not found seq:",sim["seq"],"in hashed sequences")
            sys.exit(0)
    
    if len(simulationPool) > 1:

        hashedName1 = hashedSequences[simulationPool[0]['seq']]
        hashedName2 = hashedSequences[simulationPool[1]['seq']]

        top1 = Topology(hashedName1 + '.coord', hashedName1 + '.top')
        top1.setSimId(simulationPool[0]['simId'])
        top2 = Topology(hashedName2 + '.coord', hashedName2 + '.top')
        top2.setSimId(simulationPool[1]['simId'])

        merging2File(top1, top2, topologyPath, 'none')
        for i, s in enumerate(simulationPool):
            if i < 2:
                pass
            else:
                hashedName = hashedSequences[s['seq']]

                top1 = Topology(topologyPath + '.coord', topologyPath + '.top')
                top2 = Topology(hashedName + '.coord', hashedName + '.top')

                top2.setSimId(s['simId'])

                merging2File(top1, top2, topologyPath, 'none')

    else:
        hashedName = hashedSequences[simulationPool[0]['seq']]

        top = Topology(hashedName + '.coord', hashedName + '.top')

        top.setSimId(simulationPool[0]['simId'])
        top.write(topologyPath)

import sys
import os

import random
import hashlib

import numpy as np
from scipy.spatial.transform import Rotation

import libconf

from Topology import *

createMoleculePath = os.getenv('UAMMDPATH')+"/extensions/structured/Tools/MADna/CreateMolecule/CreateMolecule.py"

######### DEFAULT OPTIONS #########

debyeFactor = 4.0
cutoffVerletDiff = 5.0

SIMULATION_OPTIONS = {

"boxSize":[10000.0,10000.0,10000.0],

"T":300.0,

"h":0.001,

"nStepsSteepestDescent":0,
"nStepsSteepestDescentProgressInterval":100,
"maxObjectiveForce":5.0,

"outPutFilePath":"output",
"outPutFormat":"lammpstrj",

"dt":0.01,
"frictionConstant":0.2,

"dielectricConstant":78.3,
"debyeLength":7.88,

"cutOffDstWCA":25.0
}

time2stepsOptions = {"simulationTime":"nSteps",
                     "infoTime":"nStepsInfoInterval",
                     "writeTime":"nStepsWriteInterval",
                     "backupTime":"nStepsBackupInterval"}

########################################################

def dict2File(inDict:dict,path:str):
    with open(path,"w") as f:
        for entry in inDict:
            if type(inDict[entry])==list:
                f.write("{:}".format(entry))
                for e in inDict[entry]:
                    if type(e)==list:
                        for ee in e:
                            f.write(" {:}".format(ee))
                    else:
                        f.write(" {:}".format(e))
                f.write("\n")
            else:
                f.write("{:} {:}\n".format(entry,inDict[entry]))

def getPairsIndex(pairBasePosition,
                  nBasisPairs):
    
    if(pairBasePosition == 0):
        print("[MADna] Invalid pair base position, it must be different than 0")
        sys.exit(0)

    nBasis = int(nBasisPairs*2)

    if(abs(pairBasePosition) > nBasisPairs ):
        print("[MADna] Invalid pair base position, its absolute value can not be larger than the total number of base pairs ({:}), but the value is: {:}".format(nBasisPairs,pairBasePosition))
        sys.exit(0)

    pairsIndex = []
    if(pairBasePosition > 0):
        pairsIndex.append(pairBasePosition-1)
        pairsIndex.append(nBasis-1-pairsIndex[0])
    else:
      pairsIndex.append(nBasisPairs+pairBasePosition)
      pairsIndex.append(nBasis-1-pairsIndex[0])

    return pairsIndex

def getAllPairsIndex(pairBasePosition,
                     top):

    allPairsIndex = []

    simIds = set()
    for s in top.propertiesLoaded['STRUCTURE']:
        simIds.add(s[5])

    for sim in simIds:
        allPairsSimIndex = []

        nBasis = set()
        for s in top.propertiesLoaded['STRUCTURE']:
            if s[5]==sim:
                nBasis.add(s[2])
        nBasisPairs = len(nBasis)//2
        res = getPairsIndex(pairBasePosition,nBasisPairs)
        for s in top.propertiesLoaded['STRUCTURE']:
            if s[5]==sim:
                if(s[2] in res):
                    allPairsSimIndex.append(s[0])
        allPairsIndex.append(allPairsSimIndex)
    return allPairsIndex

def computeCOM(indices,top):

    tM  = 0.0
    com = np.array([0.0,0.0,0.0])

    masses = {}
    for tinfo in top.propertiesLoaded['TYPES']:
        t=tinfo[0].split()[0]
        m=float(tinfo[0].split()[1])
        masses[t]=m
    for i in indices:
        pos = np.asarray([float(p) for p in top.coordLoaded[i][1].split()])
        tp   = top.propertiesLoaded['STRUCTURE'][i][1]
        tM  = tM  + masses[tp]
        com = com + pos*masses[tp]

    return com/tM

def rotateStructure(axis,angle,top):

    rot = Rotation.from_rotvec(angle*axis)

    #f = open("testRot.sp","w")
    #f.write("#\n")
    #for c in top.coordLoaded:
    #    pos = np.asarray([float(p) for p in c[1].split()])
    #    f.write("{:.4f} {:.4f} {:.4f}\n".format(pos[0],pos[1],pos[2]))
    #
    #f.write("#\n")
    #for c in top.coordLoaded:
    #    pos = np.asarray([float(p) for p in c[1].split()])
    #    posNew = rot.apply(pos)
    #    f.write("{:.4f} {:.4f} {:.4f}\n".format(posNew[0],posNew[1],posNew[2]))
    #f.close()

    for i,c in enumerate(top.coordLoaded):
        pos = np.asarray([float(p) for p in c[1].split()])
        pos = rot.apply(pos)
        pNew = [c[0],"{:.4f} {:.4f} {:.4f}".format(pos[0],pos[1],pos[2])] 
        top.coordLoaded[i]=pNew.copy()

def constraintAuto(constraintType,constInfo,top):
    if constraintType == "harmonicCOMdistance":
        bp1 = constInfo["constrainedBasePair1"]
        bp2 = constInfo["constrainedBasePair2"]
        
        bp1_indices = getAllPairsIndex(bp1,top)
        bp2_indices = getAllPairsIndex(bp2,top)

        nSim = len(bp1_indices)

        com1=np.array([0.0,0.0,0.0])
        for i in range(nSim):
            com1+=computeCOM(bp1_indices[i],top)
        com1/=nSim
        com2=np.array([0.0,0.0,0.0])
        for i in range(nSim):
            com2+=computeCOM(bp2_indices[i],top)
        com2/=nSim
        return np.linalg.norm(com2-com1)
    elif constraintType == "harmonicCOMfixed":
        fbp = constInfo["fixedBasePair1"]
        
        fbp_indices = getAllPairsIndex(fbp,top)

        nSim = len(fbp_indices)
        com=np.array([0.0,0.0,0.0])
        for i in range(nSim):
            com+=computeCOM(fbp_indices[i],top)
        com/=nSim
        return list(com)
    else:
        print("[ERROR] Auto for",constraintType,"is not implemented")
        sys.exit(0)


def genMADna(cfgFile):
    
    print("[INFO] Start")
    ################################Start################################
    
    if(type(cfgFile)==str):
        print("[INFO] Reading input from string:",cfgFile)
        with open(cfgFile,"r") as inpt:
            conf = libconf.load(inpt)
    else:
        print("[INFO] Loaded input from object:",type(cfgFile))
        conf=cfgFile
    
    ########################Generate coord and top#######################
    
    sequencesList = {}
    for protocol in conf['SEQUENCES']['generate']:
        print("[INFO] Generating protocol:",protocol)
        if protocol == "list":
            print("[INFO] Protocol \"list\"")
            dictTmp = dict(conf['SEQUENCES']['generate']['list'])
            for seq in dictTmp:
                sequencesList[seq]=dictTmp[seq]
        elif protocol == "random":
            
            seqLength  = conf['SEQUENCES']['generate'][protocol]["sequencesLength"]
            nSequences = conf['SEQUENCES']['generate'][protocol]["sequencesNumber"]
            seqCopies  = conf['SEQUENCES']['generate'][protocol]["sequencesCopies"]
            
            print("[INFO] Protocol \"random\", with sequences length:",seqLength,", number of sequences:",nSequences," and sequences copies:",seqCopies)
            
            B = ["A","T","C","G"]
            
            for n in range(nSequences):
                seq = ""
                for l in range(seqLength):
                    seq = seq + random.choice(B)
                sequencesList[seq]=seqCopies 
        else:
            print("[ERROR] Generating protocol",protocol," not implemented")
            sys.exit(0)
    
    print("[INFO] Sequence list ...")
    for s in sequencesList:
        print("      ",s,sequencesList[s])
    
    print("[INFO] Hashing sequences ...")
    sequencesHash = {}
    for s in sequencesList:
        sequencesHash[s]=hashlib.md5(s.encode('utf-8')).hexdigest()
    
    print("[INFO] Calling CreateMolecule.py (at",createMoleculePath,")")
    
    filesToRemove = []
    
    sequencesName = conf['SEQUENCES']['name']
    
    for s in sequencesList:
        command = "python {:} {:} {:} {:}".format(createMoleculePath,s,sequencesHash[s]+".coord",sequencesHash[s]+".top")
        #print("[INFO] Command:",command)
        os.system(command)
        filesToRemove.append(sequencesHash[s]+".coord")
        filesToRemove.append(sequencesHash[s]+".top")
    
    if('transformations' in conf['SEQUENCES']):
        print("[INFO] Transformations detected")
        for trans in conf['SEQUENCES']['transformations']:
            if(trans == "alingAlongAxis"):
                
                axis= conf['SEQUENCES']['transformations']['alingAlongAxis']['axis']
                bp1 = conf['SEQUENCES']['transformations']['alingAlongAxis']['basisPair1']
                bp2 = conf['SEQUENCES']['transformations']['alingAlongAxis']['basisPair2']
    
                print("[INFO] Transformation \"alingAlongAxis\", with axis:",axis,", base pair 1:",bp1," and base pair 2:",bp2)
                for s in sequencesList:
                    top = Topology(sequencesHash[s]+".coord",
                                   sequencesHash[s]+".top")
    
                    sqLen = len(s)
                    res1 = getPairsIndex(bp1,sqLen)
                    res2 = getPairsIndex(bp2,sqLen)
    
                    allIndices  = []
                    bp1Indices = []
                    bp2Indices = []
                    for struct in top.propertiesLoaded['STRUCTURE']:
                        index,tp,res,chain,mol,sim = struct  
                        allIndices.append(index)
                        if res in res1:
                            bp1Indices.append(index)
                        if res in res2:
                            bp2Indices.append(index)
    
                    #print(bp1Indices,bp2Indices)
                    comAll = computeCOM(allIndices,top)
                    #Center at COM
                    for i in allIndices:
                        pos = np.asarray([float(p) for p in top.coordLoaded[i][1].split()])
                        pos = pos - comAll
    
                        pNew = [i,"{:.4f} {:.4f} {:.4f}".format(pos[0],pos[1],pos[2])] 
                        #print(top.coordLoaded[i],pNew)
                        top.coordLoaded[i]=pNew.copy()
                    
                    com1   = computeCOM(bp1Indices,top)
                    com2   = computeCOM(bp2Indices,top)
    
                    r12   = com2-com1
                    r12   = r12/np.sqrt(np.sum(r12**2))
                    axisN = axis/np.sqrt(np.sum(np.asarray(axis)**2))
    
                    angle   = np.arccos(np.dot(r12,axisN))
                    rotAxis = np.cross(r12,axisN)
                    rotAxis = rotAxis/np.sqrt(np.sum(np.asarray(rotAxis)**2))
    
                    rotateStructure(rotAxis,angle,top)
                    
                    comAll = computeCOM(allIndices,top)
                    #print("[INFO] Final COM",comAll)
                    
                    #override old
                    top.write(sequencesHash[s])
    
            else:
                print("[ERROR] Transformation",trans," not implemented")
                sys.exit(0)
    
    
    sqlhc = []
    print("[INFO] Coping and merging sequences ...")
    for s in sequencesList:
        top = Topology(sequencesHash[s]+".coord",
                       sequencesHash[s]+".top")
        nCopies2File(top,sequencesList[s],sequencesHash[s]+str(sequencesList[s]),"simId")
        sqlhc.append(sequencesHash[s]+str(sequencesList[s]))
        filesToRemove.append(sequencesHash[s]+str(sequencesList[s])+".coord")
        filesToRemove.append(sequencesHash[s]+str(sequencesList[s])+".top")
    
    if len(sequencesList) > 1:
        top1 = Topology(sqlhc[0]+".coord",sqlhc[0]+".top")
        top2 = Topology(sqlhc[1]+".coord",sqlhc[1]+".top")
        
        merging2File(top1,top2,sequencesName,"simId")
        
        for s in sqlhc[2::]:
            top1 = Topology(sequencesName+".coord",sequencesName+".top")
            top2 = Topology(s+".coord",s+".top")
            merging2File(top1,top2,sequencesName,"simId")
    else:
        top = Topology(sqlhc[0]+".coord",sqlhc[0]+".top")
        top.write(sequencesName)
    
    print("[INFO] Removing temporal files ...")
    for f in filesToRemove:
        os.system("rm "+f)
    
    #####################################################################
    
    print("[INFO] Reading simulation options ...")
    #If dt in conf file, read it the first
    if("dt" in conf['SIMULATION']['options']):
        SIMULATION_OPTIONS['dt']=conf['SIMULATION']['options']['dt']
    for opt in conf['SIMULATION']['options']:
        if(opt in time2stepsOptions):
            if(time2stepsOptions[opt] in conf['SIMULATION']['options']):
                print("[ERROR] Both options",opt,"and",time2stepsOptions[opt],"can not be set at the same time")
                sys.exit(0)
            else:
                dt = float(SIMULATION_OPTIONS['dt'])
                SIMULATION_OPTIONS[time2stepsOptions[opt]]=int(float(conf['SIMULATION']['options'][opt])/dt)
        else:
            SIMULATION_OPTIONS[opt]=conf['SIMULATION']['options'][opt]
    
    #Update cutOffDstDH,VerletListDst
    SIMULATION_OPTIONS["cutOffDstDH"]=SIMULATION_OPTIONS["debyeLength"]*debyeFactor
    SIMULATION_OPTIONS["VerletListDst"]=max(SIMULATION_OPTIONS["cutOffDstDH"],SIMULATION_OPTIONS["cutOffDstWCA"])+cutoffVerletDiff
    
    SIMULATION_OPTIONS["inputCoordPath"]=sequencesName+".coord"
    SIMULATION_OPTIONS["inputTopologyPath"]=sequencesName+".top"
    
    ##########################Add simulation mod#########################
    
    if 'modifications' in conf['SIMULATION']:
        print("[INFO] Adding simulation modifications ...")
        for mod in conf['SIMULATION']['modifications']:
            if mod == "pulling":
                for pullingType in conf['SIMULATION']['modifications']['pulling']:
                    print("[INFO] Adding",pullingType)
                    flag = mod+(pullingType[0].capitalize()+pullingType[1::])+"Active"
                    SIMULATION_OPTIONS[flag]=""
                    options = conf['SIMULATION']['modifications']['pulling'][pullingType]
                    for opt in options:
                        SIMULATION_OPTIONS[opt]=options[opt]
            if mod == "constraints":
                for constraintType in conf['SIMULATION']['modifications']['constraints']:
                    print("[INFO] Adding",constraintType)
                    flag = "constraint"+(constraintType[0].capitalize()+constraintType[1::])+"Active"
                    SIMULATION_OPTIONS[flag]=""
                    #list of constraints of type constraintType
                    constInfo = conf['SIMULATION']['modifications']['constraints'][constraintType]
                    
                    if type(constInfo) != tuple:
                        constInfo=(constInfo,)

                    ks = (constInfo[0]).keys()
                    val = {}
                    for k in ks:
                        val[k]=[]
                        for const in constInfo:
                            if const[k] == "auto":
                                top = Topology(sequencesName+".coord",sequencesName+".top")
                                const[k]=constraintAuto(constraintType,const,top)
                            val[k].append(const[k])
                    for k in ks:
                        SIMULATION_OPTIONS[k]=val[k]
        
            if mod == "boundaries":
                for boundaryType in conf['SIMULATION']['modifications']['boundaries']:
                    print("[INFO] Adding",boundaryType)
                    flag = "boundary"+(boundaryType[0].capitalize()+boundaryType[1::])+"Active"
                    SIMULATION_OPTIONS[flag]=""
                    options = conf['SIMULATION']['modifications']['boundaries'][boundaryType]
                    if(boundaryType == "zPlates"):
                        if(options["initialPlatesSeparation"] == "auto"):
        
                            top = Topology(sequencesName+".coord",sequencesName+".top")
                            
                            maxRadius = 0
                            for tinfo in top.propertiesLoaded['TYPES']:
                                t=tinfo[0].split()[0]
                                R=float(tinfo[0].split()[2])
                                maxRadius=max(maxRadius,R)
        
                            maxHeight = 0
                            minHeight = 0
                            for c in top.coordLoaded:
                                pos = np.asarray([float(p) for p in c[1].split()])
                                z=pos[2]
                                maxHeight=max(maxHeight,z)
                                minHeight=min(minHeight,z)
        
                            maxSep = max(abs(maxHeight)+maxRadius*1.5,abs(minHeight)+maxRadius*1.5)
                            maxSep = maxSep*2.0
                            options["initialPlatesSeparation"]=maxSep
        
                        for opt in options:
                            SIMULATION_OPTIONS[opt]=options[opt]
    
    if 'measures' in conf['SIMULATION']:
        print("[INFO] Adding simulation measures ...")
        flag = "measuresActive"
        SIMULATION_OPTIONS[flag]=""
        measuresTypes = ""
        
        dt = float(SIMULATION_OPTIONS['dt'])
        SIMULATION_OPTIONS["nStepsMeasure"]=int(float(conf['SIMULATION']['measures']['measureTime'])/dt)

        for m in conf['SIMULATION']['measures']['measuresList']:
            measuresTypes+=m+" "
        SIMULATION_OPTIONS["measuresList"]=measuresTypes
    
    #########################Writeout options.dat########################
    
    print("[INFO] Writing options ...")
    dict2File(SIMULATION_OPTIONS,"./options.dat")
    print("[INFO] End")

if __name__ == "__main__":
    genMADna(sys.argv[1])

import sys

def getBaseType(info):
    bType = 'All'
    if 'pairType' in info:
        bType = info['pairType']
    return bType

def res2basePair(model, res, nBasisPairs):
    if model == "MADna":
        nBasis = int(nBasisPairs * 2)
        if res >= nBasis:
            print("[ERROR] Model:{:}.Invalid res, its value can not be larger "
                   "than the total number of basis ({:}),"
                   "but the value is: {:}".format(model,nBasis, res))
            sys.exit(0)
        if res < nBasisPairs:
            return res + 1
        return nBasis - res
    elif model == "WLC":
        if res >= nBasisBasisPairs:
            print("[ERROR] Model:{:}. Invalid res, its value can not be larger "
                   "than the total number of basis pairs({:}),"
                   "but the value is: {:}".format(model,nBasisPairs, res))
            sys.exit(0)
        return res+1
    else:
        print("[ERROR] Model",model,"is not avaible")
        sys.exit(0)

def pairBase2res(model, pairBasePosition, nBasisPairs):
    if model == "MADna":
        if pairBasePosition == 0:
            print("[ERROR] Invalid pair base position,"
                  " it must be different than 0")
            sys.exit(0)
        nBasis = int(nBasisPairs * 2)
        if abs(pairBasePosition) > nBasisPairs:
            print("[ERROR] Model:{:}, Invalid pair base position,"
                  " its absolute value can not be larger"
                  " than the total number of base pairs"
                  " ({:}), but the value is: {:}".format(model,nBasisPairs, pairBasePosition))
            sys.exit(0)
        pairsIndex = []
        if pairBasePosition > 0:
            pairsIndex.append(pairBasePosition - 1)
            pairsIndex.append(nBasis - 1 - pairsIndex[0])
        else:
            pairsIndex.append(nBasisPairs + pairBasePosition)
            pairsIndex.append(nBasis - 1 - pairsIndex[0])
        return pairsIndex
    elif model == "WLC":
        if pairBasePosition == 0:
            print("[ERROR] Invalid pair base position,"
                  " it must be different than 0")
            sys.exit(0)
        if abs(pairBasePosition) > nBasisPairs:
            print("[ERROR] Model:{:}, Invalid pair base position,"
                  " its absolute value can not be larger"
                  " than the total number of base pairs"
                  " ({:}), but the value is: {:}".format(model,nBasisPairs, pairBasePosition))
            sys.exit(0)
        pairsIndex = []
        if pairBasePosition > 0:
            pairsIndex.append(pairBasePosition - 1)
        else:
            pairsIndex.append(nBasisPairs + pairBasePosition)
        return pairsIndex

    else:
        print("[ERROR] Model",model,"is not avaible")
        sys.exit(0)

def pairBase2index(model, pairBasePosition, top, simId, bTypes):
    
    simPairsIndex = []
    nBasis = set()
    for s in top.propertiesLoaded['STRUCTURE']:
        if s[5] == simId:
            nBasis.add(s[2])
        
    if model == "MADna":
        nBasisPairs = len(nBasis) // 2
    if model == "WLC":
        nBasisPairs = len(nBasis)
    else:
        print("[ERROR] Model",model,"is not avaible")
        sys.exit(0)
    
    res = pairBase2res(model, pairBasePosition, nBasisPairs)

    for s in top.propertiesLoaded['STRUCTURE']:
        if s[5] == simId:
            if s[2] in res:
                if s[1] in bTypes:
                    simPairsIndex.append(s[0])
    
    return simPairsIndex


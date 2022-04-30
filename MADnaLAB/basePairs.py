import sys

def getBaseType(info):
    bType = 'All'
    if 'pairType' in info:
        bType = info['pairType']
    return bType

def res2basePair(model, res, nBasisPairs):
    if model == "MADna" or model == "MADnaFast":
        nBasis = int(nBasisPairs * 2)
        if res >= nBasis:
            print("[ERROR] Model:{:}.Invalid res, its value can not be larger "
                   "than the total number of basis ({:}),"
                   "but the value is: {:}".format(model,nBasis, res))
            sys.exit(1)
        if res < nBasisPairs:
            return res + 1
        return nBasis - res
    elif model == "WLC":
        if res >= nBasisBasisPairs:
            print("[ERROR] Model:{:}. Invalid res, its value can not be larger "
                   "than the total number of basis pairs({:}),"
                   "but the value is: {:}".format(model,nBasisPairs, res))
            sys.exit(1)
        return res+1
    else:
        print("[ERROR] Model",model,"is not avaible")
        sys.exit(1)

def pairBase2res(model, pairBasePosition, nBasisPairs):
    
    if type(pairBasePosition) == list:
        pass
    else:
        pairBasePosition=[pairBasePosition]

    if model == "MADna" or model == "MADnaFast":
        if 0 in pairBasePosition:
            print("[ERROR] Invalid pair base position,"
                  " it must be different than 0")
            sys.exit(1)
        nBasis = int(nBasisPairs * 2)
        for pbp in pairBasePosition:
            if abs(pbp) > nBasisPairs:
                print("[ERROR] Model:{:}, Invalid pair base position,"
                      " its absolute value can not be larger"
                      " than the total number of base pairs"
                      " ({:}), but the value is: {:}".format(model,nBasisPairs,pbp))
                sys.exit(1)
        
        pairsIndex = []
        for pbp in pairBasePosition:
            pairsIndexBuffer = []
            if pbp > 0:
                pairsIndexBuffer.append(pbp - 1)
                pairsIndexBuffer.append(nBasis - 1 - pairsIndexBuffer[0])
            else:
                pairsIndexBuffer.append(nBasisPairs + pbp)
                pairsIndexBuffer.append(nBasis - 1 - pairsIndexBuffer[0])
            pairsIndex.extend(pairsIndexBuffer)

        return pairsIndex
    elif model == "WLC":
        if 0 in pairBasePosition:
            print("[ERROR] Invalid pair base position,"
                  " it must be different than 0")
            sys.exit(1)
        for pbp in pairBasePosition:
            if abs(pbp) > nBasisPairs:
                print("[ERROR] Model:{:}, Invalid pair base position,"
                      " its absolute value can not be larger"
                      " than the total number of base pairs"
                      " ({:}), but the value is: {:}".format(model,nBasisPairs,pbp))
                sys.exit(1)

        pairsIndex = []
        for pbp in pairBasePosition:
            pairsIndexBuffer = []
            if pbp > 0:
                pairsIndexBuffer.append(pbp - 1)
            else:
                pairsIndexBuffer.append(nBasisPairs + pbp)
            pairsIndex.extend(pairsIndexBuffer)

        return pairsIndex
    else:
        print("[ERROR] Model",model,"is not avaible")
        sys.exit(1)

def pairBase2index(model, pairBasePosition, top, simId, bTypes):

    if type(pairBasePosition) == list:
        pass
    else:
        pairBasePosition=[pairBasePosition]
    
    simPairsIndex = []
    nBasis = set()
    for s in top.propertiesLoaded['STRUCTURE']:
        if s[5] == simId:
            nBasis.add(s[2])
        
    if model == "MADna" or model == "MADnaFast":
        nBasisPairs = len(nBasis) // 2
    elif model == "WLC":
        nBasisPairs = len(nBasis)
    else:
        print("[ERROR] Model",model,"is not avaible")
        sys.exit(1)
    
    res = pairBase2res(model, pairBasePosition, nBasisPairs)

    for s in top.propertiesLoaded['STRUCTURE']:
        if s[5] == simId:
            if s[2] in res:
                if s[1] in bTypes:
                    simPairsIndex.append(s[0])
    
    return simPairsIndex


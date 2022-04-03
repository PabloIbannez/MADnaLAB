import libconf

def loadCONF(cfgFile):
    
    if type(cfgFile) == str:
        print('[INFO] Reading input from string:', cfgFile)
        with open(cfgFile, 'r') as inpt:
            conf = libconf.load(inpt)
    else:
        print('[INFO] Loaded input from object:', type(cfgFile))
        conf = cfgFile

    return conf

def dict2file(inDict: dict, path: str):
    with open(path, 'w') as f:
        for entry in inDict:
            if type(inDict[entry]) == list:
                f.write('{:}'.format(entry))
                for e in inDict[entry]:
                    if type(e) == list:
                        for ee in e:
                            f.write(' {:}'.format(ee))
                    else:
                        f.write(' {:}'.format(e))
                f.write('\n')
            else:
                f.write('{:} {:}\n'.format(entry, inDict[entry]))


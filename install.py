
import os
import pathlib

export = "export PYTHONPATH=\"${PYTHONPATH}:"

HOME_path = "/home/pablo"
BIN_path =  "/home/pablo/bin"

MADnaLAB_path = pathlib.Path(__file__).parent.resolve()

UAMMD_path = "/home/pablo/Desktop/UAMMD"
UAMMD_struct_path = UAMMD_path+"/extensions/structured"

TOP_path = UAMMD_struct_path+"/Tools/TopologyUtilities"

export+=str(MADnaLAB_path)+"/src:"
export+=TOP_path
export+="\""

append=True
with open(HOME_path+"/.bashrc", 'r') as f:
    for line in f:
        if(export in line):
            append=False
if append:
    with open(HOME_path+"/.bashrc", 'a') as f:
        f.write(export+"\n")

##############################################

os.system("rm -f ./bin/MADnaLAB")
os.system("make")
os.system("cp ./bin/MADnaLAB "+BIN_path)

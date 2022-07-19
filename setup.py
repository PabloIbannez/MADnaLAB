import sys

import os
import pathlib

import subprocess

UAMMD_path= "/home/pablo/Desktop/UAMMD"

#######################################################

HOME_path = os.path.expanduser('~')

UAMMD_struct_path = UAMMD_path+"/extensions/structured"
TOP_path = UAMMD_struct_path+"/Tools/TopologyUtilities"

MADnaLAB_path = str(pathlib.Path(__file__).parent.parent.resolve())

#######################################################

export = []

export.append("export PYTHONPATH=\"${PYTHONPATH}:")

export[0]+=str(MADnaLAB_path)+":"
export[0]+=TOP_path
export[0]+="\""

#####
export.append("export UAMMDPATH="+UAMMD_path)
#####
export.append("export MADNALABPATH="+MADnaLAB_path)
#####

for ex in export:
    append=True
    with open(HOME_path+"/.bashrc", 'r') as f:
        for line in f:
            if(ex in line):
                append=False
    if append:
        with open(HOME_path+"/.bashrc", 'a') as f:
            f.write(ex+"\n")

#######################################################

subprocess.check_call([sys.executable, "-m", "pip", "install","-r","requirements.txt"])

#######################################################

#Compile
subprocess.check_call(["make", "-C", "./src"])

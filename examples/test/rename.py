import os
import json

from tqdm import tqdm

nSim = 100

#fold="$fbase/Simulations/$sequenza/f$force/sim$sim"
#${fscripts}/Compute_ext_twist_crook $fold/dump.lammpstrj $RANDOM $i0 $i1 > $fold/ext_twist_crook & 

with open('simulationPool.json') as f:
    simulationPool = json.load(f)

for sim in tqdm(simulationPool):
    old = sim["filePath"]+"/output_"+str(sim["simId"])+".lammpstrj"
    newPath = "Simulations/"+sim["alias"]+"/f"+str(int(sim["externalForceBtwCOM"][0]["force"]))+"/sim"+str((sim["simId"]%nSim)+1)
    os.makedirs(newPath, exist_ok=True)
    new = newPath+"/dump.lammpstrj"
    #os.system("cp "+old+" "+new)
    os.system("mv "+old+" "+new)

seqList = {}
for sim in simulationPool:
    if sim["alias"] not in seqList.keys():
        seqList[sim["alias"]] = sim["seq"]

for alias,seq in seqList.items():
    path = "Simulations/"+alias+"/initialization"
    os.makedirs(path, exist_ok=True)
    with open(path+"/sequence.dat","w") as f:
        f.write(seq)
    
path = "Simulations/Comparison"
os.makedirs(path, exist_ok=True)



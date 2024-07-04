import json
import numpy as np
import matplotlib.pyplot as plt

free  = "FREE"
left  = "HELIX_BOUNDARIES_LEFT"
right = "HELIX_BOUNDARIES_RIGHT"

#free  = "FREE_RANDOM"
#left  = "HELIX_BOUNDARIES_LEFT_RANDOM"
#right = "HELIX_BOUNDARIES_RIGHT_RANDOM"

skipFrames = 100

# 0 -> ste
internalTerms = [1,2,3,4,5] # bond angle dihedral dh wca
externalTerms = [6] # helix_boundaries

with open(free+"/VLMPsession.json") as f:
    session = json.load(f)

sequences = session["simulationSets"][-1][-1]

seqEnergy = { seq: {} for seq in sequences }
for seq in sequences:
    energyFree  = np.loadtxt(free+"/results/"+seq+"/energy.dat",skiprows=1)
    #For RANDOM: !!!!!
    #seqOld = seq
    #seq    = seq.split("_")[0]+seq.split("_")[1]+"_"+seq.split("_")[2]
    energyLeft  = np.loadtxt(left+"/results/"+seq+"/energy.dat",skiprows=1)
    energyRight = np.loadtxt(right+"/results/"+seq+"/energy.dat",skiprows=1)
    #seq = seqOld

    freeSteps = energyFree[:,0]
    leftSteps = energyLeft[:,0]
    rightSteps = energyRight[:,0]

    internalEnergyFree  = np.sum(energyFree[:,internalTerms],axis=1)
    internalEnergyLeft  = np.sum(energyLeft[:,internalTerms],axis=1)
    internalEnergyRight = np.sum(energyRight[:,internalTerms],axis=1)

    freeSteps = freeSteps[skipFrames:]
    leftSteps = leftSteps[skipFrames:]
    rightSteps = rightSteps[skipFrames:]

    internalEnergyFree = internalEnergyFree[skipFrames:]
    internalEnergyLeft = internalEnergyLeft[skipFrames:]
    internalEnergyRight = internalEnergyRight[skipFrames:]

    #plt.plot(freeSteps,internalEnergyFree,label="Free")
    #plt.plot(leftSteps,internalEnergyLeft,label="Left")
    #plt.plot(rightSteps,internalEnergyRight,label="Right")
    #plt.title(seq)
    #plt.legend()
    #plt.show()

    internalEnergyFreeMean  = np.mean(internalEnergyFree)
    internalEnergyLeftMean  = np.mean(internalEnergyLeft)
    internalEnergyRightMean = np.mean(internalEnergyRight)

    internalEnergyFreeErr  = np.std(internalEnergyFree)/np.sqrt(len(internalEnergyFree))
    internalEnergyLeftErr  = np.std(internalEnergyLeft)/np.sqrt(len(internalEnergyLeft))
    internalEnergyRightErr = np.std(internalEnergyRight)/np.sqrt(len(internalEnergyRight))

    seqEnergy[seq]["free"] = (internalEnergyFreeMean,internalEnergyFreeErr)
    seqEnergy[seq]["left"] = (internalEnergyLeftMean,internalEnergyLeftErr)
    seqEnergy[seq]["right"] = (internalEnergyRightMean,internalEnergyRightErr)

left = []
right = []

for seq in sequences:
    eFree = seqEnergy[seq]["free"][0]
    eLeft = seqEnergy[seq]["left"][0]
    eRight = seqEnergy[seq]["right"][0]

    errorFree = seqEnergy[seq]["free"][1]
    errorLeft = seqEnergy[seq]["left"][1]
    errorRight = seqEnergy[seq]["right"][1]

    l = eLeft - eFree
    r = eRight - eFree

    errorL = np.sqrt(errorFree**2 + errorLeft**2)
    errorR = np.sqrt(errorFree**2 + errorRight**2)

    #Convert form kcal/mol to kt (kT = 0.6 kcal/mol)

    l = l/0.6
    r = r/0.6

    errorL = errorL/0.6
    errorR = errorR/0.6

    left.append((l,errorL))
    right.append((r,errorR))

fig, ax = plt.subplots()

#xAxis -> seqName
#yAxis -> energies

xEntries = [seq[:-2] for seq in sequences]

plt.errorbar(xEntries,[l[0] for l in left],[l[1] for l in left],label="Left")
plt.errorbar(xEntries,[r[0] for r in right],[r[1] for r in right],label="Right")

plt.ylabel(r"$\Delta E_{\text{free}\rightarrow\text{nuclesosome}}$ (kT)")
plt.title(r"$\Delta E_{\text{free}\rightarrow\text{nuclesosome}}$ for different sequences")
plt.xlabel("Sequence")

plt.xticks(rotation=45)
plt.legend()
plt.show()

with open("energy.dat","w") as f:
    f.write("Sequence\tEnergyLeft\tErrorLeft\tEnergyRight\tErrorRight\n")
    for i in range(len(sequences)):
        f.write(sequences[i]+"\t"+str(left[i][0])+"\t"+str(left[i][1])+"\t"+str(right[i][0])+"\t"+str(right[i][1])+"\n")












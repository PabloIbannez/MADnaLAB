import os,sys
import json

import VLMP

from VLMP.utils.units import picosecond2KcalMol_A_time

from numpy import random

ps2AKMA = picosecond2KcalMol_A_time()

name = "FREE_RANDOM"
#Check if folder name exists
if os.path.exists(name):
    print("Folder %s already exists" % name)
    sys.exit(1)

seqLen = 147

nSimPerSeq = 1

dt       = 0.02
friction = 0.2

simulationTime = 500000.0
infoTime       = 10000.0
writeTime      = 100000.0
measureTime    = 1000.0
backupTime     = 10000.0

with open("./generation/data/randomSeq.json") as f:
#with open("./generation/data/sequences.json") as f:
    sequences_data = json.load(f)

sequences = {}
for seqAlias,seqProp in sequences_data.items():
    seqBase = seqProp[0]
    seqRep  = seqProp[1]
    sequences[seqAlias] = seqBase*seqRep
    if seqLen != -1:
        sequences[seqAlias] = sequences[seqAlias][:seqLen]
    print("Sequence %s: %s" % (seqAlias,sequences[seqAlias][:11]))

#################################################

dt       = dt*ps2AKMA
friction = friction/ps2AKMA

simulationSteps = int(simulationTime/dt)
infoSteps       = int(infoTime/dt)
writeSteps      = int(writeTime/dt)
measureSteps    = int(measureTime/dt)
backupSteps     = int(backupTime/dt)

print("Simulation steps: %d" % simulationSteps)
print("Info steps: %d" % infoSteps)
print("Write steps: %d" % writeSteps)
print("Backup steps: %d" % backupSteps)

#################################################

simulationPool = []
for seqAlias in sequences.keys():
    for i in range(nSimPerSeq):

        simName     = seqAlias + "_" + str(i)

        seq         = sequences[seqAlias]

        simulationPool.append({"system":[{"type":"simulationName","parameters":{"simulationName":simName}},
                                         {"type":"backup","parameters":{"backupIntervalStep":backupSteps}}],
                               "units":[{"type":"KcalMol_A"}],
                               "types":[{"type":"basic"}],
                               "ensemble":[{"type":"NVT","parameters":{"box":[10000.0,10000.0,10000.0],"temperature":300.0}}],
                               "integrators":[{"type":"BBK","parameters":{"timeStep":dt,"frictionConstant":friction,"integrationSteps":simulationSteps}}],
                               "models":[{"type":"MADna",
                                          "parameters":{"sequence":seq,
                                                        "debyeLength":7.8511}
                                          }],
                               "simulationSteps":[{"type":"saveState","parameters":{"intervalStep":writeSteps,
                                                                                    "outputFilePath":"output",
                                                                                    "outputFormat":"sp"}},
                                                  {"type":"thermodynamicMeasurement","parameters":{"intervalStep":measureSteps,
                                                                                     "outputFilePath":"thermo.dat"}},
                                                  {"type":"potentialEnergyMeasurement","parameters":{"intervalStep":measureSteps,
                                                                                                     "outputFilePath":"energy.dat"}},
                                                  {"type":"info","parameters":{"intervalStep":infoSteps}}]

                               })


vlmp = VLMP.VLMP()

vlmp.loadSimulationPool(simulationPool)
vlmp.distributeSimulationPool("none")
vlmp.setUpSimulation(name)


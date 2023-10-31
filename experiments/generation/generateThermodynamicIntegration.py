import os

import VLMP

from VLMP.utils.units import picosecond2KcalMol_A_time

import numpy as np
from numpy import random

ps2AKMA = picosecond2KcalMol_A_time()

outputName = "THERMODYNAMIC_INTEGRATION"

#Check if the output directory exists
if os.path.isdir(outputName):
    print("Output directory already exists. Exiting...")
    exit()

#Get file current directory
currentDirectory   = os.path.dirname(os.path.realpath(__file__))
nucleosomeFilePath = os.path.join(currentDirectory, './data/nucleosome.dat')

#Load nucleosome
nucleosomeIds = []
nucleosomePos = []
with open(nucleosomeFilePath, 'r') as f:
    for line in f:
        i,x,y,z = [float(x) for x in line.split()]
        i = int(i)
        nucleosomeIds.append(i)
        nucleosomePos.append([x,y,z])

dt       = 0.02
friction = 0.2

simulationTime = 200000.0
infoTime       = 1000.0
measureTime    = 1000.0
writeTime      = 1000.0
backupTime     = 10000.0

debyeLength    = 7.874

sequences = []

#Random sequence
seqLen = 147
bases  = ['A','T','C','G']
seq    = ''.join([random.choice(bases) for i in range(seqLen)])
sequences.append(seq)

#################################################

lambdaKs  = np.linspace(0.0, 1.0, 11)

#################################################

dt       = dt*ps2AKMA
friction = friction/ps2AKMA

simulationSteps = int(simulationTime/dt)
infoSteps       = int(infoTime/dt)
measureSteps    = int(measureTime/dt)
writeSteps      = int(writeTime/dt)
backupSteps     = int(backupTime/dt)

print("Simulation steps: %d" % simulationSteps)
print("Info steps: %d" % infoSteps)
print("Write steps: %d" % writeSteps)
print("Backup steps: %d" % backupSteps)

#################################################

simulationPool = []
for seq in sequences:
    for lambdaK in lambdaKs:

        simName     = seq + "_lambdaK_" + str(lambdaK)

        simulationPool.append({"system":[{"type":"simulationName","parameters":{"simulationName":simName}},
                                         {"type":"backup","parameters":{"backupIntervalStep":backupSteps}}],
                               "units":[{"type":"KcalMol_A"}],
                               "types":[{"type":"basic"}],
                               "ensemble":[{"type":"NVT","parameters":{"box":[10000.0,10000.0,10000.0],"temperature":300.0}}],
                               "integrators":[{"type":"BBK","parameters":{"timeStep":dt,"frictionConstant":friction,"integrationSteps":simulationSteps}}],
                               "models":[{"type":"MADna",
                                          "parameters":{"sequence":seq,
                                                        "debyeLength":debyeLength}
                                          }],
                               "modelOperations":[{"type":"setParticlePositions","parameters":{"positions":nucleosomePos,"ids":nucleosomeIds}}],
                               "modelExtensions":[{"type":"constraintParticlesPosition",
                                                   "parameters":{"K":lambdaK,
                                                                 "selection":{}}}],
                               "simulationSteps":[{"type":"saveState","parameters":{"intervalStep":writeSteps,
                                                                                    "outputFilePath":"output",
                                                                                    "outputFormat":"sp"}},
                                                  {"type":"potentialEnergyMeasurement","parameters":{"intervalStep":measureSteps,
                                                                                                     "outputFilePath":"energy.dat"}},
                                                  {"type":"info","parameters":{"intervalStep":infoSteps}}]

                               })


vlmp = VLMP.VLMP()

vlmp.loadSimulationPool(simulationPool)
vlmp.distributeSimulationPool("size", len(lambdaKs))
vlmp.setUpSimulation(outputName)


import os

import VLMP

from VLMP.utils.units import picosecond2KcalMol_A_time

import numpy as np
from numpy import random

ps2AKMA = picosecond2KcalMol_A_time()

outputName = "THERMODYNAMIC_INTEGRATION_BACKWARD"

#Check if the output directory exists
if os.path.isdir(outputName):
    print("Output directory already exists. Exiting...")
    exit()

#Get file current directory
currentDirectory   = os.path.dirname(os.path.realpath(__file__))
constraintFilePath = os.path.join(currentDirectory, './data/superHelix.dat')

#Load constraint
constraintIds = []
constraintPos = []
with open(constraintFilePath, 'r') as f:
    for line in f:
        i,x,y,z = [float(x) for x in line.split()]
        i = int(i)
        constraintIds.append(i)
        constraintPos.append([x,y,z])

dt       = 0.02
friction = 0.2

infoTime       = 1000.0
measureTime    = 1000.0
writeTime      = 1000.0
backupTime     = 1000.0

debyeLength    = 7.874

sequences = []

#Random sequence
seqLen = 147
bases  = ['A','T','C','G']
seq    = ''.join([random.choice(bases) for i in range(seqLen)])
sequences.append(seq)

#################################################

#Constraint K
K = 1.0 # Constraint potential is 0.5*lambda**2*K*(x-x0)^2+0.5*lambda**2*K*(y-y0)^2+0.5*lambda**2*K*(z-z0)^2

#Lambda
intervalLambda = 1000  # This is the interval dU/dlambda is computed, STEPS !!!!
stepLambda     = 20000 # STEPS!!!!!

lambdaValues = [0.0,0.001,0.002,0.005,0.008,0.01,0.02,0.03,0.04,0.05,0.06,0.08,0.1,0.15,0.2,0.3,0.4,0.5,1.0]

#################################################

dt       = dt*ps2AKMA
friction = friction/ps2AKMA

simulationSteps = len(lambdaValues)*stepLambda
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
    simulationPool.append({"system":[{"type":"simulationName","parameters":{"simulationName":seq}},
                                     {"type":"backup","parameters":{"backupIntervalStep":backupSteps}}],
                           "units":[{"type":"KcalMol_A"}],
                           "types":[{"type":"basic"}],
                           "ensemble":[{"type":"NVTlambda","parameters":{"box":[10000.0,10000.0,10000.0],"temperature":300.0,"lambda":lambdaValues[0]}}],
                           "integrators":[{"type":"BBK","parameters":{"timeStep":dt,"frictionConstant":friction,"integrationSteps":simulationSteps}}],
                           "models":[{"type":"MADna",
                                      "parameters":{"sequence":seq,
                                                    "debyeLength":debyeLength}
                                      }],
                           "modelExtensions":[{"type":"constraintParticlesListPositionLambda",
                                               "parameters":{"K":K,
                                                             "ids":constraintIds,
                                                             "positions":constraintPos,
                                                            }
                                              }],
                           "simulationSteps":[{"type":"saveState","parameters":{"intervalStep":writeSteps,
                                                                                "outputFilePath":"output",
                                                                                "outputFormat":"sp"}},
                                              {"type":"potentialEnergyMeasurement","parameters":{"intervalStep":measureSteps,
                                                                                                 "outputFilePath":"energy.dat"}},
                                              {"type":"thermodynamicIntegration","parameters":{"intervalStep":intervalLambda,
                                                                                               "outputFilePath":"ti.dat",
                                                                                               "stepLambda":stepLambda,
                                                                                               "lambdaValues":lambdaValues}},
                                              {"type":"info","parameters":{"intervalStep":infoSteps}}]

                           })


vlmp = VLMP.VLMP()

vlmp.loadSimulationPool(simulationPool)
vlmp.distributeSimulationPool("one")
vlmp.setUpSimulation(outputName)


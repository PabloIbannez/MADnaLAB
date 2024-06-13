import os,sys
import json

import VLMP

from VLMP.utils.units import picosecond2KcalMol_A_time

import numpy as np
from numpy import random

def helixEquation(t, pitch, radius, eps, boxZ):
    x = radius*np.cos(t)
    y = eps*radius*np.sin(t)
    z = pitch*t/(2.0*np.pi)
    z = z - boxZ/2.0
    return x,y,z

def computeTgivenContourLength(l, pitch, radius, eps):
    z = l/np.sqrt(1.0 + 4.0*np.pi**2*radius**2/pitch**2)
    t = 2.0*np.pi*z/pitch
    return t

ps2AKMA = picosecond2KcalMol_A_time()

if len(sys.argv) > 1:
    name = sys.argv[1]
else:
    name = "HELIX_CONSTRAINTS"

#Check if folder name exists
if os.path.exists(name):
    print("Folder %s already exists" % name)
    sys.exit(1)

seqLen = 147

nSimPerSeq = 1

dt       = 0.02
friction = 0.2

Ktime          = 20000.0
infoTime       = 100.0
writeTime      = 10000.0
measureTime    = 10000.0
backupTime     = 10000.0

with open("./generation/data/sequences.json") as f:
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

helixPitch       = 26.0
helixRadius      = 42.0
helixInnerRadius = 12.0

eps    = -1.0
nTurns =  125

################

KhelixConstraint = [0.01,0.1,1.0]

#################################################

dt       = dt*ps2AKMA
friction = friction/ps2AKMA

Ksteps = int(Ktime/dt)

simulationSteps = int(Ktime/dt)*len(KhelixConstraint)
infoSteps       = int(infoTime/dt)
measureSteps    = int(measureTime/dt)
writeSteps      = int(writeTime/dt)
backupSteps     = int(backupTime/dt)

print("Simulation steps: %d" % simulationSteps)
print("Info steps: %d" % infoSteps)
print("Write steps: %d" % writeSteps)
print("Backup steps: %d" % backupSteps)

#################################################

boxZ = helixPitch*nTurns

boxX = helixRadius*2.0*1.2

#Increase boxX until it is a multiple of boxZ
boxX = boxX + (boxZ - boxX%boxZ)
boxY = boxX

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
                               "ensemble":[{"type":"NVT","parameters":{"box":[boxX,
                                                                              boxY,
                                                                              boxZ],
                                                                       "temperature":300.0}}],
                               "integrators":[{"type":"BBK","parameters":{"timeStep":dt,"frictionConstant":friction,"integrationSteps":simulationSteps}}],
                               "models":[{"type":"MADna",
                                          "parameters":{"sequence":seq,
                                                        "debyeLength":7.8511}
                                          }],
                               "modelOperations":[{"type":"setCenterOfMassPosition",
                                                   "parameters":{
                                                       "position":[0.0,0.0,0.0],
                                                       "selection":{}
                                                       }
                                                  }],
                               "modelExtensions":[],
                               "simulationSteps":[{"type":"saveState","parameters":{"intervalStep":writeSteps,
                                                                                    "outputFilePath":"output",
                                                                                    "outputFormat":"sp"}},
                                                  {"type":"potentialEnergyMeasurement","parameters":{"intervalStep":measureSteps,
                                                                                                     "outputFilePath":"energy.dat"}},
                                                  {"type":"info","parameters":{"intervalStep":infoSteps}}]

                               })

        tStep = computeTgivenContourLength(32.0,helixPitch,helixRadius,eps)

        c = 1
        seqLen = len(seq)
        while c*10 < seqLen:
            t = (c-1)*tStep
            x,y,z = helixEquation(t,helixPitch,helixRadius,eps,boxZ)
            z = z + boxZ/2.0

            #Add the helix constraints to modelExtensions
            for kn,K in enumerate(KhelixConstraint):
                simulationPool[-1]["modelExtensions"].append(
                        {"name":f"helixConstraint_pairIndex{c}_K{K}",
                         "type":"constraintCenterOfMassPosition",
                         "parameters":{"position":[x,y,z],
                                       "startStep":kn*Ksteps,
                                       "endStep":kn*Ksteps + Ksteps,
                                       "selection":{"expression":{"basePairIndex":-c*10}},
                                       "r0":0.0,
                                       "K":K}
                        })
            c = c + 1

vlmp = VLMP.VLMP()

vlmp.loadSimulationPool(simulationPool)
vlmp.distributeSimulationPool("one")
vlmp.setUpSimulation(name)


import os,sys
import json

import VLMP

from VLMP.utils.units import picosecond2KcalMol_A_time

import numpy as np
from numpy import random

ps2AKMA = picosecond2KcalMol_A_time()

name = "HELIX_BOUNDARIES_RIGHT"
#Check if folder name exists
if os.path.exists(name):
    print("Folder %s already exists" % name)
    sys.exit(1)

nSimPerSeq = 1

dt       = 0.02
friction = 0.2

simulationTime = 500000.0
infoTime       = 10000.0
writeTime      = 100000.0
measureTime    = 1000.0
backupTime     = 10000.0

if len(sys.argv) != 2:
    print("Usage: python generateHelixBoundary.py <VLMPsession.json>")
    sys.exit(1)

with open(sys.argv[1]) as f:
    sequences_data = json.load(f)

sequences = {}
folderName = sequences_data["name"]
print("Session name: %s" % folderName)
for simS in sequences_data["simulationSets"]:
    if len(simS[3]) != 1:
        print("Only one sequence per simulation set is allowed")
        sys.exit(1)
    seqName  = simS[3][0].split("_")[0]
    filepath = folderName + "/" + simS[1]+"/backupEnd.json"
    sequences[seqName] = filepath
    print(seqName,filepath)

#################################################

helixPitch       = 26.0
helixRadius      = 42.0
helixInnerRadius = 12.0

eps    =  1.0
nTurns =  20

K = 1.0

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

boxZ = helixPitch*nTurns

boxX = helixRadius*2.0*1.2

#Increase boxX until it is a multiple of boxZ
boxX = boxX + (boxZ - boxX%boxZ)
boxY = boxX

nPointsHelix = 10000
nx = 500
ny = 500

dx = boxX/nx
nz = boxZ/dx

print("Box dimensions: %f %f %f" % (boxX,boxY,boxZ))
print("Number of points in the helix: %d" % nPointsHelix)
print("Grid dimensions: %d %d %d" % (nx,ny,nz))

#Check if nx,ny and nz are integers
if nx != int(nx) or ny != int(ny) or nz != int(nz):
    print("Grid dimensions are not integers. Exiting...")
    exit()
else:
    nx = int(nx)
    ny = int(ny)
    nz = int(nz)

#################################################

simulationPool = []
for seqAlias in sequences.keys():
    filepath = sequences[seqAlias]
    for i in range(nSimPerSeq):

        simName     = seqAlias + "_" + str(i)
        seq         = sequences[seqAlias]

        simulationPool.append({"system":[{"type":"simulationName","parameters":{"simulationName":simName}},
                                         {"type":"backup","parameters":{"backupIntervalStep":backupSteps}}],
                               "units":[{"type":"KcalMol_A"}],
                               "types":[{"type":"basic"}],
                               "ensemble":[{"type":"NVTlambda","parameters":{"box":[boxX,
                                                                                    boxY,
                                                                                    boxZ],
                                                                       "temperature":300.0,
                                                                       "lambda":1.0}}],
                               "integrators":[{"type":"BBK","parameters":{"timeStep":dt,
                                                                          "frictionConstant":friction,
                                                                          "integrationSteps":simulationSteps}}],
                               "models":[{"type":"FILE","parameters":{"inputFilePath":filepath,"removeInteractionsByType":["Set1"]}}],
                               "modelExtensions":[{"type":"helixBoundaries",
                                                   "parameters":{"helixPitch":helixPitch,
                                                                 "helixRadius":helixRadius,
                                                                 "helixInnerRadius":helixInnerRadius,
                                                                 "eps":eps,
                                                                 "nTurns":nTurns,
                                                                 "nPointsHelix":nPointsHelix,
                                                                 "nx":nx,"ny":ny,"nz":nz,
                                                                 "K":K}
                                                  }],
                               "simulationSteps":[{"type":"saveState","parameters":{"intervalStep":writeSteps,
                                                                                    "outputFilePath":"output",
                                                                                    "outputFormat":"sp"}},
                                                  {"type":"potentialEnergyMeasurement","parameters":{"intervalStep":measureSteps,
                                                                                                     "outputFilePath":"energy.dat"}},
                                                  {"type":"info","parameters":{"intervalStep":infoSteps}}]

                               })

vlmp = VLMP.VLMP()

vlmp.loadSimulationPool(simulationPool)
vlmp.distributeSimulationPool("none")
vlmp.setUpSimulation(name)


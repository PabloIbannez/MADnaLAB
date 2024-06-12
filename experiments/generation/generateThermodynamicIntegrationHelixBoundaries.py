import os

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

outputName = "THERMODYNAMIC_INTEGRATION_HELIX_BOUNDARIES"

#Check if the output directory exists
#if os.path.isdir(outputName):
#    print("Output directory already exists. Exiting...")
#    exit()

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

KhelixConstraint = [0.01,0.1,1.0,10.0]
Ksteps = 50000

#Lambda
intervalLambda = 1000  # This is the interval dU/dlambda is computed, STEPS !!!!
stepLambda     = 100000 # STEPS!!!!!

lambdaValues = [1.0,1.0,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.0,0.0,0.0,0.0]

#################################################

dt       = dt*ps2AKMA
friction = friction/ps2AKMA

simulationSteps = Ksteps*len(KhelixConstraint)+stepLambda*len(lambdaValues)
infoSteps       = int(infoTime/dt)
measureSteps    = int(measureTime/dt)
writeSteps      = int(writeTime/dt)
backupSteps     = int(backupTime/dt)

print("Simulation steps: %d" % simulationSteps)
print("Info steps: %d" % infoSteps)
print("Write steps: %d" % writeSteps)
print("Backup steps: %d" % backupSteps)

#################################################

helixPitch       = 26.0
helixRadius      = 42.0
helixInnerRadius = 12.0

eps    = -1.0
nTurns = 20

################

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

#Constraint K
K = 1.0


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
for seq in sequences:

    simulationPool.append({"system":[{"type":"simulationName","parameters":{"simulationName":seq}},
                                     {"type":"backup","parameters":{"backupIntervalStep":backupSteps}}],
                           "units":[{"type":"KcalMol_A"}],
                           "types":[{"type":"basic"}],
                           "ensemble":[{"type":"NVTlambda","parameters":{"box":[boxX,
                                                                                boxY,
                                                                                boxZ],
                                                                         "temperature":300.0,"lambda":lambdaValues[0]}}],
                           "integrators":[{"type":"BBK","parameters":{"timeStep":dt,"frictionConstant":friction,"integrationSteps":simulationSteps}}],
                           "models":[{"type":"MADna",
                                      "parameters":{"sequence":seq,
                                                    "debyeLength":debyeLength}
                                      }],
                           "modelOperations":[{"type":"setCenterOfMassPosition",
                                               "parameters":{
                                                   "position":[0.0,0.0,0.0],
                                                   "selection":{}
                                                   }
                                              }],
                           "modelExtensions":[
                                              {"type":"helixBoundaries",
                                               "parameters":{"startStep":Ksteps*len(KhelixConstraint),
                                                             "helixPitch":helixPitch,
                                                             "helixRadius":helixRadius,
                                                             "helixInnerRadius":helixInnerRadius,
                                                             "eps":eps,
                                                             "nTurns":nTurns,
                                                             "nPointsHelix":nPointsHelix,
                                                             "nx":nx,"ny":ny,"nz":nz,
                                                             "K":K}
                                              }
                                              ],
                           "simulationSteps":[{"type":"saveState","parameters":{"startStep":Ksteps*len(KhelixConstraint),
                                                                                "intervalStep":writeSteps,
                                                                                "outputFilePath":"output",
                                                                                "outputFormat":"sp"}},
                                              {"type":"potentialEnergyMeasurement","parameters":{"intervalStep":measureSteps,
                                                                                                 "outputFilePath":"energy.dat"}},
                                              {"type":"thermodynamicIntegration","parameters":{"startStep":Ksteps*len(KhelixConstraint),
                                                                                               "intervalStep":intervalLambda,
                                                                                               "outputFilePath":"ti.dat",
                                                                                               "stepLambda":stepLambda,
                                                                                               "lambdaValues":lambdaValues}},
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
vlmp.setUpSimulation(outputName)


import VLMP

from VLMP.utils.units import picosecond2KcalMol_A_time

from numpy import random
import hashlib

ps2AKMA = picosecond2KcalMol_A_time()

dt       = 0.02
friction = 0.2

boxZ                = 10000.0

compressionVelocity = -0.001
minPlatesSeparation =  25.0

simulationTime = 2000000.0
infoTime       = 10000.0
writeTime      = 10000.0
backupTime     = 10000.0

#
nSimPerSeq = 1

# Random sequence
sequences = []

seqLen = 300
bases  = ['A','T','C','G']

for i in range(1):
    seq    = ''.join([random.choice(bases) for i in range(seqLen)])
    sequences.append(seq)

#################################################

dt       = dt*ps2AKMA
friction = friction/ps2AKMA

simulationSteps = int(simulationTime/dt)
infoSteps       = int(infoTime/dt)
writeSteps      = int(writeTime/dt)
backupSteps     = int(backupTime/dt)

print("Simulation steps: %d" % simulationSteps)
print("Info steps: %d" % infoSteps)
print("Write steps: %d" % writeSteps)
print("Backup steps: %d" % backupSteps)

#################################################

simulationPool = []
for seq in sequences:
    for i in range(nSimPerSeq):

        seqHash = hashlib.sha1(seq.encode('utf-8')).hexdigest()
        simName  = seqHash + "_" + str(i)

        simulationPool.append({"system":[{"type":"simulationName","parameters":{"simulationName":simName}},
                                         {"type":"backup","parameters":{"backupIntervalStep":backupSteps}}],
                               "units":[{"type":"KcalMol_A"}],
                               "types":[{"type":"basic"}],
                               "ensemble":[{"type":"NVT","parameters":{"box":[1000.0,1000.0,boxZ],"temperature":300.0}}],
                               "integrators":[{"type":"BBK","parameters":{"timeStep":dt,
                                                                          "frictionConstant":friction,
                                                                          "integrationSteps":simulationSteps}}],
                               "models":[{"type":"MADna",
                                          "parameters":{"sequence":seq}
                                          }],
                               "modelOperations":[{"type":"setCenterOfMassPosition","parameters":{"position":[0.0,0.0,0.0],
                                                                                                  "selection":{}}}],
                               "modelExtensions":[{"type":"plates","parameters":{"platesSeparation":"auto",
                                                                                 "epsilon":1.0,"sigma":1.0,
                                                                                 "compressionVelocity":compressionVelocity,
                                                                                 "minPlatesSeparation":minPlatesSeparation,
                                                                                 "selection":{}}}],
                               "simulationSteps":[{"type":"saveState","parameters":{"intervalStep":writeSteps,
                                                                                    "outputFilePath":"output",
                                                                                    "outputFormat":"sp"}},
                                                  {"type":"thermodynamicMeasurement","parameters":{"intervalStep":infoSteps,
                                                                                     "outputFilePath":"thermo.dat"}},
                                                  {"type":"info","parameters":{"intervalStep":infoSteps}}]

                               })


vlmp = VLMP.VLMP()

vlmp.loadSimulationPool(simulationPool)
vlmp.distributeSimulationPool("size", nSimPerSeq)
vlmp.setUpSimulation("SURFACE")


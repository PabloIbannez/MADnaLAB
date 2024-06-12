import VLMP

from VLMP.utils.units import picosecond2KcalMol_A_time

from numpy import random

ps2AKMA = picosecond2KcalMol_A_time()

dt       = 0.02
friction = 0.2

simulationTime = 200000.0
infoTime       = 100.0
writeTime      = 10000.0
backupTime     = 10000.0

#sequences = {"AA":["CGCGAAAAAAAAAACGCG",10.8748],
#             "AC":["CGCGACACACACACCGCG",10.9725],
#             "AG":["CGCGAGAGAGAGAGCGCG",10.7969],
#             "AT":["CGCGATATATATATCGCG",10.8929],
#             "CG":["CGCGCGCGCGCGCGCGCG",10.8929],
#             "GG":["CGCGGGGGGGGGGGCGCG",10.8748]}

sequences = {"TESIS":["CGCGAATTCGCG",10.8748]}

forces = [1.0]
nSimPerSeq = 10

#################################################

dt       = dt*ps2AKMA
friction = friction/ps2AKMA

simulationSteps = int(simulationTime/dt)
infoSteps       = int(infoTime/dt)
writeSteps      = int(writeTime/dt)
backupSteps     = int(backupTime/dt)

print(simulationSteps,infoSteps,writeSteps,backupSteps)
input()

print("Simulation steps: %d" % simulationSteps)
print("Info steps: %d" % infoSteps)
print("Write steps: %d" % writeSteps)
print("Backup steps: %d" % backupSteps)

#################################################

simulationPool = []
for seqAlias in sequences.keys():
    for force in forces:
        for i in range(nSimPerSeq):

            simName     = seqAlias + "_" + str(force) + "_" + str(i)

            seq         = sequences[seqAlias][0]
            debyeLength = sequences[seqAlias][1]

            simulationPool.append({"system":[{"type":"simulationName","parameters":{"simulationName":simName}},
                                             {"type":"backup","parameters":{"backupIntervalStep":backupSteps}}],
                                   "units":[{"type":"KcalMol_A"}],
                                   "types":[{"type":"basic"}],
                                   "ensemble":[{"type":"NVT","parameters":{"box":[10000.0,10000.0,10000.0],"temperature":300.0}}],
                                   "integrators":[{"type":"GFJ","parameters":{"timeStep":dt,"frictionConstant":friction,"integrationSteps":simulationSteps}}],
                                   "models":[{"type":"MADna",
                                              "parameters":{"sequence":seq,
                                                            "debyeLength":debyeLength}
                                              }],
                                   "modelExtensions":[{"type":"constantForceBetweenCentersOfMass",
                                                       "parameters":{"force":force,
                                                                     "selection1":{"expression":{"basePairIndex":2}},
                                                                     "selection2":{"expression":{"basePairIndex":-2}}}}],
                                   "simulationSteps":[{"type":"saveState","parameters":{"intervalStep":writeSteps,
                                                                                        "outputFilePath":"output",
                                                                                        "outputFormat":"dcd"}},
                                                      {"type":"thermodynamicMeasurement","parameters":{"intervalStep":infoSteps,
                                                                                         "outputFilePath":"thermo.dat"}},
                                                      {"type":"info","parameters":{"intervalStep":infoSteps}}]

                                   })


vlmp = VLMP.VLMP()

vlmp.loadSimulationPool(simulationPool)
#vlmp.distributeSimulationPool("size", nSimPerSeq)
vlmp.distributeSimulationPool("property", ["topology","forceField","DH","parameters","debyeLength"])
vlmp.setUpSimulation("MADNA_TEST")


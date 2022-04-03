import sys

import MADnaLAB

import MADnaLAB.utils
import MADnaLAB.sequencesGeneration
import MADnaLAB.simulationPool
import MADnaLAB.simulationSets

cfg = MADnaLAB.utils.loadCONF("dev.cfg")

seqList = MADnaLAB.sequencesGeneration.generateSequencesListFromCONF(cfg)

simPool = MADnaLAB.simulationPool.generateSimulationPoolFromSequencesList(seqList)

MADnaLAB.simulationPool.loadModificationsIntoSimulationPoolFromCONF(simPool,cfg)

simPoolSet = MADnaLAB.simulationPool.splitSimulationPoolAccordingToMaxNumberOfBasePairs(simPool,int(sys.argv[1]))

MADnaLAB.simulationSets.generateSimulationSetsFromListOfSimulationPoolUsingMADnaModel(cfg,simPoolSet)


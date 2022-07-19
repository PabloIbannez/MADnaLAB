import sys,os
import pathlib

import itertools

import logging

import psutil
import signal

import multiprocessing
import subprocess

import time
import datetime

################

gpuIDList = [0]

MADnaLAB = os.environ['MADNALABPATH']+"/bin/MADnaLAB"

#################

def signalHandler(sig,frame):
    logging.info("Signal {} detected. Exiting ...".format(sig))

    mainProcess = psutil.Process(os.getpid())
    for child in mainProcess.children(recursive=True):
        child.kill()

    sys.exit(0)

def associateCPUtoGPU(gpuIDList):
    cpuName = multiprocessing.current_process().name
    try:
        cpuID = int(cpuName[cpuName.find('-') + 1:]) - 1
    except:
        cpuID = 0 
    gpuID = gpuIDList[cpuID % len(gpuIDList)]
    return gpuID

def runSimulation(options, gpuIDList):

    gpuId = associateCPUtoGPU(gpuIDList)
    cudaFlag = 'CUDA_VISIBLE_DEVICES={}'.format(gpuId)

    simulationPath = str(pathlib.Path(options).parent)+"/"
    
    fout = open(simulationPath+"out.log","w")
    ferr = open(simulationPath+"err.log","w")
    
    sim = " ".join([cudaFlag, " ".join([MADnaLAB,options])])
    simReturn = subprocess.run(sim, stdout=fout, stderr=ferr, shell=True)
    while(simReturn.returncode == 1):
        logging.warning("Simulation {} crashed. Recovering from backup.".format(sim))
        simReturn = subprocess.run(sim, stdout=fout, stderr=ferr, shell=True)

    if(simReturn.returncode == 2):
        logging.error("Simulation {} crashed. Recovery is NOT possible.".format(sim))

    fout.close()
    ferr.close()

    return simReturn

def runSimulationsSet(simulationsSet, gpuIDList):

    with multiprocessing.Pool(processes=len(gpuIDList)) as pool:
        out = list(pool.starmap(runSimulation,
                    zip(simulationsSet,
                        itertools.repeat(gpuIDList))))

    return out

if __name__ == "__main__":

    simulationsSetFile = sys.argv[1]

    logging.basicConfig(filename='launchSimulations.log',
                        filemode='w',
                        level=0,
                        format='%(asctime)s %(levelname)s: %(message)s', 
                        datefmt='%m/%d/%Y %H:%M:%S')

    logging.info("Start")
    logging.info("pid: {}".format(os.getpid()))

    signal.signal(signal.SIGINT,signalHandler)
    signal.signal(signal.SIGTERM,signalHandler)

    simulationsSet = []
    with open(simulationsSetFile,"r") as f:
        for line in f:
            simulationsSet.append(line.rstrip())

    st = time.time()
    out = runSimulationsSet(simulationsSet,gpuIDList)
    logging.info("Simulations set finished. Total time: {}".format(str(datetime.timedelta(seconds=(time.time() - st)))))

    for i in out:
        if(i.returncode!=0):
            logging.error("Something went wrong for simulation: {}".format(i))
    
    logging.info("End")

    sys.exit(0)

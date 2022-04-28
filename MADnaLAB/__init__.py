import os

CREATE_MOLECULE_PATH = os.getenv('UAMMDPATH') + '/extensions/structured/Tools/MADna/CreateMolecule/CreateMolecule.py'

SEQ_COMPULSORY_INFO     = ["alias","seq","copies"]
SIMPOOL_COMPULSORY_INFO = ["alias","seq","simId","filePath"]

AVAILABLE_MODIFICATIONS = ["externalForceBtwCOM",
                           "externalTorqueBtwCOM",
                           "externalForce",
                           "externalTorque",
                           "constraintsDistanceBtwCOM",
                           "constraintsPositionOfCOM",
                           "constraintsPositionOfBeads"]

DEBYE_FACTOR = 4.0
CUTOFF_VERLET_DIFF = 5.0

DEFAULT_SIMULATION_OPTIONS = {
        'boxSize':[10000.0, 10000.0, 10000.0], 
        'T':300.0, 
        'h':0.001, 
        'nStepsSteepestDescent':0, 
        'nStepsSteepestDescentProgressInterval':100, 
        'maxObjectiveForce':5.0, 
        'outPutFilePath':'output', 
        'outPutFormat':'lammpstrj', 
        'dt':0.01, 
        'frictionConstant':0.2, 
        'dielectricConstant':78.3, 
        'debyeLength':7.88, 
        'cutOffDstWCA':25.0
        }

TIME2STEPS_OPTIONS = {
        'simulationTime':'nSteps', 
        'infoTime':'nStepsInfoInterval', 
        'writeTime':'nStepsWriteInterval', 
        'backupTime':'nStepsBackupInterval'
        }


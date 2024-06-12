import sys
import json

if len(sys.argv) < 4:
    print("Usage: python generateContinuation.py <simulation snapshot> <additional time (ps)> <output folder>")
    sys.exit(1)

simulationSnapshot = sys.argv[1]
additionalTime     = float(sys.argv[2])
outputFolder       = sys.argv[3]

with open(simulationSnapshot, 'r') as f:
    snapshot = json.load(f)

print(snapshot["system"])

simulationPool = []


#import VLMP


import sys
import json
import random

seed = 123456
random.seed(seed)

base = ['A', 'C', 'G', 'T']

N      = int(sys.argv[1])
seqLen = int(sys.argv[2])

results = {}
for i in range(N):
    name = "randomSeq_" + str(i)
    seq = ""
    for j in range(seqLen):
        seq += random.choice(base)
    results[name] = [seq,1]

with open('randomSeq.json', 'w') as f:
    json.dump(results, f, indent=4)






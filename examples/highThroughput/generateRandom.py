import random

seqListFile = "seq.list"

seqLength  = 100
nSequences = 100

seed = 123456789

with open(seqListFile,"w") as f:
    
    random.seed(seed)            
    
    B = ["A","T","C","G"]
    
    aliasIndex=0
    for n in range(nSequences):
        seq = ""
        for l in range(seqLength):
            seq = seq + random.choice(B)
        f.write("{:} {:}\n".format("rnd_"+str(aliasIndex),seq))
        aliasIndex+=1


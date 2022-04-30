import sys,os
import gpuMultiprocessing

with open("id.log","w") as f:
    f.write(str(os.getpid())+"\n")

gpu_id_list = [0]

command_queue = []

with open(sys.argv[1],"r") as f:
    for i,line in enumerate(f):
        command_queue.append("../../bin/MADnaLAB "+line.rstrip()+" 2>stderr{}.log 1>stdout{}.log".format(i,i))
        #print(line.rstrip())

#for c in command_queue:
#    print(c)

gpuMultiprocessing.queue_runner(command_queue, gpu_id_list,
                                env_gpu_name='CUDA_VISIBLE_DEVICES',
                                processes_per_gpu=1, allowed_restarts=0)


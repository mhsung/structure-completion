#!/usr/bin/python

##########
#
#  Batch symmetry detection
#
##########

from __future__ import division
import os, string, glob, sys
import fas
import itertools
import multiprocessing as mp
import time


## Variables
dataDir = ""
expDir = ""
execType = "symmetry_detection"

def run_exec_cmd(args):
    return exec_cmd(*args)

def exec_cmd(mName, cmd, logFile):
    print("Job [" + mName + "] Started.")
    print(cmd)
    os.system("time " + cmd + " > " + logFile)
    print("Job [" + mName + "] Finished.")


## Parse arguments
if len(sys.argv) < 3:
	print("Usage: " + sys.argv[0] + " dataset_path experiment_path [options]")
	print("Options:")
	exit()
else:
    dataDir = sys.argv[1]
    expDir = sys.argv[2]
    i = 3
    while i < len(sys.argv):
        print("[ERROR] Unknown argument " + sys.argv[i])
        exit()


# Prepare files
mListTest = [os.path.basename(x) for x in glob.glob(dataDir + "/*.off")]
#print("\n".join(mListTest))
binDir = "/home/mhsung/shape2pose/code/external/gaps/bin/x86_64/"


## Batch jobs
cwd = os.getcwd()

os.chdir(expDir)
jobIDs = []
removeFiles = []
cmdList = []
count = 0

for mName in mListTest:
    mesh_name = os.path.splitext(mName)[0]

    if not os.path.isfile("output/" + mesh_name + "/" + mesh_name + "_input.pts"):
        continue

    print(mesh_name)
    
    count = count + 1
    cmd = binDir + "msh2pln" + " "
    cmd += dataDir + "/" + mName + " "
    cmd += "symmetry_detection/" + mesh_name + "/" + mesh_name + ".txt "
    cmd += "-v -input_points output/" + mesh_name + "/" + mesh_name + "_input.pts"
    
    scriptFile = "script/" + execType + "/" + mName
    
    if not os.path.isdir("script"):
        os.mkdir("script")
    if not os.path.isdir("script/" + execType):
        os.mkdir("script/" + execType)
    if not os.path.isdir("symmetry_detection"):
        os.mkdir("symmetry_detection")
    if not os.path.isdir("symmetry_detection/" + mesh_name):
        os.mkdir("symmetry_detection/" + mesh_name)
        
    logFile = scriptFile + ".out"
    f = open(scriptFile+".sh", 'w'); f.write(cmd); f.close();
    os.system("chmod 777 "+scriptFile+".sh");
    cmdList.append((mName, cmd, logFile))
    print(cmd)
    os.system(cmd)

'''
pool = mp.Pool(processes = numProcessors)
rs = pool.imap_unordered(run_exec_cmd, cmdList)
pool.close()
t = time.time()

while True:
    completed = rs._index
    if (completed == count): break
    #print "[" + execType + "]: " + expDir
    print "Waiting for", count - completed, "tasks to complete..."
    time.sleep(5)
    
elapsed = time.time() - t
print 'Elapsed time: %s' % (elapsed)
'''

os.chdir(cwd)

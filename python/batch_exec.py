#!/usr/bin/python

##########
#
#  Batch prediction
#
##########

from __future__ import division
import os, string, glob, sys
import fas
import itertools
import multiprocessing as mp
import time


## Variables
execType = ""
dataDir = ""
expDir = ""
disallowParallel = True
disallowRandomView = False
numProcessors = 4

def run_exec_cmd(args):
    return exec_cmd(*args)

def exec_cmd(mName, cmd, logFile):
    print("Job [" + mName + "] Started.")
    print(cmd)
    os.system("time " + cmd + " > " + logFile)
    print("Job [" + mName + "] Finished.")


## Parse arguments
if len(sys.argv) < 4:
	print("Usage: " + sys.argv[0] + " exec_type dataset_path experiment_path [options]")
	print("Options:")
	print("   -n")
	print("   -disallowParallel")
	print("   -disallowRandomView")
	exit()
else:
    execType = sys.argv[1]
    dataDir = sys.argv[2]
    expDir = sys.argv[3]
    i = 4
    while i < len(sys.argv):
        if (sys.argv[i] == "-n"):
            numProcessors = int(sys.argv[i+1])
            i = i + 1
        elif (sys.argv[i] == "-disallowParallel"):
            disallowParallel = True
        elif (sys.argv[i] == "-disallowRandomView"):
            disallowRandomView = True
        else:
            print("[ERROR] Unknown argument " + sys.argv[i])
            exit()
        i = i + 1


# Prepare files
mListTest = [os.path.basename(x) for x in glob.glob(dataDir + "/*.off")]
#print("\n".join(mListTest))
binDir = "../../build/OSMesaViewer/build/Build/bin/"


## Batch jobs
cwd = os.getcwd()

os.chdir(expDir)
jobIDs = []
removeFiles = []
cmdList = []
count = 0

for mName in mListTest:
	count = count + 1
	cmd = binDir + "OSMesaViewer" + " "
	cmd += "--flagfile=arguments.txt" + " "
	cmd += "--mesh_filename=" + mName + " "
	cmd += "--run_" + execType + " "

	if not disallowRandomView:
		cmd += "--occlusion_pose_filename=\"\" "
		cmd += "--random_view_seed=" + str(count) + " "

	scriptFile = "script/" + execType + "/" + mName

	if not disallowParallel:
		jobIDs += fas.ScheduleJob(cmd, mName, scriptFile)
		removeFiles += fas.TmpFilesNames(scriptFile)
	else:
		if not os.path.isdir("script"):
			os.mkdir("script")
        if not os.path.isdir("script/" + execType):
            os.mkdir("script/" + execType)

        logFile = scriptFile + ".out"
        f = open(scriptFile+".sh", 'w'); f.write(cmd); f.close();
        os.system("chmod 777 "+scriptFile+".sh");
        cmdList.append((mName, cmd, logFile))

if not disallowParallel:
    fas.WaitForJobsInArray(jobIDs)
else:
    pool = mp.Pool(processes = numProcessors)
    rs = pool.imap_unordered(run_exec_cmd, cmdList)
    pool.close()
    t = time.time()

    while True:
        completed = rs._index
        if (completed == count): break
        print "[" + execType + "]: " + expDir
        print "Waiting for", count - completed, "tasks to complete..."
        time.sleep(5)

    elapsed = time.time() - t
    print 'Elapsed time: %s' % (elapsed)

os.chdir(cwd)
#print('Generating HTML table...')
#os.system('./generate_html_table.py ' + expDir)

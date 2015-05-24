#!/usr/bin/python

##########
#
#  Batch prediction
#
##########

import os, string, glob, sys
sys.path.append(os.path.realpath("./"+os.path.dirname(sys.argv[0])))
import fas


## Variables
disallowParallel = False
meshFilename = ""
expDir = ""


## Parse arguments
if len(sys.argv) < 3:
    print("Usage: " + sys.argv[0] + "mesh_filename experiment_path [options]")
    print("Options:")
    print("   -disallowParallel")
    exit()
else:
    meshFilename = sys.argv[1]
    expDir = sys.argv[2]
    i = 3
    while i < len(sys.argv):
        if (sys.argv[i] == "-disallowParallel"):
            disallowParallel = True
        else:
            print("[ERROR] Unknown argument " + sys.argv[i])
            exit()
        i = i + 1


binDir = "../../build/OSMesaViewer/build/Build/bin/"


## Batch jobs
os.chdir(expDir)
jobIDs = []
removeFiles = []

numExperiments = 20
for count in range(numExperiments):
    cmd = binDir + "OSMesaViewer" + " "
    cmd += "--flagfile=arguments.txt" + " "
    cmd += "--mesh_filename=" + meshFilename + " "
    cmd += "--run_prediction" + " "

    scriptFile = "script/" + meshFilename

    if not disallowParallel:
        jobIDs += fas.ScheduleJob(cmd, meshFilename, scriptFile)
        removeFiles += fas.TmpFilesNames(scriptFile)
    else:
        if not os.path.isdir("script"):
            os.mkdir("script")

        print("\n====================")
        print("(%d / %d)" % (count, numExperiments))
        print("Job [" + meshFilename + "] Started.")
        os.system("time " + cmd + " > script/" + meshFilename + ".out")
        print("Job [" + meshFilename + "] Finished.")
        print("====================\n")

    os.system("mv output output_" + str(count))

fas.WaitForJobsInArray(jobIDs)

os.system("rm output/ -rf")
os.mkdir("output")
os.system("mv output_* output/ -R")

#!/usr/bin/python

##########
#
# Shape2Pose scripts library: file management, job scheduling (if running parallel), scikit-learn interfaces
#
##########
import os, string, decimal, glob, re, subprocess, shlex, sys, datetime, time, getpass

######   Job Scheduling     ######
maxJobs=10;
usrName = getpass.getuser();

def ScheduleJob(cmd, jobName, scriptFile):
	jobID = [];
	if IsParallel():
		cmdToSchedule = cmd;
		assert(not os.path.isfile(scriptFile)); # cannot schedule if file exists
		path = os.path.dirname(scriptFile);
		if not os.path.isdir(path):
			os.makedirs(path);
		f = open(scriptFile+".sh", 'w'); f.write(cmd); f.close();
		os.system("chmod 777 "+scriptFile+".sh");
		qsubParams = "-X -v display";
		qsubParams += "-l mem=8000mb,cput=4:00:00 ";	# 4 hours
		qsubParams += "-N "+jobName + " ";
		qsubParams += "-d "+os.getcwd()+" -o "+scriptFile+".log -e "+scriptFile+".err";
		cmd = "qsub "+qsubParams+" "+scriptFile+".sh ";

		if int(NumMyJobs()) >= int(maxJobs):
			while (int(NumMyJobs()) >= .9 * int(maxJobs)):
				print("Too many jobs ("+str(NumMyJobs())+" / "+str(.9 * int(maxJobs))+")... waiting 5s");
				time.sleep(5);
		
		time.sleep(.01);
		subPID = subprocess.Popen("qsub "+qsubParams+" "+scriptFile+".sh ", \
								  shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		idStr, errors = subPID.communicate();
		if errors!="":
			print(errors);
			print("[WARNING] Error executing: "+scriptFile+".sh.\n    Trying to re-run. NJobs="+str(NumMyJobs()));
			time.sleep(1);
			subPID = subprocess.Popen("qsub "+qsubParams+" "+scriptFile+".sh ", \
									  shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			idStr, errors = subPID.communicate();
			if errors!="":
				time.sleep(1);
				subPID = subprocess.Popen("qsub "+qsubParams+" "+scriptFile+".sh ", \
										  shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				idStr, errors = subPID.communicate();
				if errors!="":
					print("    [ERROR] Third run failed as well... quitting...");
					assert(False);

		jobID =[idStr.split(".")[0] +" " +cmdToSchedule];
	else:
		os.system(cmd);
	return jobID;

def TmpFilesNames(prefix):
	return [prefix+".sh", prefix+".err", prefix+".log"];

def NumMyJobs():
	subPID = subprocess.Popen("qstat | grep -c "+usrName, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	nJobs, errors = subPID.communicate();
	if nJobs=='\n' or nJobs=='':
		return 0;
	else:
		return int(nJobs);

def IsParallel():
	subPID = subprocess.Popen("type qsub", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	qsubExists, errors = subPID.communicate();
	if len(qsubExists)>=4 and qsubExists[:4]=="qsub":
		return True;	# can qsub jobs
	else:
		return False;	# cannot qsub jobs

def WriteJobsFile(jobsIDs, filename):
	if filename != "none":
		f = open(filename, "w");
		f.write("\n".join(jobsIDs));
		f.close();

def ReadJobsFile(filename):
	if (IsParallel()):
		return open(filename, "r").read().strip().split("\n");
	else:
		return [];

def WaitForAllJobs():
	if IsParallel():
		while NumMyJobs()>0:
			print("        Waiting for jobs: "+str(NumMyJobs()))
			time.sleep(5)

def WaitForJobsInFile(filename):
	if (IsParallel()):
		if (filename=="none"):
			WaitForAllJobs();
		else:
			WaitForJobsInArray(ReadJobsFile(filename));

def WaitForJobsInArray(jobsIds):
	if IsParallel():		
		# list of jobs and their status
		jobID = [];
		jobCMD = [];
		jobStatus = [];
		jobInitRunTime = [];
		lastVisit = [];
		jobRestarts = [];
		# jobID map
		jMap = dict();
		
		# initialize job data
		for job in jobsIds:
			jobSplitStr = job.split(" ");
			jobID.append(jobSplitStr[0]);
			jobCMD.append(" ".join(jobSplitStr[1:]));
			jobStatus.append("Q");
			jobInitRunTime.append(-1);
			jMap[jobSplitStr[0]] = len(jobID)-1;
			lastVisit.append(0);
			jobRestarts.append(0);
		
		# wait for jobs to be done...
		counterSleeps = 0;
		sleepTime = 5;	# in seconds
		
		visitCounter = 0;
		while len(jobID)>0:   # jobs remain
			nRunning = 0;
			nTotal = 0;
			# check if idle jobs can be udpated
			# go over running jobs - set check id
			subPID = subprocess.Popen("qstat | grep "+usrName+" ", \
									  shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			jobList, errors = subPID.communicate();
									  
			jList = jobList.split("\n");
			visitCounter = visitCounter + 1;
			for j in jList:   # check which jobs remain
			  if (j==""):
				  continue;
			  nTotal = nTotal + 1;
			  jobDescrWithEmpty = j.strip().split(" ");
			  jobDescr = [];
			  for val in jobDescrWithEmpty:
				  if val != "":
					  jobDescr.append(val);
			  jid = jobDescr[0].split(".")[0].strip();
			  jName = jobDescr[1].strip();
			  jTimeUse = jobDescr[3].strip();
			  jstat = jobDescr[4].strip();
			  if jid in jMap:
				  lid = jMap[jid];
				  lastVisit[lid] = visitCounter;
				  if jstat != "R" and jstat != "Q" and jstat != "E":
					  print("[WARNING] Invalid jstat: "+jstat);
				  if jstat != jobStatus[lid]:
					  if (jobStatus[lid]+jstat)=="QR":
						  jobStatus[lid] = jstat;
						  jobInitRunTime[lid] = 0;
					  elif (not jobStatus[lid]+jstat)=="RE":
						  print("[WARNING] Unexpected state change: "+jobStatus[lid]+" -> "+jstat);
				  if (jobStatus[lid]=="R"):
					  nRunning = nRunning + 1;

			# if some job was never visited - done. erase
			toRemove = [];
			for j in range(0, len(lastVisit)):
			  if (visitCounter!=lastVisit[j]):
				  toRemove.append(j);

			# erase all jobs that were visited
			if (len(toRemove)>0):
			  toRemove = sorted(toRemove, reverse=True);
			  for j in toRemove:
				  lastVisit.pop(j);
				  jobID.pop(j);
				  jobCMD.pop(j);
				  jobStatus.pop(j);
				  jobInitRunTime.pop(j);
				  jobRestarts.pop(j);
			  jMap = dict();
			  for j in range(0, len(jobID)):
				  jMap[jobID[j]] = j;

			counterSleeps += sleepTime;
			if (len(jMap)>0):
			  print("        Waiting for jobs: "+str(nRunning)+" / "+str(len(jMap))+" / "+str(nTotal));
			  time.sleep(sleepTime)


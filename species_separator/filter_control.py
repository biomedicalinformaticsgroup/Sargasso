#!/usr/bin/env python

# This script reads the read mapping BAM files & parallelises its filtering

# Parameters:
# 1 - Block directory filepath
# 2 - Output directory
# 3 - No. Threads
# 4 - Species 1
# 5 - Species 2

import sys
import subprocess
import os
import math
from os import listdir
from os.path import isfile, join

# Check whether output dir exists
def checkDir(dir):
        if not os.path.isdir(dir):
                print "Filepath Parameter does not exist"
		quit()
                
# Check how many instances of a process are running at the present time
def checkProcesses(processes,max):
        busy = 0
        index = 0
        i = 0
        threshold = max
        while True:
                if i<len(processes):
                        if processes[i].poll() == None:
                                busy = busy + 1
                        elif processes[i].poll() == 0:
                                del processes[i]
                                i = i - 1
                        if busy >= threshold:
                                return False
                        i = i + 1
                else:
                        break
        return True

# turn file list into dictionary; d[species1file] = species2file
def dictionary_files(files,sp1,sp2):
	fileDictionary = {}
	files = sorted(files)
	for file in files:
		sections = file.split("_")
		if sections[1] == sp1:
			sections[1] = sp2
			fileDictionary[file]=sections[0]+"_"+sections[1]+"_"+sections[2]+"_"+sections[3]
	return fileDictionary

# Cycles through previously indexed chunks of read data & assigns it to the filter subprocess script
# Chunks; for the time being I'm considering each chunk to be 10 read pairs organised as a list of dictionaries
def run_processes(params):
	chunkDir = params[1]
        outDir = params[2]
	noThreads = params[3]
	species1 = params[4]
	species2 = params[5]
	# keep track of all processes
        allProcesses = []
        processNo = -1
	maxProcessNo = int(noThreads)
	medianProcess = math.floor(maxProcessNo/2)
	chunkdex = 0
	# get block files
	blockFiles = [ f for f in listdir(chunkDir) if isfile(join(chunkDir,f)) ]
	# remove species 2 block files
	filePairs = dictionary_files(blockFiles,species1,species2)
        # cycle through chunks
        for file1 in filePairs.keys():
                processNo = processNo + 1
		# create input paths
		sp1in = chunkDir + "/"+file1
                sp2in = chunkDir + "/"+filePairs[file1]
		# create output paths
		sp1out = outDir + "/"+species1+"-block_"+str(processNo)+"-filtered"
		sp2out = outDir + "/"+species2+"-block_"+str(processNo)+"-filtered"
                commands = ["filter_reads_parallel",species1,sp1in,os.path.abspath(sp1out),species2,sp2in,os.path.abspath(sp2out)]
                proc = subprocess.Popen(commands)
                allProcesses.append(proc)
        # check all processes finished
        free = checkProcesses(allProcesses,1)
        while free == False:
		print "Waiting for Threads"
                allProcesses[0].wait()
                free = checkProcesses(allProcesses,1)
	print "Filtering Complete"
	# NEED TO CONCATENATE THE OUTPUT FILES
        # DONE

def validate_params(params):
	# check output dir exists & create if not
	#checkDir(params[1])
	#checkDir(params[2])
	#if not params[3].isnumeric():
	#	print "Number of threads specified is not a number"
	#	quit()
	#if !os.path.isfile(params[2]):
	#	print "The filepath for the first BAM file points to a non-existant file"
	#	return False
	#if !os.path.isfile(params[3]):
        #        print "The filepath for the second BAM file points to a non-existant file"
        #        return False
	return True

def filter_control():
	if validate_params(sys.argv):
		run_processes(sys.argv)

filter_control()

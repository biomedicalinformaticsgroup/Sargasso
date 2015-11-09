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
import pysam
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
		sections = split(file,"_")
		if sections[1] == sp1:
			sections[1] = sp2
			fileDictionary[file]=sections[0]+"_"+sections[1]+"_"+sections[2]"_"+sections[3]
	return fileDictionary

# Cycles through previously indexed chunks of read data & assigns it to the filter subprocess script
# Chunks; for the time being I'm considering each chunk to be 10 read pairs organised as a list of dictionaries
def run_processes(params):
	chunkDir = params[1]
        outDir = params[2]
	noThreads = params[3]
	species1 = params[4]
	species2 = params[5]
	#bam1 = pysam.AlignmentFile(params[2],"rb")
	#bam2 = pysam.AlignmentFile(params[3],"rb")
	# keep track of all processes
        allProcesses = []
        processNo = -1
	maxProcessNo = noThreads
	medianProcess = math.floor(maxProcessNo/2)
	chunkdex = 0
	#chinkSize = 20 #in lines
	# get block files
	blockFiles = [ f for f in listdir(chunkDir) if isfile(join(chunkDir,f)) ]
	# remove species 2 block files
	filePairs = dictionary_files(blockFiles)
        # cycle through chunks
        for file1 in filePairs.keys():
                processNo = processNo + 1
                commands = ["python",os.path.abspath("./filter_sample_reads"),file1,filePairs[file1],os.path.abspath(outDir)]
                proc = subprocess.Popen(commands)
                allProcesses.append(proc)
                free = checkProcesses(allProcesses,maxProcessNo)
                while free == False:
                       allProcesses[medianProcess].wait()
                        free = checkProcesses(allProcesses,maxProcessNo)
	#for thread in range(noThreads):
	#	chunkdex += 1
         #       #checkDir(outDir)
          #      #outName = outDir+"filter_out_"+str(chunkdex)
	#	sp1ChunkFile = outDir+"/Blocks/"+sample+"_"+species1+"_BLOCK_"+str(thread)
	#	sp2ChunkFile = outDir+"/Blocks/"+sample+"_"+species2+"_BLOCK_"+str(thread)
         #       processNo = processNo + 1
          #      commands = ["python",os.path.abspath("./filter_sample_reads"),sp1ChunkFile,sp2ChunkFile,sample,os.path.abspath(outDir)]
           #     proc = subprocess.Popen(commands)
            #    allProcesses.append(proc)
                #free = checkProcesses(allProcesses,maxProcessNo)
                #while free == False:
                #	allProcesses[9].wait()
                #        free = checkProcesses(allProcesses,maxProcessNo)
        # check all processes finished
        print "Finishin Up"
        free = checkProcesses(allProcesses,1)
        while free == False:
                allProcesses[0].wait()
                free = checkProcesses(allProcesses,1)
        # DONE

# read raw data 
#def read_n_index(path):
	

#def main(params):
	# read data in from BAM files & index into subprocessable chunks
	#chunks = read_n_index(params[1])
	# call the filter script with each chunk
	#run_processes(chunks,params[1])
	#DONE

def validate_params(params):
	# check output dir exists & create if not
	checkDir(params[1])
	checkDir(params[2])
	if not params[3].isnumeric():
		print "Number of threads specified is not a number"
		quit()
	#if !os.path.isfile(params[2]):
	#	print "The filepath for the first BAM file points to a non-existant file"
	#	return False
	#if !os.path.isfile(params[3]):
        #        print "The filepath for the second BAM file points to a non-existant file"
        #        return False
	#return True

if validate_params(sys.argv):
	#main(sys.argv)
	run_processes(sys.argv)

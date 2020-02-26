#!/usr/bin/env python

#######################################################
#
# constructConfigFile.py
#
# This script inputs 2 parameters: 1) a text file containing a base directory
# and a list of genome fasta and annotation gff files, followed by 2) a
# filename where output is to be written. Then, constructConfigFile.py
# constructs a configuration file for CPGwrapper.py to drive HTP 
# execution of compareGeneProfiles_main.py.  The configuration file will be
# constructed by enumerating a non-redundant list of pairwise jobs to be run.
# i.e., for genomes A, B, and C, this script will generate the config file
# for running A-B, A-C, B-C, but not B-A, C-A, or C-B.  For thoroughness in
# future, this may be modified.
# 
# Input file format example with 3 jobs:
#
#/home/zhou4/BacGenomeStudies/PAK1 
#/subdir1a/genomeFile1.fasta /subdir1b/annotationFile1.gff
#/subdir2b/genomeFile2.fasta /subdir2b/annotationFile2.gff
#/subdir3c/genomeFile3.fasta /subdir3b/annotationFile3.gff
#
# Note that there are no blank lines, the fasta file is listed first,
# and there is a single space preceeding the gff file.
#
# Programmer:  Carol L. Ecale Zhou
# Last update: 21 February 2020
#
#################################################################
'''
'''

import sys, os, re, string, copy
from subprocess import call

#### FILES

CODE_BASE_DIR = os.environ["CGP_CODE_BASE_DIR"]
CGP_OUT_DIR   = os.environ["PHATE_PIPELINE_OUTPUT_DIR"]
inFile = ""   # user designated
outFile = ""  # user designated
logFile = os.path.join(CODE_BASE_DIR, "constructConfigFile.log")
configFile = os.path.join(CODE_BASE_DIR, "CGPMwrapper.config")
LOGFILE = open(logFile,"w")
LOGFILE.write("%s\n" % ("Begin log file"))

#### CONSTANTS

ACCEPTABLE_ARG_COUNT = (2,3) # "help" or infile and outfile name expected

HELP_STRING = "This script inputs 2 parameters: 1) a text file containing a base directory and a list of genome fasta and annotation gff files, followed by 2) a filename where output is to be written. Then, constructConfigFile.py constructs a configuration file for CGPMwrapper.py to drive HTP execution of compareGeneProfiles_main.py. The configuration file will be constructed by enumerating a non-redundant list of pairwise jobs to be run. i.e., for genomes A, B, and C, this script will generate the config file for running A-B, A-C, and B-C, but not B-A, C-A, or C-B. (In theory, A-B gives the same result as B-A, though I have not tested this in multiple cases. For thoroughness in future, this should be tested and the procedure modifed if needed.)\n"

USAGE_STRING = "Usage:  constructConfigFile.py <listFile> <outFile>\nType constructConfigFile.py help for additional information.\n"

INPUT_STRING = "The input file should comprise a text file listing the common directory where all files are found, followed by a list of genome fasta and annotation gff files, separated by a single space, followed by a common directory where results are to be written.\nExample:\n/home/user/analysis\n/genomes/Ba_Ames_chromosome.fasta /annotations/Ba_Ames_chromosome.gff\n/genomes/Ba_H9401_chromosome.fasta /annotations/Ba_H9401_chromosome.gff\n/genomes/Ba_CDC_chromosome.fasta /annotations/Ba_CDC_chromosome.gff\n/home/user/analysis/results/\nSee also cgpNxN_sample.config for a sample input file to script constructConfigFile.py\n"

##### Get command-line arguments

argCount = len(sys.argv)
if argCount in ACCEPTABLE_ARG_COUNT:
    match = re.search("help", sys.argv[1].lower())
    if match:
        print (HELP_STRING)
        exit(0)
    match = re.search("input", sys.argv[1].lower())
    if match:
        print (INPUT_STRING)
        exit(0)
    match = re.search("usage", sys.argv[1].lower())
    if match:
        print (USAGE_STRING)
        exit(0)
    else:
        inFile = sys.argv[1]
        outFile = sys.argv[2]
        LOGFILE.write("%s%s%s%s\n" % ("Parameter files are, input file ", inFile, " and output (config) file ", outFile))
else:
    print (HELP_STRING)
    print (USAGE_STRING)
    print (INPUT_STRING)
    exit(0)

##### Parse input file; construct config file 

INFILE = open(inFile,"r")
CONFIG_FILE = open(outFile,"w")

baseDir = ""
fileList = [] # list of inFile dicts
genome1 = ""
genome2 = ""
annotation1 = ""
annotation2 = ""
compareDir = ""
first = True

LOGFILE.write("%s%s\n" % ("Reading input config file ",inFile))
fLines = INFILE.read().splitlines()
numLines = len(fLines)
LOGFILE.write("%s%s\n" % ("Number of lines in input file is ", numLines)) 
baseDir = fLines[0]

LOGFILE.write("%s%s\n" % ("Reading names of genome and annotation files and writing output config file",outFile))
print ("Writing output config file,",outFile) 
jobNumber = 0  # i.e., comparison number
FIRST = True
for i in xrange(1, len(fLines)):
    print ("i is", i)
    (genome1,annotation1) = fLines[i].split(' ')
    print ("genome1 is", genome1)
    for j in xrange(i+1, len(fLines)):
        jobNumber += 1 
        print ("jobNumber is", jobNumber)
        resultDir = baseDir + '/c' + str(jobNumber)
        (genome2,annotation2) = fLines[j].split(' ')
        print ("genome2 is", genome2)
        print ("resultDir is", resultDir)
        CONFIG_FILE.write("%s\n" % (baseDir))
        CONFIG_FILE.write("%s\n" % (genome1))
        CONFIG_FILE.write("%s\n" % (genome2))
        CONFIG_FILE.write("%s\n" % (annotation1))
        CONFIG_FILE.write("%s\n" % (annotation2))
        CONFIG_FILE.write("%s\n\n" % (resultDir))  
LOGFILE.write("%s%s\n" % ("Number of jobs in config file = ", jobNumber))
print ("Done! Config file is ", outFile)

##### Clean up

INFILE.close()
CONFIG_FILE.close()
LOGFILE.close()


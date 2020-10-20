#!/usr/bin/env python

################################################################
#
# Name: constructPPcgpmConfigFile.py
#
# Programmer: Carol L. Ecale Zhou
#
# Last update:  16 October 2020
#
# Description:  Constructs a config file for input to ppCGPMwrapper.py.
#    Config file is constructed by reading the directory names produced
#    by CGPMwrapper.py. These directories contain output files produced
#    by compareGeneProfiles_main.py.  This config file constructed here
#    tells ppCGPMwrapper.py where to find the .report files that need
#    to be post-processed. This code should be run in the working 
#    directory where the current analysis is being done.
#
################################################################
# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE GPL3 LICENSE. SEE INCLUDED FILE GPL-3.PDF FOR DETAILS.

import sys
import os
import re
import subprocess

dateTime = "0:0:0::0:0:0"

# Boobleans
PHATE_PROGRESS = False
PHATE_MESSAGES = False
PHATE_WARNINGS = False
PHATE_PROGRESS_STRING = os.environ["PHATE_PHATE_PROGRESS"]
PHATE_MESSAGES_STRING = os.environ["PHATE_PHATE_MESSAGES"]
PHATE_WARNINGS_STRING = os.environ["PHATE_PHATE_WARNINGS"]
if PHATE_PROGRESS_STRING.lower() == 'true':
    PHATE_PROGRESS = True
if PHATE_MESSAGES_STRING.lower() == 'true':
    PHATE_MESSAGES = True
if PHATE_WARNINGS_STRING.lower() == 'true':
    PHATE_WARNINGS = True

##### FILES
CODE_BASE_DIR = os.environ["CGP_CODE_BASE_DIR"]
inConfig  = os.path.join(CODE_BASE_DIR, "constructConfigFile.config")    # default: same file as input to constructConfigFile.py
outfile   = os.path.join(CODE_BASE_DIR, "constructPPcgpmConfigFile.out") # capture standard output from current program
outConfig = os.path.join(CODE_BASE_DIR, "ppCGPMwrapper.config")          # an outfile: constucted for input to ppCGPMwrapper.py (next process)
logfile   = os.path.join(CODE_BASE_DIR, "constructPPcgpmConfigFile.log")
LOGFILE = open(logfile,"w")
LOGFILE.write("%s%s\n" % ("Begin log file ",dateTime))

##### PATTERNS

p_config = re.compile('\.config')

##### CONSTANTS

HELP_STRING = "This code constructs a config file comprising input to ppCGPMwrapper.py.\nRun this code in the directory where output files from CGPMwrapper.py have been written.\nThis is the directory where the current analysis is being done.\nInput to this program is the config file that was passed to CGPMwrapper.py \(should be called \"CGPMwrapper.config\"\)\nType: program.py usage|input for more information\n"

USAGE_STRING = "Usage: python constructPPcgpmConfigFile.py <input_configFile>\nNote:  The input_configFile is optional. If not provided, \nthe default will be: \'./constructConfigFile.config\'"

INPUT_STRING = "Input:  None required, but you may provide a config filename if desired\n"

ACCEPTABLE_ARG_COUNT = (1,2)

#DEBUG = True
DEBUG = False

REPORT_FILE = "compareGeneProfiles_main.report"

##### VARIABLES

myDict = {
    "example" : "myStuff",
    }

##### GET INPUT PARAMETERS

argCount = len(sys.argv)
if argCount in ACCEPTABLE_ARG_COUNT:
    if argCount == 2:  # if no argument was passed then do nothing, if 2...
        match = re.search("help", sys.argv[1].lower())
        if match:
            print (HELP_STRING)
            LOGFILE.close(); exit(0)
        match = re.search("input", sys.argv[1].lower())
        if match:
            print (INPUT_STRING)
            LOGFILE.close(); exit(0)
        match = re.search("usage", sys.argv[1].lower())
        if match:
            print (USAGE_STRING)
            LOGFILE.close(); exit(0)
        inConfig = sys.argv[1]  # user has provided input config filename
else:
    print (USAGE_STRING)
    LOGFILE.write("%s\n" % ("Incorrect number of command-line arguments provided"))
    LOGFILE.close(); exit(0)

# Check files
fileError = False
match = re.search(p_config, inConfig)
if not match:
    fileError = True
if fileError:
    if PHATE_WARNINGS:
        print ("cgp_constructPPcgpmConfigFile says, ERROR: Check the formats of your input file(s):")
        print ("    config file is", inConfig)
        print (USAGE_STRING)
    LOGFILE.close(); exit(0)

# Open files
try:
    INFILE = open(inConfig,"r")
except IOError as e:
    fileError = True
    print (e)
try:
    OUTFILE = open(outfile,"w")
except IOError as e:
    fileError = True
    print (e)
try:
    OUT_CONFIG = open(outConfig,"w")
except IOError as e:
    fileError = True
    print (e)
if fileError:
    LOGFILE.close(); exit(0)


##### BEGIN MAIN 

if PHATE_PROGRESS:
    print("cgp_constructPPcgpmConfigFile says, Begin construction of PPcgpm config file.")

# Get base directory from config file
fLines = INFILE.read().splitlines()
baseDir = fLines[0]

# Get list of results directories from previous (compareGeneProfiles_main.py) calculations
command = "ls . | grep \'Results_\'"  # Get list of Results directories in current directory
proc = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
(out, err) = proc.communicate()
if err:
    print ("cgp_constructPPcgpmConfigFile says, ERROR:", err)
myResult = []
myResult = out.split('\n')

# Walk through each directory name, append the name of the report file, and write to output contig file
if PHATE_PROGRESS:
    print("cgp_constructPPcgpmConfigFile says, walking through Results directory, writing config file.")
for resultDir in myResult:
    match = re.search('\w',resultDir)
    if match:  # Process .report file in this directory
        ppInFile = resultDir + "/compareGeneProfiles.report"
        OUT_CONFIG.write("%s\n" % (ppInFile))

##### CLEAN UP

INFILE.close()
OUT_CONFIG.close()
OUTFILE.close()
dateTime = "0:0:0::0:0:9"
LOGFILE.write("%s%s\n" % ("Processing complete ",dateTime))
LOGFILE.close()
if PHATE_PROGRESS:
    print("cgp_constructPPcgpmConfigFile says, Construction of PPcgpm config file complete.")

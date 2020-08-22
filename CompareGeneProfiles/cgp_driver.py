#!/usr/bin/env python

#######################################################
#
# cgp_driver.py - command-line version 1.0.0 (/multiPhATE2/)
#
# Programmer:  Carol L. Ecale Zhou
# Last update: 21 August 2020
#
# This script inputs a config file, cgpNxN.config, which lists the
# base directory containing input fasta and annotation files for
# processing through compareGeneProfiles_main.py.
#
# cgp_driver.py reads the input config file and handles all downstream
# processing of the CompareGeneProfiles pipeline. Overall the pipeline
# processes are as follows:
# 0) Take as input, file cgpNxN.config, then...
# 1) Run code constructConfigFile.py to generate file, CGPMwrapper.config
# 2) Run code CGPMwrapper.py with input file CGPMwrapper.config
#    This executes each pairwise comparison... 
#    Generating the following output files in each pairwise comparison directory (cn):
#    - compareGeneProfiles_main.log
#    - compareGeneProfiles_main.report
#    - compareGeneProfiles_main.summary 
#    - blast result files for each binary comparison
#    Generating a set of files in each genome's subdirectory (gn):
#    - gene fasta
#    - protein fasta
#    - gene blast database files
#    - protein blast database files
# 3) Run code constructPPcgpmConfigFile.py to create ppCGPMwrapper.config
# 4) Run code ppCGPMwrapper.py with input file ppCGPMwrapper.config 
#
#################################################################
# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.PDF FOR DETAILS.

import sys, os, re, string, copy
import datetime, time
from subprocess import call

##### CONSTANTS and CONFIGURABLES
CODE_BASE_DIR = os.environ["PHATE_BASE_DIR"] + "CompareGeneProfiles/"     # Location of this and subordinate codes
CONSTRUCT_CONFIG_FILE_CODE = os.path.join(CODE_BASE_DIR, "cgp_constructConfigFile.py")
COMPARE_GENE_PROFILES_WRAPPER_CODE = os.path.join(CODE_BASE_DIR, "cgp_wrapper.py")

# Booleans
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

# Set environmental variables for CGP's dependent codes
os.environ["CGP_CODE_BASE_DIR"] = CODE_BASE_DIR

ACCEPTABLE_ARG_COUNT = (3,) #  

HELP_STRING = "This code is a wrapper for running the compareGeneProfiles pipeline over several genome comparisons.\nInput the fully qualified pathname to the project directory where the code should execute.\n"

INPUT_STRING = "Provide as input the fully qualified directory path to the users project directory.\nPrepare your config file to specify the following information:\nLine 1: The fully qualified path of the project directory\nLines 2-n: name of each genome and its corresponding annotation file, separated by a single space,\nprepended with the name of the subdirectory for this genome-annotation pair.\n"

USAGE_STRING = "Usage:  cgp_driver.py project_directory_pathname"

##### Variables

userProjectDir = os.environ["PHATE_PIPELINE_OUTPUT_DIR"] # This is where CGP output is to be written 

#### FILES

if PHATE_PROGRESS:
    print("cgp_driver says, Begin cgp_driver code")
logFile = os.path.join(CODE_BASE_DIR, "cgp_driver.log")
LOGFILE = open(logFile,"a")
LOGFILE.write("%s\n" % ("==================================="))
today = os.popen('date')
LOGFILE.write("%s%s\n" % ("Begin log file at ",today.read()))
in_configFilename    = "cgpNxN.config"
cgpm_configFilename = "cgp_wrapper.config"
cgpm_configFile = ""
cgpThreads = 0  # By default, no parallel execution of CGP code

##### Get command-line arguments
if PHATE_PROGRESS:
    print("cgp_driver says, Getting command-line arguments")
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
    # If not a help string, argument should be absolute path to input configuration file 
    in_configFile = sys.argv[1]
    if len(sys.argv) >= 3:
        cgpThreads    = sys.argv[2]
else:
    print (USAGE_STRING)
    exit(0)

##### Set absolute path/file for config files
projectDirectory = os.path.dirname(in_configFile) 
cgpm_configFile = os.path.join(projectDirectory, cgpm_configFilename)
LOGFILE.write("%s%s\n" % ("User project directory is ",projectDirectory))
LOGFILE.write("%s%s\n" % ("Input config file is ", in_configFile))
LOGFILE.write("%s%s\n" % ("Output config file is ", cgpm_configFile))

##### Parse config file; construct lists of BASE_DIRS and FILES

if PHATE_PROGRESS:
    print("cgp_driver says, Reading in_configFile",in_configFile)
IN_CONFIG_FILE = open(in_configFile,"r")
today = os.popen('date')
LOGFILE.write("%s%s\n" % ("Compare Gene Profiles pipeline execution start time is ",today.read()))
LOGFILE.write("%s%s%s%s\n" % ("Executing Step 1: Running code constructConfigFile.py with ", in_configFile, " to generate file ", cgpm_configFile))
if PHATE_PROGRESS:
    print ("cgp_driver says, Executing Step 1...")
command = "python " + CONSTRUCT_CONFIG_FILE_CODE + " " + in_configFile + " " + cgpm_configFile
os.system(command)
if PHATE_PROGRESS:
    print ("cgp_driver says, Step 1 complete.")

LOGFILE.write("%s%s\n" % ("Executing Step 2: Running code ",COMPARE_GENE_PROFILES_WRAPPER_CODE))
LOGFILE.write("%s%s\n" % ("...with parameter project directory: ",projectDirectory))
LOGFILE.write("%s\n" % ("...to generate report and accessory files..."))
today = os.popen('date')
LOGFILE.write("%s%s\n" % ("Executing Step 2 at ",today.read()))
if PHATE_PROGRESS:
    print ("cgp_driver says, Executing Step 2...")
command = "python " + COMPARE_GENE_PROFILES_WRAPPER_CODE + " " + projectDirectory + ' ' + cgpThreads
LOGFILE.write("%s%s\n" % ("Running command: ",command))
os.system(command)
if PHATE_PROGRESS:
    print ("cgp_driver says, Step 2 complete.")

today = os.popen('date')
LOGFILE.write("%s%s\n" % ("Step 2 execution complete at ",today.read()))
if PHATE_PROGRESS:
    print ("cgp_driver says, Execution of ppCGPM (steps 3 and 4) are currently manual")

##### Clean up

LOGFILE.close()



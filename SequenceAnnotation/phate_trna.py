#################################################################################
# Module: phate_trna.py
#
# Programmer:  Carol L. Ecale Zhou
#
# Most recent update: 13 January 2021
#
# Module comprising classes and data structures for predicting tRNA genes. 
#
# Classes and Methods:
#    class trna
#       setParameters
#       runTrnaScan
#       printAll2file
#       printParameters
#       printResults
#       printAll
#
#################################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE GPL3 LICENSE. SEE INCLUDED FILE GPL-3.PDF FOR DETAILS.

import re, os, copy
import subprocess
import phate_annotation as annotation

CODE_BASE_DIR  = ""
OUTPUT_DIR     = ""

# Verbosity

PHATE_WARNINGS = False
PHATE_MESSAGES = False
PHATE_PROGRESS = False

PHATE_WARNINGS_STRING = os.environ["PHATE_PHATE_WARNINGS"]
PHATE_MESSAGES_STRING = os.environ["PHATE_PHATE_MESSAGES"]
PHATE_PROGRESS_STRING = os.environ["PHATE_PHATE_PROGRESS"]

if PHATE_WARNINGS_STRING.lower() == 'true':
    PHATE_WARNINGS = True

if PHATE_MESSAGES_STRING.lower() == 'true':
    PHATE_MESSAGES = True

if PHATE_PROGRESS_STRING.lower() == 'true':
    PHATE_PROGRESS = True

# Override
PHATE_PROGRESS = True
PHATE_WARNINGS = True
PHATE_MESSAGES = True

DEBUG = False
#DEBUG = True

# Code name
TRNA_SCAN_CODE_NAME = os.environ["PHATE_TRNA_SCAN_CODE_NAME"]

# Error codes
ERROR_CODE_1 = "trnascan-se execution failed" 
SUCCESS = 0

# Locations and files
PHATE_PIPELINE_OUTPUT_DIR = os.environ["PHATE_PIPELINE_OUTPUT_DIR"] 

# Acceptable arguments for trnascan 
USE_INFERNAL    = "-I"
BOTH_COMPONENTS = "-H"
QUIET_MODE      = "-q"

# Class comparison organizes all genomic data and performs comparisons, ultimately yielding
#   homology groups for each gene/protein in a reference genome.
class trna(object):
    
    def __init__(self):

        self.codeName                    = TRNA_SCAN_CODE_NAME  # "tRNAscan-SE" or "trnascan-se", depends on the installation/OS; modify in multiPhate.py 
        self.codeVersion                 = "2"       
        self.organismType                = "-B"     # default; Options:  -B (bacteria), -E (eukaryote), -A (archea), -M (mitochondria), -O (other)
        self.useInfernal                 = ""       # set to "-I" if input to setParameters() is True; use -I for bacteriophage 
        self.inputFile                   = ""       # PipelineInput/genomeFilename set by setParameters()
        self.outputFile                  = ""       # PipelineOutput/genomeName/filename set by setParameters()
        self.trnaStructureFile           = ""       # PipelineOutput/genomeName/filename set by setParameters()
        self.statisticsFile              = ""       # PipelineOutput/genomeName/filename set by setParameters()
        self.showBothStructureComponents = ""       # Show both primary and secondary structure components to covariance model bit scores
        self.quiteMode                   = ""       # set to "-q" if input to setParameters() is True

    def setParameters(self,kvargs):
        if PHATE_MESSAGES:
            print("phate_trna says, Running setParameters.")
        if isinstance(kvargs,dict):
            if "codeVersion" in list(kvargs.keys()):
                self.codeVersion = kvargs["codeVersion"]
            if "organismType" in list(kvargs.keys()):
                organismType = kvargs["organismType"].lower()
                if re.search("archea",organismType):
                    self.organismType = "-A"
                if re.search("bacteria",organismType):
                    self.organismType = "-B"
                if re.search("eukaryot",organismType):
                    self.organismType = "-E"
                if re.search("mitochondria",organismType):
                    self.organismType = "-M"
                if re.search("other",organismType):
                    self.organismType = "-O"
            if "useInfernal" in list(kvargs.keys()):
                if kvargs["useInfernal"] == True: 
                    self.useInfernal = USE_INFERNAL 
            if "inputFile" in list(kvargs.keys()):
                self.inputFile = kvargs["inputFile"]
            if "outputFile" in list(kvargs.keys()):
                self.outputFile = kvargs["outputFile"]
            if "trnaStructureFile" in list(kvargs.keys()):
                self.trnaStructureFile = kvargs["trnaStructureFile"]
            if "statisticsFile" in list(kvargs.keys()):
                self.statisticsFile = kvargs["statisticsFile"]
            if "showBothStructureComponents" in list(kvargs.keys()):
                if kvargs["showBothStructureComponents"] == True:
                    self.showBothStructureComponents = BOTH_COMPONENTS
            if "quietMode" in list(kvargs.keys()):
                if kvargs["quietMode"] == True:
                    self.quietMode = QUIET_MODE
        return

    def runTrnaScan(self):
        if PHATE_MESSAGES:
            print("phate_trna says, Running runTrnaScan.")
        # Construct command
        command = self.codeName + ' ' 
        if self.inputFile != "":
            command += self.inputFile + ' ' + self.organismType + ' '
        if self.outputFile != "":
            command += "-o " + self.outputFile + ' '
        if self.trnaStructureFile != "":
            command += "-f " + self.trnaStructureFile + ' '
        if self.statisticsFile != "":
            command += "-m " + self.statisticsFile + ' '
        if self.useInfernal != "":
            command += self.useInfernal + ' '
        if self.showBothStructureComponents != "":
            command += self.showBothStructureComponents + ' '
        if self.quietMode != "":
            command += self.quietMode + ' '
        # Run command
        try:
            if PHATE_PROGRESS:
                print ("phate_trna says, Running trnascan-se")
            # Clean out old files, if exist
            if os.path.exists(self.outputFile):
                rm_command = 'rm ' + self.outputFile
                os.system(rm_command)
            if os.path.exists(self.trnaStructureFile):
                rm_command = 'rm ' + self.trnaStructureFile
                os.system(rm_command)
            if os.path.exists(self.statisticsFile):
                rm_command = 'rm ' + self.statisticsFile
                os.system(rm_command)
            # Run transcan-se
            os.system(command)     
        except:
            print ("phate_trna says, ERROR: trnascan-se execution failed.")
            return ERROR_CODE_1
        return SUCCESS

    def printAll2file(self,FILE_H):
        if PHATE_MESSAGES:
            print("phate_trna says, Running printAll2file.")
        FILE_H.write("%s\n" % ("=========tRNAs============="))
        return

    def printParameters(self):
        if PHATE_MESSAGES:
            print("phate_trna says, Running printParameters.")
        print("=========parameters========")
        print("codeName:",self.codeName,", version",self.codeVersion)
        print("organism type:",self.organismType)
        print("useInfernal:",self.useInfernal)
        print("inputFile:",self.inputFile)
        print("outputFile:",self.outputFile)
        print("trnaStructureFile:",self.trnaStructureFile)
        print("statisticsFile:",self.statisticsFile)
        print("showBothStructureComponents",self.showBothStructureComponents)
        print("quietMode:",self.quietMode)
        return

    def printResults(self):
        if PHATE_MESSAGES:
            print("phate_trna says, Running printResults.")
        print("=========tRNA Genes========")
        print("nothing to show just yet")
        return

    def printAll(self):
        if PHATE_MESSAGES:
            print("phate_trna says, Running printAll.")
        self.printParameters()
        self.printResults()

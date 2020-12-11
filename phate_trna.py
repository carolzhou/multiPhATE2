#################################################################################
# Module: phate_trna.py
#
# Programmer:  Carol L. Ecale Zhou
#
# Most recent update: 27 October 2020
#
# Module comprising classes and data structures for comparing genomes
#
# Classes and Methods:
#    class trna
#       setParameters
#       runTrnaScan
#       write2gffFile
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

CODE_BASE_DIR    = ""
OUTPUT_DIR       = ""

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

# Error codes
ERROR_CODE_1 = "trnascan-se execution failed" 
SUCCESS = 0

# Locations and files
PHATE_PIPELINE_OUTPUT_DIR       = os.environ["PHATE_PIPELINE_OUTPUT_DIR"] 

# Acceptable arguments
USE_INFERNAL    = "-I"
BOTH_COMPONENTS = "-H"
QUIET_MODE      = "-q"

# Class comparison organizes all genomic data and performs comparisons, ultimately yielding
#   homology groups for each gene/protein in a reference genome.
class trna(object):
    
    def __init__(self):

        self.codeName                    = "trnascan-se"
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
                print ("phate_trna says, Running trnascan-se; command is ",command)
                if os.path.exists(self.outputFile):
                    rm_command = 'rm ' + self.outputFile
                    os.system(rm_command)
                if os.path.exists(self.trnaStructureFile):
                    rm_command = 'rm ' + self.trnaStructureFile
                    os.system(rm_command)
                if os.path.exists(self.statisticsFile):
                    rm_command = 'rm ' + self.statisticsFile
                    os.system(rm_command)
            os.system(command)     
        except:
            print ("phate_trna says, ERROR: trnascan-se execution failed.")
            return ERROR_CODE_1
        return SUCCESS

    def write2gffFile(self,infile):   # FILE_H should be abs_path/phate_sequenceAnnotation_main.gff
        # INFILE contains sequence annotations and exists populated or was created anew by phate_runPipeline.py
        contig = ""; previousContig = ""; header = ""
        if os.path.exists(self.outputFile):
            fileSize = os.path.getsize(self.outputFile)
            if fileSize == 0:
                if PHATE_PROGRESS:
                    print("phate_trna says, trnaOutFile is empty: no tRNA genes to add to GFF output")
            else:
                # Open Files; capture data
                TRNA_OUTFILE_H = open(self.outputFile,'r')
                tLines = TRNA_OUTFILE_H.read().splitlines()
                TRNA_OUTFILE_H.close()
                INFILE_H = open(infile,'a') # File might already hold annotations; append to existing data
                # If starting with an empty gff file, then write header
                fileSize = os.path.getsize(infile)
                if fileSize == 0:
                    header = "##gff-version 3"
                    INFILE_H.write("%s\n" % (header))
                # Parse data in trnascan results file
                START_DATA = False
                for tLine in tLines:
                    if START_DATA:
                        # Parse data
                        fields = tLine.split('\t')
                        sequenceName = fields[0]
                        tRNAnumber   = fields[1]
                        start        = fields[2]
                        end          = fields[3]
                        tRNAtype     = fields[4]
                        antiCodon    = fields[5]
                        intronStart  = fields[6]
                        intronEnd    = fields[7]
                        infScore     = fields[8]
                        note         = fields[9]
                        # Reformulate data and load into gff file
                        sequenceName = sequenceName.rstrip()
                        contig = sequenceName
                        tRNAname = 'ID=' + sequenceName + '_' + tRNAnumber + ';'
                        if intronStart == 0:
                            intronStart = '.'
                        if intronEnd == 0:
                            intronEnd = '.'
                        gffLine = sequenceName + '\t' + self.codeName + '\ttRNA\t' + start + '\t' + end + '\t' + infScore + '\t' + intronStart + '\t' + intronEnd + '\t' + tRNAname + note
                        # Print header, if required
                        if contig != previousContig:
                            header = "##sequence-region;contig " + contig + " 0 0"
                            INFILE_H.write("%s\n" % (header))
                        INFILE_H.write("%s\n" % (gffLine))
                        previousContig = contig
                    else:
                        if re.search('^-',tLine):
                            START_DATA = True
                INFILE_H.close()
        return

    def printAll2file(self,FILE_H):
        FILE_H.write("%s\n" % ("=========tRNAs============="))
        return

    def printParameters(self):
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
        print("=========tRNA Genes========")
        print("nothing to show just yet")
        return

    def printAll(self):
        self.printParameters()
        self.printResults()

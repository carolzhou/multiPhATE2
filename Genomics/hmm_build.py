#################################################################################
# Module: hmm_build.py
#
# Programmer:  Carol L. Ecale Zhou
#
# Most recent update: 19 November 2020
#
# Module comprising classes and data structures for building HMM profiles 
#
# Classes and Methods:
#    class build1hmm
#       setParameters
#       runHmmbuild
#       build1hmm
#       printParameters
#       printResults
#       printAll
#       printAll2file
#    class build 
#       performBuild
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
#PHATE_PROGRESS = True 
#PHATE_WARNINGS = True
#PHATE_MESSAGES = True

DEBUG = False
#DEBUG = True

# Error codes
ERROR_CODE_1 = "" 
ERROR_CODE_2 = "" 
SUCCESS = 0

# Acceptable Clustal Omega output formats that are recognized by Hmmbuild 
#LUSTALO_FORMAT = "fasta"      #
#CLUSTALO_FORMAT = "clustal"    # 
#CLUSTALO_FORMAT = "phylip"     # This one doesn't work on server, fnt input; complains that couldn't read format
#CLUSTALO_FORMAT = "selex"      # This one doesn't work on laptop or server, fnt or faa input; complains about illegal characters
#CLUSTALO_FORMAT = "stockholm"  #
CLUSTALO_FORMAT = "vienna"     # This one doesn't work on server, fnt input; complains that couldn't read format

# Locations and files
PHATE_PIPELINE_OUTPUT_DIR = os.environ["PHATE_PIPELINE_OUTPUT_DIR"] 

# templates 

# Class comparison organizes all genomic data and performs comparisons, ultimately yielding
#   homology groups for each gene/protein in a reference genome.
class build1hmm(object):
    
    def __init__(self):

        # parameters set by setParameters() via build class
        self.codeName                    = "hmmbuild" # hmmbuild (default)
        self.codeVersion                 = "2"        # 2 (default at time of code creation)
        self.alignmentCodeName           = "clustalo" # Clustal Omega (default)
        self.alignmentCodeVersion        = ""         # unknown 
        self.groupFile                   = ""         # PipelineOutput/GENOMICS/... fastaFile
        self.alignmentFile               = ""         # Output from alignment code
        self.hmmFile                     = ""         # results from hmmbuild
        self.outFile                     = ""         # captures screen output from hmmbuild
        self.sequenceType                = ""         # "nucleotide" or "protein"
        self.alignmentFormat             = ""         # set by build object

    def setParameters(self,kvargs):
        if isinstance(kvargs,dict):
            if "codeName" in list(kvargs.keys()):
                self.codeName = kvargs["codeName"]
            if "codeVersion" in list(kvargs.keys()):
                self.codeVersion = kvargs["codeVersion"]
            if "alignmentCodeName" in list(kvargs.keys()):
                self.alignmentCodeName = kvargs["alignmentCodeName"]
            if "alignmentCodeVersion" in list(kvargs.keys()):
                self.alignmentCodeVersion = kvargs["alignmentCodeVersion"]
            if "groupFile" in list(kvargs.keys()):
                self.groupFile = kvargs["groupFile"]
            if "alignmentFile" in list(kvargs.keys()):
                self.alignmentFile = kvargs["alignmentFile"]
            if "hmmFile" in list(kvargs.keys()):
                self.hmmFile = kvargs["hmmFile"]
            if "outFile" in list(kvargs.keys()):
                self.outFile = kvargs["outFile"]
            if "sequenceType" in list(kvargs.keys()):
                self.sequenceType = kvargs["sequenceType"]
            if "alignmentFormat" in list(kvargs.keys()):
                self.alignmentFormat = kvargs["alignmentFormat"]
        return

    def runHmmAnalysis(self):
        self.runAlignment()
        self.runHmmbuild()

    def runAlignment(self):  # Currently coded for running clustal omega
        # Construct command
        command = ''

        # If parameters have been set, construct alignment command
        if self.groupFile != '' and self.alignmentFile != '' and self.alignmentFormat != '':
            command = self.alignmentCodeName + ' -i ' + self.groupFile + ' -o ' + self.alignmentFile + ' --outfmt=' + self.alignmentFormat
        else:
            return ERROR_CODE_1

        # Run command
        try:
            if os.path.exists(self.alignmentFile):
                if PHATE_MESSAGES:
                    print ("hmm_build says, Clearing previous results...removing ",self.alignmentFile)
                rm_command = 'rm ' + self.alignmentFile
                os.system(rm_command)
            if PHATE_MESSAGES:
                print ("hmm_build says, Running clustal omega; command is ",command)
            os.system(command)
        except:
            print("hmm_build says, ERROR: alignment execution failed for alignment file ",self.alignmentFile)
            return ERROR_CODE_2

        return SUCCESS

    def runHmmbuild(self):  # Run method runAlignment before running this method 

        #Help out hmmbuild by determining sequence type
        if self.sequenceType.lower() == 'nucleotide' or self.sequenceType.lower() == 'dna' or self.sequenceType.lower() == 'nt':
            seqType = '--dna'
        elif self.sequenceType.lower() == 'rna':
            seqType = '--rna'
        elif self.sequenceType.lower() == 'protein' or self.sequenceType.lower() == 'peptide' or self.sequenceType.lower() == 'aa':
            seqType = '--amino'

        # Construct command
        command = ''
        if self.hmmFile != '' and self.alignmentFile != '' and self.outFile != '' and self.groupFile != '':
            command = self.codeName + ' ' + seqType + ' -o ' + self.outFile + ' ' + self.hmmFile + ' ' + self.alignmentFile
        else:
            return ERROR_CODE_1

        # Run command
        try:
            if os.path.exists(self.hmmFile):
                if PHATE_MESSAGES:
                    print ("multiPhATE_build says, Clearing out previous results...")
                rm_command = 'rm ' + self.hmmFile
                os.system(rm_command)
            if PHATE_MESSAGES:
                print ("hmm_build says, Running hmmbuild; command is ",command)
            os.system(command)     
        except:
            print ("hmm_build says, ERROR: hmmbuild execution failed for group ", self.groupFile)
            return ERROR_CODE_2

        return SUCCESS

    def printParameters(self):
        print("=========parameters========")
        print("codeName:",self.codeName,", version",self.codeVersion)
        print("alignmentCodeName:",self.alignmentCodeName,", version",self.alignmentCodeVersion)
        print("groupFile:",self.groupFile)
        print("alignmentFile:",self.alignmentFile)
        print("hmmFile:",self.hmmFile)
        print("outFile:",self.outFile)
        print("alignmentFormat:",self.alignmentFormat)
        return

    def printParameters2file(self,FILE_H):
        FILE_H.write("%s\n" % ("=========parameters========"))
        FILE_H.write("%s%s%s%s\n" % ("codeName:",self.codeName,", version",self.codeVersion))
        FILE_H.write("%s%s%s%s\n" % ("alignmentCodeName:",self.alignmentCodeName,", version",self.alignmentCodeVersion))
        FILE_H.write("%s%s\n" % ("groupFile:",self.groupFile))
        FILE_H.write("%s%s\n" % ("alignmentFile:",self.alignmentFile))
        FILE_H.write("%s%s\n" % ("hmmFile:",self.hmmFile))
        FILE_H.write("%s%s\n" % ("outFile:",self.outFile))
        FILE_H.write("%s%s\n" % ("alignmentFormat:",self.alignmentFormat))
        return

    def printResults(self):
        print("=========HMM Profiles========")
        print("nothing to show just yet")
        return

    def printResults2file(self,FILE_H):
        FILE_H.write("%s\n" % ("=========HMM Profiles========"))
        FILE_H.write("%s\n" % ("nothing to show just yet"))
        return

    def printAll(self):
        self.printParameters()
        self.printResults()
        return

    def printAll2file(self,FILE_H):
        self.printParameters2file(FILE_H)
        self.printResults2file(FILE_H)
        return

class build(object):

    def __init__(self):
        self.codeName                    = "hmmbuild" # hmmbuild by default
        self.codeVersion                 = "2"        # 2 at time of code creation
        self.alignmentCodeName           = "clustalo" # Clustal Omega by default
        self.alignmentCodeVersion        = ""         # 
        self.buildDirectory              = ""         # set by calling code (ie, genomics_driver.py)
        self.buildSubDirectory           = ""         # computed by self.setParameters() 
        self.fileList_fnt                = []         # list of nucleotide fasta group files for hmm build
        self.fileList_faa                = []         # list of peptide fasta group files for hmm build
        self.alignmentFormat             = CLUSTALO_FORMAT  # Default is defined at top of this file

    def setParameters(self,kvargs):
        if isinstance(kvargs,dict):
            if "codeName" in list(kvargs.keys()):
                self.codeName = kvargs["codeName"]
            if "codeVersion" in list(kvargs.keys()):
                self.codeVersion = kvargs["codeVersion"]
            if "alignmentCodeName" in list(kvargs.keys()):
                self.alignmentCodeName = kvargs["alignmentCodeName"]
            if "alignmentCodeVersion" in list(kvargs.keys()):
                self.alignmentCodeVersion = kvargs["alignmentCodeVersion"]
            if "buildDirectory" in list(kvargs.keys()):
                self.buildDirectory = kvargs["buildDirectory"]
                self.buildSubDirectory = os.path.join(self.buildDirectory,"HOMOLOGY_GROUPS")
            if "alignmentFormat" in list(kvargs.keys()):
                self.alignmentFormat = kvargs["alignmentFormat"]
        return

    def getListsOfGroups(self):

        if PHATE_PROGRESS:
            print("hmm_build says, Getting lists of groups")

        if os.path.isdir(self.buildSubDirectory):
            command = 'ls ' + self.buildSubDirectory
            proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True) 
            (rawresult, err) = proc.communicate()
            result = rawresult.decode('utf-8')
            fileList = str(result).split('\n')   # Python 3
            for filename in fileList:
                if filename.endswith('fnt'):
                    self.fileList_fnt.append(filename)
                elif filename.endswith('faa'):
                    self.fileList_faa.append(filename)
        else:
            print("hmm_build says, ERROR: Cannot open directory ",self.buildSubDirectory)
            return
        return

    def performBuild(self):

        # Query directory for list of group filenames
        self.getListsOfGroups()

        # Create HMM build object
        myHmm = build1hmm()

        if PHATE_PROGRESS:
            print("hmm_build says, Building HMMs for nucleotide files")

        # Determine file extension appropriate for alignmentFormat
        alignmentExtension = ''
        if self.alignmentFormat.lower() == 'fasta':
            alignmentExtension = '.afa'  #???
        elif self.alignmentFormat.lower() == 'clustal':
            alignmentExtension = '.clu'
        elif self.alignmentFormat.lower() == 'phylip':
            alignmentExtension = '.phy'
        elif self.alignmentFormat.lower() == 'selex':
            alignmentExtension = '.sel'
        elif self.alignmentFormat.lower() == 'stockholm':
            alignmentExtension = '.sto'
        elif self.alignmentFormat.lower() == 'vienna':
            alignmentExtension = '.vie'
        else:
            alignmentExtension = '.aln'

        # Create HMMs for nucleotide groups
        for group in self.fileList_fnt:
            # load parameters into standard build1hmm object
            groupFile = os.path.join(self.buildSubDirectory,group)
            alignmentFile = groupFile + alignmentExtension  
            hmmFile       = groupFile + '.hmm'
            outFile       = groupFile + '.out'
            parameterSet = {
                    "codeName"          : self.codeName,
                    "codeVersion"       : self.codeVersion,
                    "alignmentCodeName" : self.alignmentCodeName,
                    "groupFile"         : groupFile,
                    "alignmentFile"     : alignmentFile,
                    "hmmFile"           : hmmFile,
                    "outFile"           : outFile,
                    "sequenceType"      : "nucleotide",
                    "alignmentFormat"   : self.alignmentFormat,
                    }
            myHmm.setParameters(parameterSet)
            if PHATE_MESSAGES:
                print("hmm_build says, creating HMM for group ",group)
            myHmm.runHmmAnalysis()

        # Create HMMs for peptide groups

        if PHATE_PROGRESS:
            print("hmm_build says, Building HMMs for protein files")

        for group in self.fileList_faa:
            # Load parameters into standard build1hmm object
            groupFile = os.path.join(self.buildSubDirectory,group)
            alignmentFile = groupFile + alignmentExtension
            hmmFile       = groupFile + '.hmm'
            outFile       = groupFile + '.out'
            parameterSet = {
                    "codeName"          : self.codeName,
                    "codeVersion"       : self.codeVersion,
                    "alignmentCodeName" : self.alignmentCodeName,
                    "groupFile"         : groupFile,
                    "alignmentFile"     : alignmentFile,
                    "hmmFile"           : hmmFile,
                    "outFile"           : outFile,
                    "sequenceType"      : "protein",
                    "alignmentFormat"   : self.alignmentFormat,
                    }
            myHmm.setParameters(parameterSet)
            if PHATE_MESSAGES:
                print("hmm_build says, creating HMM for group ",group)
            myHmm.runHmmAnalysis()

        return

    def printParameters(self):
        print("Build Parameters:")
        print("codeName: ",self.codeName," codeVersion: ",self.codeVersion)
        print("alginmentCodeName: ",self.alignmentCodeName," alignmentCodeVersion: ",self.alignmentCodeVersion)
        print("buildDirectory: ",self.buildDirectory)
        return

    def printParameters2file(self,FILE_H):
        FILE_H.write("%s\n" % ("Build Parameters:"))
        FILE_H.write("%s%s%s%s\n" % ("codeName: ",self.codeName," codeVersion: ",self.codeVersion))
        FILE_H.write("%s%s%s%s\n" % ("alginmentCodeName: ",self.alignmentCodeName," alignmentCodeVersion: ",self.alignmentCodeVersion))
        FILE_H.write("%s%s\n" % ("buildDirectory: ",self.buildDirectory))
        return

    def printAll(self):
        self.printParameters()
        print("Files for hmm build:")
        print("Nucleotide files:")
        print(self.fileList_fnt)
        print("Protein files:")
        print(self.fileList_faa)
        return

    def printAll2file(self,FILE_H):
        self.printParameters2file(FILE_H)
        FILE_H.write("%s\n" % ("Files for hmm build:"))
        FILE_H.write("%s\n" % (self.fileList_fnt)) 
        FILE_H.write("%s\n" % (self.fileList_faa)) 
        return

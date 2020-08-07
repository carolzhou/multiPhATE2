#!usr/bin/env python


###################################################################
#
# compareGeneProfiles_main.py
#
# Programmer:  Carol L. Ecale Zhou
# Last update: 05 June 2020
#
# This program compares the gene calls from 2 genomes and identifies genes that 
# match, are similar, or are unique in each genome. This code is a re-write of
# compareGeneProfiles.py, which was written using beginner's knowledge of Python.
# This new program implements the same functionality,  but is written
# as object-oriented, modular code.
#
# Input parameters, format specifications, and default values:
# Files:
#    genome sequence files (in multi-fasta format; -g1 and -g2)
#    annotation files (in RAST .gff or .gff3 format; -a1 and -a2)
# Match/similarity detection (integer 0-100) with default value:
#    %identity gene-match cutoff (default 95)
#    %coverage gene-match cutoff (default 95)
#    %identity gene-similarity cutoff (default 60)
#    %coverage gene-similarity cutoff (default 75)
#    %identity domain-match cutoff (default 95)
#    %coverage domain-match cutoff (default 45)
#    %identity domain-similarity cutoff (default 60)
#    %coverage domain-similarity cutoff (default 45)
#    %identity paralog-match cutoff (default 95)
#    %coverage paralog-match cutoff (default 95)
#    %identity paralog-domain similarity cutoff (default 60)
#    %coverage paralog-domain similarity cutoff (default 45) 
#    %identity protein-match cutoff (default 95)
#    %coverage protein-match cutoff (default 95)
#    %identity protein-similarity cutoff (default 60)
#    %coverage protein-similarity cutoff (default 75)
#    %identity protein-domain-match cutoff (default 95)
#    %coverage protein-domain-match cutoff (default 45)
#    %identity protein-domain-similarity cutoff (default 60)
#    %coverage protein-domain-similarity cutoff (default 45)
#    %identity protein-paralog-match cutoff (default 95)
#    %coverage protein-paralog-match cutoff (default 95)
#    %identity protein-paralog-domain similarity cutoff (default 60)
#    %coverage protein-paralog-domain similarity cutoff (default 45) 
#
###################################################################
# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.PDF FOR DETAILS.

#SERVER = False  # True if this version of the code is running on the Django server
SERVER = True

BACTERIAL_CODE = 11

import sys, os, re, string, copy
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from time import strftime
from time import gmtime
import datetime
import subprocess
from subprocess import call
import cgp_fastaSequence as fastaSequence
import cgp_annotation as annotation
import cgp_genomeSequence as genomeSequence
import cgp_blastAnalysis as blastAnalysis

# Set messaging booleans

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

# More environmental variables

PHATE_EMBOSS_HOME = os.environ["PHATE_EMBOSS_HOME"]

# The following switches are for the benefit of code development, primarily.

DEBUG = False
#DEBUG = True

BLAST_ON = True    # controls whether blast is performed (used during development to skip lengthy blast for testing)
#BLAST_ON = False

PROTEIN = True     # controls whether protein sequences are blasted 
#PROTEIN = False

REPORT_STATS = True  # controls whether reports are written
#REPORT_STATS = False

REPORT_ON = True  # controls whether reports are written
#REPORT_ON = False

#Constants
MAX_TARGET_SEQS = 1
PARALOG_MAX     = 5 # max number of paralogs to be detected

try:
    maketrans = ''.maketrans
except:
    from string import maketrans

##### FILE names

##### FILES and DIRECTORIES
CODE_BASE_DIR       = os.environ["CGP_CODE_BASE_DIR"]
PIPELINE_OUTPUT_DIR = os.environ["PHATE_PIPELINE_OUTPUT_DIR"]
BASE_DIR = ""
EMAIL_DIR = BASE_DIR  # Move files to this dir for email results to user
OUT_DIR          = ""   # Default
GENOME_FILE1     = ""   # Input parameters... 
GENOME_FILE2     = ""   #  
ANNOTATION_FILE1 = ""   # 
ANNOTATION_FILE2 = ""   # 
OUT_FILE         = "compareGeneProfiles_main.out"
REPORT_FILE      = "compareGeneProfiles_main.report"
SUMMARY_FILE     = "compareGeneProfiles_main.summary"
PARALOG_FILE     = "compareGeneProfiles_main.paralogs"
LOG_FILE         = "compareGeneProfiles_main.log"
LOG_COPY         = "compareGeneProfiles_main.log.copy"
OUT              = ""   # File handles...
REPORT           = ""   #
SUMMARY          = ""   #
LOG              = ""   #
ERROR_LOG        = ""   #

##########################################################################
##### PATTERNS

p_fasta = re.compile('\.(fasta)|(fna)|(fas)|(fnt)|(fa)')
p_gff   = re.compile('\.gff*')

########################################################################
##### DECLARATIONS 

##### CONSTANTS

HELP_STRING = "This script inputs 2 genome fasta files and 2 RAST annotation files and returns a file (compareGeneProfiles_main.report containing a report specifying a comparison between the gene profiles of the 2 genomes, including matching genes, gene homologs and paralogs, and putative domain-fusion events. Optionally you may provide percent identity and coverage cutoffs for gene matching at the nucleotide level. By default these values for matching genes are 95% identity and 95% coverage. (All identity and coverage values for gene/domain/paralog match and similarity detection are configurable. For complete information on input parameters, type: compareGeneProfiles_main -input)"

USAGE_STRING = "Usage:  compareGeneProfiles_main.py -g1 <genome1.fasta> -a1 <genome1.gff> -g2 <genome2.fasta> -a2 <genome2.gff>\nYou may enter parameters interactively by typing: compareGeneProfiles_main.py -interactive" 

INPUT_STRING = "You may enter parameters interactively (via prompt) by typing: compareGeneProfiles_main.py -interactive\nYou will then be promted for the following required command-line parameters: \n   genome1.fasta\n   genome2.fasta\n   annotation1.gff\n   annotation2.gff\nOptionally, you may provide any/all of the following cutoff values as integers ranging from 0 to 100, with restrictions, by entering them at the command line prompt, or you may accept the defaults:\n   gene match sequence identity\n   gene match minumum coverage\n   gene similarity identity\n   gene similarity coverage\n   domain match identity\n   domain match coverage\n   domain similarity identity\n   domain similarity coverage\n   paralog match identity\n   paralog match coverage\n   paralog domain similarity identity\n   paralog domain similarity coverage\nType: compareGeneProfiles_main.py -interactive"   

ACCEPTABLE_ARG_COUNT = (2,3,9,11) # 2 if "help", "input", or "interactive"; 9 if 4 files provided with command-line labels (ie, '-g1')
                               # 3 if "config" mode; 11 if also including projectDirectory

##### BLAST default parameters

GENETIC_CODE    = 11 # bacterial, by default
BLAST_IDENTITY  = 60 
SCORE_EDGE      = 0.05
MAX_TARGET_SEQS = 1
OVERHANG        = 0.25
OUTPUT_FORMAT   = 7  # 7 is tabbed list output format for blast
EVALUE          = 1 

##### VARIABLES: command-line arguments

argCount         = 0 # number of command-line arguments provided
argument         = ""
genomeFile1      = "" # g1
genomeFile2      = "" # g2
rootFilename     = ""
annotationFile1  = "" # a1
annotationFile2  = "" # a2
projectDirectory = "" # d

#### VARIABLES: files

outFile      = ""     # Contains misc output
reportFile   = ""     # Contains binary genome comparison data 
summaryFile  = ""     # Contains summary data pertaining to binary comparison
paralogFile  = ""     # Contains list of paralogs detected
logFile      = ""     # script log

errorLog = os.path.join(CODE_BASE_DIR, "compareGeneProfiles_main.err") # error log (records prior to script main body)

##### VARIABLES: other

fileError       = False    # assume all is ok until proven otherwise
extensionString = ""       # for capturing file extension string

##### DATA STRUCTURES

parameters = { # data structure for input to class methods; defaults as indicated
"geneMatchIdentity"                      : 95, 
"geneMatchCoverage"                      : 95,
"geneSimilarityIdentity"                 : 60,
"geneSimilarityCoverage"                 : 75,
"domainMatchIdentity"                    : 95,
"domainMatchCoverage"                    : 45,
"domainSimilarityIdentity"               : 60,
"domainSimilarityCoverage"               : 45,
"paralogMatchIdentity"                   : 60,  
"paralogMatchCoverage"                   : 60,  
"paralogSimilarityIdentity"              : 60,
"paralogSimilarityCoverage"              : 75,
"paralogDomainMatchIdentity"             : 95,
"paralogDomainMatchCoverage"             : 45,
"paralogDomainSimilarityIdentity"        : 60,
"paralogDomainSimilarityCoverage"        : 45,
"proteinMatchIdentity"                   : 95, 
"proteinMatchCoverage"                   : 95,
"proteinSimilarityIdentity"              : 60,
"proteinSimilarityCoverage"              : 75,
"proteinDomainMatchIdentity"             : 95,
"proteinDomainMatchCoverage"             : 45,
"proteinDomainSimilarityIdentity"        : 60,
"proteinDomainSimilarityCoverage"        : 45,
"proteinParalogMatchIdentity"            : 60,  
"proteinParalogMatchCoverage"            : 60,  
"proteinParalogSimilarityIdentity"       : 60,
"proteinParalogSimilarityCoverage"       : 75,
"proteinParalogDomainMatchIdentity"      : 95,
"proteinParalogDomainMatchCoverage"      : 45,
"proteinParalogDomainSimilarityIdentity" : 60,
"proteinParalogDomainSimilarityCoverage" : 45,
"analysis"                               : None, # will fill in as "cross comparison" or "paralog"
"type"                                   : None # will fill in as "gene" or "protein"
}

files = { # User's input files and other generated files go here
"genomeFile1"          : "",  # user input
"genomeFile2"          : "",  # user input
"annotationFile1"      : "",  # user input
"annotationFile2"      : "",  # user input
"geneFile1"            : "",  # user input
"geneFile2"            : "",  # user input
"proteinFile1"         : "",  # user input
"proteinFile2"         : "",  # user input
"baseDir"              : "",  # user input
"geneFile1_root"       : "",  # file name sans all directory prefixing
"geneFile2_root"       : "",
"protFile1_root"       : "",
"protFile2_root"       : "",
"projectDirectory"     : "",  # if SERVER, input from wrapper code
}

##### FUNCTIONS ##############################################################

# For inserting "gene" or "prot" into genome filename, for new gene, protein fasta files
# This method is deprecated
def ConstructFilename(inFile,infix):  # infix is typically "gene" or "prot"
    newFile = ""
    stringList = re.findall('\.\w*', inFile)
    extensionString = stringList[-1] # Take last element in list
    rootFilename = inFile.rstrip(extensionString)
    newFile = rootFilename + "_" + infix + extensionString
    if PHATE_MESSAGES:
        print ("newFile name: ", newFile)
    return newFile

# For inserting "gene" or "prot" into genome filename, for new gene, protein fasta files
def ConstructNewFilename(inFile,infix,extension):
    newFile = ""
    (myPath,myFile) = os.path.split(inFile)
    (fileRoot,exten) = myFile.split('.')
    filename = fileRoot + '_' + infix + '.' + extension
    newFile = os.path.join(PIPELINE_OUTPUT_DIR, filename)  
    return newFile

# FUNCTION GetRootFile - Captures file name from directory-path/filename string
def GetRootFile(inFile):   # strip the inFile of all directory prefixes
    if PHATE_MESSAGES:
        print ("inFile is", inFile)
    rootFile = ""
    stringList = re.findall('\/[\w\d\.]*', inFile)
    if stringList:
        rootFile = stringList[-1]  # Take last element in list
    else:
        rootFile = inFile
    rootFile = rootFile.lstrip('/')
    return rootFile 

# FUNCTION GetConfig - parses a configuration string and assigns file names
def GetConfig(configString):
    stringList = []
    configString = configString.lstrip('\"')
    configString = configString.rstrip('\"')
    stringList = configString.split('#')
    BASE_DIR                 = stringList[0]
    files["genomeFile1"]     = stringList[1]
    files["genomeFile2"]     = stringList[2]
    files["annotationFile1"] = stringList[3]
    files["annotationFile2"] = stringList[4]
    (pathRoot1, fileName1)   = os.path.split(files["annotationFile1"])
    (pathRoot2, fileName2)   = os.path.split(files["annotationFile2"])
    (pathHead1,pathTail1)    = os.path.split(pathRoot1)
    (pathHead2,pathTail2)    = os.path.split(pathRoot2)
    joinName1 = pathTail1 + "_cgp_gene.fnt"
    joinName2 = pathTail2 + "_cgp_gene.fnt"
    files["geneFile1"]       = os.path.join(pathRoot1, joinName1)
    files["geneFile2"]       = os.path.join(pathRoot2, joinName2)
    joinName1 = pathTail1 + "_cgp_protein.faa"
    joinName2 = pathTail2 + "_cgp_protein.faa"
    files["proteinFile1"]    = os.path.join(pathRoot1, joinName1)
    files["proteinFile2"]    = os.path.join(pathRoot2, joinName2)

# FUNCTION GetArguments - gets all input parameters from user interactively 
def GetArguments(parameters,files):
    # get Files
    sys.stdout.write("Please enter the base directory where your files are located: ")
    files["baseDir"] = sys.stdin.readline().rstrip()
    #files["baseDir"] = files["baseDir"].rstrip('/')  # remove terminal slash if exsits
    #files["baseDir"] = files["baseDir"] + '/'        # and put it there for sure
    sys.stdout.write("Please enter the genome #1 filename:")
    files["genomeFile1"] = sys.stdin.readline().rstrip()
    sys.stdout.write("Genome #2 filename:")
    files["genomeFile2"] = sys.stdin.readline().rstrip()
    sys.stdout.write("Genome #1 annotation filename:")
    files["annotationFile1"] = sys.stdin.readline().rstrip()
    sys.stdout.write("Genome #2 annotation filename:")
    files["annotationFile2"] = sys.stdin.readline().rstrip()
    GET_PARAMS = False
    sys.stdout.write("Do you want to accept all defaults? enter y/n  ")
    decision = sys.stdin.readline().rstrip()
    if decision.lower() == "n" or decision.lower() == "no":
        GET_PARAMS = True 
    if GET_PARAMS:
        #GENE###########################################################
        sys.stdout.write("gene match sequence identity (default=")
        sys.stdout.write(str(parameters["geneMatchIdentity"]))
        sys.stdout.write("): ")
        new = sys.stdin.readline().rstrip()
        if new:
            if int(new) >= 0 and int(new) <= 100:
                parameters["geneMatchIdentity"] = new

        sys.stdout.write("gene match minimum coverage (default=")
        sys.stdout.write(str(parameters["geneMatchCoverage"]))
        sys.stdout.write("): ")
        new = sys.stdin.readline().rstrip()
        if new:
            if int(new) >= 0 and int(new) <= 100:
                parameters["geneMatchCoverage"] = new

        sys.stdout.write("gene similarity sequence identity (default=")
        sys.stdout.write(str(parameters["geneSimilarityIdentity"]))
        sys.stdout.write("): ")
        new = sys.stdin.readline().rstrip()
        if new:
            if int(new) >= 0 and int(new) <= 100:
                parameters["geneSimilarityIdentity"] = new

        sys.stdout.write("gene similarity coverage (default=")
        sys.stdout.write(str(parameters["geneSimilarityCoverage"]))
        sys.stdout.write("): ")
        new = sys.stdin.readline().rstrip()
        if new:
            if int(new) >= 0 and int(new) <= 100:
                parameters["geneSimilarityCoverage"] = new

        #PROTEIN###########################################################
        sys.stdout.write("protein match sequence identity (default=")
        sys.stdout.write(str(parameters["proteinMatchIdentity"]))
        sys.stdout.write("): ")
        new = sys.stdin.readline().rstrip()
        if new:
            if int(new) >= 0 and int(new) <= 100:
                parameters["proteinMatchIdentity"] = new

        sys.stdout.write("protein match minimum coverage (default=")
        sys.stdout.write(str(parameters["proteinMatchCoverage"]))
        sys.stdout.write("): ")
        new = sys.stdin.readline().rstrip()
        if new:
            if int(new) >= 0 and int(new) <= 100:
                parameters["proteinMatchCoverage"] = new

        sys.stdout.write("protein similarity sequence identity (default=")
        sys.stdout.write(str(parameters["proteinSimilarityIdentity"]))
        sys.stdout.write("): ")
        new = sys.stdin.readline().rstrip()
        if new:
            if int(new) >= 0 and int(new) <= 100:
                parameters["proteinSimilarityIdentity"] = new

        sys.stdout.write("protein similarity coverage (default=")
        sys.stdout.write(str(parameters["proteinSimilarityCoverage"]))
        sys.stdout.write("): ")
        new = sys.stdin.readline().rstrip()
        if new:
            if int(new) >= 0 and int(new) <= 100:
                parameters["proteinSimilarityCoverage"] = new

        return(0) 

# FUNCTION Translate2protein inputs gene sequence and metadata and translates gene sequence to protein
def Translate2protein(kvargs):
    translation = ""
    geneSequence = ""
    geneticCode = BACTERIAL_CODE
    tempGeneFile        = "./tempGeneSequence.txt"
    tempTranslationFile = "./tempTranslation.txt"
    if isinstance(kvargs,dict):
        if "geneSequence" in list(kvargs.keys()):
            geneSequence = kvargs["geneSequence"]
        if "geneticCode" in list(kvargs.keys()):
            geneticCode = kvargs["geneticCode"]
    try:
        GENE_H = open(tempGeneFile,"w")
        GENE_H.write("%s\n" % (geneSequence))
        GENE_H.close()
    except:
        print("cgp_compareGeneProfiles_main says, ERROR: Cannot open temporary gene file,",tempGeneFile)
        return "ERROR" 
    command = PHATE_EMBOSS_HOME + "transeq" + " -sequence " + tempGeneFile + " -outseq " + tempTranslationFile + " -table " + str(geneticCode) 
    result = os.system(command)
    PROT_H = open(tempTranslationFile,"r")
    fLines = PROT_H.read().splitlines()
    for i in range(1,(len(fLines))):
        translation += fLines[i]
    return translation

########################################################################################################
# BEGIN MAIN
########################################################################################################
# Get command-line arguments; open user input/output files
########################################################################################################

today = os.popen('date')
ERROR_LOG = open(errorLog,"a")
ERROR_LOG.write("%s%s\n" % ("Reading command-line input at ",today.read()))

argCount = len(sys.argv)
if PHATE_PROGRESS:
    print ("cgp_compareGeneProfiles_main says, Reading input arguments")
if PHATE_MESSAGES:
    print ("  Number of command-line arguments:", argCount)
if argCount in ACCEPTABLE_ARG_COUNT:
    match = re.search("help", sys.argv[1].lower())
    if match:
        print (HELP_STRING)
        print (USAGE_STRING)
        ERROR_LOG.close(); exit(0)

    match = re.search("input", sys.argv[1].lower())
    if match:
        print (INPUT_STRING)
        ERROR_LOG.close(); exit(0)

    match = re.search("interactive", sys.argv[1].lower())
    if match:
        GetArguments(parameters,files)

    match = re.search("config", sys.argv[1].lower())
    if match:
        configString = sys.argv[2]
        GetConfig(configString)

    if argCount == 9 or argCount == 11:
        if PHATE_PROGRESS:
            print ("cgp_compareGeneProfiles_main says, Reading input parameters")
        for i in range(argCount):
            if sys.argv[i] == "-g1":
                files["genomeFile1"] = sys.argv[i+1]
            if sys.argv[i] == "-g2":
                files["genomeFile2"] = sys.argv[i+1]
            if sys.argv[i] == "-a1":
                files["annotationFile1"] = sys.argv[i+1]
            if sys.argv[i] == "-a2":
                files["annotationFile2"] = sys.argv[i+1]
            if sys.argv[i] == "-d":
                files["projectDirectory"] = sys.argv[i+1]
                BASE_DIR = files["projectDirectory"]
        if PHATE_PROGRESS:
            print ("  genomeFile1 is", files["genomeFile1"])
            print ("  genomeFile2 is", files["genomeFile2"])

        # compute gene files
        (pathRoot1, fileName1)   = os.path.split(files["annotationFile1"])
        (pathRoot2, fileName2)   = os.path.split(files["annotationFile2"])
        (pathHead1,pathTail1)    = os.path.split(pathRoot1)
        (pathHead2,pathTail2)    = os.path.split(pathRoot2)
        joinName1 = pathTail1 + "_cgp_gene.fnt"
        joinName2 = pathTail2 + "_cgp_gene.fnt"
        files["geneFile1"]       = os.path.join(pathRoot1, joinName1)
        files["geneFile2"]       = os.path.join(pathRoot2, joinName2)
        joinName1 = pathTail1 + "_cgp_protein.fnt"
        joinName2 = pathTail2 + "_cgp_protein.fnt"
        files["proteinFile1"]    = os.path.join(pathRoot1, joinName1)
        files["proteinFile2"]    = os.path.join(pathRoot2, joinName2)

        if PHATE_MESSAGES:
            print ("  geneFile1 is",       files["geneFile1"])
            print ("  geneFile2 is",       files["geneFile2"])
            print ("  proteinFile1 is",    files["proteinFile1"])
            print ("  proteinFile2 is",    files["proteinFile2"])
            print ("  annotationFile1 is", files["annotationFile1"])
            print ("  annotationFile2 is", files["annotationFile2"])
            print ("  BASE_DIR is",        files["projectDirectory"])
            print ("  Done reading input parameters")
    else:
        print (USAGE_STRING)
        ERROR_LOG.close();  exit(0)
else:
    print (USAGE_STRING)
    ERROR_LOG.write("%s\n" % ("Incorrect number of command-line arguments provided"))
    ERROR_LOG.close(); exit(0)

# Check files
match = re.search(p_fasta, files["genomeFile1"])
if not match:
    fileError = True
match = re.search(p_fasta, files["genomeFile2"])
if not match:
    fileError = True
match = re.search(p_gff, files["annotationFile1"])
if not match:
    fileError = True
match = re.search(p_gff, files["annotationFile2"])
if not match:
    fileError = True
if fileError:
    print ("cgp_compareGeneProfiles says, ERROR: Check the formats of your input files:")  
    print ("  Genome file #1:         ", files["genomeFile1"])
    print ("  Annotation file #1:     ", files["annotationFile1"])
    print ("  Genome file #2:         ", files["genomeFile2"])
    print ("  Annotation file #2:     ", files["annotationFile2"])
    print ("  Project directory:      ", files["projectDirectory"])
    print (USAGE_STRING)
    ERROR_LOG.write("%s\n" % ("ERROR in input file parameter"))
    ERROR_LOG.close(); exit(0)
else:
    ERROR_LOG.close()

# Prepend filenames with projectDirectory

if SERVER:
    OUT_DIR      = os.path.join(files["projectDirectory"], OUT_DIR)
    reportFile   = os.path.join(OUT_DIR, "compareGeneProfiles_main.report")
    profilesFile = os.path.join(OUT_DIR, "compareGeneProfiles_main.out")
    paralogFile  = os.path.join(OUT_DIR, "compareGeneProfiles_main.paralogs")
    logFile      = os.path.join(OUT_DIR, LOG_FILE)

if PHATE_PROGRESS:
    print ("cgp_compareGeneProfiles_main says...")
    print ("  1: OUT_DIR is", OUT_DIR)
    print ("  2: SERVER is", SERVER)
    print ("  baseDir is", files["baseDir"])

# Prepend file names with baseDir; ""/no change if not provided:
files["genomeFile1"] = os.path.join(files["baseDir"], files["genomeFile1"])
files["genomeFile2"] = os.path.join(files["baseDir"], files["genomeFile2"])
files["annotationFile1"] = os.path.join(files["baseDir"], files["annotationFile1"])
files["annotationFile2"] = os.path.join(files["baseDir"], files["annotationFile2"])

##### Open files

try:
    GENOME_FILE1 = open(files["genomeFile1"],"r")
except IOError as e:
    fileError = True
    print (e)
try:
    GENOME_FILE2 = open(files["genomeFile2"],"r")
except IOError as e:
    fileError = True
    print (e)
try:
    ANNOTATION_FILE1 = open(files["annotationFile1"],"r")
except IOError as e:
    fileError = True
    print (e)
try:
    ANNOTATION_FILE2 = open(files["annotationFile2"],"r")
except IOError as e:
    fileError = True
    print (e)
if fileError:
    ERROR_LOG.write("%s\n" % ("ERROR with input filename"))
    ERROR_LOG.close(); exit(0)

dateTime = str(datetime.datetime.now().time()) # need time down to sub-seconds 
today    = os.popen('date')
if SERVER:
    OUT_DIR = os.path.join(files["projectDirectory"], "Results_" + dateTime + "/")  #
else:
    OUT_DIR = "./"

if PHATE_MESSAGES:
    print ("cgp_compareGeneProfiles_main says, Making directory", OUT_DIR)
command = "mkdir " + OUT_DIR
os.system(command)
outFile      = os.path.join(OUT_DIR, OUT_FILE)     # Contains misc output
reportFile   = os.path.join(OUT_DIR, REPORT_FILE)  # Contains binary genome comparison data
summaryFile  = os.path.join(OUT_DIR, SUMMARY_FILE) # Contains summary data pertaining to binary comparison
paralogFile  = os.path.join(OUT_DIR, PARALOG_FILE) # Contains list of paralogs detected
logFile      = os.path.join(OUT_DIR, LOG_FILE)

LOG          = open(logFile,"w")
LOG.write("%s%s\n" % ("Begin script log ",today.read()))
LOG.write("%s%s\n" % ("Result directory suffix is ",dateTime))

# Record to log and keep in touch with user 

LOG.write("%s\n" % ("Parameters are: "))
LOG.write("%s%s\n" % ("genome file #1: ",files["genomeFile1"]))
LOG.write("%s%s\n" % ("genome file #2: ",files["genomeFile2"]))
LOG.write("%s%s\n" % ("annotation file #1: ",files["annotationFile1"]))
LOG.write("%s%s\n" % ("annotation file #2: ",files["annotationFile2"]))
LOG.write("%s%s\n" % ("user directory is: ",files["projectDirectory"]))
LOG.write("%s%s\n" % ("outFile is: ",outFile))

if PHATE_MESSAGES:
    print ("cgp_compareGeneProfiles_main says, Parameters are:")
    print ("  Genome file #1:", files["genomeFile1"])
    print ("  Genome file #2:", files["genomeFile2"])
    print ("  Annotation file #1:", files["annotationFile1"])
    print ("  Annotation file #2:", files["annotationFile2"])
    print ("  User directory:", files["projectDirectory"])
    print ("  Out file:", outFile)
    keys = list(parameters); keys.sort()
    for key in keys:
        print (key, "is", parameters[key])

#########################################################################################################
# Construct genome data structures
#########################################################################################################

fLines = GENOME_FILE1.read().splitlines() # read lines into list, removing newlines
genome1 = genomeSequence.genome()
genome1.filename = files["genomeFile1"]
genome1.addContigs(fLines)
GENOME_FILE1.close()

fLines = GENOME_FILE2.read().splitlines()
genome2 = genomeSequence.genome()
genome2.filename = files["genomeFile2"]
genome2.addContigs(fLines)
GENOME_FILE2.close()

#########################################################################################################
# Extract gene sequences from genome files based on gene calls from annotation files;
# Translate gene sequences to protein.
# Construct gene and protein fastas and add to multi-fasta objects for each.
# Create gene and protein blast databases. 
# Save gene and protein fasta sequences to files.
#########################################################################################################

if PROTEIN: 
    if PHATE_PROGRESS:
        print ("cgp_compareGeneProfiles_main says...")
        print ("  Blast databases for genome 1 genes & proteins will be:")
        print ('  ',files["geneFile1"], "and")
        print ('  ',files["proteinFile1"])
        print ("  Blast databases for genome 2 genes & proteins will be:")
        print ('  ',files["geneFile2"], "and")
        print ('  ',files["proteinFile2"])
else:
    if PHATE_PROGRESS:
        print ("cgp_compareGeneProfiles_main says...")
        print ("  Blast databases for genome 1 genes will be:")
        print ('  ',files["geneFile1"])
        print ("  Blast databases for genome 2 genes will be:")
        print ('  ',files["geneFile2"])

################################################################

#*** Move this function up to where it belongs with other functions.
# Extract gene calls from GFF file (lines), then create and store new gene and protein objects.
def ExtractGeneCalls(genomeX,lines): 
    geneCount    = 0
    geneCountCDS = 0
    geneTemplate       = fastaSequence.fasta()  # generic gene/protein object
    proteinTemplate    = fastaSequence.fasta()  # generic gene/protein object
    annotationTemplate = annotation.annotationRecord()  # generic annotation object  
    gff = {  # for passing data to annotation class
        "contig"       : "unknown",
        "method"       : "PhATE",
        "type"         : "gene",
        "start"        : 0,
        "end"          : 0,
        "score"        : 0,
        "strand"       : 0,
        "readingFrame" : 0,
        "annotation"   : "",
        "source"       : "multiPhATE",  #***
        }  
    gene = {  # for passing data to fasta class
        "contig"         : "",
        "header"         : "",
        "name"           : "",
        "sequence"       : "",
        "type"           : "nt",
        "parentSequence" : "",
        "order"          : 0,
        "start"          : 0,
        "end"            : 0,
        "annotationList" : [],
        }
    protein = {  # for passing data to fasta class
        "contig"         : "",
        "header"         : "",
        "name"           : "",
        "sequence"       : "",
        "type"           : "nt",
        "parentSequence" : "",
        "order"          : 0,
        "start"          : 0,
        "end"            : 0,
        "annotationList" : [],
        }
    translationArgs = {  # for passing to translate2protein function
        "geneSequence"   : "",
        "geneticCode"    : BACTERIAL_CODE,
        }

    if PHATE_MESSAGES:
        print("cgp_compareGeneProfiles, extractGeneCalls() says, removing this line:",lines[0])
    lines.remove(lines[0]) # skip header line 

    # Extract gene fastas; create fasta object for gene; add to genome's geneSet
    # Then, create corresponding protein object, and add that to genome's proteinSet
    for line in lines: 
        match_comment = re.search('^#',line)
        if match_comment:
            continue 
        geneCount += 1  # Don't need this, but diagnostic
        fields = line.split('\t')
        geneCallType = fields[2]
        if geneCallType.lower() == "cds":  # ie, skip '*RNA', 'gene', and other entries (for now)
            geneCountCDS += 1

            # Extract data from input line
            gff["contig"]          = fields[0]
            gff["method"]          = fields[1]       # should be "PhATE"
            gff["type"]            = fields[2]       # "gene" or "CDS"
            gff["start"]           = int(fields[3])
            gff["end"]             = int(fields[4])
            gff["score"]           = fields[5]       # not used in this code 
            gff["strand"]          = fields[6]
            gff["readingFrame"]    = fields[7]
            gff["annotation"]      = fields[8]

            # Transfer data to gene dict
            gene["order"]          = geneCountCDS
            gene["start"]          = gff["start"]
            gene["end"]            = gff["end"]
            gene["name"]           = "cds" + str(geneCountCDS)
            gene["parentSequence"] = fields[0]  # the contig this gene is on (a shortHeader if from RAST)
            gene["contig"]         = fields[0]
            gene["sequence"] = genomeX.getSubsequence(gff["start"],gff["end"],gene["parentSequence"])
            if gff["strand"] == '-':  # Reverse complement sequence if on reverse strand
                #reverseComplement = gene["sequence"].translate(complements)[::-1]
                #gene["sequence"] = reverseComplement
                tempGene = Seq(gene["sequence"])
                gene["sequence"] = tempGene.reverse_complement()
            header = gene["name"] + "/" + gff["strand"] + "/" + str(gff["start"]) + "/" + str(gff["end"]) + "/"
            gene["header"] = header
            newGene = copy.deepcopy(geneTemplate) # dynamically allocate memory for next gene object
            newGene.enterGeneData(gene)
            newAnnotation = copy.deepcopy(annotationTemplate)
            newAnnotation.enterGFFdata(gff)
            newGene.addAnnotation(newAnnotation)

            # Transfer data to protein dict
            protein["contig"]               = gene["contig"]
            protein["order"]                = gene["order"]
            protein["start"]                = gene["start"]
            protein["end"]                  = gene["end"]
            protein["name"]                 = gene["name"]
            protein["parentSequence"]       = gene["parentSequence"]
            protein["header"]               = gene["header"]
            translationArgs["geneSequence"] = gene["sequence"]
            protein["sequence"]             = Translate2protein(translationArgs)
            newProtein = copy.deepcopy(proteinTemplate)
            newProtein.enterProteinData(protein)
            newAnnotation = copy.deepcopy(annotationTemplate)
            newAnnotation.enterGFFdata(gff)
            newProtein.addAnnotation(newAnnotation)

            # Add new gene and new protein to genome
            genomeX.addGene(newGene)
            genomeX.addProtein(newProtein)
            genomeX.cleanUpAfterEMBOSS()

#########################################################################################
# Extract gene sequences from genomes based on annotation file gene calls...
# ...construct gene fasta sequence objects and insert into multiFasta objects;
# ...translate gene sequences to protein, then create and insert protein multiFasta objects.
# Write gene and protein fastas to files.
#########################################################################################

##### Genome #1 genes and proteins

# Extract and process gene calls from GFF file (also creates protein fasta list object) 
fLines = ANNOTATION_FILE1.read().splitlines()
ExtractGeneCalls(genome1,fLines)
ANNOTATION_FILE1.close()

# Print gene sequences to file
printFastas2fileArgs = {  # for passing parameters to class genome/
    "mtype"      : "gene",      # "gene" (default), "protein", or "contig"
    "headerType" : "short",     # "full", "short" (default), "truncated", "compound" # compound header is header + contig
    "filename"   : files["geneFile1"],
    }
if PHATE_PROGRESS:
    print("cgp_compareGeneProfiles says, Printing genome1 gene sequences to file,",printFastas2fileArgs["filename"])
success = genome1.printFastas2file(printFastas2fileArgs)

# Print protein sequences to file
printFastas2fileArgs = {  # for passing parameters to class genome/
    "mtype"      : "protein",   
    "headerType" : "short",    
    "filename"   : files["proteinFile1"],
    }
if PHATE_PROGRESS:
    print("cgp_compareGeneProfiles says, Printing genome1 protein sequences to file,",printFastas2fileArgs["filename"])
success = genome1.printFastas2file(printFastas2fileArgs)

##### Genome #2 genes

# Extract and process gene calls from GFF file
fLines = ANNOTATION_FILE2.read().splitlines()
ExtractGeneCalls(genome2,fLines)
ANNOTATION_FILE2.close()

# Print gene sequences to file
printFastas2fileArgs = {  # for passing parameters to class genome/
    "mtype"      : "gene",   
    "headerType" : "short",  
    "filename"   : files["geneFile2"],
    }
if PHATE_PROGRESS:
    print("cgp_compareGeneProfiles says, Printing genome2 gene sequences to file,",printFastas2fileArgs["filename"])
success = genome2.printFastas2file(printFastas2fileArgs)

# Print protein sequences to file
printFastas2fileArgs = {  # for passing parameters to class genome/
    "mtype"      : "protein",   
    "headerType" : "short",    
    "filename"   : files["proteinFile2"],
    }
if PHATE_PROGRESS:
    print("cgp_compareGeneProfiles says, Printing genome2 protein sequences to file,",printFastas2fileArgs["filename"])
success = genome2.printFastas2file(printFastas2fileArgs)

#########################################################################################################
# Create blast databases for gene and protein multi-fasta files 
#########################################################################################################

blastDBargs = {  # for passing parameters to class genomeSequence/makeBlastDBs method
    "dbType"   : "nucl", # default
    "filename" : ""
    }

if PHATE_PROGRESS:
    if PROTEIN:
        print ("Creating blast databases for genome 1 and genome 2 gene and protein sets...")
    else:
        print ("Creating blast databases for genome 1 and genome 2 gene sets...")

myBlast = blastAnalysis.blast()                # Make blast object

blastDBargs["dbType"] = "nucl"         # Set up for nucleotide blast
blastDBargs["filename"] = files["geneFile1"]
myBlast.makeBlastDB(blastDBargs)
blastDBargs["filename"] = files["geneFile2"]
myBlast.makeBlastDB(blastDBargs)

if PROTEIN:
    blastDBargs["dbType"] = "prot"         # Set up for protein blast
    blastDBargs["filename"] = files["proteinFile1"]
    myBlast.makeBlastDB(blastDBargs)
    blastDBargs["filename"] = files["proteinFile2"]
    myBlast.makeBlastDB(blastDBargs)

if PHATE_PROGRESS:
    print ("cgp_compareGeneProfiles_main says, Blast database creation complete.")

#####################################################################################
# Perform blast of gene and protein sets between genomes and against self
#####################################################################################

# Prepare parameters for blast

blastArgs = {
    "query"         : "",
    "subject"       : "",
    "mtype"         : "nucl",
    "evalue"        : EVALUE, 
    "identity"      : BLAST_IDENTITY,
    "scoreEdge"     : SCORE_EDGE,
    "maxTargetSeqs" : MAX_TARGET_SEQS,
    "overhang"      : OVERHANG,
    "outputFormat"  : OUTPUT_FORMAT, 
    "outfile"       : ""
}

if PHATE_PROGRESS:
    print ("cgp_compareGeneProfiles_main says, Blasting gene sets...")
blastArgs["mtype"]   = "nucl"

# Need the raw filename (without directory information) to construct blast output filenames
files["geneFile1_root"] = GetRootFile(files["geneFile1"])
files["geneFile2_root"] = GetRootFile(files["geneFile2"])
files["protFile1_root"] = GetRootFile(files["proteinFile1"])
files["protFile2_root"] = GetRootFile(files["proteinFile2"])

### Genome 1 - Genome 2

if PHATE_PROGRESS:
    print ("Genome 1 genes against genome 2 genes...")
blastArgs["query"]   = files["geneFile1"]
blastArgs["subject"] = files["geneFile2"]
blastArgs["maxTargetSeqs"] = 1 
outfile = OUT_DIR + files["geneFile1_root"] + "_" + files["geneFile2_root"] + "_" + "blastn_" +\
    str(blastArgs["evalue"]) + "_" + str(blastArgs["identity"]) + ".out"
if PHATE_PROGRESS:
    print("cpg_compareGeneProfiles_main says, genome1-genome2 outfile is",outfile)
blastArgs["outfile"] = outfile
files["g1_g2_blastn"] = outfile
if BLAST_ON:
    result = myBlast.performBlast(blastArgs)
    if PHATE_MESSAGES:
        print ("cgp_compareGeneProfiles_main says, Result of blasting genes genome1-genome2:", result)

### Genome 2 - Genome 1

if PHATE_PROGRESS:
    print ("cgp_compareGeneProfiles_main says, Genome 2 genes against genome 1 genes...")
blastArgs["query"]   = files["geneFile2"]
blastArgs["subject"] = files["geneFile1"]
blastArgs["maxTargetSeqs"] = MAX_TARGET_SEQS 
outfile = OUT_DIR + files["geneFile2_root"] + "_" + files["geneFile1_root"] + "_" + "blastn_" +\
    str(blastArgs["evalue"]) + "_" + str(blastArgs["identity"]) + ".out"
if PHATE_MESSAGES:
    print("cpg_compareGeneProfiles_main says, genome2-genome1 outfile is",outfile)
blastArgs["outfile"] = outfile
files["g2_g1_blastn"] = outfile
if BLAST_ON:
    result = myBlast.performBlast(blastArgs)
    if PHATE_PROGRESS:
        print ("cgp_compareGeneProfiles_main says, Result of blasting genes genome2-genome1:", result)

### Genome 1 - Genome 1

if PHATE_PROGRESS:
    print ("cgp_compareGeneProfiles_main says, Genome 1 genes against self...")
blastArgs["query"]   = files["geneFile1"]
blastArgs["subject"] = files["geneFile1"]
blastArgs["maxTargetSeqs"] = PARALOG_MAX 
outfile = OUT_DIR + files["geneFile1_root"] + "_" + files["geneFile1_root"] + "_" + "blastn_" +\
    str(blastArgs["evalue"]) + "_" + str(blastArgs["identity"]) + ".out"
blastArgs["outfile"] = outfile
files["g1_g1_blastn"] = outfile
if BLAST_ON:
    result = myBlast.performBlast(blastArgs)
    if PHATE_PROGRESS:
        print ("cgp_compareGeneProfiles_main says, Result of blasting genes genome1-genome1:", result)

### Genome 2 - Genome 2 

if PHATE_PROGRESS:
    print ("cgp_compareGeneProfiles_main says, Genome 2 genes against self...")
blastArgs["query"]   = files["geneFile2"]
blastArgs["subject"] = files["geneFile2"]
blastArgs["maxTargetSeqs"] = PARALOG_MAX 
outfile = OUT_DIR + files["geneFile2_root"] + "_" + files["geneFile2_root"] + "_" + "blastn_" +\
    str(blastArgs["evalue"]) + "_" + str(blastArgs["identity"]) + ".out"
blastArgs["outfile"] = outfile
files["g2_g2_blastn"] = outfile
if BLAST_ON:
    result = myBlast.performBlast(blastArgs)
    if PHATE_PROGRESS:
        print ("cgp_compareGeneProfiles_main says, Result of blasting genes genome2-genome1:", result)

if PHATE_PROGRESS:
    print ("cgp_compareGeneProfiles_main says, Blasting of gene sets complete.")

### Proteome 1 - Proteome 2
if PROTEIN:
    if PHATE_PROGRESS:
        print ("cgp_compareGeneProfiles_main says, Blasting protein sets...")
    blastArgs["mtype"]   = "prot"

    if PHATE_PROGRESS:
        print ("cgp_compareGeneProfiles_main says, Genome 1 proteins against genome 2 proteins...")
    blastArgs["query"]   = files["proteinFile1"]
    blastArgs["subject"] = files["proteinFile2"]
    blastArgs["maxTargetSeqs"] = MAX_TARGET_SEQS 
    outfile = OUT_DIR + files["protFile1_root"] + "_" + files["protFile2_root"] + "_" + "blastn_" +\
        str(blastArgs["evalue"]) + "_" + str(blastArgs["identity"]) + ".out"
    blastArgs["outfile"] = outfile
    files["g1_g2_blastp"] = outfile
    if BLAST_ON:
        result = myBlast.performBlast(blastArgs)
        if PHATE_PROGRESS:
            print ("cgp_compareGeneProfiles_main says, Result of blasting proteins genome1-genome2:", result)

    ### Proteome 2 - Proteome 1

    if PHATE_PROGRESS:
        print ("cgp_compareGeneProfiles says, Genome 2 proteins against genome 1 proteins...")
    blastArgs["query"]   = files["proteinFile2"]
    blastArgs["subject"] = files["proteinFile1"]
    blastArgs["maxTargetSeqs"] = MAX_TARGET_SEQS 
    outfile = OUT_DIR + files["protFile2_root"] + "_" + files["protFile1_root"] + "_" + "blastn_" +\
        str(blastArgs["evalue"]) + "_" + str(blastArgs["identity"]) + ".out"
    blastArgs["outfile"] = outfile
    files["g2_g1_blastp"] = outfile
    if BLAST_ON:
        result = myBlast.performBlast(blastArgs)
        if PHATE_PROGRESS:
            print ("cgp_compareGeneProfiles says, Result of blasting proteins genome2-genome1:", result)

    ### Proteome 1 - Proteome 1
 
    if PHATE_PROGRESS:
        print ("cgp_compareGeneProfiles_main says, Genome 1 proteins against self...")
    blastArgs["query"]   = files["proteinFile1"]
    blastArgs["subject"] = files["proteinFile1"]
    blastArgs["maxTargetSeqs"] = PARALOG_MAX 
    outfile = OUT_DIR + files["protFile1_root"] + "_" + files["protFile1_root"] + "_" + "blastn_" +\
        str(blastArgs["evalue"]) + "_" + str(blastArgs["identity"]) + ".out"
    blastArgs["outfile"] = outfile
    files["g1_g1_blastp"] = outfile
    if BLAST_ON:
        result = myBlast.performBlast(blastArgs)
        if PHATE_PROGRESS:
            print ("cgp_compareGeneProfiles_main says, Result of blasting proteins genome1-genome1:", result)
 
    ### Proteome 2 - Proteome 2

    if PHATE_PROGRESS:
        print ("cgp_compareGeneProfiles_main says, Genome 2 proteins against self...")
    blastArgs["query"]   = files["proteinFile2"]
    blastArgs["subject"] = files["proteinFile2"]
    blastArgs["maxTargetSeqs"] = PARALOG_MAX 
    outfile = OUT_DIR + files["protFile2_root"] + "_" + files["protFile2_root"] + "_" + "blastn_" +\
        str(blastArgs["evalue"]) + "_" + str(blastArgs["identity"]) + ".out"
    blastArgs["outfile"] = outfile
    files["g2_g2_blastp"] = outfile
    if BLAST_ON:
        result = myBlast.performBlast(blastArgs)
        if PHATE_PROGRESS:
            print ("cgp_compareGeneProfiles_main says, Result of blasting proteins genome2-genome2:", result)

    if PHATE_PROGRESS:
        print ("cgp_compareGeneProfiles_main says, Blasting of protein sets complete.")

#########################################################################################################
# Parse blast output files; write hits to data structures 
#########################################################################################################

if PHATE_PROGRESS:
    print ("cgp_compareGeneProfiles_main says, Recording gene hits...",)

gene12hitList = myBlast.recordHits(files["g1_g2_blastn"])
gene21hitList = myBlast.recordHits(files["g2_g1_blastn"])
gene11hitList = myBlast.recordHits(files["g1_g1_blastn"])
gene22hitList = myBlast.recordHits(files["g2_g2_blastn"])

if PHATE_PROGRESS:
    print ("cgp_compareGeneProfiles_main says, Recording of gene hits complete.")

if PROTEIN:

    if PHATE_PROGRESS:
        print ("cgp_compareGeneProfiles_main says, Recording protein hits...",)

    prot12hitList = myBlast.recordHits(files["g1_g2_blastp"])
    prot21hitList = myBlast.recordHits(files["g2_g1_blastp"])
    prot11hitList = myBlast.recordHits(files["g1_g1_blastp"])
    prot22hitList = myBlast.recordHits(files["g2_g2_blastp"])

    if PHATE_PROGRESS:
        print ("cgp_compareGeneProfiles_main says, Recording of protein hits complete.")

#########################################################################################################
# Identify mutual best hits, singluar best hits, loners, homologs, and paralogs
##########################################################################################################

lonerArgs = {
    "seqList1"     : None,
    "seqList2"     : None,
    "comparedHits" : None  # homology object returned by blast.compareHits method
}
#**************************************************************

if PHATE_PROGRESS:
    print ("cgp_compareGeneProfiles_main says, Identifying gene mutual best hits and singular hits...",)
parameters["type"] = "gene"
geneComparison = myBlast.compareHits(gene12hitList,gene21hitList,parameters)
if PHATE_PROGRESS:
    print ("cgp_compareGeneProfiles_main says, Identificaton of mutual-best and singular hits compelte.")

if PHATE_PROGRESS:
    print ("cgp_compareGeneProfiles_main says, Identifying gene loners...",)
lonerArgs["seqList1"] = genome1.geneSet
lonerArgs["seqList2"] = genome2.geneSet
lonerArgs["comparedHits"] = geneComparison
geneComparison = myBlast.identifyLoners(lonerArgs)
if PHATE_PROGRESS:
    print ("cgp_compareGeneProfiles_main says, Identification of loners complete.")

if PHATE_PROGRESS:
    print ("cgp_compareGeneProfiles_main says, Identifying gene paralogs...",)
myBlast.identifyParalogs(genome1.geneSet,gene11hitList,parameters)
myBlast.identifyParalogs(genome2.geneSet,gene22hitList,parameters)
if PHATE_PROGRESS:
    print ("cgp_compareGeneProfiles_main says, Identification of paralogs complete.")

#**************************************************************

if PROTEIN:
    if PHATE_PROGRESS:
        print ("cgp_compareGeneProfiles_main says, Identifying protein mutual best hits and singular hits...",)
    parameters["type"] = "protein"
    proteinComparison = myBlast.compareHits(prot12hitList,prot21hitList,parameters)
    if PHATE_PROGRESS:
        print ("cgp_compareGeneProfiles_main says, Identificaton of mutual-best and singular hits complete.")

    if PHATE_PROGRESS:
        print ("cgp_compareGeneProfiles_main says, Identifying protein loners...",)
    lonerArgs["seqList1"] = genome1.proteinSet
    lonerArgs["seqList2"] = genome2.proteinSet
    lonerArgs["comparedHits"] = proteinComparison
    proteinComparison = myBlast.identifyLoners(lonerArgs)
    if PHATE_PROGRESS:
        print ("cgp_compareGeneProfiles_main says, Identification of loners complete.")

    if PHATE_PROGRESS:
        print ("cgp_compareGeneProfiles_main says, Identifying protein paralogs...",)
    myBlast.identifyParalogs(genome1.proteinSet,prot11hitList,parameters)
    myBlast.identifyParalogs(genome2.proteinSet,prot22hitList,parameters)
    if PHATE_PROGRESS:
        print ("cgp_compareGeneProfiles_main says, Identification of paralogs complete.")

#######################################################################################################
# Tally and report genome1/2 gene/protein comparison statistics  
#######################################################################################################
#*** Testing:  are summary statistics correct?

if PHATE_MESSAGES:
    print("cgp_compareGeneProfiles_main says, Tallying genome1/2 gene comparison statistics.")
geneComp_stats = geneComparison.reportStats()  # reminder: a homology object
genome1_gene_stats = genome1.geneSet.reportStats() # reminder: a multi-fasta object
genome2_gene_stats = genome2.geneSet.reportStats()

if PROTEIN:
    if PHATE_MESSAGES:
        print("cgp_compareGeneProfiles_main says, Tallying genome1/2 protein comparison statistics.")
    protComp_stats = proteinComparison.reportStats()
    genome1_prot_stats = genome1.proteinSet.reportStats()
    genome2_prot_stats = genome2.proteinSet.reportStats()

if PHATE_MESSAGES:
    print("cgp_compareGeneProfiles_main says, Reporting genome1/2 gene and protein statistics.")

SUMMARY = open(summaryFile,"w")
SUMMARY.write("\n%s\n" % ("************Next Report************"))
SUMMARY.write("%s%s\n" % ("Genome: ",files["genomeFile1"]))
SUMMARY.write("%s%s\n" % ("Genome: ",files["genomeFile2"]))

for stat in geneComp_stats:
    SUMMARY.write("%s\n" % (stat))

if PROTEIN:
    for stat in protComp_stats:
        SUMMARY.write("%s\n" % (stat))

for stat in genome1_gene_stats:
    SUMMARY.write("%s\n" % (stat))
for stat in genome2_gene_stats:
        SUMMARY.write("%s\n" % (stat))

if PROTEIN:
    for stat in genome1_prot_stats:
        SUMMARY.write("%s\n" % (stat))
    for stat in genome2_prot_stats:
        SUMMARY.write("%s\n" % (stat))

SUMMARY.close()

if REPORT_STATS:
    geneComparison.reportStats()
    genome1.geneSet.reportStats()
    genome2.geneSet.reportStats()
    if PROTEIN:
        proteinComparison.reportStats()
        genome1.proteinSet.reportStats()
        genome2.proteinSet.reportStats()

if PHATE_PROGRESS:
    print("cgp_compareGeneProfiles_main says, genome1/2 gene and protein statistics complete.")

#######################################################################################################
# Combine mutual and singular best hits and unique genes/proteins in order wrt genome1
#######################################################################################################

# Perform data line merging
if PHATE_PROGRESS:
    print("cgp_compareGeneProfiles_main says, Performing data line merging.")

geneMergeList = geneComparison.mergeAll(genome1.geneSet,genome2.geneSet)
if PROTEIN:
    proteinMergeList = proteinComparison.mergeAll(genome1.proteinSet,genome2.proteinSet)

# Write merged data lines for gene analysis to out file
if PHATE_PROGRESS:
    print("cgp_compareGeneProfiles_main says, Recording comparison results to file", outFile)
OUT = open(outFile,"w")
OUT.write("%s" % ("GENE COMPARISON RESULTS:\n"))
for hitItem in geneMergeList:
    OUT.write("%s%s" % (hitItem,"\n"))

# Write merged data lines for protein analysis to out file
if PROTEIN:
    OUT.write("%s" % ("PROTEIN COMPARISON RESULTS:\n"))  #*** Skipping protein for now
    for hitItem in proteinMergeList:
        OUT.write("%s%s" % (hitItem,"\n"))

if PHATE_PROGRESS:
    print("cgp_compareGeneProfiles_main says, Ordering of mutual, singular-best, and unique genes complete.")
OUT.close()

#######################################################################################################
# Create report
#######################################################################################################

# Create gene report

if PHATE_MESSAGES:
    print("cgp_compareGeneProfiles_main says, Creating report.")
REPORT = open(reportFile,"w")
REPORT.write("%s" % ("#GENE HITS\n"))
for line in geneMergeList:
    for item in line:
        REPORT.write("%s%s" % (item,"\t"))
    REPORT.write("%s" % ("\n"))

# Create protein report
if PROTEIN:
    REPORT.write("%s" % ("#PROTEIN HITS\n"))   #*** Skipping protein for now
    for line in proteinMergeList:
        for item in line:
            REPORT.write("%s%s" % (item,"\t"))
        REPORT.write("%s" % ("\n"))
REPORT.close()
if PHATE_PROGRESS:
    print("cgp_compareGeneProfiles_main says, Report(s) created.")

#######################################################################################################
# Create paralog reports
#######################################################################################################

if PHATE_MESSAGES:
    print("cgp_compareGeneProfiles_main says, Writing paralog file.")
PARALOG = open(paralogFile,"w")

# genome1
PARALOG.write("%s%s\n" % ("# PARALOGS for genome ",genome1.filename))
PARALOG.write("%s\n" % ("# Gene Paralogs"))
for fasta_obj in genome1.geneSet.fastaList:
    for paralog_obj in fasta_obj.paralogList:
        paralog_obj.printAll2file_tab(PARALOG)
PARALOG.write("%s\n" % ("# Protein Paralogs"))
for fasta_obj in genome1.proteinSet.fastaList:
    for paralog_obj in fasta_obj.paralogList:
        paralog_obj.printAll2file_tab(PARALOG)
# genome2
PARALOG.write("%s%s\n" % ("# PARALOGS for genome ",genome2.filename))
PARALOG.write("%s\n" % ("# Gene Paralogs"))
for fasta_obj in genome2.geneSet.fastaList:
    for paralog_obj in fasta_obj.paralogList:
        paralog_obj.printAll2file_tab(PARALOG)
PARALOG.write("%s\n" % ("# Protein Paralogs"))
for fasta_obj in genome2.proteinSet.fastaList:
    for paralog_obj in fasta_obj.paralogList:
        paralog_obj.printAll2file_tab(PARALOG)

PARALOG.close()
if PHATE_PROGRESS:
    print("cgp_compareGeneProfiles_main says, gene paralog file created.")

#########################################################################################################
# Clean Up
#########################################################################################################

LOG.close()

##### Copy result files out from Results_xxx directory into base user directory
#***  For pairwise comparison, that means: out, summary, report files
#***  For NxN comparisons, that means overall results files (not yet generated) 

EMAIL_DIR    = BASE_DIR 
logCopy      = os.path.join(EMAIL_DIR, LOG_COPY)
outCopy      = os.path.join(EMAIL_DIR, OUT_FILE)
reportCopy   = os.path.join(EMAIL_DIR, REPORT_FILE)
summaryCopy  = os.path.join(EMAIL_DIR, SUMMARY_FILE)
paralogCopy  = os.path.join(EMAIL_DIR, PARALOG_FILE)
call(["cp", logFile,      logCopy])
call(["cp", outFile,      outCopy])
call(["cp", reportFile,   reportCopy])
call(["cp", summaryFile,  summaryCopy])
call(["cp", paralogFile,  paralogCopy])

if PHATE_PROGRESS:
    print ("cpg_compareGeneProfiles_main says, CompareGeneProfiles completed.")

#######################################################################################################

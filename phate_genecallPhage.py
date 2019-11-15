################################################################
#
# phate_genecallPhage.py
#
# Programmers: Jeff Kimbrel, Carol Zhou
#
# Description: Single command to run PHANOTATE, Prodigal, Glimmer and GeneMarkS on a fasta file
#
################################################################

# This code was originally developed by Jeff Kimbrel and adapted by Carol Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.PDF FOR DETAILS.

import os
import sys
import re

DEBUG = False 
#DEBUG = True

# booleans control which gene finder(s) to call
GENEMARKS_CALLS = False   
PRODIGAL_CALLS  = False 
GLIMMER_CALLS   = False 
PHANOTATE_CALLS = False 

# patterns
p_genemarks = re.compile('[gG][eE][nN][eE][mM][aA][rR][kK]')
p_glimmer   = re.compile('[gG][lL][iI][mM][mM][eE][rR]')
p_prodigal  = re.compile('[pP][rR][oO][dD][iI][gG][aA][lL]')
p_phanotate = re.compile('[pP][hH][aA][nN][oO][tT][aA][tT][eE]')

########## HOUSEKEEPING ##########

#paths
prodigalPath  = os.environ["PRODIGAL_PATH"]
glimmerPath   = os.environ["GLIMMER_PATH"] 
geneMarkSPath = os.environ["GENEMARKS_PATH"]
phanotatePath = os.environ["PHANOTATE_PATH"]
cgcPath       = os.environ["CGC_PATH"]

# Verbosity
PHATE_MESSAGES = os.environ["PHATE_MESSAGES"]
PHATE_WARNINGS = os.environ["PHATE_WARNINGS"]
PHATE_PROGRESS = os.environ["PHATE_PROGRESS"]

# Data Structures

files = {
    'Raw Prodigal GFF'    : '',
    'Raw Glimmer Output'  : '',
    'Raw GeneMarkS GFF'   : '',
    }

# Input Parameters

if PHATE_MESSAGES == 'True':
    print("There are", len(sys.argv), "input parameters:", sys.argv)

if len(sys.argv) == 1:
    if PHATE_WARNINGS == 'True':
        print("Usage: /usr/local/bin/python3.4 annotatePhage.py fastaFile.fa outFolder", "Exiting now.")
    exit(0)

fastaFileName = sys.argv[1]
outputFolder = sys.argv[2] + "/"
cgcLog = outputFolder + "cgc.log"
cgcGff = outputFolder + "cgc.gff"

# booleans to control gene finding
if len(sys.argv) == 4:
    # get instructions for which gene finders to run; typically, running in the automated pipeline here
    # argument is a string containing the names of gene callers to be run
    genecallParams = sys.argv[3]  # a string listing gene callers to use 
    match_genemarks = re.search(p_genemarks,genecallParams)
    match_prodigal  = re.search(p_prodigal, genecallParams)
    match_glimmer   = re.search(p_glimmer,  genecallParams)
    match_phanotate = re.search(p_phanotate,genecallParams)
    if match_genemarks:
        GENEMARKS_CALLS = True 
    if match_prodigal:
        PRODIGAL_CALLS = True
    if match_glimmer:
        GLIMMER_CALLS = True 
    if match_phanotate:
        PHANOTATE_CALLS = True 

#logfile = open("./phate_genecallPhage.log","w")
logfilefullpath = outputFolder + "phate_genecallPhage.log"
logfile = open(logfilefullpath,"w")
logfile.write("%s%s\n" % ("Input parameters are:",sys.argv))
workingFolder = os.getcwd()
if PHATE_MESSAGES == 'True':
    print("workingFolder is", workingFolder)
if not os.path.exists(outputFolder):
    os.mkdir(outputFolder)
logfile.write("%s%s\n" % ("output folder is ", outputFolder))
logfile.write("%s%s\n" % ("working folder is ", workingFolder))
resultsFile = outputFolder + "results.txt"
results = open(resultsFile,"w")
files = {'Results File' : resultsFile}  #*** CHECK THIS
logfile.write("%s%s\n" % ("results file is ",resultsFile))
logfile.write("%s%s\n" % ("GENEMARKS_CALLS is ",GENEMARKS_CALLS))
logfile.write("%s%s\n" % ("PRODIGAL_CALLS is ",PRODIGAL_CALLS))
logfile.write("%s%s\n" % ("GLIMMER_CALLS is ",GLIMMER_CALLS))
logfile.write("%s%s\n" % ("PHANOTATE_CALLS is ",PHANOTATE_CALLS))
logfile.write("%s%s\n" % ("DEBUG is ",DEBUG))
logfile.write("%s%s\n" % ("CGC log file is ",cgcLog))

callCounts = {'prodigal' : 0, 'glimmer' : 0, 'genemarks' : 0, 'phanotate' : 0}

iterateGlimmer = False
#iterateGlimmer = True 
runCGC = False  # Turn on (below) if there are at least 2 gene callers being used 

########## CLASSES ##########
class geneCall:
    geneCallList = []

    def __init__(self, ID, method, contig, start, stop, strand, score):
        geneCall.geneCallList.append(self)
        self.ID = ID
        self.method = method
        self.contig = contig
        self.start = start
        self.stop = stop
        self.strand = strand
        self.score = score

    def printall(self):
        print(self.ID, self.method, self.contig, self.start, self.stop, self.strand, self.score)
        for gene in geneCallList:
            print(gene)

########## FUNCTIONS ##########

def systemCall(command):
    #print("\nSYSTEM CALL: "+command)
    if PHATE_MESSAGES == 'True':
        print("\nSYSTEM CALL: ", command)
    logfile.write("%s%s\n" % ("command is ",command))
    os.system(command)

def getProdigalId(x):
    colonSplit = x.split(";")
    equalSplit = colonSplit[0].split("=")
    return(equalSplit[1])

def getGeneMarkSId(x):
    colonSplit = x.split(",")
    equalSplit = colonSplit[0].split("=")
    return(equalSplit[1])

def processProdigal(line):
    if not line.startswith("#"):
        lineSplit = line.split("\t")

        contig = lineSplit[0]
        method = lineSplit[1]
        start = lineSplit[3]
        stop = lineSplit[4]
        score = lineSplit[5]
        strand = lineSplit[6]

        ID = getProdigalId(lineSplit[8])

        if PHATE_MESSAGES == 'True':
            print("    ID=",ID,"method=",method,"contig=",contig,"start/stop=",start,'/',stop,"strand=",strand,"score=",score)
        geneCall(ID, method, contig, start, stop, strand, score)
        callCounts['prodigal'] += 1

def processGlimmer(line,contig):
    contigSplit = contig.split(" ")
    contig = contigSplit[0]

    lineSplit = line.split()
    method = "glimmer3"
    ID = lineSplit[0]

    start = lineSplit[1]
    stop = lineSplit[2]
    strand = lineSplit[3][0]
    if strand == "-":
        start = lineSplit[2]
        stop = lineSplit[1]

    score = lineSplit[4]
    geneCall(ID, method, contig, start, stop, strand, score)

    callCounts['glimmer'] += 1

def processGeneMarkS(line):
    if not line.startswith("#"):
        lineSplit = line.split("\t")

        if len(lineSplit) > 1:

            contigSplit = lineSplit[0].split(" ")
            contig = contigSplit[0]

            method = lineSplit[1]
            start  = lineSplit[3]
            stop   = lineSplit[4]
            score  = lineSplit[5]
            strand = lineSplit[6]

            ID = getGeneMarkSId(lineSplit[8])

            geneCall(ID, method, contig, start, stop, strand, score)
            callCounts['genemarks'] += 1

def processPhanotate(line):
    if not line.startswith('#'):

        split = line.split("\t")
        geneCall("NA", "PHANOTATE", split[3], split[0], split[1], split[2], "NA")
        #geneCall("NA", "PHANOTATE", "NA", split[0], split[1], split[2], "NA")
        callCounts['phanotate'] += 1

def Convert_cgc2gff(cgcFile,gffFile):
    p_caller   = re.compile('^#\s([\w\d\.\-]+)\sgene\scalls')  # caller is names in first line of file
    p_dataLine = re.compile('^\d')   # data lines begin with a digit

    # Initialize gff fields and caller
    seqid = '.'; source = '.'; type = 'cds'
    start = '0'; end = '0'; strand = '.'
    phase = '.'; attributes = '.'; score = '.'
    caller = 'unknown'

    # Open input/output files
    CGC_H = open(cgcFile,"r")
    GFF_H = open(gffFile,"w")
    GFF_H.write("%s\n" % ("##gff-version 3"))

    # Walk through input file (cgc-formatted gene-calls), and write to gff output file
    cLines = CGC_H.read().splitlines()
    for cLine in cLines:
        match_caller   = re.search(p_caller,cLine)
        match_dataLine = re.search(p_dataLine,cLine) 
        if match_dataLine:
            (geneNo,strand,leftEnd,rightEnd,length,contig,protein) = cLine.split('\t')
            seqid  = contig
            source = caller
            if strand == '+':           # GFF lists start/end, rather than left/right
                start = str(leftEnd)
                end   = str(rightEnd)
            elif strand == '-':
                start = str(rightEnd)
                end   = str(leftEnd) 
            else:
                if PHATE_WARNINGS == 'True':
                    print("WARNING: phate_geneCall says, Unexpected strand:", strand)
            GFF_H.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (seqid,source,type,start,end,strand,phase,attributes))
        elif match_caller:
            caller = match_caller.group(1)

    # Clean up
    CGC_H.close()
    GFF_H.close()

########## PRODIGAL ##########

# system call to prodigal
## the '-p' option will currently give an error that there is not enough sequence info to train.
## This will be addressed in version 3 (https://github.com/hyattpd/Prodigal/issues/11)

if PRODIGAL_CALLS:
    if PHATE_PROGRESS == 'True':
        print("Genecall module says: running Prodigal.")
    if PHATE_MESSAGES == 'True':
        print("\n########## Prodigal ##########")
    logfile.write("%s\n" % ("Processing Prodigal"))

    command = prodigalPath + 'prodigal -q -i ' + fastaFileName + ' -o ' + outputFolder + 'prodigal.gff -f gff -p meta'
    systemCall(command)

    prodigalPeptideFile = outputFolder + "prodigal.proteins.faa"
    prodigalPotentialFile = outputFolder + "prodigal.genes.potential"

    command = prodigalPath + 'prodigal -i ' + fastaFileName + ' -o ' + outputFolder + 'prodigal.genes.sco -f sco -p meta -d ' + prodigalPeptideFile + ' -s ' + prodigalPotentialFile
    systemCall(command)
    files['Raw Prodigal GFF'] = outputFolder + 'prodigal.gff'

    prodigalFile = open(files['Raw Prodigal GFF'], 'r')
    lines = prodigalFile.read().splitlines()
    for line in lines:
        processProdigal(line)
    logfile.write("%s\n" % ("Prodigal processing complete."))

else:
    logfile.write("%s\n" % ("Not running Prodigal gene calling"))

########## GLIMMER ##########

## This runs the standard version of glimmer3.02b ("from scratch")

if DEBUG:
    logfile.write("%s\n" % ("Preparing to process Glimmer calls"))
    logfile.write("%s%s\n" % ("GLIMMER_CALLS is ",GLIMMER_CALLS))
if GLIMMER_CALLS:
    if PHATE_PROGRESS == 'True':
        print("Genecall module says: running Glimmer.")
    if PHATE_MESSAGES == 'True':
        print("\n########## Glimmer ##########")
    logfile.write("%s\n" % ("Processing Glimmer"))
    logfile.write("%s%s\n" % ("glimmerPath is ",glimmerPath))
    logfile.write("%s%s\n" % ("fastaFileName is ",fastaFileName))
    logfile.write("%s%s\n" % ("outputFolder is ",outputFolder))

    systemCall(glimmerPath + 'long-orfs -n -t 1.15 ' + fastaFileName + ' ' + outputFolder + 'glimmer.longorfs')
    systemCall(glimmerPath + 'extract -t '           + fastaFileName + ' ' + outputFolder + 'glimmer.longorfs > ' + outputFolder + 'glimmer.train')
    systemCall(glimmerPath + 'build-icm -r ' + outputFolder + 'glimmer.icm < ' + outputFolder + 'glimmer.train')
    systemCall(glimmerPath + 'glimmer3 -o50 -g110 -t30 ' + fastaFileName + ' ' + outputFolder + 'glimmer.icm ' + outputFolder + 'glimmer')
    systemCall('tail -n +2 ' + outputFolder + 'glimmer.predict > ' + outputFolder + 'glimmer.coords') 

    glimmerOutputHandle = "glimmer"
    
    if iterateGlimmer == True:
        systemCall('tail -n +2 ' + outputFolder + 'glimmer.predict > ' + outputFolder + 'long-orfs.coords')    
        systemCall(glimmerPath + '/scripts/upstream-coords.awk 25 0 ' + outputFolder + 'long-orfs.coords | ' + glimmerPath + '/bin/extract ' + fastaFileName + ' - > ' + outputFolder + 'glimmer.upstream')
        systemCall('/data/data1/softwares/ELPH/bin/Linux-i386/elph ' + outputFolder + 'glimmer.upstream LEN=6 | ' + glimmerPath + '/scripts/get-motif-counts.awk > ' + outputFolder + 'glimmer.motif')
        systemCall('set startuse = \'' + glimmerPath + '/bin/start-codon-distrib -3 ' + fastaFileName + ' ' + outputFolder + '/long-orfs.coords\'')
        systemCall(glimmerPath + '/bin/glimmer3 -o50 -g110 -t30 -b ' + outputFolder + 'glimmer.motif -P $startuse ' + fastaFileName + ' ' + outputFolder + 'glimmer.icm ' + outputFolder + 'glimmerIterative')
        glimmerOutputHandle = "glimmerIterative"


    files['Raw Glimmer Output'] = outputFolder + '' + glimmerOutputHandle + '.predict'
    logfile.write("%s%s\n" % ("Raw Glimmer Output is ", files['Raw Glimmer Output']))

    currentGlimmerContig = ""
    for line in open(outputFolder + '' + glimmerOutputHandle + '.predict', 'rt'):
        line = line.rstrip()
        if line.startswith(">"):
            currentGlimmerContig = line[1:]
        else:
            processGlimmer(line,currentGlimmerContig)
    logfile.write("%s\n" % ("Processing Glimmer complete."))
else:
    logfile.write("%s\n" % ("Not running Glimmer gene calling"))

########## GENEMARKS ##########

if DEBUG:
    logfile.write("%s\n" % ("Preparing to process Genemarks calls"))
    logfile.write("%s%s\n" % ("GENEMARKS_CALLS is ",GENEMARKS_CALLS))
if GENEMARKS_CALLS:
    if PHATE_PROGRESS == 'True':
        print("Genecall module says: running GeneMarkS.")
    if PHATE_MESSAGES == 'True':
        print("\n########## GeneMarkS ##########")
    logfile.write("%s\n" % ("Processing GeneMarkS"))

    systemCall(geneMarkSPath + 'gmhmmp -m ' + geneMarkSPath + '/heu_11.mod ' + fastaFileName + ' -o ' + outputFolder + 'geneMarkS.lst -r')
    systemCall(geneMarkSPath + 'gmhmmp -m ' + geneMarkSPath + '/heu_11.mod ' + fastaFileName + ' -o ' + outputFolder + 'geneMarkS.gff -r -f G')
    files['Raw GeneMarkS GFF'] = outputFolder + 'geneMarkS.gff'

    for line in open(outputFolder + 'geneMarkS.gff', 'rt'):
        line = line.rstrip()
        processGeneMarkS(line)
    logfile.write("%s\n" % ("Processing Genemarks complete."))
else:
    logfile.write("%s\n" % ("Not running GeneMarkS gene calling"))

########## PHANOTATE ##########

if DEBUG:
    logfile.write("%s\n" % ("Preparing to process PHANOTATE calls"))
    logfile.write("%s%s\n" % ("PHANOTATE_CALLS is ",PHANOTATE_CALLS))

if PHANOTATE_CALLS:
    if PHATE_PROGRESS == 'True':
        print("Genecall module says: running Phanotate.")
    if PHATE_MESSAGES == 'True':
        print("\n########## PHANOTATE ##########")
    logfile.write("%s\n" % ("Processing PHANOTATE"))

    os.chdir(phanotatePath)
    systemCall('python phanotate.py ' + fastaFileName + ' > ' + outputFolder + 'phanotateOutput.txt' )
    os.chdir(workingFolder)

    for line in open(outputFolder + 'phanotateOutput.txt', 'rt'):
        line = line.rstrip()
        processPhanotate(line)
    logfile.write("%s\n" % ("Processing PHANOTATE complete."))
else:
    logfile.write("%s\n" % ("Not running PHANOTATE gene calling"))

########## RESULTS ##########

print("\n########## RESULTS ##########")
logfile.write("%s\n" % ("Preparing results"))
results.write("ID\tSTART\tSTOP\tSTRAND\tSCORE\tMETHOD\tCONTIG\n")

if DEBUG:
    logfile.write("%s\n" % ("GENE_CALL_TESTING"))
    print("Printing all gene calls:")
    #geneCall.printall()

logfile.write("%s\n" % ("Writing genecall data to results file..."))
for gene in geneCall.geneCallList:
    results.write(gene.ID + "\t" + gene.start + "\t" + gene.stop + "\t" + gene.strand + "\t" + gene.score + "\t" + gene.method + "\t" + gene.contig + "\n")
results.close

logfile.write("%s\n" % ("Printing tallied genecall method call counts..."))
for method in callCounts:
    print((method + ": " + str(callCounts[method])))

print()

logfile.write("%s\n" % ("Printing file info..."))
for fileType in files:
    #print(fileType,files[fileType],sep=": ")
    print(fileType, ": ", files[fileType], end=' ')
print()

logfile.write("%s\n" % ("Parsing genecall files into CGC format..."))
callerCount = 0
if GENEMARKS_CALLS:
    callerCount += 1
    systemCall('python ' + cgcPath + '/CGC_parser.py Genemarks ' + outputFolder + 'geneMarkS.gff ' + outputFolder + 'genemark.cgc')
if PRODIGAL_CALLS:
    callerCount += 1
    systemCall('python ' + cgcPath + '/CGC_parser.py Prodigal ' + outputFolder + 'prodigal.genes.sco ' + outputFolder + 'prodigal.cgc')
if GLIMMER_CALLS:
    callerCount += 1
    #systemCall('python ' + cgcPath + '/CGC_parser.py Glimmer ' + outputFolder + 'glimmer.coords ' + outputFolder + 'glimmer.cgc')
    systemCall('python ' + cgcPath + '/CGC_parser.py Glimmer ' + outputFolder + 'glimmer.predict ' + outputFolder + 'glimmer.cgc')
if PHANOTATE_CALLS:
    callerCount += 1
    systemCall('python ' + cgcPath + '/CGC_parser.py PHANOTATE ' + outputFolder + 'phanotateOutput.txt ' + outputFolder + 'phanotate.cgc')

logfile.write("%s%s\n" % ("callerCount is ",callerCount))
if callerCount >= 2:
    runCGC = True

if DEBUG:
    print("DEBUG: next, running CGC_main.py")
if runCGC:
    #systemCall('python ' + cgcPath + 'CGC_main.py ' + outputFolder + '*.cgc > ' + outputFolder + 'CGC_results.txt')
    if DEBUG:
        print("DEBUG: calling CGC, cgcLog is", cgcLog)
    #systemCall('python ' + cgcPath + 'CGC_main.py log=' + cgcLog + ' '                         + outputFolder + '*.cgc > ' + outputFolder + 'CGC_results.txt')
    systemCall('python ' + cgcPath + 'CGC_main.py log=' + cgcLog + ' ' + 'gff=' + cgcGff + ' ' + outputFolder + '*.cgc > ' + outputFolder + 'CGC_results.txt')
    
else:
    logfile.write("%s\n" % ("Not running CGC code: too few gene callers to meaningfully compare"))

# Lastly, convert cgc-formatted gene-call files to gff and write to output
# For each gene caller, open and process its caller.cgc file, format to gff, then write to cgc_caller.gff
if GENEMARKS_CALLS:
    cgcFile = outputFolder + 'genemark'       + '.cgc'
    gffFile = outputFolder + 'phate_geneMarkS' + '.gff'
    Convert_cgc2gff(cgcFile,gffFile)
if PRODIGAL_CALLS:
    cgcFile = outputFolder + 'prodigal'        + '.cgc'
    gffFile = outputFolder + 'phate_prodigal'  + '.gff'
    Convert_cgc2gff(cgcFile,gffFile)

if GLIMMER_CALLS:
    cgcFile = outputFolder + 'glimmer'         + '.cgc'
    gffFile = outputFolder + 'phate_glimmer'   + '.gff'
    Convert_cgc2gff(cgcFile,gffFile)

if PHANOTATE_CALLS:
    cgcFile = outputFolder + 'phanotate'       + '.cgc'
    gffFile = outputFolder + 'phate_phanotate' + '.gff'
    Convert_cgc2gff(cgcFile,gffFile)

logfile.close()

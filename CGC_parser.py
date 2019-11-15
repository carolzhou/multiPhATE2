#!/usr/bin/env python

################################################################
#
# CGC_parser.py
#
# Programmer:  Carol Zhou
#
# Description:  This code inputs the name of a gene caller plus
#    the gene-caller's output file, and outputs a properly formatted
#    file for input to CGC_main.py.  The purpose of this code is to
#    normalize all gene caller outputs to a standard format, which is
#    recognized by CGC_main.py (a code that compares gene calls
#    among gene callers.) Note that regardless of how the gene caller
#    designates the start/end of a gene call, this code sets the 
#    "leftEnd" to the smaller position number, and "rightEnd" to the
#    larger position number.  The strand ('+' or '-') determines 
#    whether these numbers represent start or end positions in the
#    script output file. 
#
# Input Files:  For GeneMarkS, use the XXX.fasta.lst file; for
#    Glimmer2, use the XXX.g2.coord file; and for Prodigal, use the
#    XXX.genes.sco file, for Glimmer3, use the run3.coords file. 
#
################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.pdf FOR DETAILS.

import datetime
import sys
import os
import re
import string
import copy
from subprocess import call

##### Verbosity

CGC_WARNINGS  = 'True'
CGC_MESSAGES  = 'True'
CGC_PROGRESS  = 'True'

#os.environ["CGC_WARNINGS"] = CGC_WARNINGS
#os.environ["CGC_MESSAGES"] = CGC_MESSAGES
#os.environ["CGC_PROGRESS"] = CGC_PROGRESS


##### CONFIGURABLE

# GLIMMER3 bool controls which version of glimmer was run
# You get many more gene call matches when you adjust (by 3) the glimmer2 coordinates
# Glimmer3 changed this, so it should match better to other gene calls 
# You may adjust the Prodigal setting for the Prodigal out file you are using (sco vs. gff)
# RAST needs further testing; we are not running RAST for phage genomes

GLIMMER3     = True   # if False => glimmer2
PRODIGAL_sco = True   # using the XXX.genes.sco file
PRODIGAL_gff = False  # using the XXX.genes.gff file
RAST_GFF3    = True   # using RAST gff3 file; other RAST formats not yet supported
GFF3         = True   # using a properly formatted Genbank .faa file (format as GFF3)

##### FILES

CODE_BASE      = "CGC_parser"
CODE           = CODE_BASE + ".py"
logfile        = "./" + CODE_BASE + ".log"  # Generic log; this should be converted to a global pipeline log
runlogFilename = CODE_BASE + ".log"         # Log recording to user's results directory
runlogFile     = ""                         # Will be path/filename to user's results directory
outFilename    = CODE_BASE + ".out"  # Standard out file
tmpFilename    = CODE_BASE + ".tmp"  # temporary out file for un-sorted gene calls
tmpFile        = ""                  # will be path/filename to user's results directory
userOutdir     = ""   # In pipeline mode, extract user's output directory from input parameter specifying output file
infile         = ""  # user provided
USER_OUT_PROVIDED = False 

LOGFILE = open(logfile,"w")

##### PATTERNS

p_comment   = re.compile('^#')
p_empty     = re.compile('^$')  # blank line
p_prodigal  = re.compile('[Pp][Rr][Oo][Dd][Ii][Gg][Aa][Ll]')
p_genemark  = re.compile('[Gg][Ee][Nn][Ee][Mm][Aa][Rr][Kk]')
p_glimmer   = re.compile('[Gg][Ll][Ii][Mm]+[Ee][Rr]')
p_rast      = re.compile('[Rr][Aa][Ss][Tt]')
p_thea      = re.compile('[Tt][Hh][Ee][Aa]')
p_phate     = re.compile('[Pp][Hh][Aa][Tt][Ee]')  # PHANOTATE, actually
p_phanotate = re.compile('[Pp][Hh][Aa][Nn][Oo][Tt][Aa][Tt][Ee]')
p_gff3      = re.compile('[Gg][Ff][Ff]')  # take gff or gff3
p_genbank   = re.compile('[Gg][Ee][Nn][Bb][Aa][Nn][Kk]')

##### IDIOMS

CHATTY = True
#CHATTY = False

##### CONSTANTS

HELP_STRING = "Script " + CODE + " inputs the name of a gene caller plus the output file arising \nfrom that gene caller. Then, the script converts the data to a format that is acceptable as input to \nscript CGC_main.py, which compares gene calls among a set of gene caller outputs.\nType: python" + CODE + " usage|input for more information\n"

USAGE_STRING = "Usage:  python " + CODE + " <geneCaller_name> <geneCall_filename> (optional)<output_filename>\n"

INPUT_STRING = "You may enter the name of a gene caller (e.g., Prodigal, GeneMark, Glimmer, RAST, PHANOTATE), followed by the gene-call file that the program produced. For Prodigal, use the Name.genes.sco file. For GeneMarkS, use the Name.fasta.lst file. For Glimmer2, use the Name.g2.coord file, but for Glimmer3 use the run3.coords file. For RAST, use gff3 output. For PHANOTATE... TBD.\n"

ACCEPTABLE_ARG_COUNT = (2,3,4)  # 2 if 'help'|'usage'|'input', or 3 if gene-caller and gene-caller.out, 4 if optional output file

##### GET INPUT PARAMETERS

geneCaller    = ""
geneCallerOut = ""
userOutfile   = ""
userOutdir    = ""
RUNLOGOPEN    = False

argCount = len(sys.argv)
if argCount in ACCEPTABLE_ARG_COUNT:
    match = re.search("help", sys.argv[1].lower())
    if match:
        print(HELP_STRING)
        LOGFILE.close()
        exit(0)
    match = re.search("input", sys.argv[1].lower())
    if match:
        print(INPUT_STRING)
        LOGFILE.close()
        exit(0)
    match = re.search("usage", sys.argv[1].lower())
    if match:
        print(USAGE_STRING)
        LOGFILE.close()
        exit(0)

    # Capture name of gene caller and its output file
    if argCount == 3 or argCount == 4:
        geneCaller    = sys.argv[1].lower()  # case insensitive
        geneCallerOut = sys.argv[2]
    else:
        print(USAGE_STRING)
        LOGFILE.write("%s\n" % ("Incorrect number of command-line arguments provided"))
        LOGFILE.close()
        exit(0)
    if argCount == 4:
        USER_OUT_PROVIDED = True
        userOutfile   = sys.argv[3]
        userOutdir = os.path.dirname(userOutfile)
else:
    print(USAGE_STRING)
    LOGFILE.write("%s\n" % ("Incorrect number of command-line arguments provided"))
    LOGFILE.close()
    exit(0)

# Open file

fileError = False

if userOutdir != "":
    try:
        runlogfile = userOutdir + '/' + runlogFilename
        RUNLOGFILE = open(runlogfile,"w")
        RUNLOGOPEN = True
    except IOError as e:
        fileError = True
        print(e)

if userOutdir != "":
    try:
        tmpFile = userOutdir + '/' + tmpFilename
        TMPFILE = open(tmpFile,"w")
        TMPOPEN = True
    except IOError as e:
        fileError = True
        print(e)

try:
    INFILE = open(geneCallerOut,"r")
except IOError as e:
    fileError = True
    print(e)

try:
    if userOutdir == "":
        outfile = './' + outFilename
    else:
        outfile = userOutdir + '/' + outFilename
    OUTFILE = open(outfile,"w")
except IOError as e:
    fileError = True
    print(e)

if fileError:
    LOGFILE.write("%s%s%s\n" % ("ERROR: problem with input file:",geneCallerOut,e))
    LOGFILE.close(); exit(0)
    if RUNLOGOPEN:
        RUNLOGFILE.write("%s%s%s\n" % ("ERROR: problem with input file:",geneCallerOut,e))
        RUNLOGFILE.close(); exit(0)

if USER_OUT_PROVIDED:
    try:
        USER_OUT = open(userOutfile,"w")
    except IOError as e:
        fileError = True
        print(e)

##### FUNCTIONS

def ProcessGenemark(fLines,OUT):
    geneNo = 0; contig = ''; strand = ''; leftEnd = ''; rightEnd = ''; length = 0; count = 0; cclass = ''; protein = 'unknown' 
    p_dataLine   = re.compile('\s+(\d+)\s+([+-])\s+([\d\>\<]+)\s+(\d+)\s+(\d+)\s+\d+')
    p_contigLine = re.compile('FASTA\sdefinition\sline:\s(.*)\s+length=\d+')
    for line in fLines:
        match_comment    = re.search(p_comment,line)
        match_contigLine = re.search(p_contigLine,line)
        match_dataLine   = re.search(p_dataLine,line)
        if match_comment:
            continue 
        elif match_contigLine:
            contig = match_contigLine.group(1)
        elif match_dataLine:
            geneNo   = int(match_dataLine.group(1))
            strand   =     match_dataLine.group(2)     
            leftEnd  =     match_dataLine.group(3)   # Note: left/right end could have '<' or '>' symbol
            rightEnd =     match_dataLine.group(4)
            length   = int(match_dataLine.group(5))
            #cclass   =     match_dataLine.group(6)  # skip for now
            # Note: left/right end could have '>' or '<' symbol, if gene call spanned across contigs
            # I am removing the symbol, so if start=1 this may imply a gene that spans 
            # For a gene spanning 2 contigs or wrapping around, genemark estimates length based on
            # raw start/stop without symbol
            if re.search('^>',leftEnd) or re.search('^<',leftEnd):
                leftEnd = leftEnd[1:] 
            if re.search('^>',rightEnd) or re.search('^<',rightEnd):
                rightEnd = rightEnd[1:] 
            if strand != '+' and strand != '-':
                if RUNLOGOPEN:
                    RUNLOGFILE.write("%s%s\n" % ("ERROR: unknown strand designator, ",strand))
                OUT.write("%s\n" ("ERROR encountered: unknown strand designator\n"))
                if USER_OUT_PROVIDED:
                    USER_OUT.write("%s\n" ("ERROR encountered: unknown strand designator\n"))
                if CGC_WARNINGS == 'True':
                    print("ERROR in CGC_parser module: unexpected strand designator,", strand)
                return
            if contig == '':
                contig = 'unknown'  # Contig name may be absent in input file
            #TMPFILE.write("%s\t%s\t%s\t%s\t%s\t%s\t%si\t%s\n" % (geneNo,strand,leftEnd,rightEnd,length,contig,protein))
            OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%si\t%s\n" % (geneNo,strand,leftEnd,rightEnd,length,contig,protein))
            if USER_OUT_PROVIDED:
                USER_OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (geneNo,strand,leftEnd,rightEnd,length,contig,protein))
    return
            
def ProcessGlimmer(fLines,OUT):
    # NOTE:  Glimmer2 appears to truncate the initial met (sometimes?): test with other data sets (e.g., + strand calls)
    geneNo = 0; contig = ''; strand = ''; leftEnd = ''; rightEnd = ''; length = 0; count = 0; protein = 'unknown'
    #p_contigLine = re.compile('^>(.*)\s+length=\d+\s+numreads=\d+') # Matches only for glimmer3
    p_contigLine = re.compile('^>(.*)')

    if GLIMMER3:
        #p_dataLine   = re.compile('orf(\d+)\s+(\d+)\s+(\d+)\s+([+-])\d\s+([\d\.]+)')
        p_dataLine   = re.compile('orf(\d+)\s+(\d+)\s+(\d+)\s+([+-])\d\s+([\d\.]+)')
    else:
        p_dataLine = re.compile('\s+(\d+)\s+(\d+)\s+(\d+)\s+\[([+-])\d\sL=\s*(\d+)\sr=.*\]') 

    for line in fLines:
        match_comment    = re.search(p_comment,line)
        match_contigLine = re.search(p_contigLine,line)
        match_dataLine   = re.search(p_dataLine,line)
        if match_comment:
            continue 
        elif match_contigLine:
            contig = match_contigLine.group(1)
        elif match_dataLine:
            fields   = line.split('\t')
            if len(fields) >= 6:
                contig = fields[5]
            geneNo   = int(match_dataLine.group(1))
            left     = int(match_dataLine.group(2))
            right    = int(match_dataLine.group(3))
            strand   =     match_dataLine.group(4)
            if strand == '+':
                leftEnd  = left 
                if GLIMMER3:
                    rightEnd = right
                else:
                    rightEnd = right + 3
            elif strand == '-':  
                if GLIMMER3:
                    leftEnd = right
                else:
                    leftEnd  = right - 3    
                rightEnd = left  
            else:
                if RUNLOGOPEN:
                    RUNLOGFILE.write("%s%s\n" % ("ERROR: unknown strand designator, ",strand))
                OUT.write("%s\n" ("ERROR encountered: unknown strand designator\n"))
                if USER_OUT_PROVIDED:
                    USER_OUT.write("%s\n" ("ERROR encountered: unknown strand designator\n"))
                if CGC_WARNINGS == 'True':
                    print("ERROR in CGC_parser module: unexpected strand designator,", strand)
                return

            length = rightEnd - leftEnd + 1
            if contig == '':    # contig name may be left out of input file
                contig = 'unknown'
            count += 1; geneNo = count  # re-assign gene number over that assigned by Glimmer
            OUT.write("%s\t%s\t%s\t%s\t%s\t%si\t%s\n" % (geneNo,strand,str(leftEnd),str(rightEnd),str(length),contig,protein))

            if USER_OUT_PROVIDED:
                USER_OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (geneNo,strand,str(leftEnd),str(rightEnd),str(length),contig,protein))
    return

def ProcessRAST(fLines,OUT):
    geneNo = 0; contig = ""; strand = ''; leftEnd = ''; rightEnd = ''; length = 0; count = 0; 
    if RAST_GFF3:
        p_dataLine = re.compile('(.*)\t.*\tCDS\t(\d+)\t(\d+)\t\.\t([+-])\t\d\t.*')
        for line in fLines:
            match_comment = re.search(p_comment,line)
            match_dataLine = re.search(p_dataLine,line) 
            if match_comment:
                continue
            elif match_dataLine:
                count += 1
                geneNo = count 
                contig    =     match_dataLine.group(1)
                leftEnd   = int(match_dataLine.group(2))
                rightEnd  = int(match_dataLine.group(3))
                strand    =     match_dataLine.group(4) 
                length = rightEnd - leftEnd + 1
                OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (geneNo,strand,str(leftEnd),str(rightEnd),str(length),contig,protein))
                if USER_OUT_PROVIDED:
                    USER_OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (geneNo,strand,str(leftEnd),str(rightEnd),str(length),contig,protein))
    else:
        pass # Not using other RAST format (for now)
    return
 
def ProcessGFF3(fLines,OUT):
    geneNo = 0; contig = ''; source = ''; feature = ''; strand = ''; leftEnd = ''; rightEnd = ''; length = 0; count = 0; protein = 'unknown'
    fields = []
    if GFF3:
        #p_dataLine = re.compile('(.*)\t.*\tCDS\t(\d+)\t(\d+)\t\.\t([+-])\t\d\t(.*)')
        for line in fLines:
            match_comment = re.search(p_comment,line)
            match_empty   = re.search(p_empty,line)
            #match_dataLine = re.search(p_dataLine,line) 
            if match_comment or match_empty:
                continue
            else:
                fields = line.split('\t')
                if fields:
                    count += 1
                    geneNo = count
                    if fields[0] == '':
                        contig = 'unknown'
                    else:
                        contig = fields[0]
                    source    =     fields[1]  #*** not using this currently
                    feature   =     fields[2]
                    leftEnd   = int(fields[3])
                    rightEnd  = int(fields[4])
                    strand    =     fields[6]
                    if fields[8] != '':
                        protein = fields[8]
                    if feature == 'CDS':
                        length = rightEnd - leftEnd + 1
                        OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (geneNo,strand,str(leftEnd),str(rightEnd),str(length),contig,protein))
                        if USER_OUT_PROVIDED:
                            USER_OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (geneNo,strand,str(leftEnd),str(rightEnd),str(length),contig,protein))
    else:
        pass  # do nothing if GFF3 not enabled
    return
 
def ProcessProdigal(fLines,OUT):
    geneNo = 0; contig = "unknown"; strand = ''; leftEnd = ''; rightEnd = ''; length = 0; count = 0; protein = 'unknown' 

    if PRODIGAL_sco:  # using XXX.genes.sco file
        p_header   = re.compile('seqhdr=\"(.*)\"')
        p_dataLine = re.compile('^>(\d+)_(\d+)_(\d+)_([+-])')
        for line in fLines:
            match_comment  = re.search(p_comment, line)
            match_header   = re.search(p_header,  line)
            match_dataLine = re.search(p_dataLine,line)
            if match_header:
                contig = match_header.group(1)
            elif match_comment:
                continue 
            elif match_dataLine:
                geneNo   = match_dataLine.group(1)
                leftEnd  = int(match_dataLine.group(2))
                rightEnd = int(match_dataLine.group(3))
                strand   =     match_dataLine.group(4)
                length   = int(rightEnd) - int(leftEnd) + 1 
                count += 1; geneNo = count  # Re-number gene number assigned by Prodigal
                print("Writing next prodigal gene", geneNo, strand, leftEnd, rightEnd, length, contig, protein)
                OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (geneNo,strand,leftEnd,rightEnd,length,contig,protein))
                if USER_OUT_PROVIDED:
                    USER_OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (geneNo,strand,leftEnd,rightEnd,length,contig,protein))

    else: # using XXX.genes file
        p_header   = re.compile('seqhdr=\"(.*)\"')
        p_dataLine = re.compile('(.*)\t.*\tCDS\t(\d+)\t(\d+)\t([\d\.]+)\t([+-])\t.\t(.*)') 
        for line in fLines:
            match_header   = re.search(p_header,  line)
            match_comment = re.search(p_comment,line)
            match_dataLine = re.search(p_dataLine,line)
            if match_header:
                contig = match_header.group(1)
            elif match_comment:
                continue 
            elif match_dataLine:
                count += 1
                geneNo   = count 
                contig   =     match_dataLine.group(1)
                leftEnd  = int(match_dataLine.group(2))
                rightEnd = int(match_dataLine.group(3))
                strand   =     match_dataLine.group(5)
                length   = int(rightEnd) - int(leftEnd) + 1
            OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (geneNo,strand,leftEnd,rightEnd,length,contig,protein))
            if USER_OUT_PROVIDED:
                USER_OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (geneNo,strand,leftEnd,rightEnd,length,contig,protein))
    return

def ProcessPHANOTATE(fLines,OUT):      # SDSU code
    p_comment  = re.compile('^#')
    p_dataLine = re.compile('^[\d><]')  # Look for initical character = digit, '>', or '<'
    geneNo = 0; contig = 'unknown'; strand = '?'; leftEnd = 0; rightEnd = 0; length = 0; count = 0; protein = 'unknown'; temp = ''
    for fLine in fLines:
        fields = []
        match_comment  = re.search(p_comment,fLine)
        match_dataLine = re.search(p_dataLine,fLine)
        if match_comment:
            continue
        if match_dataLine:
            count += 1
            fields   = fLine.split('\t')
            geneNo   = count
            left     = fields[0]
            right    = fields[1]
            strand   = fields[2]
            contig   = fields[3]
            matchLeft = re.search('\d+',left)   # take only digits; ie, omit '>' or '<'
            if matchLeft:
                leftEnd = matchLeft.group(0)
            matchRight = re.search('\d+',right) 
            if matchRight:
                rightEnd = matchRight.group(0)
            if int(rightEnd) < int(leftEnd):
                temp = leftEnd
                leftEnd = rightEnd
                rightEnd = temp 
            length   = int(rightEnd) - int(leftEnd) + 1
            OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (geneNo,strand,leftEnd,rightEnd,str(length),contig,protein))
            if USER_OUT_PROVIDED:
                USER_OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (geneNo,strand,leftEnd,rightEnd,str(length),contig,protein))
    return

##### BEGIN MAIN 

if CGC_PROGRESS == 'True':
    print("CGC_parser: Parsing gene call outputs.")

# First, determine which gene caller was used

match_glimmer     = re.search(p_glimmer,geneCaller)
match_genemark    = re.search(p_genemark,geneCaller)
match_prodigal    = re.search(p_prodigal,geneCaller)
match_rast        = re.search(p_rast,geneCaller)
match_thea        = re.search(p_thea,geneCaller)
match_phanotate   = re.search(p_phanotate,geneCaller)
match_gff3        = re.search(p_gff3,geneCaller)
match_genbank     = re.search(p_genbank,geneCaller)
match_phate       = re.search(p_phate,geneCaller)

OUTFILE.write("%s%s%s%s%s\n" % ('# ',geneCaller, " gene calls",", taken from file ",geneCallerOut))
OUTFILE.write("%s\n" % ("Gene No.\tStrand\tLeftEnd\tRightEnd\tLength\tContig\tProtein"))

if USER_OUT_PROVIDED:
    USER_OUT.write("%s%s%s%s%s\n" % ('# ',geneCaller, " gene calls",", taken from file ",geneCallerOut))
    USER_OUT.write("%s\n" % ("Gene No.\tStrand\tLeftEnd\tRightEnd\tLength\tContig\tProtein"))

fileLines = INFILE.read().splitlines()

if match_genemark:
    #ProcessGenemark(fileLines,OUTFILE)  #*** use GFF3 for now; need to troubleshoot
    ProcessGFF3(fileLines,OUTFILE)
elif match_glimmer:
    ProcessGlimmer(fileLines,OUTFILE)
elif match_prodigal:
    ProcessProdigal(fileLines,OUTFILE)
elif match_rast:
    ProcessRAST(fileLines,OUTFILE)
elif match_phanotate:
    ProcessPHANOTATE(fileLines,OUTFILE)
elif match_phate:
    ProcessPHANOTATE(fileLines,OUTFILE)
elif match_phanotate:
    ProcessPHANOTATE(fileLines,OUTFILE)
elif match_gff3:
    ProcessGFF3(fileLines,OUTFILE)
elif match_genbank:
    ProcessGFF3(fileLines,OUTFILE)  # Need to format genbank's protein fasta file as GFF3
else:
    if RUNLOGOPEN:
        RUNLOGFILE.write("%s%s\n" % ("ERROR: Cannot process unknown gene caller output file:",geneCaller))

if USER_OUT_PROVIDED:
    USER_OUT.write("%s\n" % ("# END"))
OUTFILE.write("%s\n" % ("# END"))

##### CLEAN UP

INFILE.close()
OUTFILE.close()
if USER_OUT_PROVIDED:
    USER_OUT.close()
if RUNLOGOPEN:
    RUNLOGFILE.write("%s\n" % ("Processing complete"))
    RUNLOGFILE.close()
LOGFILE.write("%s%s\n" % ("Processing complete at ", datetime.datetime.now()))
LOGFILE.close()

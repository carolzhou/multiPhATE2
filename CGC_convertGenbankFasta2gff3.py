#!/usr/bin/env python

################################################################
#
# CGC_convertGenbankFasta2gff3.py  # 
#
# Programmer: Carol Zhou
#
# Description:  Accepts a genbank protein fasta file and outputs a GFF3-formatted list of genes 
#    Use this script when you need to input genbank's proteins into CGC to compare with PHANOTATE's (or other).
#    This script will extract information from genbank's header, and construct a GFF3-formatted list of genes.
#    You need to input:
#         genbank protein fasta file
#         contig name
#
################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.pdf FOR DETAILS.

import sys
import os
import re
import string
import copy
from subprocess import call

##### Verbosity

CGC_WARNINGS = os.environ["CGC_WARNINGS"]
CGC_MESSAGES = os.environ["CGC_MESSAGES"]
CGC_PROGRESS = os.environ["CGC_PROGRESS"]

##### FILES

CODE_BASE = "./CGC_convertGenbankFasta2gff3"
CODE_FILE = CODE_BASE + ".py"
LOG_FILE  = CODE_BASE + ".log"
OUT_FILE  = CODE_BASE + ".out"

infile = ""
OUT = open(OUT_FILE,"w") 
LOG = open(LOG_FILE,"w")

##### PATTERNS

p_comment  = re.compile('^#')
p_order    = re.compile('Order')

##### PRINT CONTROL 

#DEBUG = True    # Print even more!
DEBUG = False

##### CONSTANTS

HELP_STRING = "This code inputs a genbank protein fasta file and a contig name, and outputs a GFF3-formatted file for input to CGC.  Type: python " + CODE_FILE + " usage|input for more information\n"

USAGE_STRING = "Usage:  python " + CODE_FILE + " <genbankProteinFasta> <outfilename>n\n"

INPUT_STRING = "Input for " + CODE_FILE + " comprises a genbank protein fasta path/filename followed by a contig name\nExample:  python " + CODE_FILE + " /mydir/NC_000001.faa myContig\n"

##### GET INPUT PARAMETERS

fileSet = []
contig = 'unknown'
argCount = len(sys.argv)
if argCount == 2 or argCount == 3:
    match = re.search("help", sys.argv[1].lower())
    if match:
        print(HELP_STRING)
        LOG.close(); exit(0)
    match = re.search("input", sys.argv[1].lower())
    if match:
        print(INPUT_STRING)
        LOG.close(); exit(0)
    match = re.search("usage", sys.argv[1].lower())
    if match:
        print(USAGE_STRING)
        LOG.close(); exit(0)
    else:
        infile = sys.argv[1]  # skip 0th element = name of code
        contig = sys.argv[2]
else:
    LOG.write("%s\n" % ("Incorrect number of command-line arguments provided"))
    print(USAGE_STRING)
    LOG.close(); exit(0)

# Open input file

IN = open(infile,"r")

# Construct patterns

p_header         = re.compile('^>')
p_contig         = re.compile('\|(NC\_\d\d\d\d\d\d\.\d)\_')
p_locusTag       = re.compile('locus_tag=([\w\d_\-]+)\]')
p_protein        = re.compile('protein=([\w\d\_\-\+\'\.\:\/\(\)\s]+)\]')
p_location_plus  = re.compile('location=(\d+)\.\.(\d+)\]')
p_location_minus = re.compile('location=complement*\((\d+)\.\.(\d+)\)\]')
p_location_join  = re.compile('join\((\d+)\.\.(\d+),(\d+)\.\.(\d+)\)\]')

##### BEGIN MAIN 

count = 0

# For each user-provided gene call file, create a call set and add to list of call sets

if CGC_MESSAGES == 'True':
    print("CGC convertGenbank module says: Input parameters are: genbank protein fasta file is", infile, "and contig name is", contig)

iLines = IN.read().splitlines()
OUT.write("%s\n" % ("#gff-version3"))
OUT.write("%s%s\n" % ("#",contig))
for iLine in iLines:
    contig   = 'unknown'
    locusTag = 'unknown'
    protein  = 'unknown'
    left     = ''
    right    = ''
    strand   = '?'
    outLine  = ''
    match_header = re.search(p_header,iLine)
    if match_header:
        count += 1
        match_contig         = re.search(p_contig,iLine)
        match_locusTag       = re.search(p_locusTag,iLine)
        match_protein        = re.search(p_protein,iLine)
        match_location_plus  = re.search(p_location_plus,iLine)
        match_location_minus = re.search(p_location_minus,iLine)
        match_location_join  = re.search(p_location_join,iLine)
        if match_contig:
            contig = match_contig.group(1)
        if match_locusTag:
            locusTag = match_locusTag.group(1)
        if match_protein:
            protein = match_protein.group(1)
        if match_location_plus:
            left  = match_location_plus.group(1)
            right = match_location_plus.group(2)
            strand = '+'
        if match_location_minus:
            left  = match_location_minus.group(1)
            right = match_location_minus.group(2)
            strand = '-'
        if match_location_join:                   #*** Note: I have not seen an ex. of join/complement, so strand might be wrong in that case
            left  = match_location_join.group(1)
            right = match_location_join.group(4)
            strand = '+'
            protein = protein + ' JOINED'
        outLine = contig + '\tgenbank\tCDS\t' + left + '\t' + right + '\t.\t' + strand + '\t.\t' + locusTag + ' ; ' + protein
        OUT.write("%s\n" % (outLine))

##### CLEAN UP

IN.close()
OUT.close()
LOG.close()

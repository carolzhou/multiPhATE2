#!/usr/bin/env python

#####################################################################
#
# phate_getLast3percent.py
#
# Description: This code captures the remaining 3% of pVOG sequences
#    from the NR database, which could not be retrieved in the original
#    capture.  Recall that a subset (query on anything that looked like
#    phage) of NR was used locally to pull sequences. Pulling either
#    from a downloaded complete NR or via web query would have been
#    very slow. However, the query subset was missing about 3% of the
#    pVOG sequences. Hense, this code to retrieve the remaining 3%
#    using directy query from a local instance of the complete NR.
#    Should be fast enough to do just 3% locally.
#
# Usage:  python phate_getLast3percent.py  (IMPORTANT: READ BELOW)
#    NOTE:  First run phate_createPvogFastaFile.py. This will search the
#       NCBI phage subset database and tag the fasta headers with the
#       pVOG information. Then, take the output from that code, 
#       pVOGs_missing.lst, and use that as input to this (current) code.
#       The output will be a fasta file containing the missing fasta
#       sequences. Then, you need to add those sequences to the NR
#       phage fasta subset file, and RE-RUN phate_createPvogFastaFile.py
#       (or run separately and add when they are done).
#       Then, the output from that code (i.e., pVOGs.faa) should 
#       contain 100% of the pVOG sequences, pulled from NR.
#
# Programmers Notes:
#    NOTE:  As currently written, this code finds only the 1st ncbiID in
#       the nr database, then 'skips' through the rest. Need to fix that
#       (ie, figure out how the SeqIO is actually working.) Nevertheless,
#       querying the 1st ncbiID required 1/2 hour, which means that the
#       set of 5000 will require 3 months on a single mpath processor.
#       !!! I need to find a better way.
#    1) Filenames are coded as constants
#    2) Only need to run this code one time (until update)
#    3) This code only pulls the missing fasta sequences. The user must
#       than add these sequences to the NR subset and re-run the original
#       code, phate_createPvogFastaFile.py (see NOTE above).
#
# Programmer:  CEZhou
#
####################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.PDF FOR DETAILS.

import sys
import os
import re
import string
import copy
import time
import datetime
from subprocess import call 
from Bio import SeqIO

# Control
CHATTY = True
DEBUG  = False

CODE_BASE = "phate_getLast3percent"

# Paths to files; when running on mpath

PVOG_BASE_DIR = "/home/zhou4/PhATE/Databases/pVOGs/"

# Files

NR_FASTA_DB           = "/data/data1/sandbox/BLAST/FASTA/nr"
PVOG_MISSING_PROTEINS = PVOG_BASE_DIR + "pVOGs_missing.lst"
PVOG_MISSING_FASTAS   = PVOG_BASE_DIR + "pVOGs_missing_fastas.faa"
LOG                   = PVOG_BASE_DIR + CODE_BASE + ".log"

# Dat# Environment variables, which are global to any instance of this code's execution
MPATH = True
SANDBOX = False

# Import PhATE Modules


# Open Files

PVOG_MISSING_PROTEINS_H = open(PVOG_MISSING_PROTEINS, 'r') 
PVOG_MISSING_FASTAS_H   = open(PVOG_MISSING_FASTAS, 'w')
LOG_H                   = open(LOG, 'w')
LOG_H.write("%s%s\n" % ("Processing began at ",datetime.datetime.now()))
print("LOG file is", LOG)

##### Find fasta files for pVOGs in "missing" list

# First, collect legitimate pVOG identifiers from "missing" list.
# Note, some members of this list are "predicted ORFn". This is some
# bogus data that the pVOGs folk included in their data sets.

# Read in lines from the "missing" file
vLines = PVOG_MISSING_PROTEINS_H.read().splitlines()

# Patterns
p_bogus = re.compile('predicted')
p_ncbi  = re.compile('[\w\d\_\.]')

# Collect only those that have an NCBI identifier
fields = []
pVOGwantedList = []
for vLine in vLines:
    fields = vLine.split(' ')
    match_bogus = re.search(p_bogus,fields[5])
    if match_bogus:
        if DEBUG:
            LOG_H.write("%s%s\n" % ("The following line contains bogus identifier: ",vLine))
    else:
        match_ncbi = re.search(p_ncbi,fields[5])
        if match_ncbi:
            if DEBUG:
                LOG_H.write("%s%s\n" % ("Found an NCBI identifier: ",fields[5]))
            pVOGwantedList.append(fields[5])
        else:
            if DEBUG:
                LOG_H.write("%s%s\n" % ("This string is unaccounted for: ",vLine))

LOG_H.write("%s%s\n" % ("pVOG wanted list has this many members: ",len(pVOGwantedList)))

# For each ncbi identifier, pull the fasta from NR and write to out file
count = 0
if CHATTY:
    print(pVOGwantedList)

fastaSequences = SeqIO.parse(open(NR_FASTA_DB),'fasta')
for ncbiID in pVOGwantedList:
    count += 1
    print(count, "Processing ncibID", ncbiID)
    for seq in fastaSequences:
        if seq.id == ncbiID:
            print("Writing fasta sequence for", ncbiID)
            SeqIO.write([seq], PVOG_MISSING_FASTAS_H, "fasta") 

print(pVOGwantedList)

# Clean up

PVOG_MISSING_PROTEINS_H.close()
PVOG_MISSING_FASTAS_H.close()

print("Processing complete!")
LOG_H.write("%s%s\n" % ("Processing completed at ",datetime.datetime.now()))
LOG_H.close()

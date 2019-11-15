#!/usr/bin/env python

#####################################################################
#
# Program:  dbPrep_createPvogFastaFile.py
#
# Description: This code reads input files, AllFamilyProteinList.tab and NCBI_Phage.faa,
#    and creates a file containing pVOG-tagged protein fastas.
#
# Usage:  python createPvogFastaFile.py
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

# Control
CHATTY = True
#CHATTY = False 
#DEBUG = True
DEBUG = False

CODE_BASE = "createPvogFastaFile"

# Paths to files; when running on mpath

#PVOG_BASE_DIR         = "/home/zhou4/PhATE/Databases/pVOGs/"
PVOG_BASE_DIR         = "/Users/carolzhou/DEV/PhATE/Databases/pVOGs/"
NCBI_BASE_DIR         = "/Users/carolzhou/DEV/PhATE/Databases/NCBI/Virus_Protein/"

# Files

PROTEIN_LIST          = PVOG_BASE_DIR + "AllFamilyProteinList.tab"
#NCBI_PHAGE_FAA        = PVOG_BASE_DIR + "NCBI_phage.faa"
#NCBI_PHAGE_FAA        = NCBI_BASE_DIR + "viral.nonredundant_protein.1.protein.faa"
NCBI_PHAGE_FAA        = NCBI_BASE_DIR + "viral.protein.faa"
PVOG_TAGGED_PROTEINS  = PVOG_BASE_DIR + "pVOGs.faa"
PVOG_MISSING_PROTEINS = PVOG_BASE_DIR + "pVOGs_missing.lst"
LOG                   = PVOG_BASE_DIR + CODE_BASE + ".log"
PVOG_LIB_DIR          = PVOG_BASE_DIR + "Allvogtables/"

# Dat# Environment variables, which are global to any instance of this code's execution
#MPATH = True
#SANDBOX = False

# Import PhATE Modules

import dbPrep_annotation
import dbPrep_fastaSequence
import dbPrep_pVOG

# Open Files

PROTEIN_LIST_H          = open(PROTEIN_LIST, 'r')
NCBI_PHAGE_FAA_H        = open(NCBI_PHAGE_FAA, 'r')
PVOG_TAGGED_PROTEINS_H  = open(PVOG_TAGGED_PROTEINS, 'w')
PVOG_MISSING_PROTEINS_H = open(PVOG_MISSING_PROTEINS, 'w') 
LOG_H                   = open(LOG, 'w')
LOG_H.write("%s%s\n" % ("Processing began at ",datetime.datetime.now()))

# Capture fastas from ncbi phage fasta file

if CHATTY:
    print("Reading ncbi phage fasta file")
ncbiPhageSeqs = dbPrep_fastaSequence.multiFasta()
fastaLines = NCBI_PHAGE_FAA_H.read().splitlines()
ncbiPhageSeqs.addFastas(fastaLines,'aa')
if CHATTY:
    print("done!")

# Read in the pVOG identifiers and their associated accession numbers

if CHATTY:
    print("Reading pVOGs from pVOG file")
pVOGdb = dbPrep_pVOG.pVOGs()
pVOGlines = PROTEIN_LIST_H.read().splitlines()
#pVOGdb.addPvogs_old(pVOGlines,LOG_H)
pVOGdb.addPvogs(PVOG_LIB_DIR,LOG_H)
if CHATTY:
    print("pVOGs have been recorded")
LOG_H.write("%s%s%s\n" % ("There are ",len(pVOGdb.pVOGlist)," pVOGs"))
accessionCount = pVOGdb.getAccessionCount()
LOG_H.write("%s%s\n" % ("The total number of accessions is ",accessionCount))

# Visit each pVOG and find the fasta sequence that corresponds to each member accession
# Modify the header of each identified fasta to reflect its membership in the pVOG cluster

if CHATTY:
    print("Searching sequences for each pVOG-associated accession")

# Create a fasta object (to be replicated as needed)
nextFasta = dbPrep_fastaSequence.fasta()

# For each pVOG, get its associated peptide accessions, then
# Find that fasta in the ncbi database subset, and 
# Tag the fasta header with the pVOG information
accnCount = 0
for pVOG in pVOGdb.pVOGlist:
    foundCount = 0; missingCount = 0; missingList = []
    for accession in pVOG.accessionList:                            # For each accession (approx 200k of them, members of pVOG groups)
        accnCount += 1
        if CHATTY:
            print("Processing pVOG", pVOG.pVOGid, "and accession", accession)
        LOG_H.write("%s%s%s%s\n" % ("Processing pVOG ",pVOG.pVOGid," and accession ",accession))
        nextFasta = ncbiPhageSeqs.findStringInHeader(accession)     # Search sequence headers in phage/virus database; find this accn, if exists
        if nextFasta:
            nextFasta.pVOGassociationList.append(pVOG.pVOGid)       # Record this pVOG association
            pVOGstring = ""
            for pVOGidentifier in nextFasta.pVOGassociationList:
                pVOGstring += pVOGidentifier + '|'                  # Construct new header containing each associated pVOG
            nextFasta.customHeader = pVOGstring + nextFasta.header  # Update the custom header to reflect multiple pVOG membership
            if len(nextFasta.pVOGassociationList) > 1:
                LOG_H.write("%s%s%s%s\n" % ("NOTE: Fasta associated with multiple pVOGs: ", len(nextFasta.pVOGassociationList), "pVOGs:", nextFasta.pVOGassociationList))
                LOG_H.write("%s%s\n" % ("   header is: ", nextFasta.customHeader))
            foundCount += 1
        else:
            if accession not in missingList:
                missingList.append(accession)
                missingCount += 1
                PVOG_MISSING_PROTEINS_H.write("%s%s\n" % ("No pVOG found for accession ", accession))
if CHATTY:
    print("pVOG accession sequences have been tagged")
# Write the pVOG-tagged fastas to file
if CHATTY:
    print("Writing pVOG-tagged fasta sequences to file")
ncbiPhageSeqs.printMultiFasta2file_custom(PVOG_TAGGED_PROTEINS_H)
if CHATTY:
    print("done!")

LOG_H.write("%s%s\n" % ("Number accessions = ", accnCount))
LOG_H.write("%s%s\n" % ("Number of fasta sequences found = ", foundCount))
LOG_H.write("%s%s\n" % ("Number missing = ", missingCount))

# Clean up

PROTEIN_LIST_H.close()
NCBI_PHAGE_FAA_H.close()
PVOG_TAGGED_PROTEINS_H.close()
PVOG_MISSING_PROTEINS_H.close()
print("Processing complete!")
LOG_H.write("%s%s\n" % ("Processing completed at ",datetime.datetime.now()))
LOG_H.close()

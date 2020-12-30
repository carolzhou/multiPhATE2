#!/usr/bin/env python3

#####################################################################
#
# Program:  dbPrep_createVogFastaFile.py
#
# Description: This code reads input files downloaded from the VOG website,
#    (http://...),
#    and creates files containing VOG-tagged gene and protein fasta sequences.
#    Note that unlike the pVOGs, the VOG data files provide all of the 
#    fasta sequences (gene and protein), so these do not need to be queried 
#    using a large local protein database or at NCBI for slow download.
#    These files are: vog.genes.all.fa, and vog.proteins.all.fa.
#    Thus, this code uses file vog.members.tsv, which is a
#    mapping of VOGid to sequence identifier (header), and files 
#    (VOG fasta files) to re-write the fastas to include the VOG 
#    identifier(s) at the front of the fasta header for each sequence.
#
# Last Update: 21 December 2020
#
# Usage:  python3 createVogFastaFile.py
#
# Programmer:  CEZhou
#
####################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE GPL-3 LICENSE. SEE INCLUDED FILE GPL-3.PDF FOR DETAILS.

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

CODE_BASE = "createVogFastaFile"

# Paths to files; when running on mpath

VOG_BASE_DIR          = "/Users/carolzhou/DEV/PhATE/Databases/VOGs/FASTA/"

# Files

VOG_FASTA_MAP         = PVOG_BASE_DIR + "vog.members.tsv"
GENE_FNT              = PVOG_BASE_DIR + "vog.genes.all.fa"
PROTEIN_FAA           = PVOG_BASE_DIR + "vog.proteins.all.fa"
VOG_TAGGED_GENES      = PVOG_BASE_DIR + "VOGs.fnt"
VOG_TAGGED_PROTEINS   = PVOG_BASE_DIR + "VOGs.faa"
LOG                   = PVOG_BASE_DIR + CODE_BASE + ".log"
PVOG_LIB_DIR          = PVOG_BASE_DIR + "Allvogtables/"

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

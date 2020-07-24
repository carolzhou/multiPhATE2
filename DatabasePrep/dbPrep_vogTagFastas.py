#!/usr/bin/env python

###########################################################################
#
# Program: dbPrep_vogTagFastas.py
#
# Description: Processes VOG data files. Inserts VOG identifiers in header
#    string of each fasta file.  
#
# Instructions: Run this code within the /DatabasePrep/ folder
#
# Programmer:  Carol L. Ecale Zhou
#
# Most recent update: 23 July 2020 
#
###########################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.PDF FOR DETAILS.

import os, sys, re, copy, subprocess
import dbPrep_vog

# CONFIGURABLE
DATABASE_NAME = "VOGs"
DOWNLOAD_DATE = "July 2020"
VERSION       = "vog99"

##### MAIN ################################################################

if len(sys.argv) != 2:
    print ("Usage: python dbPrep_vogTagFastas.py VOGsDir")
    exit(0)
VOG_DIR = sys.argv[1]

# FILES
VOG_DIR = "../Databases/VOGs/"  # Run this code within the /DatabasePrep/ folder
#VOG_BLASTDB_DIR            = os.path.join(VOG_DIR,"BlastDBs/")
VOG_MAP_FILE               = os.path.join(VOG_DIR,"vog.members.tsv")
VOG_ANNOTATION_FILE        = os.path.join(VOG_DIR,"vog.annotations.tsv")
VOG_GENE_FASTA_FILE        = os.path.join(VOG_DIR,"vog.genes.all.fa")
VOG_PROTEIN_FASTA_FILE     = os.path.join(VOG_DIR,"vog.proteins.all.fa")
VOG_GENE_FASTA_OUT_FILE    = os.path.join(VOG_DIR,"vog.genes.tagged.all.fa")
VOG_PROTEIN_FASTA_OUT_FILE = os.path.join(VOG_DIR,"vog.proteins.tagged.all.fa")

paramSet = {
    "databaseName"           : DATABASE_NAME,
    "downloadData"           : DOWNLOAD_DATE,
    "version"                : VERSION,
    "VOGmapFile"             : VOG_MAP_FILE,
    "VOGannotationFile"      : VOG_ANNOTATION_FILE,
    "VOGgeneFastaFile"       : VOG_GENE_FASTA_FILE,
    "VOGproteinFastaFile"    : VOG_PROTEIN_FASTA_FILE,
    "VOGgeneFastaOutFile"    : VOG_GENE_FASTA_OUT_FILE,
    "VOGproteinFastaOutFile" : VOG_PROTEIN_FASTA_OUT_FILE,
}

# Create VOGs object, load data from files, create vog-tagged fasta files.
print("dbPrep_vogTagFastas says, Tagging gene and protein fasta sequences with VOG identifiers.")
cwd = os.getcwd()

myVog = dbPrep_vog.VOGs()
myVog.tagVogFastas(paramSet)

print("dbPrep_vogTagFastas says, VOG-tagging is complete.")
os.chdir(cwd)

##### END ################################################################

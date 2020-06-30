#!/usr/bin/env python

###########################################################################
#
# Program: dbPrep_vogTagFastas.py
#
# Description: Processes VOG data files. Inserts VOG identifiers in header
#    string of each fasta file.
#
# Programmer:  Carol L. Ecale Zhou
#
# Most recent update: 10 June 2020 
#
###########################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.PDF FOR DETAILS.

import os, re, copy, subprocess
import dbPrep_vog

# CONFIGURABLE
DATABASE_NAME = "VOGs"
DOWNLOAD_DATE = "June 2020"
VERSION       = "vog98"

# FILES
VOG_DIR = "/Users/zhou4/DEV/PhATE/multiPhATE2/Databases/VOGs/"
VOG_FASTA_DIR              = os.path.join(VOG_DIR,"FASTA/")
#VOG_MAP_FILE               = os.path.join(VOG_DIR,"vog.members.tsv")
VOG_MAP_FILE               = os.path.join(VOG_DIR,"vog.members.tsv")
VOG_ANNOTATION_FILE        = os.path.join(VOG_DIR,"vog.annotations.tsv")
VOG_GENE_FASTA_FILE        = os.path.join(VOG_FASTA_DIR,"vog.genes.all.fa")
VOG_PROTEIN_FASTA_FILE     = os.path.join(VOG_FASTA_DIR,"vog.proteins.all.fa")
VOG_GENE_FASTA_OUT_FILE    = os.path.join(VOG_FASTA_DIR,"vog.genes.tagged.all.fa")
VOG_PROTEIN_FASTA_OUT_FILE = os.path.join(VOG_FASTA_DIR,"vog.proteins.tagged.all.fa")

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


##### MAIN ################################################################

# Create VOGs object, load data from files, create vog-tagged fasta files.
print("dbPrep_vogTagFastas says, Tagging gene and protein fasta sequences with VOG identifiers.")
myVog = dbPrep_vog.VOGs()
myVog.tagVogFastas(paramSet)
#myVog.printAll()

print("dbPrep_vogTagFastas says, VOG-tagging is complete.")
#!/usr/bin/env python

#######################################################################################
#
# genomics_driver.py
#
# Description: This driver code runs the Genomics module, which identifies gene correspondences
#    among genomes that have been compared using CompareGeneProfiles (CGP), creates fasta
#    clusters of homologus gene and protein sequences, and builds hmms accordingly.
#
# Programmer:  C. E. Zhou
#
# Latest update:  05 May 2020
# multiPhATE2 v. 1.5
#
########################################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.PDF FOR DETAILS.

import sys, os, re, string, copy
import time, datetime
from subprocess import call
import genomics_compareGenomes

DEBUG = True
#DEBUG = False

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
PHATE_PROGRESS = True
PHATE_MESSAGES = True
PHATE_WARNINGS = True

# Environmental Variables
GENOMICS_RESULTS_DIR = os.environ["PHATE_GENOMICS_RESULTS_DIR"]

###########################################################################################
# BEGIN MAIN
#

if PHATE_PROGRESS:
    print("genomics_driver says, Performing comparisons among genomes.")

# Create genome comparison object
genomeComparison = genomics_compareGenomes.comparison()

# Create output directory
try:
    os.stat(GENOMICS_RESULTS_DIR)
except:
    os.mkdir(GENOMICS_RESULTS_DIR)

# Perform genome comparisons
genomeComparison.performComparison()

# Clean up
if PHATE_PROGRESS:
    print("genomics_driver says, Genome comparison complete.")


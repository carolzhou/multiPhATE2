#!/usr/bin/env python

#######################################################################################
#
# Name: genomics_driver.py
#
# Programmer:  C. E. Zhou
#
# Latest update:  07 December 2020
#
# Description: This driver code runs the Genomics module, which identifies gene correspondences
#    among genomes that have been compared using CompareGeneProfiles (CGP), creates fasta
#    clusters of homologus gene and protein sequences, and builds hmms accordingly.
#
########################################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory
# THIS CODE IS COVERED BY THE GPL3 LICENSE. SEE INCLUDED FILE GPL-3.PDF FOR DETAILS.

import sys, os, re, string, copy
import time, datetime
from subprocess import call
import genomics_compareGenomes
import hmm_build

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
#PHATE_PROGRESS = True
#PHATE_MESSAGES = True
#PHATE_WARNINGS = True

# Environmental Variables
GENOMICS_RESULTS_DIR = os.environ["PHATE_GENOMICS_RESULTS_DIR"]

###########################################################################################
# BEGIN MAIN
#

# Collect input parameter
referenceGenome = "" 
if len(sys.argv) == 2:
    referenceGenome = sys.argv[1]
if referenceGenome == "":
    print("genomics_driver says, ERROR: reference genome not set.")

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
paramSet = {
        "referenceGenome" : referenceGenome,
        }
genomeComparison.setParameters(paramSet)
genomeComparison.performComparison()

if PHATE_PROGRESS:
    print("genomics_driver says, Creating hmm profiles.")

# Create hmm build object
hmmBuild = hmm_build.build()
buildDirectory = ""
buildParameters = {
        "codeName"          : "hmmbuild",
        "codeVersion"       : '2',
        "alignmentCodeName" : "clustalo",
        "buildDirectory"    : GENOMICS_RESULTS_DIR,
        }
hmmBuild.setParameters(buildParameters)
if PHATE_MESSAGES:
    hmmBuild.printParameters()

# Build hmms for gene/protein groups
if PHATE_MESSAGES:
    print("genomics_driver says, fnt files for hmm build:")
    print(hmmBuild.fileList_fnt)
    print("genomics_driver says, faa files for hmm build:")
    print(hmmBuild.fileList_faa)
hmmBuild.performBuild()

# Clean up
if PHATE_PROGRESS:
    print("genomics_driver says, Genome comparison complete.")


#############################################
# Name: cleanUpFasta.py
#
# Programmer:  Carol L. Ecale Zhou
#
# Last update: 16 October 2020
#
# Description:
# This script inputs a multi-fasta file and outputs a
# multi-fasta file identical to the input but pesky
# characters are removed from the header and sequence.
# Spaces in the header are converted to '_' and 
# spaces in the sequence are removed.
#
#############################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE GPL3 LICENSE. SEE INCLUDED FILE GPL-3.PDF FOR DETAILS.

import sys
import os
import re

import cpg_fastaSequence as fastaSequence

##### FILES

LOGFILE = open ("./cleanUpFasta.log","w")  # log file
today = os.popen('date')
LOGFILE.write("%s\n" % (today.read()))
INFILE = ""
OUTFILE = ""

infile = sys.argv[1]
if infile == "-help" or infile == "help":
    print ("This script inputs a multi-fasta file and outputs a multi-fasta file identical to the input\nexcept that pesky characters are removed from the header and sequence.\nSpaces and specials in the header are converted to '_' and spaces in the sequence are removed.")
    exit(0)
elif infile == '-usage' or infile == "usage":
    print ("Usage:  cleanUpFasta.py multi.fasta <nt|aa> (optional, default is nt)")
    exit(0)
INFILE = open(infile,"r")
outfile = infile + ".new"
OUTFILE = open(outfile,"w")
sequenceType = '' 
if len(sys.argv) == 3:
    sequenceType = sys.argv[2]
else:
    sequenceType = 'nt'
if sequenceType != 'nt' and sequenceType != 'aa':
    print ("Incorrect sequence type:", sequenceType)
    exit(0)

fLines = INFILE.read().splitlines()
myFastaList = fastaSequence.multiFasta()
myFastaList.addFastas(fLines,sequenceType)
myFastaList.printMultiFasta()
myFastaList.printMultiFasta2file_case(OUTFILE,"upper")

INFILE.close()
OUTFILE.close()


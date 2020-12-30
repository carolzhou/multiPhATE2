#!/usr/bin/env python3

#######################################################################################
#
# script cleanHeader.py replaces with '_' (underscore) any characters is a fasta header
#    that is not alphanumeric or '_'.
#
# Programmer:  Carol Zhou
#
# Last Update: 21 December 2020
#
#######################################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE GPL3 LICENSE. SEE INCLUDED FILE GPL-3.pdf FOR DETAILS.

import re, sys, os

if len(sys.argv) != 2:
    print("usage: python3 cleanHeader.py <fastaFile>")
    exit(0)
else:
    INFILE_H = open(sys.argv[1],'r')
    OUTFILE_H = open("./clean.fasta",'w')
    for line in INFILE_H.read().splitlines():
        match_header = re.search('^>',line)
        if match_header:
            cleanHeader = re.sub('[^>a-zA-Z\d_]','_',line)
            OUTFILE_H.write("%s\n" % (cleanHeader))
            print(cleanHeader)
        else:
            OUTFILE_H.write("%s\n" % (line))
            print(line)
    INFILE_H.close()
    OUTFILE_H.close()

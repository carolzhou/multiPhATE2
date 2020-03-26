#!/usr/bin/env python

#################
#
# script cleanHeader.py replaces with '_' (underscore) any characters is a fasta header
#    that is not alphanumeric or '_'.
#
#################

import re, sys, os

if len(sys.argv) != 2:
    print("usage: python cleanHeader.py <fastaFile>")
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

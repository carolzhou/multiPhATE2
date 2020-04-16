#!/usr/bin/env python

#############################################################################
#
# program: convertNCBIcalls2cgc.py
#
# programmer:  C. E. Zhou
#
# Summary:  This script inputs an NCBI gene fasta file and writes a CGC-formatted file (single-line tabbed).
#
# Most recent update:  13 Jan 2020 
#
##############################################################################

import os, sys, re, time, datetime

##############################################################################

##### CONSTANTS, CONFIGURABLES

geneNo = 0
strand = '+'
start  = 0
end    = 0
length = 0
contig = "Ecoli_K12"
caller = "custom"
score  = '.'
protein = "unknown"

inStrings = []

INFILE_H = open("./PipelineOutput/ForMatt/EcoliK12.ncbi.genes.fnt",'r')
OUTFILE_H = open("custom.gff",'w')
OUTFILE_H.write("%s\n" % ("##gff-version 3"))
OUTFILE_H.write("%s\n" % ("# Sequence Data: seqnum=1;seqlen=999999;seqhdr=\"EcoliK12\""))

fLines = INFILE_H.read().splitlines()
for fLine in fLines:
    match = re.search('>',fLine)
    if match:
        print("Found a match:", fLine)
        inStrings = fLine.split(';')
        print (inStrings[0])
        fields = inStrings[0].split('#')
        print (fields)
        startCol  = fields[1]
        endCol    = fields[2]
        strandCol = fields[3]
        geneNo += 1
        if re.search('-',strandCol):
            strand = '-'
        else:
            strand = '+'
        matchNo = re.search('\d+',startCol)
        start = matchNo.group(0)
        matchNo = re.search('\d+',endCol)
        end = matchNo.group(0)
        length = abs(int(end) - int(start)) + 1
        OUTFILE_H.write("%s\t%s\t%s\t%s\t%s\t%c\t%c\t%c\t%c\n" % (contig,caller,"cds",start,end,score,strand,'.','.'))

OUTFILE_H.close()
INFILE_H.close()

#############################################################################

print ("Done!")

##############################################################################
##############################################################################

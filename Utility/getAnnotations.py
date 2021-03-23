#!/usr/bin/env python3

#######################################################################################
#
# getAnnotations.py
#
# Summary: Inputs a gene identifier (see format below) and output the annotation from
# the phate_sequenceAnnotation_main.gff file)
#
# Programmer:  Carol Zhou
#
# Last Update: 23 March 2021
#
#######################################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE GPL3 LICENSE. SEE INCLUDED FILE GPL-3.pdf FOR DETAILS.

import re, sys, os

# Help strings
HELP_STRING = """\nThis script pulls the annotation field from a genome's PhATE phate_sequenceAnnotation_main.gff file.\nInput the gene or protein identifier, followed by the full path to the genome's PhATE output subdirectory.\nNote: Use this script for printing the annotation of a single gene/protein.\nIf processing a list of genes/proteins, use getAnnotations4list.py.\nFor more information, type: python getAnnotations.py usage, or format.\n"""
FORMAT_STRING = """\nA properly formatted gene or protein identifier has the following format:\ngenomeName:contig:cds123/+/456/789\nwhere 123 is the gene number, + or - is the strand, and 456 is left end, 789 is right end.\n"""
USAGE_STRING = """\nUsage: getAnnotations.py <geneID> <PhATE subdirectory>\n  Example: getAnnotations.py A1234_c:A1234_contig_1:cds1/-/2/181 /users/me/multiPhATE2/PipelineOutput/A1234_c/\n"""

# Gather parameters. User should input a gene identifier and the PhATE annotation path/filename
if len(sys.argv) == 1:
    print(HELP_STRING)
    print(USAGE_STRING)
    exit(0)
if len(sys.argv) == 2:
    parameter = sys.argv[1]
    if re.search('help',parameter.lower()):
        print(HELP_STRING)
        exit(0)
    elif re.search('format',parameter.lower()):
        print(FORMAT_STRING)
        exit(0)
    elif re.search('usage',parameter.lower()):
        print(USAGE_STRING)
        exit(0)
if len(sys.argv) == 3:
    geneIDstring = sys.argv[1]  # Example:  A1234_c:A1234:cds235/-/5572/5897
    directory    = sys.argv[2]  # Example:  /users/home/multiPhATE2/PipelineOutput/MyGenome/"
    geneID = geneIDstring.lstrip('>') # Remove '>', just in case user entered the fasta header
    
# Check that gff file is where it needs to be
annotFile = os.path.join(directory,"phate_sequenceAnnotation_main.gff")
try:
    ANNOT_H = open(annotFile,'r') 
except(RuntimeError):
    print("Annotation file,",annotFile,"cannot be opened. Please check your input parameters.")
    exit(0)

# Check for identifier format
idGenome = ''; idContig = ''; idGeneString = ''
idString = ''; idStrand = ''; idLeftEnd = ''; idRightEnd = ''; cds = '' 
# Peel off genome name and contig
try:
    (idGenome, idContig, idGeneString) = geneID.split(':')
except(RuntimeError):
    print("ERROR: Input gene/protein identifier has unexpected format:", geneID)
# Separate gene call fields
try:
    (idString, idStrand, idLeftEnd, idRightEnd) = idGeneString.split('/')
except(RuntimeError):
    print("ERROR: Input gene/protein identifier has unexpected format:", geneID)
    exit(0)
# Get gene number
match_cds = re.search('cds(\d+)',idString)
if match_cds:
    cds = match_cds.group(1) 
else:
    print("ERROR: unexpected cds string format:",idString)
    exit(0)

# Scan through annotation file to find annotations for the input gene/protein
contig = ''; source = ''; seqType = ''; leftEnd = ''; rightEnd = ''; score = ''; strand = ''; phase = ''; annot = ''
annotationString = ''

# Identify data line(s) corresponding to the input gene identifier
aLines = ANNOT_H.read().splitlines()
for aLine in aLines:
    fields = aLine.split('\t')
    if len(fields) == 9:
        contig    = fields[0]
        source    = fields[1]
        seqType   = fields[2]
        leftEnd   = fields[3]
        rightEnd  = fields[4]
        score     = fields[5]
        strand    = fields[6]
        phase     = fields[7]
        annot     = fields[8]
    else:
        continue

    # Determine if the current data line corresponds to input gene identifier
    if idContig == contig and idStrand == strand and idLeftEnd == leftEnd and idRightEnd == rightEnd: 
        if seqType.lower() == 'gene':
            annotationString += 'gene: ' + annot         
        if seqType.lower() == 'cds': 
            if annotationString != '':
                annotationString += ' | '
            annotationString += 'protein: ' + annot
# Report
#print("Annotation(s) for gene identifier,",geneID,": ")
#print(annotationString)
print(geneID,'\t',annotationString)

# Clean up
ANNOT_H.close()


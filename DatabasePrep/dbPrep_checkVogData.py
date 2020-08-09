#!/usr/bin/env python

#######################################
#
# checkVogData.py
#
# Date: 20 July 2020
#
# Summary: This script compares the vog.genes/proteins.all/fa files to vog.genes/proteins.tagged.all.fa files
#   to determine why the latter files are smaller than the former.
#
######################################

import os, sys, re

# Load fasta headers from header files

geneHeaders          = []
proteinHeaders       = []
geneTaggedHeaders    = []
proteinTaggedHeaders = []

gDuplicates          = []
tDuplicates          = []

# Gene files

GENE_H     = open("Databases/VOGs/gene.headers",'r')
GENE_TAG_H = open("Databases/VOGs/gene_tagged.headers",'r')

gLines = GENE_H.read().splitlines()
tLines = GENE_TAG_H.read().splitlines()

for gLine in gLines:
    if gLine not in geneHeaders:
        geneHeaders.append(gLine)
    else:
        gDuplicates.append(gLine)

for tLine in tLines:
    if tLine not in geneTaggedHeaders:
        geneTaggedHeaders.append(tLine)
    else:
        tDuplicates.append(tLine)

print("Number of unique gene headers:    ",len(geneHeaders))
print("Number of duplicate gene Headers: ",len(gDuplicates))
print("Number of tagged gene headers:    ",len(geneTaggedHeaders))
print("Number of tagged duplicates:      ",len(tDuplicates))

# Remove VOG identifier(s) from tagged headers; reload into a list
geneCleanHeaders = []
p_header = re.compile(">.*\|(.*)")
duplicateHeaderCount = 0
for header in geneTaggedHeaders:
    match = re.search(p_header,header)
    cleanHeader = '>' + match.group(1)
    if cleanHeader in geneCleanHeaders:
        duplicateHeaderCount += 1
        print ("Duplicate header: ",cleanHeader)
    else:
        geneCleanHeaders.append(cleanHeader)

# Clean up; Reset (free memory)
geneTaggedHeaders = []

# Compare gene headers to gene clean headers: which headers are missing?
missingHeaders = []
missingCount = 0
matchCount = 0
for header in geneHeaders:
    if header in geneCleanHeaders:
        matchCount += 1
    else:
        missingHeaders.append(header)
        missingCount += 1
print("gene vs. gene_clean")
print("matchCount:   ",matchCount)
print("missingCount: ",missingCount)

OUT_H = open ("./missingHeaders.lst",'w')
for header in missingHeaders:
    OUT_H.write("%s\n" % (header))
"""
missingCount = 0
matchCount = 0
for header in geneCleanHeaders:
    if header in geneHeaders:
        matchCount += 1
    else:
        missingCount += 1

print("gene_clean vs. gene")
print("matchCount:   ",matchCount)
print("missingCount: ",missingCount)

# Clean up; Reset (free memory)
geneHeaders = []
missingHeaders = []
"""
#####
"""
PROTEIN_H    = open("Databases/VOGs/protein.headers",'r')
PROTEIN_TAG_H = open("Databases/VOGs/protein_tagged.headers",'r')

pLines = PROTEIN_H.read().splitlines() 
tLine  = PROTEIN_TAG_H.read().splitlines()

for pLine in pLines:
    if pLine not in proteinHeaders:
        proteinHeaders.append(pLine)
    else:
        pDuplicates.append(pLine)

tDuplicates = []  # reset
for tLine in tLines:
    if tLine not in proteinTaggedHeaders:
        proteinTaggedHeaders.append(tLine)
    else:
        tDuplicates.append(tLine)
print("*****")
print("Number of unique protein headers:    ",len(proteinHeaders))
print("Number of duplicate protein headers: ",len(pDuplicates))
print("Number of tagged protein headers:    ",len(geneTaggedHeaders))
print("Number of tagged duplicates:         ",len(tDuplicates))
"""
############################################################################

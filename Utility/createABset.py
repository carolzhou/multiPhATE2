#!/usr/bin/env python3

#######################################################################################
#
# script createABset.py reads in a compareGeneProfiles.report file and creates a "match" 
#   table between  genes and proteins of 2 genomes. 
#
# Programmer:  Carol Zhou
#
# Last Update: 03 March 2021
#
#######################################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE GPL3 LICENSE. SEE INCLUDED FILE GPL-3.pdf FOR DETAILS.

import re, sys, os

# Gather parameters
if len(sys.argv) != 2:
    print("usage: python3 createABset.py compareGeneProfiles_report_file")
    exit(0)
else:
    CGP_FILE_H = open(sys.argv[1],'r')
    OUTFILE_H  = open("./createABset.out",'w')  

# Open out files
CGP_LIST_GENE_H    = open("createABlist.gene.out",'w')
CGP_LIST_PROTEIN_H = open("createABlist.protein.out",'w')

# Extract gene/protein matches from input file.
cgpList_gene = []; cgpList_protein = []
cgpLines = CGP_FILE_H.read().splitlines()
GENES = True; PROTEINS = False
for cgpLine in cgpLines:
    if re.search('PROTEIN HITS',cgpLine):
        GENES = False; PROTEINS = True 
    if re.search('^#',cgpLine):
        continue
    if re.search('loner',cgpLine):
        continue
    fields = cgpLine.split('\t')
    if re.search('mutual',cgpLine):
        matchString = fields[8] + '|' + fields[9] + '|' + fields[10] + '|' + fields[15] + '|' + fields[16] + '|' + fields[17]
        if GENES:
            cgpList_gene.append(matchString)
            CGP_LIST_GENE_H.write("%s\n" % (matchString))
        elif PROTEINS:
            cgpList_protein.append(matchString) 
            CGP_LIST_PROTEIN_H.write("%s\n" % (matchString))

CGP_LIST_GENE_H.close()
CGP_LIST_PROTEIN_H.close()

for item in cgpList_gene:
    print(item)

print("Number of binary matches at the gene level:",len(cgpList_gene))
print("Number of binary matches at the protein level:",len(cgpList_protein))

OUTFILE_H.close()



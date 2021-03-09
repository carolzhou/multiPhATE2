#!/usr/bin/env python3

#######################################################################################
#
# script computeVenn.py calculates the portions within a Venn diagram for a 4x4 comparison
#
# Programmer:  Carol Zhou
#
# Last Update: 05 March 2021
#
#######################################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE GPL3 LICENSE. SEE INCLUDED FILE GPL-3.pdf FOR DETAILS.

import re, sys, os

CORE_GENOME     = False   # If True, then user inputs genomics_coreGenome.out
CORRESPONDENCES = False   # If True, then user inputs genomics_correspondences.out
GENE_SETS       = False   # Processing gene data
PROTEIN_SETS    = False   # Processing protein data

#CORE_GENOME     = True    # If True, then user inputs genomics_coreGenome.out
CORRESPONDENCES = True    # If True, then user inputs genomics_correspondences.out
GENE_SETS       = True    # Processing gene data
#PROTEIN_SETS    = True    # Processing protein data

abcList    = []

# Gather parameters

if len(sys.argv) != 8:
    print("usage: python3 computeVenn.py coreGenomeFile ABfile order1 ACfile order2 BCfile order3")
    exit(0)
else:
    CORE_H    = open(sys.argv[1],'r')  # {A,B,C,D}
    AB_H      = open(sys.argv[2],'r')  # {A,B}
    AB_order  = sys.argv[3]
    AC_H      = open(sys.argv[4],'r')  # {A,C}
    AC_order  = sys.argv[5]
    BC_H      = open(sys.argv[6],'r')  # {B,C}
    BC_order  = sys.argv[7]
    OUTFILE_H = open("./computeVenn.out",'w')  # {A,B,C}

print("TESTING: parameters are:",sys.argv)

# Create gene lists from binary-match files

A1list = []; B1list = []
A2list = []; C1list = []
B2list = []; C2list = []

A1LIST_OUT_H = open("A1list.out",'w')
A2LIST_OUT_H = open("A2list.out",'w')
B1LIST_OUT_H = open("B1list.out",'w')
B2LIST_OUT_H = open("B2list.out",'w')
C1LIST_OUT_H = open("C1list.out",'w')
C2LIST_OUT_H = open("C2list.out",'w')

abLines = AB_H.read().splitlines()
for abLine in abLines:
    fields  = abLine.split('|')
    if AB_order == '1':
        query   = fields[2] + '_' + fields[1]
        subject = fields[5] + '_' + fields[4]
    else:
        query   = fields[5] + '_' + fields[4]
        subject = fields[2] + '_' + fields[1]
    if query not in A1list:
        A1list.append(query)
    if subject not in B1list:
        B1list.append(subject)
print("Number of genes in A1list:",len(A1list))
print("Number of genes in B1list:",len(B1list))
print()
for gene in A1list:
    A1LIST_OUT_H.write("%s\n" % (gene))
for gene in B1list:
    B1LIST_OUT_H.write("%s\n" % (gene))

acLines = AC_H.read().splitlines()
for acLine in acLines:
    fields  = acLine.split('|')
    if AC_order == '1':
        query   = fields[2] + '_' + fields[1]
        subject = fields[5] + '_' + fields[4]
    else:
        query   = fields[5] + '_' + fields[4]
        subject = fields[2] + '_' + fields[1]
    if query not in A2list:
        A2list.append(query)
    if subject not in C1list:
        C1list.append(subject)
print("Number of genes in A2list:",len(A2list))
print("Number of genes in C1list:",len(C1list))
print()
for gene in A2list:
    A2LIST_OUT_H.write("%s\n" % (gene))
for gene in C1list:
    C1LIST_OUT_H.write("%s\n" % (gene))

bcLines = BC_H.read().splitlines()
for bcLine in bcLines:
    fields  = bcLine.split('|')
    if BC_order == '1':
        query   = fields[2] + '_' + fields[1]
        subject = fields[5] + '_' + fields[4]
    else:
        query   = fields[5] + '_' + fields[4]
        subject = fields[2] + '_' + fields[1]
    if query not in B2list:
        B2list.append(query)
    if subject not in C2list:
        C2list.append(subject)
print("Number of genes in B2list:",len(B2list))
print("Number of genes in C2list:",len(C2list))
print()
for gene in B2list:
    B2LIST_OUT_H.write("%s\n" % (gene))
for gene in C2list:
    C2LIST_OUT_H.write("%s\n" % (gene))

A1LIST_OUT_H.close()
A2LIST_OUT_H.close()
B1LIST_OUT_H.close()
B2LIST_OUT_H.close()
C1LIST_OUT_H.close()
C2LIST_OUT_H.close()

# Computer intersections between corresponding lists

Alist = []; Blist = []; Clist = []

for gene in A1list:
    if gene not in Alist:
        if gene in A2list:
            Alist.append(gene)
for gene in A2list:
    if gene not in Alist:
        if gene in A1list:
            Alist.append(gene)
print("Number of genes in Alist:",len(Alist))

for gene in B1list:
    if gene not in Blist:
        if gene in B2list:
            Blist.append(gene)
for gene in B2list:
    if gene not in Blist:
        if gene in B1list:
            Blist.append(gene)
print("Number of genes in Blist:",len(Blist))

for gene in C1list:
    if gene not in Clist:
        if gene in C2list:
            Clist.append(gene)
for gene in C2list:
    if gene not in Clist:
        if gene in C1list:
            Clist.append(gene)
print("Number of genes in Clist:",len(Clist))
print()

ALIST_OUT_H = open("Alist.out",'w')
BLIST_OUT_H = open("Blist.out",'w')
CLIST_OUT_H = open("Clist.out",'w')

for gene in Alist:
    ALIST_OUT_H.write("%s\n" % (gene))
for gene in Blist:
    BLIST_OUT_H.write("%s\n" % (gene))
for gene in Clist:
    CLIST_OUT_H.write("%s\n" % (gene))

ALIST_OUT_H.close()
BLIST_OUT_H.close()
CLIST_OUT_H.close()

# Make list of contig's genes in the core set

coreList_gene = []; coreList_protein = []; dataCount = 0
geneMatches = 0; proteinMatches = 0
referenceGeneCount = 0; referenceProteinCount = 0 # Number of reference genes/proteins that match to 3 other genomes 

if CORE_GENOME:
    PROTEIN = False
    for cLine in CORE_H.read().splitlines():
        match_protein  = re.search('CORE GENOME: PROTEIN',cLine)
        if match_protein:
            PROTEIN = True
        match_skip1 = re.search('^\*',cLine)
        match_skip2 = re.search('Set',cLine)
        if match_skip1 or match_skip2:
            continue
        dataCount += 1
        fields = cLine.split(':')
        gene = fields[1] + '_' + fields[2]
        if PROTEIN:
            if gene not in coreList_protein:
                coreList_protein.append(gene)
        else:
            if gene not in coreList_gene:
                coreList_gene.append(gene)
    CORE_LIST_H = open("coreList.out",'w')
    CORE_LIST_H.write("%s\n" % ("Core Genome Genes:"))
    for gene in coreList_gene:
        CORE_LIST_H.write("%s\n" % (gene))
    CORE_LIST_H.write("%s\n" % ("Core Genome Proteins:"))
    for protein in coreList_protein:
        CORE_LIST_H.write("%s\n" % (protein))
    CORE_LIST_H.close()
    CORE_H.close()

elif CORRESPONDENCES:
    PROTEIN = False
    coreList_gene = []; coreList_protein = []
    geneCount = 0; proteinCount = 0; # Number of genes/proteins from other genomes that match to reference gene/protein.
    GENE2 = True; PROTEIN2 = False
    for cLine in CORE_H.read().splitlines():
        nextData = ''; nextGene = ''; nextProtein = ''
        match_geneMatch = re.search('Genes corresponding to',cLine)
        match_protein = re.search('PROTEIN CORRESPONDENCES',cLine)
        match_gene    = re.search('GENE CORRESPONDENCES',cLine)
        if match_geneMatch:
            geneMatches += 1
        match_proteinMatch = re.search('Proteins corresponding to',cLine)
        if match_proteinMatch:
            proteinMatches += 1
        match_skip    = re.search('^#',cLine)
        if match_skip:
            continue
        if match_protein:
            PROTEIN = True
            continue
        if match_gene:
            continue

        if PROTEIN:
            match_text = re.search('Proteins corresponding to(.*)',cLine)
        else:
            match_text = re.search('Genes corresponding to(.*)',cLine)

        if match_text:
            nextData = match_text.group(1)   # Must remove text before data
        else:
            nextData = cLine
            nextData = nextData.lstrip(' ')  # Remove spaces from front of data line

        #print("TESTING: nextData is,",nextData)
        dataCount += 1
        fields = nextData.split(':')
        gene = fields[1] + '_' + fields[2]  

        if PROTEIN:
            if gene not in coreList_protein:
                coreList_protein.append(gene) 
        else:
            if gene not in coreList_gene:
                coreList_gene.append(gene)

        # Count number of reference genes/proteins that have 3 corresponding genes/proteins

        match_referenceGene = re.search('Genes corresponding',cLine)
        if match_referenceGene:
            GENE2 = True
            if geneCount == 3:
                referenceGeneCount += 1
            geneCount = 0 # reset 
        else:
            geneCount += 1
 
        match_referenceProtein = re.search('Proteins corresponding',cLine) 
        if match_referenceProtein or PROTEIN2:
            if GENE2:  # Record final match set
                if geneCount == 3:
                    referenceGeneCount += 1
                GENE2 = False
                PROTEIN2 = True
            if proteinCount == 3:
                referenceProteinCount += 1
                proteinCount = 0 # reset
        else:
            if PROTEIN2:
                proteinCount += 1
    # Count final set
    if proteinCount == 3:
        referenceProteinCount += 1

    CORE_LIST_H = open("coreList.out",'w')
    CORE_LIST_H.write("%s\n" % ("Core Genome Genes:"))
    for gene in coreList_gene:
        CORE_LIST_H.write("%s\n" % (gene))
    CORE_LIST_H.write("%s\n" % ("Core Genome Proteins:"))
    for protein in coreList_protein:
        CORE_LIST_H.write("%s\n" % (protein))
    CORE_LIST_H.close()
    CORE_H.close()

print("Number of reference genes with matches to 3 genomes:  ",referenceGeneCount)
print("Number of reference proteins with matches to 3 genomes:",referenceProteinCount)

print("Number of total matches in {A,B,C,D}:",dataCount)
print("Number of gene matches to reference:   ",geneMatches)
print("Number of protein matches to reference:",proteinMatches)

# Remove genes from A, B, and C lists that are already accounted for in coreList

vennListA = []
for gene in Alist:
    if GENE_SETS:
        if gene not in coreList_gene:
            vennListA.append(gene)
    elif PROTEIN_SETS:
        if gene not in coreList_protein:
            vennListA.append(gene)
vennCount = len(vennListA)
print("Number of genes in {abc}, based on A genes, is",vennCount)

vennListB = []
for gene in Blist:
    if GENE_SETS:
        if gene not in coreList_gene:
            vennListB.append(gene)
    elif PROTEIN_SETS:
        if gene not in coreList_protein:
            vennListB.append(gene)
vennCount = len(vennListB)
print("Number of genes in {abc}, based on B genes, is",vennCount)

vennListC = []
for gene in Clist:
    if GENE_SETS:
        if gene not in coreList_gene:
            vennListC.append(gene)
    elif PROTEIN_SETS:
        if gene not in coreList_protein:
            vennListC.append(gene)
vennCount = len(vennListC)
print("Number of genes in {abc}, based on C genes, is",vennCount)

# Record in files

OUTFILE_H.write("%s\n" % ("Venn List A:"))
for gene in vennListA:
    OUTFILE_H.write("%s\n" % (gene))

OUTFILE_H.write("%s\n" % ("Venn List B:"))
for gene in vennListA:
    OUTFILE_H.write("%s\n" % (gene))

OUTFILE_H.write("%s\n" % ("Venn List C:"))
for gene in vennListA:
    OUTFILE_H.write("%s\n" % (gene))

OUTFILE_H.close()


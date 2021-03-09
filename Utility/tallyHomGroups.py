#!/usr/bin/env python3

#######################################################################################
#
# script tallyHomGroups.py counts the numbers of homology groups with differenct numbers of genes 
#
# Programmer:  Carol Zhou
#
# Last Update: 07 March 2021
#
#######################################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE GPL3 LICENSE. SEE INCLUDED FILE GPL-3.pdf FOR DETAILS.

import re, sys, os
import subprocess

# Gather parameters

if len(sys.argv) != 2:
    print("usage: python3 tallyHomGroups directory")
    exit(0)
else:
    homologyDir = sys.argv[1]  # directory holding homology group files 

print("TESTING: parameters are:",sys.argv)

# Read files from directory

print("Reading files from directory",homologyDir)
p_fileName_nt = re.compile('homologyGroup_(\d+)\.fnt$')
p_fileName_aa = re.compile('homologyGroup_(\d+)\.faa$')
fileList = []
command = "ls " + homologyDir
proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
(rawresult, err) = proc.communicate()
result = rawresult.decode('utf-8')
fileList = result.split('\n')

# Count numbers of gene and protein homology groups
print("Counting numbers of gene and protein homology groups")
geneFileCount = 0; proteinFileCount = 0   # fnt and faa files, respectively
geneCoreCount = 0; proteinCoreCount = 0   # 
A0462_A0463_A0374_count_nt = 0
A0462_A0463_A0035_count_nt = 0
A0462_A0374_A0035_count_nt = 0
A0463_A0374_A0035_count_nt = 0
A0462_A0463_A0374_count_aa = 0
A0462_A0463_A0035_count_aa = 0
A0462_A0374_A0035_count_aa = 0
A0463_A0374_A0035_count_aa = 0
A0462_A0463_count_nt = 0
A0462_A0374_count_nt = 0
A0462_A0035_count_nt = 0
A0463_A0374_count_nt = 0
A0463_A0035_count_nt = 0
A0374_A0035_count_nt = 0
A0462_A0463_count_aa = 0
A0462_A0374_count_aa = 0
A0462_A0035_count_aa = 0
A0463_A0374_count_aa = 0
A0463_A0035_count_aa = 0
A0374_A0035_count_aa = 0
A0462_count_nt = 0
A0463_count_nt = 0
A0374_count_nt = 0
A0035_count_nt = 0
A0462_count_aa = 0
A0463_count_aa = 0
A0374_count_aa = 0
A0035_count_aa = 0

for fileName in fileList:
    if re.search(p_fileName_nt,fileName):
        A0462 = False; A0463 = False; A0374 = False; A0035 = False
        geneFileCount += 1
        command = "grep \'>\' " + os.path.join(homologyDir,fileName)
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        (rawresult, err) = proc.communicate()
        result = rawresult.decode('utf-8')
        if re.search('A0462',result):
            A0462 = True
        if re.search('A0463',result):
            A0463 = True
        if re.search('A0374',result):
            A0374 = True
        if re.search('A0035',result):
            A0035 = True
        if len(re.findall('A0462_c:',result)) > 1 or len(re.findall('A0463_c:',result)) > 1 or len(re.findall('A0374_c:',result)) > 1 or len(re.findall('A0035_c:',result)) > 1:
            print('TESTING: Detected a paralog inclusion for',result)
             
        # Quads 
        if A0462 and A0463 and A0374 and A0035:
            geneCoreCount += 1

        # Triplets 
        if A0462 and A0463 and A0374 and not A0035:
            A0462_A0463_A0374_count_nt += 1
        if A0462 and A0463 and A0035 and not A0374:
            A0462_A0463_A0035_count_nt += 1
        if A0462 and A0374 and A0035 and not A0463:
            A0462_A0374_A0035_count_nt += 1
        if A0463 and A0374 and A0035 and not A0462:
            A0463_A0374_A0035_count_nt += 1

        # Doublets
        if A0462 and A0463 and not A0374 and not A0035:
            A0462_A0463_count_nt += 1
        if A0462 and A0374 and not A0463 and not A0035:
            A0462_A0374_count_nt += 1
        if A0462 and A0035 and not A0463 and not A0374:
            A0462_A0035_count_nt += 1
        if A0463 and A0374 and not A0462 and not A0035:
            A0463_A0374_count_nt += 1
        if A0463 and A0035 and not A0462 and not A0374:
            A0463_A0035_count_nt += 1
        if A0374 and A0035 and not A0463 and not A0462:
            A0374_A0035_count_nt += 1

        # Singlets
        if A0462 and not (A0462 or A0374 or A0035):
            A0462_count_nt += 1
        if A0463 and not (A0462 or A0374 or A0035):
            A0463_count_nt += 1
        if A0374 and not (A0462 or A0463 or A0035):
            A0374_count_nt += 1
        if A0035 and not (A0462 or A0463 or A0374):
            A0035_count_nt += 1

    elif re.search(p_fileName_aa,fileName):
        A0462 = False; A0463 = False; A0374 = False; A0035 = False
        proteinFileCount += 1
        command = "grep \'>\' " + os.path.join(homologyDir,fileName)
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        (rawresult, err) = proc.communicate()
        result = rawresult.decode('utf-8')
        if re.search('A0462',result):
            A0462 = True
        if re.search('A0463',result):
            A0463 = True
        if re.search('A0374',result):
            A0374 = True
        if re.search('A0035',result):
            A0035 = True
        if len(re.findall('A0462_c:',result)) > 1 or len(re.findall('A0463_c:',result)) > 1 or len(re.findall('A0374_c:',result)) > 1 or len(re.findall('A0035_c:',result)) > 1:
            print('TESTING: Detected a paralog inclusion for',result)

        # Quads
        if A0462 and A0463 and A0374 and A0035:
            proteinCoreCount += 1

        # Triplets
        if A0462 and A0463 and A0374 and not A0035:
            A0462_A0463_A0374_count_aa += 1
        if A0462 and A0463 and A0035 and not A0374:
            A0462_A0463_A0035_count_aa += 1
        if A0462 and A0374 and A0035 and not A0463:
            A0462_A0374_A0035_count_aa += 1
        if A0463 and A0374 and A0035 and not A0462:
            A0463_A0374_A0035_count_aa += 1

        # Doublets
        if A0462 and A0463 and not A0374 and not A0035:
            A0462_A0463_count_nt += 1
        if A0462 and A0374 and not A0463 and not A0035:
            A0462_A0374_count_nt += 1
        if A0462 and A0035 and not A0463 and not A0374:
            A0462_A0035_count_nt += 1
        if A0463 and A0374 and not A0462 and not A0035:
            A0463_A0374_count_nt += 1
        if A0463 and A0035 and not A0462 and not A0374:
            A0463_A0035_count_nt += 1
        if A0374 and A0035 and not A0463 and not A0462:
            A0374_A0035_count_nt += 1

        # Singlets
        if A0462 and not (A0462 or A0374 or A0035):
            A0462_count_nt += 1
        if A0463 and not (A0462 or A0374 or A0035):
            A0463_count_nt += 1
        if A0374 and not (A0462 or A0463 or A0035):
            A0374_count_nt += 1
        if A0035 and not (A0462 or A0463 or A0374):
            A0035_count_nt += 1

print("Data presented below are computed based on the input PipelineOutput data sets.")
print("Record the data that corresponds to the reference genome only. For other genomes,")
print("re-run multiPhate2.py with each genome as the reference. Then, re-run this script.")
print()
print("Number of gene files:",geneFileCount,"and number of protein files:",proteinFileCount)
print("Core genome:")
print("{A0462,A0463,A0374,A0035} @ gene level:",geneCoreCount)
print("{A0462,A0463,A0374,A0035} @ protein level:",proteinCoreCount)
print()
print("Number of core genes:",geneCoreCount,"and number of core proteins:",proteinCoreCount)
print("Numbers of triplets at nucleotide level (not accounted for in higher set):")
print("{A0462,A0463,A0374}:",A0462_A0463_A0374_count_nt)
print("{A0462,A0463,A0035}:",A0462_A0463_A0035_count_nt)
print("{A0462,A0374,A0035}:",A0462_A0374_A0035_count_nt)
print("{A0463,A0374,A0035}:",A0463_A0374_A0035_count_nt)
print("Numbers of triplets at protein level (not accounted for in higher sets):")
print("{A0462,A0463,A0374}:",A0462_A0463_A0374_count_aa)
print("{A0462,A0463,A0035}:",A0462_A0463_A0035_count_aa)
print("{A0462,A0374,A0035}:",A0462_A0374_A0035_count_aa)
print("{A0463,A0374,A0035}:",A0463_A0374_A0035_count_aa)
print("Numbers of doublets at nucleotide level (not accounted for in higher sets:")
print("{A0462,A0463}:",A0462_A0463_count_nt)
print("{A0462,A0374}:",A0462_A0374_count_nt)
print("{A0462,A0035}:",A0462_A0035_count_nt)
print("{A0463,A0374}:",A0463_A0374_count_nt)
print("{A0463,A0035}:",A0463_A0035_count_nt)
print("{A0374,A0035}:",A0374_A0035_count_nt)
print("Numbers of doublets at protein level (not accounted for in higher sets:")
print("{A0462,A0463}:",A0462_A0463_count_aa)
print("{A0462,A0374}:",A0462_A0374_count_aa)
print("{A0462,A0035}:",A0462_A0035_count_aa)
print("{A0463,A0374}:",A0463_A0374_count_aa)
print("{A0463,A0035}:",A0463_A0035_count_aa)
print("{A0374,A0035}:",A0374_A0035_count_aa)

A0462_SUM = geneCoreCount \
          + A0462_A0463_A0374_count_nt + A0462_A0463_A0035_count_nt  + A0462_A0374_A0035_count_nt \
          + A0462_A0463_count_nt + A0462_A0374_count_nt + A0462_A0035_count_nt + A0462_count_nt
A0463_SUM = geneCoreCount \
          + A0462_A0463_A0374_count_nt + A0462_A0463_A0035_count_nt + A0463_A0374_A0035_count_nt \
          + A0462_A0463_count_nt + A0463_A0374_count_nt + A0463_A0035_count_nt + A0463_count_nt
A0374_SUM = geneCoreCount \
          + A0462_A0463_A0374_count_nt + A0462_A0374_A0035_count_nt + A0463_A0374_A0035_count_nt \
          + A0462_A0374_count_nt + A0463_A0374_count_nt + A0374_A0035_count_nt + A0374_count_nt
A0035_SUM = geneCoreCount \
          + A0462_A0463_A0035_count_nt + A0462_A0374_A0035_count_nt + A0463_A0374_A0035_count_nt \
          + A0462_A0035_count_nt + A0463_A0035_count_nt + A0374_A0035_count_nt + A0035_count_nt

print("Sums of Venn segments, nucleotide level")
print("Sums of A0462 quadruplets, triplets, doublets, and singlets = ",A0462_SUM)
print("Sums of A0463 quadruplets, triplets, doublets, and singlets = ",A0463_SUM)
print("Sums of A0374 quadruplets, triplets, doublets, and singlets = ",A0374_SUM)
print("Sums of A0035 quadruplets, triplets, doublets, and singlets = ",A0035_SUM)

A0462_SUM = proteinCoreCount \
          + A0462_A0463_A0374_count_aa + A0462_A0463_A0035_count_aa  + A0462_A0374_A0035_count_aa \
          + A0462_A0463_count_aa + A0462_A0374_count_aa + A0462_A0035_count_aa + A0462_count_aa
A0463_SUM = proteinCoreCount \
          + A0462_A0463_A0374_count_aa + A0462_A0463_A0035_count_aa + A0463_A0374_A0035_count_aa \
          + A0462_A0463_count_aa + A0463_A0374_count_aa + A0463_A0035_count_aa + A0463_count_aa
A0374_SUM = proteinCoreCount \
          + A0462_A0463_A0374_count_aa + A0462_A0374_A0035_count_aa + A0463_A0374_A0035_count_aa \
          + A0462_A0374_count_aa + A0463_A0374_count_aa + A0374_A0035_count_aa + A0374_count_aa
A0035_SUM = proteinCoreCount \
          + A0462_A0463_A0035_count_aa + A0462_A0374_A0035_count_aa + A0463_A0374_A0035_count_aa \
          + A0462_A0035_count_aa + A0463_A0035_count_aa + A0374_A0035_count_aa + A0035_count_aa

print("Sums of Venn segments, amino-acid level")
print("Sums of A0462 quadruplets, triplets, doublets, and singlets = ",A0462_SUM)
print("Sums of A0463 quadruplets, triplets, doublets, and singlets = ",A0463_SUM)
print("Sums of A0374 quadruplets, triplets, doublets, and singlets = ",A0374_SUM)
print("Sums of A0035 quadruplets, triplets, doublets, and singlets = ",A0035_SUM)


#!/usr/bin/env python3

#######################################################################################
#
# script tallyHomGroups.py takes as input the full path to the multiPhATE2 directory
#   containing homology groups. This is normally located under the PipelineOutput/
#   directory. For example:  $HOME/multiPhATE2/PipelineOutput/GENOMICS/HOMOLOGY_GROUPS/ 
#
# Programmer:  Carol Zhou
#
# Last Update: 25 March 2021
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
    #homologyDir = sys.argv[1]  # directory holding homology group files 
    parameter = sys.argv[1]  # should be directory holding homology group files 
if parameter.lower() == 'help':
    print("""This code takes as input the full path to the multiPhATE2 directory containing\n homology groups and writes homologous genes/proteins to files.""")
    print("Type: python3 tallyHomGroups.py <path_to_multiPhATE2_homology_groups_dir>")
    exit(0)

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

##### Count numbers of gene and protein homology groups
print("Counting numbers of gene and protein homology groups")

# Count variables
geneFileCount = 0; proteinFileCount = 0   # fnt and faa files, respectively
geneCoreCount = 0; proteinCoreCount = 0   # 
A_B_C_count_nt = 0
A_B_D_count_nt = 0
A_C_D_count_nt = 0
B_C_D_count_nt = 0
A_B_C_count_aa = 0
A_B_D_count_aa = 0
A_C_D_count_aa = 0
B_C_D_count_aa = 0
A_B_count_nt = 0
A_C_count_nt = 0
A_D_count_nt = 0
B_C_count_nt = 0
B_D_count_nt = 0
C_D_count_nt = 0
A_B_count_aa = 0
A_C_count_aa = 0
A_D_count_aa = 0
B_C_count_aa = 0
B_D_count_aa = 0
C_D_count_aa = 0
A_count_nt = 0
B_count_nt = 0
C_count_nt = 0
D_count_nt = 0
A_count_aa = 0
B_count_aa = 0
C_count_aa = 0
D_count_aa = 0

# Files for writing subsets
# Genes
GENE_CORE_FILE_H = open ("ABCD.lst",'w')
ABC_H = open("ABC.lst",'w')
ABD_H = open("ABD.lst",'w')
ACD_H = open("ACD.lst",'w')
BCD_H = open("BCD.lst",'w')
AB_H  = open("AB.lst", 'w')
AC_H  = open("AC.lst", 'w')
AD_H  = open("AD.lst", 'w')
BC_H  = open("BC.lst", 'w')
BD_H  = open("BD.lst", 'w')
CD_H  = open("CD.lst", 'w')
A_H   = open("A.lst",  'w')
B_H   = open("B.lst",  'w')
C_H   = open("C.lst",  'w')
D_H   = open("D.lst",  'w')
# Proteins
PROT_CORE_FILE_H = open ("ABCDp.lst",'w')
ABCp_H = open("ABCp.lst",'w')
ABDp_H = open("ABDp.lst",'w')
ACDp_H = open("ACDp.lst",'w')
BCDp_H = open("BCDp.lst",'w')
ABp_H  = open("ABp.lst", 'w')
ACp_H  = open("ACp.lst", 'w')
ADp_H  = open("ADp.lst", 'w')
BCp_H  = open("BCp.lst", 'w')
BDp_H  = open("BDp.lst", 'w')
CDp_H  = open("CDp.lst", 'w')
Ap_H   = open("Ap.lst",  'w')
Bp_H   = open("Bp.lst",  'w')
Cp_H   = open("Cp.lst",  'w')
Dp_H   = open("Dp.lst",  'w')
resultCount = 0

##### COUNT

for fileName in fileList:
    # GENES
    if re.search(p_fileName_nt,fileName):
        A = False; B = False; C = False; D = False
        geneFileCount += 1
        command = "grep \'>\' " + os.path.join(homologyDir,fileName)
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        (rawresult, err) = proc.communicate()
        result = rawresult.decode('utf-8')
        if re.search('A_c:',result):
            A = True
        if re.search('B_c:',result):
            B = True
        if re.search('C_c:',result):
            C = True
        if re.search('D_c:',result):
            D = True
        if len(re.findall('A_c:',result)) > 1 or len(re.findall('B_c:',result)) > 1 or len(re.findall('C_c:',result)) > 1 or len(re.findall('D_c:',result)) > 1:
            pass; #print('TESTING: Detected a paralog inclusion for',result)
             
        # Quads 
        if A and B and C and D:
            geneCoreCount += 1
            resultCount += 1
            GENE_CORE_FILE_H.write("%s\n%s\n" % (resultCount,result))

        # Triplets 
        if A and B and C and not D:
            A_B_C_count_nt += 1
            ABC_H.write("%s\n%s\n" % (A_B_C_count_nt,result))
        if A and B and D and not C:
            A_B_D_count_nt += 1
            ABD_H.write("%s\n%s\n" % (A_B_D_count_nt,result))
        if A and C and D and not B:
            A_C_D_count_nt += 1
            ACD_H.write("%s\n%s\n" % (A_C_D_count_nt,result))
        if B and C and D and not A:
            B_C_D_count_nt += 1
            BCD_H.write("%s\n%s\n" % (B_C_D_count_nt,result))

        # Doublets
        if A and B and not C and not D:
            A_B_count_nt += 1
            AB_H.write("%s\n%s\n" % (A_B_count_nt,result))
        if A and C and not B and not D:
            A_C_count_nt += 1
            AC_H.write("%s\n%s\n" % (A_C_count_nt,result))
        if A and D and not B and not C:
            A_D_count_nt += 1
            AD_H.write("%s\n%s\n" % (A_D_count_nt,result))
        if B and C and not A and not D:
            B_C_count_nt += 1
            BC_H.write("%s\n%s\n" % (B_C_count_nt,result))
        if B and D and not A and not C:
            B_D_count_nt += 1
            BD_H.write("%s\n%s\n" % (B_D_count_nt,result))
        if C and D and not A and not B:
            C_D_count_nt += 1
            CD_H.write("%s\n%s\n" % (C_D_count_nt,result))

        # Singlets
        if A and not (A or C or D):
            A_count_nt += 1
            A_H.write("%s\n%s\n" % (A_count_nt,result))
        if B and not (A or C or D):
            B_count_nt += 1
            B_H.write("%s\n%s\n" % (B_count_nt,result))
        if C and not (A or B or D):
            C_count_nt += 1
            C_H.write("%s\n%s\n" % (C_count_nt,result))
        if D and not (A or B or C):
            D_count_nt += 1
            D_H.write("%s\n%s\n" % (D_count_nt,result))

    # PROTEINS
    elif re.search(p_fileName_aa,fileName):
        A = False; B = False; C = False; D = False
        proteinFileCount += 1
        command = "grep \'>\' " + os.path.join(homologyDir,fileName)
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        (rawresult, err) = proc.communicate()
        result = rawresult.decode('utf-8')
        if re.search('A_c:',result):
            A = True
        if re.search('B_c:',result):
            B = True
        if re.search('C_c:',result):
            C = True
        if re.search('D_c:',result):
            D = True
        if len(re.findall('A_c:',result)) > 1 or len(re.findall('B_c:',result)) > 1 or len(re.findall('C_c:',result)) > 1 or len(re.findall('D_c:',result)) > 1:
            pass; #print('TESTING: Detected a paralog inclusion for',result)
            PROT_CORE_FILE_H.write("%s\n%s\n" % (A_B_C_D_count_aa,result))

        # Quads
        if A and B and C and D:
            proteinCoreCount += 1
            PROT_CORE_FILE_H.write("%s\n%s\n" % (resultCount,result))

        # Triplets
        if A and B and C and not D:
            A_B_C_count_aa += 1
            ABCp_H.write("%s\n%s\n" % (A_B_C_count_aa,result))
        if A and B and D and not C:
            A_B_D_count_aa += 1
            ABDp_H.write("%s\n%s\n" % (A_B_D_count_aa,result))
        if A and C and D and not B:
            A_C_D_count_aa += 1
            ACDp_H.write("%s\n%s\n" % (A_C_D_count_aa,result))
        if B and C and D and not A:
            B_C_D_count_aa += 1
            BCDp_H.write("%s\n%s\n" % (B_C_D_count_aa,result))

        # Doublets
        if A and B and not C and not D:
            A_B_count_aa += 1
            ABp_H.write("%s\n%s\n" % (A_B_count_aa,result))
        if A and C and not B and not D:
            A_C_count_aa += 1
            ACp_H.write("%s\n%s\n" % (A_C_count_aa,result))
        if A and D and not B and not C:
            A_D_count_aa += 1
            ADp_H.write("%s\n%s\n" % (A_D_count_aa,result))
        if B and C and not A and not D:
            B_C_count_aa += 1
            BCp_H.write("%s\n%s\n" % (B_C_count_aa,result))
        if B and D and not A and not C:
            B_D_count_aa += 1
            BDp_H.write("%s\n%s\n" % (B_D_count_aa,result))
        if C and D and not A and not B:
            C_D_count_aa += 1
            CDp_H.write("%s\n%s\n" % (C_D_count_aa,result))

        # Singlets
        if A and not (A or C or D):
            A_count_aa += 1
            Ap_H.write("%s\n%s\n" % (A_count_aa,result))
        if B and not (A or C or D):
            B_count_aa += 1
            Bp_H.write("%s\n%s\n" % (B_count_aa,result))
        if C and not (A or B or D):
            C_count_aa += 1
            Cp_H.write("%s\n%s\n" % (C_count_aa,result))
        if D and not (A or B or C):
            D_count_aa += 1
            Dp_H.write("%s\n%s\n" % (D_count_aa,result))

##### WRITE TO FILES

GENE_CORE_FILE_H.close()
PROT_CORE_FILE_H.close()
ABC_H.close()
ABD_H.close()
ACD_H.close()
BCD_H.close()
AB_H.close()
AC_H.close()
AD_H.close()
BC_H.close()
BD_H.close()
CD_H.close()
A_H.close()
B_H.close()
C_H.close()
D_H.close()
ABCp_H.close()
ABDp_H.close()
ACDp_H.close()
BCDp_H.close()
ABp_H.close()
ACp_H.close()
ADp_H.close()
BCp_H.close()
BDp_H.close()
CDp_H.close()
Ap_H.close()
Bp_H.close()
Cp_H.close()
Dp_H.close()

##### REPORT

print("Data presented below are computed based on the input PipelineOutput data sets.")
print("Record the data that corresponds to the reference genome only. For other genomes,")
print("re-run multiPhate2.py with each genome as the reference. Then, re-run this script.")
print()
print("Number of gene files:",geneFileCount,"and number of protein files:",proteinFileCount)
print("Core genome:")
print("{A,B,C,D} @ gene level:",geneCoreCount)
print("{A,B,C,D} @ protein level:",proteinCoreCount)
print()
print("Number of core genes:",geneCoreCount,"and number of core proteins:",proteinCoreCount)
print("Numbers of triplets at nucleotide level (not accounted for in higher set):")
print("{A,B,C}:",A_B_C_count_nt)
print("{A,B,D}:",A_B_D_count_nt)
print("{A,C,D}:",A_C_D_count_nt)
print("{B,C,D}:",B_C_D_count_nt)
print("Numbers of triplets at protein level (not accounted for in higher sets):")
print("{A,B,C}:",A_B_C_count_aa)
print("{A,B,D}:",A_B_D_count_aa)
print("{A,C,D}:",A_C_D_count_aa)
print("{B,C,D}:",B_C_D_count_aa)
print("Numbers of doublets at nucleotide level (not accounted for in higher sets:")
print("{A,B}:",A_B_count_nt)
print("{A,C}:",A_C_count_nt)
print("{A,D}:",A_D_count_nt)
print("{B,C}:",B_C_count_nt)
print("{B,D}:",B_D_count_nt)
print("{C,D}:",C_D_count_nt)
print("Numbers of doublets at protein level (not accounted for in higher sets:")
print("{A,B}:",A_B_count_aa)
print("{A,C}:",A_C_count_aa)
print("{A,D}:",A_D_count_aa)
print("{B,C}:",B_C_count_aa)
print("{B,D}:",B_D_count_aa)
print("{C,D}:",C_D_count_aa)

A_SUM = geneCoreCount \
          + A_B_C_count_nt + A_B_D_count_nt  + A_C_D_count_nt \
          + A_B_count_nt + A_C_count_nt + A_D_count_nt + A_count_nt
B_SUM = geneCoreCount \
          + A_B_C_count_nt + A_B_D_count_nt + B_C_D_count_nt \
          + A_B_count_nt + B_C_count_nt + B_D_count_nt + B_count_nt
C_SUM = geneCoreCount \
          + A_B_C_count_nt + A_C_D_count_nt + B_C_D_count_nt \
          + A_C_count_nt + B_C_count_nt + C_D_count_nt + C_count_nt
D_SUM = geneCoreCount \
          + A_B_D_count_nt + A_C_D_count_nt + B_C_D_count_nt \
          + A_D_count_nt + B_D_count_nt + C_D_count_nt + D_count_nt

print("Sums of Venn segments, nucleotide level")
print("Sums of A quadruplets, triplets, doublets, and singlets = ",A_SUM)
print("Sums of B quadruplets, triplets, doublets, and singlets = ",B_SUM)
print("Sums of C quadruplets, triplets, doublets, and singlets = ",C_SUM)
print("Sums of D quadruplets, triplets, doublets, and singlets = ",D_SUM)

A_SUM = proteinCoreCount \
          + A_B_C_count_aa + A_B_D_count_aa  + A_C_D_count_aa \
          + A_B_count_aa + A_C_count_aa + A_D_count_aa + A_count_aa
B_SUM = proteinCoreCount \
          + A_B_C_count_aa + A_B_D_count_aa + B_C_D_count_aa \
          + A_B_count_aa + B_C_count_aa + B_D_count_aa + B_count_aa
C_SUM = proteinCoreCount \
          + A_B_C_count_aa + A_C_D_count_aa + B_C_D_count_aa \
          + A_C_count_aa + B_C_count_aa + C_D_count_aa + C_count_aa
D_SUM = proteinCoreCount \
          + A_B_D_count_aa + A_C_D_count_aa + B_C_D_count_aa \
          + A_D_count_aa + B_D_count_aa + C_D_count_aa + D_count_aa

print("Sums of Venn segments, amino-acid level")
print("Sums of A quadruplets, triplets, doublets, and singlets = ",A_SUM)
print("Sums of B quadruplets, triplets, doublets, and singlets = ",B_SUM)
print("Sums of C quadruplets, triplets, doublets, and singlets = ",C_SUM)
print("Sums of D quadruplets, triplets, doublets, and singlets = ",D_SUM)


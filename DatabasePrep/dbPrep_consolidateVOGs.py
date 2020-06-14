##########################################################################
#
# module dbPrep_consolidateVOGs.py 
#
# Description:  Lists vog.hmm directory, concatenates all .hmm files into one. 
# Note:  Mac command-line utilities cannot handle the size/length of arguments required to do this simply.
#
# Programmer:  C. E. Zhou
#
#########################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.PDF FOR DETAILS.

#DEBUG = True
DEBUG = False 
SYSTEM = 'MAC'     # Controls which method used for system call
#SYSTEM = 'LINUX'

import re
import os
import copy
import subprocess

#hmmDir  = "../Databases/VOGhmms/"
hmmDir  = "/Users/zhou4/DEV/PhATE/multiPhATE2/Databases/VOGhmms/"
hmmSubdir = os.path.join(hmmDir,"vog.hmm/")
consolidationFile = os.path.join(hmmDir,"vogs.hmm")
try:
    subprocess.call(['touch', consolidationFile])
except:
    print("Cannot touch file ",consolidationFile)
command = 'ls ' + hmmSubdir
proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
(rawresult, err) = proc.communicate()
result = rawresult.decode('utf-8')
fileList = str(result).split('\n')

print("hmmSubdir is ",hmmSubdir)
print("consolidationFile is ",consolidationFile)
print("command is ",command)

count = 0
for filename in fileList:
    count += 1
    print("Concatenating ",filename)
    file2cat = os.path.join(hmmSubdir,filename)
    command = 'cat ' + file2cat + ' >> ' + consolidationFile 
    print("Concatenation command is ",command)
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    (rawresult, err) = proc.communicate()
    result = rawresult.decode('utf-8')
    print(result)
   
print("Concatenated hmms are in file ",consolidationFile)

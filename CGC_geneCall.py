################################################################################################
#
# Module:  CGC_geneCall.py
#
# Programmer:  Carol Zhou
#
# Description:  Module containing classes and methods for handling data output from a gene caller 
#
# Classes and Methods:
#    GeneCall()
#        AssignGeneCall(<input parameters>)
#        PrintAll()
#        PrintAll_brief()
#    GeneCallSet()
#        AddGeneCall(newGeneCall)
#        AddGeneCalls(GENE_FILE_HANDLE)
#        IsLesser(gene1,gene2)
#        UpdateGeneCount()
#        GetGeneCalls()
#        SortGeneCalls()
#        PrintAll()
#        PrintAll_brief()
#
#################################################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.pdf FOR DETAILS.

import re
import copy
import os

PHATE_PIPELINE = True  # Running this code within the PhATE pipeline. Set this to False if running code independently
#PHATE_PIPELINE = False

##### Verbosity

#if PHATE_PIPELINE:
#    CGC_WARNINGS = os.environ["CGC_WARNINGS"]
#    CGC_MESSAGES = os.environ["CGC_MESSAGES"]
#    CGC_PROGRESS = os.environ["CGC_PROGRESS"]
#else:
#    CGC_WARNINGS = 'True'
#    CGC_MESSAGES = 'True'
#    CGC_PROGRESS = 'True'

CGC_WARNINGS = 'True'
CGC_MESSAGES = 'True'
CGC_PROGRESS = 'True'
#DEBUG = True
DEBUG = False

p_comment    = re.compile('^#')
p_caller     = re.compile('([\w\d]+)\sgene\scalls')
p_callerName = re.compile('[Gg][Ee][Nn][Ee][Mm][Aa][Rr][Kk]|[Gg][Ll][Ii][Mm][Mm][Ee][Rr]|[Pp][Rr][Oo][Dd][Ii][Gg][Aa][Ll]|[Rr][Aa][Ss][Tt]|[Tt][Hh][Ee][Aa]|[Pp][Hh][Aa][Nn][Oo][Tt][Aa][Tt][Ee]|[Gg][Ff][Ff][3]')
#p_callerName = re.compile('[Gg][Ee][Nn][Ee][Mm][Aa][Rr][Kk]|[Gg][Ll][Ii][Mm][Mm][Ee][Rr]|[Pp][Rr][Oo][Dd][Ii][Gg][Aa][Ll]|[Rr][Aa][Ss][Tt]|[Pp][Hh][Aa][Tt][Ee]')
p_dataLine   = re.compile('^(\d+)\t([+-])\t(\d+)\t(\d+)\t(\d+)\t([\d\w\.\-\_]+)')

class GeneCall(object):
    
    def __init__(self):
        self.geneName    = "unknown"
        self.geneCaller  = "unknown"
        self.geneNumber  = 0
        self.strand      = 'x'  # Typically '+' or '-'; could be '?'; 'x' indicates NULL
        self.leftEnd     = 0    # may be start or stop, depending on strand (orientation)
        self.rightEnd    = 0
        self.geneLength  = 0
        self.contig      = "unknown"
        self.proteinName = "unknown"
        self.score1      = 1.0  # Gene-call score, based on consensus among gene callers; calculated in CGC_compare/Score().
        self.score2      = 0.0  # Annotation score, based broadly on similarity to entries in the databases.

    def AssignGeneCall(self,geneName,geneCaller,geneNumber,strand,leftEnd,rightEnd,geneLength,contig="unknown",protein="unknown"):
        self.geneName    = geneName
        self.geneCaller  = geneCaller
        self.geneNumber  = geneNumber
        self.strand      = strand 
        self.leftEnd     = leftEnd 
        self.rightEnd    = rightEnd 
        self.geneLength  = geneLength
        self.contig      = contig
        self.proteinName = protein
        return

    def PrintAll(self):
        print("\ngeneName =", self.geneName)
        print("geneCaller =", self.geneCaller)
        print("geneNumber =", self.geneNumber)
        print("leftEnd =",    self.leftEnd)
        print("rightEnd =",   self.rightEnd)
        print("strand =",     self.strand)
        print("length =",     self.geneLength)
        print("contig =",     self.contig)
        print("protein =",    self.proteinName)
        print("score1 =",     self.score1)
        print("score2 =",     self.score2)
        return

    def PrintAll2file(self,FILE_H):
        FILE_H.write("%s%s\n" % ("\ngeneName = ", self.geneName))
        FILE_H.write("%s%s\n" % ("geneCaller = ", self.geneCaller))
        FILE_H.write("%s%s\n" % ("geneNumber = ", self.geneNumber))
        FILE_H.write("%s%s\n" % ("leftEnd = ",    self.leftEnd))
        FILE_H.write("%s%s\n" % ("rightEnd = ",   self.rightEnd))
        FILE_H.write("%s%s\n" % ("strand = ",     self.strand))
        FILE_H.write("%s%s\n" % ("length = ",     self.geneLength))
        FILE_H.write("%s%s\n" % ("contig = ",     self.contig))
        FILE_H.write("%s%s\n" % ("protein = ",    self.proteinName))
        FILE_H.write("%s%s\n" % ("score1 = ",     self.score1))
        FILE_H.write("%s%s\n" % ("score2 = ",     self.score2))
        return

    def PrintAll_brief(self):
        print("Gene No.", self.geneNumber, "gene caller: ", self.geneCaller, ", leftEnd:", self.leftEnd, ", rightEnd:", self.rightEnd, ", strand:", self.strand, ", length:", self.geneLength, ", contig:", self.contig, ", protein:", self.proteinName, ", score1:", self.score1, ", score2:", self.score2)
        return

    def PrintAll_brief_2file(self, File_H):
        File_H.write("%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s\n" % ("Gene No. ",self.geneNumber,", gene caller: ",self.geneCaller,", leftEnd: ", self.leftEnd,", rightEnd: ",self.rightEnd,", strand: ",self.strand,", length: ",self.geneLength,", contig: ",self.contig,", protein: ",self.proteinName,", score1: ", self.score1,", score2: ", self.score2))
        return

class GeneCallSet(object):

    def __init__(self):
        self.geneCaller     = ""  # Typically, 'GeneMark', 'Glimmer', 'Prodigal', 'RAST', 'PHANOTATE'
        self.geneCount      = 0
        self.geneCallList   = []  # list of GeneCall objects 
        self.geneCall_obj   = GeneCall()

    def UpdateGeneCount(self):

        self.geneCount = len(self.geneCallList)
        return self.geneCount

    def GetGeneCalls(self,fLines,geneCaller):

        for line in fLines:
            match_data = re.search(p_dataLine,line)
            if match_data:
                geneNumber = match_data.group(1)
                strand     = match_data.group(2)
                leftEnd    = match_data.group(3)
                rightEnd   = match_data.group(4)
                geneLength = match_data.group(5) 
                contig     = match_data.group(6)
                geneName   = self.geneCaller + '_' + geneNumber
                newGeneCall = copy.deepcopy(self.geneCall_obj)
                newGeneCall.AssignGeneCall(geneName,geneCaller,geneNumber,strand,leftEnd,rightEnd,geneLength,contig)
                self.AddGeneCall(newGeneCall)
        return

    def AddGeneCall(self,newGeneCall):

        self.geneCallList.append(newGeneCall)
        self.UpdateGeneCount()
        return

    def AddGeneCalls(self,GENE_FILE_HANDLE):

        fLines = GENE_FILE_HANDLE.read().splitlines()
        for line in fLines:
            match_caller = re.search(p_caller,line)
            if match_caller:
                caller = match_caller.group(1).lower()
                match_callerName = re.search(p_callerName,caller)
                if match_callerName:
                    self.geneCaller = caller 
                    self.GetGeneCalls(fLines,caller)
                else:
                    if CGC_WARNINGS == 'True':
                        print("ERROR in CGC_geneCall: gene caller not recognized in geneCall.GeneCallSet,", caller, line)
        return

    # Determine which of 2 gene calls occurs first along the sequence (left to right, regardless of orientation) 
    def IsLesser(self,gene1,gene2):  # Input is 2 geneCall objects

        # Sort on left end position
        if (int(gene1.leftEnd) < int(gene2.leftEnd)):
            return True

        # Sort on end position if start positions are equal
        if (int(gene1.leftEnd) == int(gene2.leftEnd)) and (int(gene1.rightEnd) < int(gene2.rightEnd)):
            return True 

        # If you're still here, then gene1 > gene2
        return False 
 
    def SortGeneCalls(self):

        # First, sort the gene call list by contig name
        self.geneCallList.sort(key=lambda x: x.contig) 

        # Second, make a list of the contigs on which genes were called; then sort (alphabetically)
        contigList = []
        for index in range(0,len(self.geneCallList)):
            contig = self.geneCallList[index].contig
            if contig not in contigList:
                contigList.append(contig) 
        contigList.sort()

        # Now, walk through the list of gene call objects, for each contig group, order by gene coordinates 
        for contig in contigList:  # Order each contig's gene calls
            startIndex = 0; endIndex = 0   # initialize; where a given contig's data starts/ends within self.geneCallList
            index1 = 0; index2 = 0         # for looping through self.geneCallList
            temp = []                      # temp list for sorting gene-call data for a given contig
            END_FOUND = False              # controls search for given contig's data start/end within self.geneCallList

            # Find position where this contig's data begins
            for index1 in range(0,len(self.geneCallList)):
                if self.geneCallList[index1].contig == contig:  # walk through list until position where contig starts is found
                    startIndex = index1
                    break
            # Find position where this contig's data ends 
            for index2 in range(index1,len(self.geneCallList)):
                if self.geneCallList[index2].contig != contig:
                    endIndex = index2
                    END_FOUND = True
                    break
            if not END_FOUND:  # For the last contig, a change in contig will not be found, so endIndex is the end of the list
                endIndex = len(self.geneCallList)

            # Order the gene call objects for the current contig only
            for i in range(startIndex,endIndex):
                temp.append(self.geneCallList[i])

            # First sort by gene caller, then by second coordinate, then by first
            temp.sort(key=lambda x: x.geneCaller)
            temp.sort(key=lambda x: int(x.rightEnd))
            temp.sort(key=lambda x: int(x.leftEnd))

            # Replace the gene call list segment corresponding to current contig with the sorted temp
            last = len(self.geneCallList) + 1
            self.geneCallList = self.geneCallList[0:startIndex] + temp + self.geneCallList[endIndex:last]

        return

    def Swap(self, geneCallObj1, geneCallObj2):
        temp = GeneCall()
        temp = geneCallObj1
        geneCallObj1 = geneCallObj2
        geneCallObj2 = temp
        return 

    def PrintAll(self):
        print("Gene Caller: ",self.geneCaller)
        for gene in self.geneCallList:
            gene.PrintAll()
        return

    def PrintAll2file(self, FILE_H):
        FILE_H.write("%s%s\n" % ("Gene Caller: ",self.geneCaller))
        for gene in self.geneCallList:
            gene.PrintAll2file(FILE_H)
        return

    def PrintAll_brief(self):
        print("Gene Caller: ",self.geneCaller)
        for gene in self.geneCallList:
            gene.PrintAll_brief()
        return

    def PrintAll_brief_2file(self,File_H):
        File_H.write("%s%s\n" % ("Gene Caller: ",self.geneCaller))
        for gene in self.geneCallList:
            gene.PrintAll_brief_2file(File_H)
        return


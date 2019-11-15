###################################################################################################
#
# Module:  CGC_compare.py
#
# Programmer:  Carol Zhou
#
# Description:  Module containing classes and methods for comparing results from different gene
#    callers.
#
# Updates:
#    12 October 2018
#
# Programmer's Notes:
#
# Classes and Methods:
#    Comparison
#        IdentifyCallers()
#        IdentifyCommonCore()
#        IsLesser(gene1,gene2)
#        Merge(nextGeneSet)
#        Compare()
#        PrintMergeList()
#        PrintUniqueList()
#        PrintCommonCore()
#        PrintCallerList()
#        PrintGenecallGrid()
#        PrintReport()
#        PrintStats()
#        PrintAll()
#        PrintAll_verbose()
#
###################################################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.pdf FOR DETAILS.

import os, re
import copy
import CGC_geneCall
import math
from collections import defaultdict

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
CGC_MESSAGES = 'False'
CGC_PROGRESS = 'True'
#DEBUG = True  # Controls verbosity in this code only
DEBUG = False

p_comment   = re.compile('^#')

class Comparison(object):
    
    def __init__(self):
        self.commonCore      = []  # list of lists of identical CGC_geneCall object calls (ie, different callers, same call) 
        self.mergeList       = []  # combined list of CGC_geneCall objects, merged by self.Merge()
        self.uniqueList      = []  # list of lists of unique gene calls over all callers; each item in list is a common gene call (>=1 gene caller)
        self.callerList      = []  # non-redundant list of callers 
        self.geneCall        = CGC_geneCall.GeneCall()  # a geneCall object
        self.averageScores   = defaultdict(dict)  # Holds average gene-call score for each caller 

    ##### IDENTIFY CALLERS ##### 
    # Create a non-redundant list of gene callers
    #################
    def IdentifyCallers(self):
        if self.mergeList:
            for gene in self.mergeList:
                if gene.geneCaller not in self.callerList:
                    self.callerList.append(gene.geneCaller)
            self.callerList.sort()
            if DEBUG:
                print("DEBUG: CGC_compare/IdentifyCallers(): self.callerList is", self.callerList)
            return len(self.callerList)
        else:
            if CGC_WARNINGS == 'True':
                print("WARNING in CGC_compare module: IdentifyCallers(): No callers to extract: call method Merge() to establish mergeList before calling this method") 
            return 0


    ##### IDENTIFY COMMON CORE #####
    # Run Merge() and Compare() before running this method   
    #################
    def IdentifyCommonCore(self):
        # First, determine the callers used  
        if self.uniqueList:  # Must have previously called self.Compare() to fill this list 
            if self.mergeList:
                callerCount = self.IdentifyCallers() 
                if callerCount > 0:
                    for commonCalls in self.uniqueList:
                        count = 1
                        if len(commonCalls) == callerCount:
                            newCommonCoreCall = copy.deepcopy(self.geneCall)
                            geneName   = "CommonCoreGene_" + str(count)
                            strand     = commonCalls[0].strand
                            leftEnd    = commonCalls[0].leftEnd
                            rightEnd   = commonCalls[0].rightEnd
                            geneLength = commonCalls[0].geneLength
                            newCommonCoreCall.AssignGeneCall(geneName,"All_callers",count,strand,leftEnd,rightEnd,geneLength)
                            self.commonCore.append(newCommonCoreCall)
                            count += 1
                else:
                    if CGC_WARNINGS == 'True':
                        print("WARNING in CGC_compare: IdentifyCommonCore(): callerCount is zero! cannot process")
            else:
                if CGC_WARNINGS == 'True':
                    print("WARNING in CGC_compare module: IdentifyCommonCore(): MergeList is empty:  need to run self.Merge()")
        else:
            if CGC_WARNINGS == 'True':
                print("WARNING in CGC_compare module: IdentifyCommonCore(): No data available to identify common core")
        return 

 
    ##### IS LESSER ? #####
    # Determine which gene call occurs first along the sequence
    #################
    def IsLesser(self,gene1,gene2):  # input is 2 geneCall objects
        if (int(gene1.leftEnd) < int(gene2.leftEnd)):
            return True
        if (int(gene1.leftEnd) == int(gene2.leftEnd)) and (int(gene1.rightEnd) < int(gene2.rightEnd)):
            return True
        return False


    ##### MERGE #####
    # Merge a list of gene call objects with self.mergeList
    # Call this method once for each caller's output (i.e., loop over the set of gene caller outputs) 
    #################
    def Merge(self,nextGeneSet):  # Merge a list of gene call objects with self.mergeList
        contigList = []           # non-redundant list of contigs from *both* lists (self.mergeList and incoming nextGeneSet)

        # Add new contigs to the merge list
        for geneCall in nextGeneSet:
            newGeneCall = copy.deepcopy(geneCall)
            self.mergeList.append(newGeneCall)

        # For bookkeeping, compile a non-redundant list of the contigs upon which genes were called; then sort the list
        for index in range(0,len(self.mergeList)):
            contig = self.mergeList[index].contig
            if contig not in contigList:
                contigList.append(contig)
        contigList.sort()

        # Next, sort the merge list by contig name; this groups all the gene calls together that are on the same contig
        self.mergeList.sort(key=lambda x: x.contig)
        
        # Now walk through the list of gene call objects, for each contig group, and order the gene calls by coordinates 
        for contig in contigList:
            startIndex = 0; endIndex = 0   # Initialize; where a given contig's data starts/ends within self.mergeList
            index1 = 0; index2 = 0         # for looping through self.mergeList
            temp = []                      # temp list for sorting gene-call data for a given contig
            END_FOUND = False              # catches last contig group at end of list

            # Find position where this contig's data begins
            for index1 in range(0,len(self.mergeList)):
                if self.mergeList[index1].contig == contig:   # walk through list until position where contig starts is found
                    startIndex = index1
                    break
            # Find position where this contig's data ends
            for index2 in range(index1,len(self.mergeList)):
                if self.mergeList[index2].contig != contig:  # continue walking through list until contig changes
                    endIndex = index2 
                    END_FOUND = True
                    break
            if not END_FOUND:  # For the last contig, a change in contig will not be found, so endIndex is the end of the list
                endIndex = len(self.mergeList)

            # Order the gene call objects for the current contig only
            for i in range(startIndex,endIndex):   
                temp.append(self.mergeList[i])  # copy out these gene call objects

            # Sort by gene caller, then by second coordinate, then by first
            temp.sort(key=lambda x: x.geneCaller)
            temp.sort(key=lambda x: int(x.rightEnd)) 
            temp.sort(key=lambda x: int(x.leftEnd)) 

            # Replace merge list segment corresponding to current contig with the sorted temp 
            last = len(self.mergeList) + 1
            self.mergeList = self.mergeList[0:startIndex] + temp + self.mergeList[endIndex:last]  

        return

    ##### COMPARE #####
    # Compares the genes in self.mergeList; creates a list of ordered, unique gene calls 
    # Run this method after having merged all of your gene call sets into self.mergeList
    #################
    def Compare(self):  
        identityList = []
        if self.mergeList:
            callCount = len(self.mergeList)               # number of total gene calls, all callers
            nextCall  = copy.deepcopy(self.mergeList[0])  # capture 1st gene call
            identityList.append(nextCall)                 # each addition to list is identical to existing in list
            for i in range(1,callCount):                 # start with 2nd gene call
                nextCall = copy.deepcopy(self.mergeList[i])
                if self.mergeList[i].strand   != self.mergeList[i-1].strand  or \
                   self.mergeList[i].leftEnd  != self.mergeList[i-1].leftEnd or \
                   self.mergeList[i].rightEnd != self.mergeList[i-1].rightEnd or \
                   self.mergeList[i].contig   != self.mergeList[i-1].contig:
                    self.uniqueList.append(identityList)  # all identicals for this gene call are identified
                    identityList = []                     # reset
                identityList.append(nextCall)
            if identityList:
                self.uniqueList.append(identityList)
        else:
            if CGC_MESSAGES == 'True':
                print("CGC_compare says: Compare(): Nothing to Compare")
        return

    ##### SCORE ##### 
    # Scores each gene call based on its consensus with other callers
    # Run this method last. This method must be run AFTER merge list is complete,
    #   and the self.Compare method has been run, generating the self.uniqueList,
    #   which is a superset of the gene calls called by each of the callers.
    # Unique gene calls are scored 0.0
    # Unanimous gene calls are scored 1.0
    # Gene calls that are in common with one or more other gene callers are scored based on how many
    #   gene callers were in agreement. If there is disagreement in the start coordinate only, then additional
    #   (partial) score is awarded.
    #################
    def Score(self):
        self.IdentifyCallers()  # First, make sure that we have a complete, non-redundant list of callers
        geneCallerCount = len(self.callerList)  # And we know how many callers there were
        geneScore_a = 0.0                       # Gene-call score based on identical calls
        geneScore_b = 0.0                       # Gene-call score "boost" based on same stop coordinate
        contig_strand_start_stop = ""           # Functions to tally number of identical calls 
        contig_strand_stop       = ""           # Functions to tally number of calls with same stop coordinate
        callTally = defaultdict(dict)           # Hash for tallying identical and similar calls; defaultdict is like a Perl hash
        start = 0; stop = 0                     # start and stop coordinates

        if DEBUG:
            self.PrintUniqueList()

        # Compute gene-call scores
        if geneCallerCount > 0:  # No need to score if only one gene caller was invoked
        
            # First, tally up the number of instances of identical and same-stop (similar) gene calls, on each contig
            for identityList in self.uniqueList:
                for geneCall in identityList:
                    start = 0; stop = 0
                    if geneCall.strand == '+':
                        start = geneCall.leftEnd
                        stop  = geneCall.rightEnd
                    elif geneCall.strand == '-':
                        start = geneCall.rightEnd
                        stop  = geneCall.leftEnd
                    else:
                        if CGC_WARNINGS == 'True':
                            print("WARNING: CGC_compare says, Unrecognized strand:", geneCall.strand) 
                    contig_strand_start_stop = geneCall.contig + '|' + geneCall.strand + '|' + str(start) + '|' + str(stop) 
                    contig_strand_stop       = geneCall.contig + '|' + geneCall.strand + '|'                    + str(stop)

                    if callTally[contig_strand_start_stop]:      # If returns 'none', then this key does not yet exist in "hash"
                        callTally[contig_strand_start_stop] += 1
                    else:
                        callTally[contig_strand_start_stop] = 1  

                    if callTally[contig_strand_stop]:
                        callTally[contig_strand_stop] += 1   
                    else:
                        callTally[contig_strand_stop] = 1  
            if DEBUG:
                print("gene caller count:", geneCallerCount)
                print("callTally:", callTally)

            # Next, compute score for each gene call; this is done is two steps: identity, then similarity
            for identityList in self.uniqueList:
                for geneCall in identityList:
                    # initialize
                    geneScore_a = 0.0; geneScore_b = 0.0 
                    identicals = 0; similars = 0        
                    start = 0; stop = 0

                    # Formulate keys for callTally hash
                    if geneCall.strand == '+':
                        start = geneCall.leftEnd
                        stop  = geneCall.rightEnd
                    elif geneCall.strand == '-':
                        start = geneCall.rightEnd
                        stop  = geneCall.leftEnd
                    else:
                        if DEBUG:
                            print("DEBUG: CGC_compare says, Unrecognized strand:", geneCall.strand) 
                    contig_strand_start_stop = geneCall.contig + '|' + geneCall.strand + '|' + str(start) + '|' + str(stop) 
                    contig_strand_stop       = geneCall.contig + '|' + geneCall.strand + '|'                    + str(stop)

                    # Compute number of identicals and similars wrt current gene call
                    if callTally[contig_strand_start_stop]:
                        identicals = callTally[contig_strand_start_stop]
                    if callTally[contig_strand_stop]:
                        similars = callTally[contig_strand_stop] 

                    # Compute gene score part a  # score based on identical calls
                    geneScore_a = len(identityList)/float(geneCallerCount)

                    # Compute gene score part b  # score based on similar calls: stop coordinates same, start coordinates differ
                    # Only need to compute geneScore_b if at least 1 similar call exists (different start)
                    if (geneScore_a < 1.0) or (similars - identicals) > 0:  # geneScore_a == 1.0 ==> unanimous among all callers
                        if callTally[contig_strand_stop] > 1:               # must be at least 2 similars 
                            geneScore_b = 0.5 / float(geneCallerCount)      # at a fractional score to "boost" due to similarity
                    if DEBUG:
                        print("geneScore_a is", geneScore_a, "for gene call", contig_strand_start_stop) 
                        print("geneScore_b is", geneScore_b, "for gene call", contig_strand_stop)

                    # Compute and record final gene-call score
                    geneCall.score1 = geneScore_a + geneScore_b

            # Lastly, compute the average gene-call scores for each gene caller
            self.AverageGeneCallScores()

        return

    ##### COMPUTE SUMMARY GENE SCORE
    # Averages the gene-call scores for each gene caller
    #################
    def AverageGeneCallScores(self):
        cumulativeScores = defaultdict(dict)
        callCounts       = defaultdict(dict)

        for identityList in self.uniqueList:
            for geneCall in identityList:
                if cumulativeScores[geneCall.geneCaller]:   # First check if this key exists yet
                    cumulativeScores[geneCall.geneCaller] += geneCall.score1 
                    callCounts[geneCall.geneCaller]       += 1
                else:
                    cumulativeScores[geneCall.geneCaller] = geneCall.score1  # Initialize this key and record 1st value
                    callCounts[geneCall.geneCaller]       = 1

        for caller in self.callerList:
            self.averageScores[caller] = cumulativeScores[caller] / callCounts[caller]

        return

    ##### IS IDENTICAL #####
    # Determines whether two gene calls are identical
    #################
    def IsIdentical(self, geneCall1, geneCall2):
        if geneCall1.geneCaller == geneCall2.geneCaller and \
            geneCall1.contig    == geneCall2.contig     and \
            geneCall1.leftEnd   == geneCall2.leftEnd    and \
            geneCall1.rightEnd  == geneCall2.rightEnd   and \
            geneCall1.strand    == geneCall2.strand:
            return True
        else:
            return False 

    ##### SAME SAVE START COORDINATES #####
    # Determines whether two gene calls are identical except for the start coordinates
    ################
    def SameSaveStartCoordinates(self, geneCall1, geneCall2):
        if geneCall1.geneCaller == geneCall2.geneCaller and \
            geneCall1.contig    == geneCall2.contig     and \
            geneCall1.strand    == geneCall2.strand:
            if geneCall1.strand == '+':
                if geneCall1.leftEnd   != geneCall2.leftEnd and \
                    geneCall1.rightEnd == geneCall2.rightEnd: 
                    return True
            elif geneCall1.strang == '-':
                if geneCall1.leftEnd   == geneCall2.leftEnd and \
                    geneCall1.rightEnd != geneCall2.rightEnd:
                    return True
            else:
                return False 

    ##### PRINT METHODS ##########################################################################################

    def PrintGenecalls2file_gff(self,FILE_H,dataSet):
        seqid = '.'; source = '.'; type = 'cds'
        start = '0'; end = '0'; strand = '.' 
        phase = '.'; attributes = '.'; score = '.'
        FILE_H.write("%s\n" % ("##gff-version 3"))

        if dataSet.lower() == 'superset':     # all gene calls called by all callers
            for callList in self.uniqueList:
                seqid  = callList[0].contig
                source = callList[0].geneCaller
                if callList[0].strand == '+':
                    strand = callList[0].strand
                    start  = str(callList[0].leftEnd)
                    end    = str(callList[0].rightEnd)
                elif callList[0].strand == '-':
                    strand = callList[0].strand
                    start  = str(callList[0].rightEnd)
                    end    = str(callList[0].leftEnd)
                else:
                    if CGC_WARNINGS == 'True':
                        print("WARNING: CGC_compare says, Unexpected value for strand:", callList[0].strand)
                if len(callList) == len(self.callerList):      # all callers agreed
                    attributes = 'unanimous call'
                    source     = 'all callers'
                elif len(callList) == 1:                       # only 1 caller called this one
                    attributes = 'unique call'
                elif len(self.callerList) - len(callList) > 0: # multiple callers
                    source = ''
                    for geneCall in callList:
                        source += geneCall.geneCaller + '_'
                    source = source.rstrip('_')
                    attributes = 'multiple callers'
                else:
                    attributes = '.'
                FILE_H.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (seqid,source,type,start,end,strand,phase,attributes))
        elif dataSet.lower() == 'commoncore': # gene calls that were identical among all gene callers
            for callList in self.commonCore:
                seqid  = callList[0].contig
                source = callList[0].geneCaller
                if callList[0].strand == '+':
                    strand = callList[0].strand
                    start  = str(callList[0].leftEnd)
                    end    = str(callList[0].rightEnd)
                elif callList[0].strand == '-':
                    strand = callList[0].strand
                    start  = str(callList[0].rightEnd)
                    end    = str(callList[0].leftEnd)
                else:
                    if CGC_WARNINGS == 'True':
                        print("WARNING: CGC_compare says, Unexpected value for strand:", callList[0].strand)
                FILE_H.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (seqid,source,type,start,end,strand,phase,attributes))
        else:
            if CGC_WARNINGS == 'True':
                print("WARNING: CGC_compare says, dataSet not recognized:", dataSet) 
        return

    def PrintMergeList(self):
        print("\n***************Merge List")
        count = 1 
        for gene in self.mergeList:
            print(count, end=' ')
            gene.PrintAll_brief()
            count += 1
        return
  
    def PrintMergeList2file(self,FILE_H):
        FILE_H.write("%s\n" % ("\n***************Merge List"))
        count = 1 
        for gene in self.mergeList:
            FILE_H.write("%s\t" % (count))
            gene.PrintAll_brief_2file(FILE_H)
            count += 1
        return

    def PrintUniqueList(self):
        print("\n***************Unique List")
        count = 1 
        for list in self.uniqueList:
            print("Unique gene call", count)
            for gene in list:
                print("  --  ", end=' ') 
                gene.PrintAll_brief()
            count += 1
        return

    def PrintUniqueList2file(self,FILE_H):
        FILE_H.write("%s\n" % ("\n***************Unique List"))
        count = 1 
        for list in self.uniqueList:
            FILE_H.write("%s%s\n" % ("Unique gene call ", count))
            for gene in list:
                FILE_H.write("%s" % ("  --  ")) 
                gene.PrintAll_brief_2file(FILE_H)
            count += 1
        return

    def PrintCommonCore(self):
        print("\n***************Common Core")
        count = 1 
        for gene in self.commonCore:
            print(count, end=' ')
            gene.PrintAll_brief()
            count += 1
        return

    def PrintCommonCore2file(self,FILE_H):
        FILE_H.write("%s\n" % ("\n***************Common Core"))
        count = 1 
        for gene in self.commonCore:
            FILE_H.write("%s%s\t" % ("CC-gene ",count))
            gene.PrintAll_brief_2file(FILE_H)
            count += 1
        return

    def PrintCallerList(self):
        print("\n***************List of Gene Callers")
        for caller in self.callerList:
            print(caller) 

    def PrintCallerList2file(self,FILE_H):
        FILE_H.write("%s\n" % ("\n***************List of Gene Callers"))
        for caller in self.callerList:
            FILE_H.write("%s\t" % (caller)) 

    def PrintConsensusScores(self):
        print("\nGene-call Consensus Scores:")
        for caller in self.callerList:
            print("Caller", caller, "gene-call consensus score: %0.2f" % self.averageScores[caller])
        return

    def PrintConsensusScores2file(self,FILE_H):
        FILE_H.write("%s\n" % ("\nGene-call Consensus Scores:"))
        for caller in self.callerList:
            FILE_H.write("%s%s%s%0.2f\n" % ("Caller ",caller," average gene-call score: ",self.averageScores[caller]))
        return

    # Formats the unique calls list and prints to standard out; this is the final comparison data set
    def PrintGenecallGrid(self): # Prints an ordered, complete list of gene calls, each caller's in a column, identical calls in same row
        if self.callerList:
            if self.uniqueList: # Recall, uniqueList is list of unique gene calls, many of which were called by >1 caller
                count = 1

                # Print column headers, for as many gene callers as we have
                print("\nGene-call Table:")
                print("count\t", end=' ')
                for i in range(0,len(self.callerList)):
                    print("caller\tstrand\tleftEnd\trightEnd\tlength\tcontig\t", end=' ')
                print() 

                # Format each gene call as a single line of output, arranging gene callers in order left to right
                for geneList in self.uniqueList:

                    # Create an empty array for printing identical gene calls in order going across by gene caller
                    printArray = [] # initialize
                    for i in range(0,len(self.callerList)): # Make the array as big as it needs to be for current gene call 
                        printArray.append('')

                    # Fill the print array in order by gene caller # will be at least 1 gene caller's call 
                    for i in range(0,len(geneList)):
                        currentCaller = geneList[i].geneCaller
                        printColumn = self.callerList.index(currentCaller) # capture index of this gene caller in self.callerList
                        printArray[printColumn] = geneList[i].geneCaller + '\t' + geneList[i].strand   + '\t' \
                                                + geneList[i].leftEnd    + '\t' + geneList[i].rightEnd + '\t' \
                                                + geneList[i].geneLength + '\t' + geneList[i].contig + '\t'

                    # Print the current row: horizontal list of identical gene calls 
                    print(count, '\t', end=' ')
                    for geneCallString in printArray:
                        if geneCallString == '':
                            print("\t\t\t\t\t\t", end=' ')
                        else:
                            print(geneCallString, end=' ') 
                    print() 
                    count += 1
            else:
                print("PrintGenecallGrid(): uniqueList is empty")
        else:
            print("PrintGenecallGrid(): callerList is empty")
        return

    # Formats the unique calls list and prints to a file; this is the final comparison data set
    def PrintGenecallGrid2file(self,FILE_H): # Prints an ordered, complete list of gene calls, each caller's in a column, identical calls in same row
        if self.callerList:
            if self.uniqueList: # Recall, uniqueList is list of unique gene calls, many of which were called by >1 caller
                count = 1

                # Print column headers, for as many gene callers as we have
                FILE_H.write("%s\n" % ("Gene-call Table:"))
                FILE_H.write("%s" % ("count\t"))
                for i in range(0,len(self.callerList)):
                    FILE_H.write("%s\n" % ("caller\tstrand\tleftEnd\trightEnd\tlength\tcontig\t"))

                # Format each gene call as a single line of output, arranging gene callers in order left to right
                for geneList in self.uniqueList:

                    # Create an empty array for printing identical gene calls in order going across by gene caller
                    printArray = [] # initialize
                    for i in range(0,len(self.callerList)): # Make the array as big as it needs to be for current gene call 
                        printArray.append('')

                    # Fill the print array in order by gene caller # will be at least 1 gene caller's call 
                    for i in range(0,len(geneList)):
                        currentCaller = geneList[i].geneCaller
                        printColumn = self.callerList.index(currentCaller) # capture index of this gene caller in self.callerList
                        printArray[printColumn] = geneList[i].geneCaller + '\t' + geneList[i].strand   + '\t' \
                                                + geneList[i].leftEnd    + '\t' + geneList[i].rightEnd + '\t' \
                                                + geneList[i].geneLength + '\t' + geneList[i].contig + '\t'

                    # Print the current row: horizontal list of identical gene calls 
                    FILE_H.write("%s\t" % (count))
                    for geneCallString in printArray:
                        if geneCallString == '':
                            FILE_H.write("%s" % ("\t\t\t\t\t\t"))
                        else:
                            FILE_H.write("%s" % (geneCallString)) 
                    FILE_H.write("\n")
                    count += 1
            else:
                FILE_H.write("%s\n" % ("PrintGenecallGrid(): uniqueList is empty"))
        else:
            FILE_H.write("%s\n" % ("PrintGenecallGrid(): callerList is empty"))
        return

    def PrintReport(self):  # Final output
        self.PrintStats()
        self.PrintGenecallGrid()
        self.PrintGenecallScores()
        self.PrintConsensusScores()
        return

    def PrintReport2file(self,FILE_H):  # Final output
        self.PrintStats2file(FILE_H)
        self.PrintGenecallGrid2file(FILE_H)
        self.PrintGenecallScores2file(FILE_H)
        return

    def PrintGenecallScores(self):
        print()
        print("Gene-call Scores:")
        print('Caller\t', 'Contig\t', 'Gene No.\t', 'Left End\t', 'Right End\t', 'Strand\t', 'gcScore')
        for geneCallList in self.uniqueList:
            for geneCall in geneCallList:
                print(geneCall.geneCaller, '\t', geneCall.contig, '\t', geneCall.geneNumber, '\t', geneCall.leftEnd, '\t', geneCall.rightEnd, '\t', geneCall.strand, '\t', '%1.2f' % geneCall.score1)
        return

    def PrintGenecallScores2file(self,FILE_H):
        FILE_H.write("\n")
        FILE_H.write("%s\n" % ("Gene-call Scores:"))
        FILE_H.write("%s%s%s%s%s%s%s\n" % ('Caller\t','Contig\t','Gene No.\t','Left End\t','Right End\t','Strand\t','gcScore'))
        for geneCallList in self.uniqueList:
            for geneCall in geneCallList:
                FILE_H.write("%s\t%s\t%s\t%s\t%s\t%s\t%0.2f\n" % (geneCall.geneCaller,geneCall.contig,geneCall.geneNumber,geneCall.leftEnd,geneCall.rightEnd,geneCall.strand,geneCall.score1))
        return

    def PrintStats(self):

        # Print a list of the callers 
        print("The following gene callers were considered:", end=' ')
        for caller in self.callerList:
            print(',', caller, end=' ')
        print()
        print("The number of distinct gene calls over all gene callers is", len(self.uniqueList))
        print("The number of gene calls in common among all callers is", len(self.commonCore)) 

        # Calculate number of gene calls that are not shared between any 2 gene callers
        loneCallCount = 0
        for callerList in self.uniqueList:
            if len(callerList) == 1:
                loneCallCount += 1
        print("The number of unique (non-matching) gene calls is", loneCallCount) 

        # For each gene caller, calculate the number of calls it made 
        for caller in self.callerList:
            callCount        = 0
            cumulativeLength = 0 
            maxLength        = 0
            minLength        = 1000000 
            aveLength        = 0
            for call in self.mergeList:
                if (call.geneCaller == caller):
                    callCount += 1
                    intLength = int(call.geneLength)
                    cumulativeLength += intLength
                    if maxLength < intLength:
                        maxLength = intLength
                    if minLength > intLength:
                        minLength = intLength
            aveLength = cumulativeLength / callCount
            print("Caller", caller, "produced", callCount, "gene calls.")
            print("Caller", caller, "gene-call length stats:  min:", minLength, ", max:", maxLength, ", ave:", aveLength)

    def PrintStats2file(self,FILE_H):

        # Print a list of the callers 
        FILE_H.write("%s" % ("The following gene callers were considered: "))
        for caller in self.callerList:
            FILE_H.write("%s%s" % (caller, ', '))
        FILE_H.write("%s%s\n" % ("The number of distinct gene calls over all gene callers is ", len(self.uniqueList)))
        FILE_H.write("%s%s\n" % ("The number of gene calls in common among all callers is ", len(self.commonCore)))

        # Calculate number of gene calls that are not shared between any 2 gene callers
        loneCallCount = 0
        for callerList in self.uniqueList:
            if len(callerList) == 1:
                loneCallCount += 1
        FILE_H.write("%s%s\n" % ("The number of unique (non-matching) gene calls is", loneCallCount)) 

        # For each gene caller, calculate the number of calls it made 
        for caller in self.callerList:
            callCount        = 0
            cumulativeLength = 0 
            maxLength        = 0
            minLength        = 1000000 
            aveLength        = 0
            for call in self.mergeList:
                if (call.geneCaller == caller):
                    callCount += 1
                    intLength = int(call.geneLength)
                    cumulativeLength += intLength
                    if maxLength < intLength:
                        maxLength = intLength
                    if minLength > intLength:
                        minLength = intLength
            aveLength = cumulativeLength / callCount
            FILE_H.write("%s%s%s%s%s\n" % ("Caller ", caller, " produced ", callCount, " gene calls."))
            FILE_H.write("%s%s%s%s%s%s%s%s\n" % ("Caller ", caller, " gene-call length stats:  min: ", minLength, ", max: ", maxLength, ", ave: ", aveLength))

    def PrintAll(self):  # Print a dump of everything (debug/diagnostic) 
        self.PrintCallerList()
        self.PrintMergeList()
        self.PrintUniqueList()
        self.PrintCommonCore()
        self.PrintGenecallGrid()
        return

    def PrintAll2file(self,FILE_H):  # Print a dump of everything (debug/diagnostic) 
        self.PrintCallerList2file(FILE_H)
        self.PrintMergeList2file(FILE_H)
        self.PrintUniqueList2file(FILE_H)
        self.PrintCommonCore2file(FILE_H)
        self.PrintGenecallGrid2file(FILE_H)
        return

    def PrintAll_verbose(self):
        self.PrintCallerList()
        print("\n***************Merge List")
        for gene in self.mergeList:
            gene.PrintAll()
        print("\n***************Unique List")
        for list in self.uniqueList:
            for gene in list:
                gene.PrintAll()
        print("\n***************Common Core")
        for gene in self.commonCore:
            gene.PrintAll()
        print("\n***************GeneCall Grid")
        self.PrintGenecallGrid()
        return



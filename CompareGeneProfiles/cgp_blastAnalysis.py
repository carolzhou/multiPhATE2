##############################################
# Module: blastAnalysis.py
# Programmer:  Carol L. Ecale Zhou
# Date of last update:  16 May 2020
#
# Module comprising data structures and methods for blasting the genes and proteins
#    of two genome objects and comparing the gene profiles and the nt and aa levels.
#
# Classes and methods: 
#     hit 
#         printAll
#         printAll2file
#     hitList
#         append
#         printAll
#         printAll2file
#     homology
#         reportStats
#         printAll
#         printAll2file
#         mergeAll
#     paralog
#         printAll
#         printAll2file
#     blast
#         identifyLoners
#         identifyParalogs
#         compareHits
#         recordHits
#         printHits
#         printHits2file
#         makeBlastDB
#         performBlast
#

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.PDF FOR DETAILS.

import string
import cgp_fastaSequence as fastaSequence
import cgp_genomeSequence as genomeSequence
import cgp_annotation as annotation
import re, os, copy
import operator
from itertools import groupby   # used in mergeAll method # changed sort method 16july2013

# Boolean control: verbosity

PHATE_PROGRESS = False
PHATE_MESSAGES = False
PHATE_WARNINGS = False

PHATE_PROGRESS_STRING = os.environ["PHATE_PHATE_PROGRESS"] 
PHATE_MESSAGES_STRING = os.environ["PHATE_PHATE_MESSAGES"] 
PHATE_WARNINGS_STRING = os.environ["PHATE_PHATE_WARNINGS"] 

if PHATE_PROGRESS_STRING.lower() == 'true':
    PHATE_PROGRESS = True 
if PHATE_MESSAGES_STRING.lower() == 'true':
    PHATE_MESSAGES = True 
if PHATE_WARNINGS_STRING.lower() == 'true':
    PHATE_WARNINGS = True 

DEBUG = True
#DEBUG = False

# Default values used if user- or program-defined defaults not provided
GENE_MATCH_IDENTITY_DEFAULT = 95
GENE_MATCH_COVERAGE_DEFAULT = 60 
GENE_SIMILARITY_IDENTITY_DEFAULT = 60
GENE_SIMILARITY_COVERAGE_DEFAULT = 75
DOMAIN_MATCH_IDENTITY_DEFAULT = 95
DOMAIN_MATCH_COVERAGE_DEFAULT = 45
DOMAIN_SIMILARITY_IDENTITY_DEFAULT = 60
DOMAIN_SIMILARITY_COVERAGE_DEFAULT = 45
PARALOG_MATCH_IDENTITY_DEFAULT = 10             #
PARALOG_MATCH_COVERAGE_DEFAULT = 10             #
PARALOG_SIMILARITY_IDENTITY_DEFAULT = 60
PARALOG_SIMILARITY_COVERAGE_DEFAULT = 75
PARALOG_DOMAIN_MATCH_IDENTITY_DEFAULT = 95
PARALOG_DOMAIN_MATCH_COVERAGE_DEFAULT = 45
PARALOG_DOMAIN_SIMILARITY_IDENTITY_DEFAULT = 60
PARALOG_DOMAIN_SIMILARITY_COVERAGE_DEFAULT = 45
PROTEIN_MATCH_IDENTITY_DEFAULT = 95
PROTEIN_MATCH_COVERAGE_DEFAULT = 95
PROTEIN_SIMILARITY_IDENTITY_DEFAULT = 60
PROTEIN_SIMILARITY_COVERAGE_DEFAULT = 75
PROTEIN_DOMAIN_MATCH_IDENTITY_DEFAULT = 95
PROTEIN_DOMAIN_MATCH_COVERAGE_DEFAULT = 45
PROTEIN_DOMAIN_SIMILARITY_IDENTITY_DEFAULT = 60
PROTEIN_DOMAIN_SIMILARITY_COVERAGE_DEFAULT = 45
PROTEIN_PARALOG_MATCH_IDENTITY_DEFAULT = 10      #
PROTEIN_PARALOG_MATCH_COVERAGE_DEFAULT = 10      # 
PROTEIN_PARALOG_SIMILARITY_IDENTITY_DEFAULT = 60 
PROTEIN_PARALOG_SIMILARITY_COVERAGE_DEFAULT = 75
PROTEIN_PARALOG_DOMAIN_MATCH_IDENTITY_DEFAULT = 95
PROTEIN_PARALOG_DOMAIN_MATCH_COVERAGE_DEFAULT = 45
PROTEIN_PARALOG_DOMAIN_SIMILARITY_IDENTITY_DEFAULT = 60
PROTEIN_PARALOG_DOMAIN_SIMILARITY_COVERAGE_DEFAULT = 45

p_comment = re.compile('^#')

##############################################################################################
class hit(object):
   
    def __init__(self):
        self.queryHeader     = ""
        self.subjectHeader   = ""
        self.identity        = ""
        self.alignmentLength = 0
        self.mismatches      = None
        self.gapopens        = None
        self.queryStart      = 0
        self.queryEnd        = 0
        self.subjectStart    = 0
        self.subjectEnd      = 0
        self.evalue          = 0
        self.bitscore        = 0
        self.queryCoverage   = 0  # calculate later, when we get q/s lengths
        self.subjectCoverage = 0
        self.hitType         = "unknown"  # fill in with "mutual best", "singular best"
        self.queryLength     = 0
        self.subjectLength   = 0

    def computeCoverage(self,queryLength,subjectLength):
        errorCode = [] # for debug
        if queryLength == 0 or subjectLength == 0:
            errorCode.append(1) 
            return errorCode
        qLen = float(queryLength)
        sLen = float(subjectLength)
        qSpan = float(int(self.queryEnd) - int(self.queryStart) + 1)
        sSpan = float(int(self.subjectEnd) - int(self.subjectStart) + 1)
        self.queryCoverage = int(round(100*qSpan/qLen))
        self.subjectCoverage = int(round(100*sSpan/sLen))
        return (self.queryCoverage,self.subjectCoverage) 

    def printAll(self):
        print ("query:", self.queryHeader)
        print ("subject:", self.subjectHeader)
        print ("hit type:", self.hitType)
        print ("identity:", self.identity)
        print ("alignment length:", self.alignmentLength)
        print ("mismatches:", self.mismatches)
        print ("gapopens:", self.gapopens)
        print ("query start/end:", self.queryStart, "/", self.queryEnd)
        print ("subject start/end:", self.subjectStart, "/", self.subjectEnd)
        print ("evalue:", self.evalue)
        print ("bitscore:", self.bitscore)
        print ("query/subject coverage:", self.queryCoverage, "/", self.subjectCoverage)

    def printAll2file(self,FILE_HANDLE):
        FILE_HANDLE.write("%s%s%s" % ("query:",self.queryHeader,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("subject:",self.subjectHeader,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("hit type:",self.hitType,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("identity:",self.identity,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("alignment length:",self.alignmentLength,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("mismatches:",self.mismatches,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("gapopens:",self.gapopens,"\n"))
        FILE_HANDLE.write("%s%s%s%s%s" % ("query start/end:",self.queryStart,"/",self.queryEnd,"\n"))
        FILE_HANDLE.write("%s%s%s%s%s" % ("subject start/end:",self.subjectStart,"/",self.subjectEnd,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("evalue:",self.evalue,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("bitscore:",self.bitscore,"\n"))
        FILE_HANDLE.write("%s%s%s%s%s" % ("query/subject coverage:",self.queryCoverage,"/",self.subjectCoverage,"\n"))

    def printAll2file_tab(self,FILE_HANDLE):
        FILE_HANDLE.write("%s%s%s" % ("query:",self.queryHeader,"\t"))
        FILE_HANDLE.write("%s%s%s" % ("subject:",self.subjectHeader,"\t"))
        FILE_HANDLE.write("%s%s%s" % ("hit type:",self.hitType,"\t"))
        FILE_HANDLE.write("%s%s%s" % ("identity:",self.identity,"\t"))
        FILE_HANDLE.write("%s%s%s" % ("alignment length:",self.alignmentLength,"\t"))
        FILE_HANDLE.write("%s%s%s" % ("mismatches:",self.mismatches,"\t"))
        FILE_HANDLE.write("%s%s%s" % ("gapopens:",self.gapopens,"\t"))
        FILE_HANDLE.write("%s%s%s%s%s" % ("query start/end:",self.queryStart,"/",self.queryEnd,"\t"))
        FILE_HANDLE.write("%s%s%s%s%s" % ("subject start/end:",self.subjectStart,"/",self.subjectEnd,"\t"))
        FILE_HANDLE.write("%s%s%s" % ("evalue:",self.evalue,"\t"))
        FILE_HANDLE.write("%s%s%s" % ("bitscore:",self.bitscore,"\t"))
        FILE_HANDLE.write("%s%s%s%s%s" % ("query/subject coverage:",self.queryCoverage,"/",self.subjectCoverage,"\n"))

###########################################################################################################
class hitList(object): # Holds output from a blast run plus meta-data about the run

    def __init__(self):
        self.blastHits = []  # list of hit objects
        self.type = "unknown"  # "mutual best", "singular", "loner", "paralog"
        self.identityCutoff = "unknown" 
        self.evalueCutoff   = "unknown" 

    def append(self,newHit):
        self.blastHits.append(newHit)

    def printAll(self):
        for hit in self.blastHits:
            hit.printAll()

    def printAll2file(self,FILE_HANDLE):
        for hit in self.blastHits:
            hit.printAll2file(FILE_HANDLE)

###########################################################################################################
class homology(object):  # holds comparative information between 2 gene/protein sets
    
    def __init__(self):
        self.mutualBestHits = { #*** CHECK WHETHER THE VALUES SHOULD BE TYPE "None" (can Python accept user-defined type?) 
            "set1" : [],   # NOTE: using literal list of hit objects here; python doesn't like my user-defined type
            "set2" : []    # type list of hit objects 
        }
        self.singularHits = {
            "set1" : [],   # type list of hit objects 
            "set2" : []    # type list of hit objects 
        } 
        self.loners = {
            "set1" : [],   # type FastaSequence
            "set2" : []    # type FastaSequence
        }

    def computeCoverage(self,seqList1,seqList2):  # NEEDS TESTING 
        # This method fills in fields hit.queryCoverage and hit.subjectCoverage
        seqLength1 = {}
        seqLength2 = {}
        for seq in seqList1:
            seqLength1[seqList1.header] = len(seqList1.sequence)
        for seq in seqList2:
            seqLength2[seqList2.header] = len(seqList2.sequence)
        for hit in self.mutualBestHits["set1"]:
            hit.computeCoverage(seqLength1[hit.queryHeader],seqLength2[hit.subjectHeader])
        for hit in self.mutualBestHits["set2"]:
            hit.computeCoverage(seqLength2[hit.queryHeader],seqLength1[hit.subjectHeader])
        for hit in self.singularHits["set1"]:
            hit.computeCoverage(seqLength1[hit.queryHeader],seqLength2[hit.subjectHeader])
        for hit in self.singluarHits["set2"]:
            hit.computeCoverage(seqLength2[hit.queryHeader],seqLength1[hit.subjectHeader])

    def mergeAll(self,seqList1,seqList2):  
        # This method merges all of the hit types (plus loners) for both genomes.
        # First, a list is constructed 
        # comprising the reference genome's genes and their mutual best and singular
        # hits to the other genome's genes; this list includes genome 1's loners.
        # Then, genome 2's singular hits are added, with reference location based
        # on the reference genome--thus, according to the start positions of the
        # subject (gene in reference genome that was hit by genome 2's gene). Then
        # this combined list is sorted by genome 1's gene start positions.
        # Genome 2's loners are then merged in by fitting them in between the closest
        # (mutual or singular) hits (i.e., subject start) to reference genome genes.
        # Note that there may not be a single uniquely best location to place genome 
        # 2's loners, but a reasonably good position is found using this approach.

        errorCode      = [] # for debugging
        mergeList      = [] # array/list of hits and loners for gene/protein set1 (i.e., seqList1)
        hitLine        = [] # array if items to be appended to mergeList

        #########################################################################################
        # Fields for hitLine (0-19):
        # 0: sortPostion                               - used below to sort all hitLines relative to genome 1
        # 1: genome_hitType - "genome1"|"genome2" + "_" + "mutual"|"singular"|"loner"
        # 2-5: qStart,qEnd,sStart,sEnd                 - query/subject start/end values from blast output 
        # 6-7: identity,evalue                         - from blast output
        # 8: "query"|"subject"|"loner"                 - role in blast hit, or if no hit was found 
        # 9: g1header                                  - header of query sequence
        # 10: g1 query sequence contig (ie, parent sequence of gene)  #*** new
        # 11: g1annotations                             - annotations attached to query sequence
        # 12-14: g1Start_onGenome,g1End_onGenome,g1Strand - gene/protein start/end on its genome
        # 15: "query"|"subject"|"loner"
        # 16: g2header                                 - header of subject sequence
        # 17: g2 subject sequence contig (ie, parent sequence of gene) #*** new
        # 18: g2annotations                            - annotations attached to subject sequence
        # 19-21: g2Start_onGenome,g2sEnd_onGenome,g2Strand - subject gene/protein start/end on its genome
        # 22: G1coverage                               - genome1's gene/protein coverage in blast hit
        # 23: G2coverage                               - genome2's gene/protein coverage in blast hit
        # 24: alignmentLength                          - reported by blast
        # 25: gapopens                                 - reported by blast
        # 26: g1segment                                - span within blast alignment
        # 27: g1length                                 - length of g1 gene
        # 28: g2segment                                - span within blast alignment
        # 29: g2length                                 - length of g2 gene
        #########################################################################################

        # Fields for hitLine (0-32):
        SORT_POSITION      = 0
        GENOME_HIT_TYPE    = 1
        Q_START            = 2
        Q_END              = 3
        S_START            = 4
        S_END              = 5
        IDENTITY           = 6
        EVALUE             = 7
        G1_HIT_TYPE        = 8
        G1_HEADER          = 9
        G1_CONTIG          = 10
        G1_ANNOTATIONS     = 11
        G1_START_ON_GENOME = 12
        G1_END_ON_GENOME   = 13
        G1_STRAND          = 14
        G2_HIT_TYPE        = 15
        G2_HEADER          = 16
        G2_CONTIG          = 17
        G2_ANNOTATIONS     = 18
        G2_START_ON_GENOME = 19
        G2_END_ON_GENOME   = 20
        G2_STRAND          = 21 
        G1_COVERAGE        = 22
        G2_COVERAGE        = 23
        ALIGNMENT_LENGTH   = 24
        GAP_OPENS          = 25
        G1_SPAN            = 26
        G1_LENGTH          = 27
        G2_SPAN            = 28
        G2_LENGTH          = 29
        ALERT_LIST         = 30  # list of detected items of note and anomalies
        G1_SEQUENCE        = 31  # fill this in if possible alternate start site(s)
        G2_SEQUENCE        = 32  # ditto
 
        #########################################################################################

        ##### Set up data structures for quick access to annotations and sequence lengths
        #seqAnnot1  = {}   # key = header, value = annotations, for genome1 genes/proteins
        #seqAnnot2  = {}   # ditto, for genome2
        #seqLength1 = {}   # key = header, value = length of sequence for genome1 genes/proteins
        #seqLength2 = {}   # ditto, for genome2
        #seqContig1 = {}   # key = header, value = parent sequence (ie, contig name) for genome1 gene

        ##### Set up data structures for quick access to annotations and sequence lengths
        seqAnnot1  = {}   # key = header, value = annotations, for genome1 genes/proteins
        seqAnnot2  = {}   # ditto, for genome2
        seqLength1 = {}   # key = header, value = length of sequence for genome1 genes/proteins
        seqLength2 = {}   # ditto, for genome2
        seqContig1 = {}   # key = header, value = parent sequence (ie, contig name) for genome1 gene
        seqContig2 = {}   # ditto, for genome2

        annotations = []
        for seq in seqList1.fastaList:
            aList = list(seq.annotationList)
            for annot in aList:
                annotations.append(list(annot.annotationList))
            seqAnnot1[seq.header] = list(annotations) 
            annotations = []  # reset
            seqLength1[seq.header] = len(seq.sequence)
            seqContig1[seq.header] = seq.parentSequence
        for seq in seqList2.fastaList:
            aList = list(seq.annotationList)
            for annot in aList:
                annotations.append(list(annot.annotationList))
            seqAnnot2[seq.header] = list(annotations)
            annotations = []  # reset 
            seqLength2[seq.header] = len(seq.sequence)
            seqContig2[seq.header] = seq.parentSequence

        # Create lists to hold different categories of hits and loners 
        quickMutual1   = []
        quickSingular1 = []
        quickLoner1    = []
        quickSingular2 = []
        quickLoner2    = [] 

        # Capture mutual Best Hits for Genome 1   #*** CONTINUE HERE
        for hit in self.mutualBestHits["set1"]:  
            (qLocusTag,qStrand,qStart,qEnd,junk) = hit.queryHeader.split('/') #*** should split to list then pull data by index
            (sLocusTag,sStrand,sStart,sEnd,junk) = hit.subjectHeader.split('/')
            try:
                qAnnotations = seqAnnot1[hit.queryHeader]
                sAnnotations = seqAnnot2[hit.subjectHeader]
            except:
                print("cgp_blastAnalysis says, WARNING: error in retrieving qAnnotations or sAnnotations")
                print("  hit.queryHeader is",hit.queryHeader)
                print("  hit.subjectHeader is",hit.subjectHeader)
                continue 
            g1segment    = int(hit.queryEnd)   - int(hit.queryStart) + 1
            g2segment    = int(hit.subjectEnd) - int(hit.subjectStart) + 1
            g1length     = seqLength1[hit.queryHeader]
            g2length     = seqLength2[hit.subjectHeader]
            g1Contig     = seqContig1[hit.queryHeader]
            g2Contig     = seqContig2[hit.subjectHeader]
            (g1coverage,g2coverage) = hit.computeCoverage(g1length,g2length)
            sortPosition = qStart   # start position of query gene/protein, relative to genome1
            hitLine = [sortPosition,"genome1_mutual",hit.queryStart,hit.queryEnd,hit.subjectStart,hit.subjectEnd,hit.identity,hit.evalue,"query",hit.queryHeader,g1Contig,qAnnotations,qStart,qEnd,qStrand,"subject",hit.subjectHeader,g2Contig,sAnnotations,sStart,sEnd,sStrand,g1coverage,g2coverage,hit.alignmentLength,hit.gapopens,g1segment,g1length,g2segment,g2length]
            quickMutual1.append(hitLine)

        # Capture singular Best Hits for Genome 1
        for hit in self.singularHits["set1"]:
            (qLocusTag,qStrand,qStart,qEnd,junk) = hit.queryHeader.split('/')
            (sLocusTag,sStrand,sStart,sEnd,junk) = hit.subjectHeader.split('/')
            try:
                qAnnotations = seqAnnot1[hit.queryHeader]
                sAnnotations = seqAnnot2[hit.subjectHeader]
            except:
                if PHATE_WARNINGS:
                    print("cgp_blastAnalysis says, WARNING: error in retrieving qAnnotations or sAnnotations")
                    print("  hit.queryHeader is",hit.queryHeader)
                    print("  hit.subjectHeader is",hit.subjectHeader)
                continue 
            g1segment    = int(hit.queryEnd) - int(hit.queryStart) + 1
            g2segment    = int(hit.subjectEnd) - int(hit.subjectStart) + 1
            g1length     = seqLength1[hit.queryHeader]
            g2length     = seqLength2[hit.subjectHeader]
            g1Contig     = seqContig1[hit.queryHeader]
            g2Contig     = seqContig2[hit.subjectHeader]
            (g1coverage,g2coverage) = hit.computeCoverage(g1length,g2length)
            sortPosition = qStart     # start position of query gene/protein, relative to genome1
            hitLine = [sortPosition,"genome1_singular",hit.queryStart,hit.queryEnd,hit.subjectStart,hit.subjectEnd,hit.identity,hit.evalue,"query",hit.queryHeader,g1Contig,qAnnotations,qStart,qEnd,qStrand,"subject",hit.subjectHeader,g2Contig,sAnnotations,sStart,sEnd,sStrand,g1coverage,g2coverage,hit.alignmentLength,hit.gapopens,g1segment,g1length,g2segment,g2length]
            quickSingular1.append(hitLine)

        # Capture loners for Genome 1
        for seq in self.loners["set1"]:
            annotations = seqAnnot1[seq.header] 
            g1Contig    = seqContig1[seq.header]
            (locusTag,strand,start,end,junk) = seq.header.split('/')
            sortPosition = start             # no hit, actually; reference position is gene/protein start on genome 1
            hitLine = [sortPosition,"genome1_loner","","","","","","","loner",seq.header,g1Contig,annotations,start,end,strand,"","","","","","","","","","","","","","",""] 
            quickLoner1.append(hitLine)

        # Load hits & loners for genome 1 (i.e., reference genome's genes/proteins) into merge list 
        for hitLine in quickMutual1:
            mergeList.append(hitLine)
        for hitLine in quickSingular1:
            mergeList.append(hitLine)
        for loner in quickLoner1:
            mergeList.append(loner)

        # Prepare to insert genome2's singular hits and loners to mergeList
        # First, capture genome2's singular hits and loners in 'hitList' format 

        # Capture singular hits for genome2
        for hit in self.singularHits["set2"]:
            try:
                qAnnotations = seqAnnot1[hit.queryHeader]
                sAnnotations = seqAnnot2[hit.subjectHeader]
            except:
                if PHATE_WARNINGS:
                    print("cgp_blastAnalysis says, WARNING: error in retrieving qAnnotations or sAnnotations")
                    print("  hit.queryHeader is",hit.queryHeader)
                    print("  hit.subjectHeader is",hit.subjectHeader)
                continue 
            g1segment    = int(hit.subjectEnd) - int(hit.subjectStart) + 1
            g2segment    = int(hit.queryEnd)   - int(hit.queryStart) + 1
            g1length     = seqLength1[hit.subjectHeader]
            g2length     = seqLength2[hit.queryHeader]
            g1Contig     = seqContig1[hit.subjectHeader]
            g2Contig     = seqContig2[hit.queryHeader]
            (g2coverage,g1coverage) = hit.computeCoverage(g2length,g1length)
            (qLocusTag,qStrand,qStart,qEnd,junk) = hit.queryHeader.split('/')
            (sLocusTag,sStrand,sStart,sEnd,junk) = hit.subjectHeader.split('/')
            sortPosition = sStart  # relative to genome1 (the subject)
            hitLine = [sortPosition,"genome2_singular",hit.queryStart,hit.queryEnd,hit.subjectStart,hit.subjectEnd,hit.identity,hit.evalue,"subject",hit.subjectHeader,g1Contig,sAnnotations,sStart,sEnd,sStrand,"query",hit.queryHeader,g2Contig,qAnnotations,qStart,qEnd,qStrand,g1coverage,g2coverage,hit.alignmentLength,hit.gapopens,g1segment,g1length,g2segment,g2length]
            quickSingular2.append(hitLine) 

        # Capture loners for genome 2, but note that there is no hit, so no reference to genome 1
        for seq in self.loners["set2"]:
            annotations = seqAnnot2[seq.header]
            g2Contig    = seqContig2[seq.header]
            (locusTag,strand,start,end,junk) = seq.header.split('/')
            hitLine = ["0","genome2_loner","","","","","","","","","","","","","","loner",seq.header,g2Contig,annotations,start,end,strand,"","","","","","","",""]
            quickLoner2.append(hitLine) 

        # Append genome2's singular hits to mergeList and sort mergeList
        for hitLine in quickSingular2:
            mergeList.append(hitLine)
        g1contig_position = lambda x: (x[G1_CONTIG], int(x[SORT_POSITION]))
        mergeList.sort(key=g1contig_position)

        # Find a suitable location for each genome2 loner (place according to cds number along genome 2) #*** Try, but maybe better to lump at top
        for loner in quickLoner2:
            #lonerIndex = 0  # will place loner at front of mergeList unless a better location is found (see below)
            lonerIndex = len(mergeList)  # will place loner at end of mergeList unless a better location is found (see below)
            lonerContig     = loner[G2_CONTIG]
            lonerHeader     = loner[G2_HEADER]
            lonerHeaderList = lonerHeader.split('/')
            lonerCDSlist    = lonerHeaderList[0].split('cds') 
            lonerCDSnumber  = lonerCDSlist[1] 
            # Identify location for genome2's loner after line containing that genome's previous cds number
            #for hitLine in mergeList: 
            #    if lonerContig == hitLine[G2_CONTIG]: # look for location within same contig on genome 2
            #        if hitLine[G2_HEADER]:  #*** Check this
            #            g2Header     = hitLine[G2_HEADER]
            #            g2HeaderList = g2Header.split('/')
            #            g2CDSlist    = g2HeaderList[0].split('cds') 
            #            g2CDSnumber  = g2CDSlist[1] 
            #            # CDSs should be (pretty much) in order, based on having been sorted (above)
            #            if int(lonerCDSnumber) > int(g2CDSnumber):  
            #                # find index of this g2 data line in mergeList...
            #                lonerIndex = mergeList.index(hitLine) + 1   #  place loner just after 
            #                break  # found where to insert
            # Lastly, add genome2's loner to mergeList at position lonerIndex
            mergeList.insert(lonerIndex,loner)

        # Create a header for the mergeList; Make deep copy (new memory) of mergeList and return it
        headerLine = ["#sortPos","genome_type","qStart","qEnd","sStart","sEnd","identity","evalue","Q|S|Loner","G1header","G1contig","G1annotations","G1geneStart","G1geneEnd","G1strand","Q|S|Loner","G2header","G2contig","G2annotations","G2geneStart","G2geneEnd","G2strand","G1coverage","G2coverage","alignmentLength","gapOpens","g1span","g1length","g2span","g2length"]
        mergeList.insert(0,headerLine)
        newMergeList = copy.deepcopy(mergeList)  #*** ? Pretty sure dynamic memory is needed
        return newMergeList  

    def reportStats(self):
        statsList = []
        set1mutuals = len(self.mutualBestHits["set1"])
        set2mutuals = len(self.mutualBestHits["set2"])
        set1singulars = len(self.singularHits["set1"])
        set2singulars = len(self.singularHits["set2"])
        set1loners = len(self.loners["set1"])
        set2loners = len(self.loners["set2"])

        if PHATE_MESSAGES:
            print ("Set 1 mutual best hits:", set1mutuals)
            print ("Set 2 mutual best hits:", set2mutuals)
            print ("Set 1 singular hits:", set1singulars)
            print ("Set 2 singular hits:", set2singulars)
            print ("Set 1 loners:", set1loners)
            print ("Set 2 loners:", set2loners) 

        statsList.append("Set 1 mutual best hits: " + str(set1mutuals))
        statsList.append("Set 2 mutual best hits: " + str(set2mutuals))
        statsList.append("Set 1 singular hits: " + str(set1singulars))
        statsList.append("Set 2 singular hits:" + str(set2singulars))
        statsList.append("Set 1 loners: " + str(set1loners))
        statsList.append("Set 2 loners: " + str(set2loners))
        return statsList

    def printAll(self):
        print ("List of mutual best hits for set 1:")
        for hit in self.mutualBestHits["set1"]:   #*** Test if it works to say, "self.mutualBestHits["set1"].printAll()"
            hit.printAll()
        print ("End list of mutual best hits for set 1.")
        print ("List of mutual best hits for set 2:")
        for hit in self.mutualBestHits["set2"]:
            hit.printAll()
        print ("End of list of mutual best hits for set 2.")
        print ("List of singular best hits for set 1:")
        for hit in self.singularHits["set1"]:
            hit.printAll()
        print ("End of list of singular best hits for set 1.")
        print ("List of singular best hits for set 2:")
        for hit in self.singularHits["set2"]:
            hit.printAll()
        print ("End of list of singular best hits for set 2.")
        print ("List of loners for set 1:")
        for seq in self.loners["set1"]:
            seq.printAll()
        print ("End of list of loners for set 1.")
        print ("List of loners for set 2:")
        for seq in self.loners["set2"]:
            seq.printAll()
        print ("End of list of loners for set 2.")

    def printAll2file(self,FILE_HANDLE): 
        FILE_HANDLE.write("%s" % ("List of mutual best hits for set 1:\n"))
        for hit in self.mutualBestHits["set1"]:   
            hit.printAll2file(FILE_HANDLE)
        FILE_HANDLE.write("%s" % ("End list of mutual best hits for set 1.\n"))
        FILE_HANDLE.write("%s" % ("List of mutual best hits for set 2:\n"))
        for hit in self.mutualBestHits["set2"]:
            hit.printAll2file(FILE_HANDLE)
        FILE_HANDLE.write("%s" % ("End of list of mutual best hits for set 2.\n"))
        FILE_HANDLE.write("%s" % ("List of singular best hits for set 1:\n"))
        for hit in self.singularHits["set1"]:
            hit.printAll2file(FILE_HANDLE)
        FILE_HANDLE.write("%s" % ("End of list of singular best hits for set 1.\n"))
        FILE_HANDLE.write("%s" % ("List of singular best hits for set 2:\n"))
        for hit in self.singularHits["set2"]:
            hit.printAll2file(FILE_HANDLE)
        FILE_HANDLE.write("%s" % ("End of list of singular best hits for set 2.\n"))
        FILE_HANDLE.write("%s" % ("List of loners for set 1:\n"))
        for seq in self.loners["set1"]:
            seq.printAll2file(FILE_HANDLE)
        FILE_HANDLE.write("%s" % ("End of list of loners for set 1.\n"))
        FILE_HANDLE.write("%s" % ("List of loners for set 2:\n"))
        for seq in self.loners["set2"]:
            seq.printAll2file(FILE_HANDLE)
        FILE_HANDLE.write("%s" % ("End of list of loners for set 2.\n"))

###########################################################################################################
class paralog(object):  #
    def __init__(self):
        self.header = ""      # header of paralog's fasta sequence
        self.blastHit = None  # hit object that links this paralog
        self.coverage = 0.0

    def printAll(self):
        print ("Paralogs information:")
        print ("header:", self.header)
        print ("blast hit:")
        self.blastHit.printAll()
        print ("coverage:",self.coverage)

    def printAll2file(self,FILE_HANDLE):
        FILE_HANDLE.write("%s\n" % ("Paralogs information:"))
        FILE_HANDLE.write("%s%s%s\n" % ("header:",self.header))
        FILE_HANDLE.write("%s\n" % ("blast hit:"))
        self.blastHit.printAll2file(FILE_HANDLE)
        FILE_HANDLE.write("%s%s\n" % ("coverage:",self.coverage))
        
    def printAll2file_tab(self,FILE_HANDLE):
        FILE_HANDLE.write("%s%s\n" % ("header:",self.header))
        FILE_HANDLE.write("%s\n" % ("blast hit:"))
        self.blastHit.printAll2file_tab(FILE_HANDLE)
        FILE_HANDLE.write("%s%s\n" % ("coverage:",self.coverage))
        
###########################################################################################################
class blast(object):

    # Class blast creates blast databases and performs blast between 2 fasta sets and compares results

    def __init__(self):
        self.blastHits = hitList()  # List of hit objects w/associated data 
        self.blastHit  = hit()
        self.paralogT  = paralog()   # paralog template
        self.homologs  = homology()

    def identifyLoners(self,kvargs):  # kvargs contains 2 sequence lists and a homology object
        if "seqList1" in kvargs.keys():       # a multi-fasta object
            seqList1 = kvargs["seqList1"]
        else:
            return False
        if "seqList2" in kvargs.keys():       # a multi-fasta object
            seqList2 = kvargs["seqList2"]
        else:
            return False
        if "comparedHits" in kvargs.keys():   # a homology object
            comparedHits = kvargs["comparedHits"]
        else:
            return False

        mutualHitList1   = comparedHits.mutualBestHits["set1"]
        mutualHitList2   = comparedHits.mutualBestHits["set2"]
        singularHitList1 = comparedHits.singularHits["set1"]
        singularHitList2 = comparedHits.singularHits["set2"]
        lonerList1       = comparedHits.loners["set1"]
        lonerList2       = comparedHits.loners["set2"]

        hadahit = False
        for seq in seqList1.fastaList:
            for hit in mutualHitList1:  
                if seq.header == hit.queryHeader:
                    hadahit = True
            for hit in singularHitList1:
                if seq.header == hit.queryHeader:
                    hadahit = True 
            if not hadahit:
                lonerList1.append(seq)
            hadahit = False  # reset

        hadahit = False
        for seq in seqList2.fastaList:
            for hit in mutualHitList2:
                if seq.header == hit.queryHeader:
                    hadahit = True
            for hit in singularHitList2:
                if seq.header == hit.queryHeader:
                    hadahit = True
            if not hadahit:
                lonerList2.append(seq)
            hadahit = False  # reset
        return comparedHits   # returns a homology object
 
    def identifyParalogs(self,inList1,inList2,kvargs):
        paralogCount = 0
        if "type" in kvargs.keys():
            mType = kvargs["type"]
        if mType.lower() == "gene":
            if "paralogMatchIdentity" in kvargs.keys():
                identity = kvargs["paralogMatchIdentity"]
            else:
                identity = PARALOG_MATCH_IDENTITY_DEFAULT
            if "paralogMatchCoverage" in kvargs.keys():
                coverage = kvargs["paralogMatchCoverage"]
            else:
                coverage = PARALOG_MATCH_IDENTITY_DEFAULT
        elif mType.lower() == "protein":
            if "proteinParalogMatchIdentity" in kvargs.keys():
                identity = kvargs["proteinParalogMatchIdentity"]
            else:
                identity = PROTEIN_PARALOG_MATCH_IDENTITY_DEFAULT
            if "proteinParalogMatchCoverage" in kvargs.keys():
                coverage = kvargs["proteinParalogMatchCoverage"]
            else:
                coverage = PROTEIN_PARALOG_MATCH_IDENTITY_DEFAULT
        else:
            return False

        # For each sequence, record any qualifying hits to other sequences in the genome 
        for seq in inList1.fastaList:
            qLength = abs(int(seq.start) - int(seq.end))
            for nextHit in inList2.blastHits:  # check if it's a hit of seq against non-self seq
                qSpan = abs(int(nextHit.queryStart) - int(nextHit.queryEnd))
                try:
                    seqCoverage = 100 * (float(qSpan) / float(qLength))
                except:
                    seqCoverage = 0.0
                if seq.header == nextHit.queryHeader and nextHit.queryHeader != nextHit.subjectHeader:
                    if float(nextHit.identity) >= float(identity) and seqCoverage >= coverage: # check hit quality
                        newParalog = copy.deepcopy(self.paralogT) # replicate the paralog template
                        newParalog.header = nextHit.subjectHeader
                        newParalog.coverage = seqCoverage
                        newParalog.blastHit = nextHit
                        seq.paralogList.append(newParalog) # add to list of paralogs for this sequence object
                        paralogCount += 1 
        return paralogCount 

    def compareHits(self,inList1,inList2,kvargs): # Compare hits between 2 gene/protein sets
        newComparison = copy.deepcopy(self.homologs)
        errorCode = []  # for debug
        if "type" in kvargs.keys():
            mType = kvargs["type"]
        if mType.lower() == "gene":
            if "geneMatchIdentity" in kvargs.keys():
                identity = kvargs["geneMatchIdentity"]
            else:
                identity = GENE_MATCH_IDENTITY_DEFAULT
            if "geneMatchCoverage" in kvargs.keys():
                coverage = kvargs["geneMatchCoverage"]
            else:
                coverage = GENE_MATCH_COVERAGE_DEFAULT
        elif mType.lower() == "protein":
            if "proteinMatchIdentity" in kvargs.keys():
                identity = kvargs["proteinMatchIdentity"]
            else:
                identity = PROTEIN_MATCH_IDENTITY_DEFAULT
            if "proteinMatchCoverage" in kvargs.keys():
                coverage = kvargs["proteinMatchCoverage"]
            else:
                coverage = PROTEIN_MATCH_COVERAGE_DEFAULT
        else:
            errorCode.append(1)

        # Compare hits between 2 input hitList objects 
        for hit1 in inList1.blastHits:
            for hit2 in inList2.blastHits: 
                if hit1.subjectHeader == hit2.queryHeader and hit2.subjectHeader == hit1.queryHeader:
                    newComparison.mutualBestHits["set1"].append(hit1)
                    newComparison.mutualBestHits["set2"].append(hit2)
                    hit1.hitType = "mutual best"
                    hit2.hitType = "mutual best"
                elif hit1.subjectHeader == hit2.queryHeader:  # hit1's top (but not mutual)
                    newComparison.singularHits["set1"].append(hit1)
                    hit1.hitType = "singular best"
                elif hit2.subjectHeader == hit1.queryHeader:  # hit2's top (but not mutual)
                    newComparison.singularHits["set2"].append(hit2)
                    hit2.hitType = "singlular best"
        if errorCode:
            print ("cgp_blastAnalysis says, ERROR: compareHits() errorCode:", errorCode)
        return newComparison

    def recordHits(self,filename):
        HIT_FILE = open(filename,"r")
        newHitList = copy.deepcopy(self.blastHits) 
        fLines = HIT_FILE.read().splitlines()  # read lines into list, removing newlines
        for line in fLines:
            match = re.search(p_comment, line)
            if not match:
                fields = line.split('\t')
                newHit = copy.deepcopy(self.blastHit)
                newHit.queryHeader     = fields[0]
                newHit.subjectHeader   = fields[1]
                newHit.identity        = fields[2]
                newHit.alignmentLength = fields[3]
                newHit.mismatches      = fields[4]
                newHit.gapopens        = fields[5]
                newHit.queryStart      = fields[6]
                newHit.queryEnd        = fields[7]
                newHit.subjectStart    = fields[8]
                newHit.subjectEnd      = fields[9]
                newHit.evalue          = fields[10]
                newHit.bitscore        = fields[11]
                newHitList.append(newHit)
        HIT_FILE.close()
        return(newHitList)

    def printHits(self,hitList):
        hitList.printAll()

    def printHits2file(self,hitList,FILE_HANDLE):
        hitList.printAll2file(FILE_HANDLE)

    def makeBlastDB(self,kvargs):  # Create blastDBs for contigs, genes, proteins
        if "dbType" in kvargs.keys():
            databaseType = kvargs["dbType"].lower()
        else:
            databaseType = "nucl"
        if "filename" in kvargs.keys():
            filename = kvargs["filename"]
        else:
            return False
        if databaseType == "nucl" or databaseType == "nt" or databaseType == "dna":
            command = "makeblastdb -in " + filename + " -dbtype nucl -logfile " + filename + ".blastdb1_nucl_log"
        elif databaseType == "prot" or databaseType == "aa" or databaseType == "protein":
            command = "makeblastdb -in " + filename + " -dbtype prot -logfile " + filename + ".blastdb1_prot_log"
        else:
            return False
        myResult = os.system(command)
        return myResult

    def performBlast(self,kvargs):
        command = ""
        errorList = []

        # Gather parameters for blast
        if "query" in kvargs.keys():        # The sequences being blasted
            query = kvargs["query"]
        else:
            errorList.append(1) 
        if "subject" in kvargs.keys():      # The database that is being blasted against
            subject = kvargs["subject"]
        else:
            errorList.append(2)
        if "mtype" in kvargs.keys():
            mtype = kvargs["mtype"]
        else:
            errorList.append(3) 
        if "evalue" in kvargs.keys():
            evalue = str(kvargs["evalue"])
        else:
            evalue = str(0.01)              # default
        if "identity" in kvargs.keys():
            identity = str(kvargs["identity"])
        else:
            identity = str(50)              # default
        if "scoreEdge" in kvargs.keys():    # Best hit score edge
            scoreEdge = str(kvargs["scoreEdge"])
        else:
            scoreEdge = str(0.05)           # default
        if "maxTargetSeqs" in kvargs.keys():
            maxTargetSeqs = str(kvargs["maxTargetSeqs"])
        else:
            maxTargetSeqs = str(1)          # default
        if "overhang" in kvargs.keys():     # Blast hit overhang
            overhang = str(kvargs["overhang"])
        else:
            overhang = str(0.25)            # default
        if "outputFormat" in kvargs.keys():
            outputFormat = str(kvargs["outputFormat"])
        else:
            outputFormat = str(7)           # tabbed list by default
        if "outfile" in kvargs.keys():
            outfile = kvargs["outfile"]
        else:
            errorList.append(4)

        # Blast
        if mtype == "nucl" or mtype == "nt" or mtype == "gene" or mtype == "nucleotide":
            command = "blastn -query " + query + " -out " + outfile + " -task blastn -db " + subject + \
                " -evalue " + evalue + " -best_hit_score_edge " + scoreEdge + " -best_hit_overhang " + \
                overhang + " -outfmt " + outputFormat + " -perc_identity " + identity + \
                " -max_target_seqs " + maxTargetSeqs + " -max_hsps 1"
        elif mtype == "prot" or mtype == "aa" or mtype == "protein" or mtype == "peptide":
            command = "blastp -query " + query + " -out " + outfile + " -task blastp -db " + subject + \
                " -evalue " + evalue + " -best_hit_score_edge " + scoreEdge + " -best_hit_overhang " + \
                overhang + " -outfmt " + outputFormat + " -max_target_seqs " + maxTargetSeqs + \
                " -max_hsps 1" 
        else:
            errorList.append(10)

        result = os.system(command)
        errorList.append(result)
 
        return errorList 
  

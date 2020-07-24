##########################################################################
#
# module dbPrep_vog
#
# Description:  Handles data and processing for VOG data type
#
# Last update:  23 July 2020
#
# Programmer:  C. E. Zhou
#
#########################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.PDF FOR DETAILS.

import os
import re
import copy
import subprocess
import timeit
import datetime
import dbPrep_fastaSequence

#DEBUG = True
DEBUG = False 

VERBOSE = False
if "dbPrep_VERBOSE" in os.environ.keys():
    if os.environ["dbPrep_VERBOSE"] == 'True':
        VERBOSE = True

DO_GENE    = True
DO_PROTEIN = True

fastaObj = dbPrep_fastaSequence.fasta()

class VOG(object):
    
    def __init__(self):
        self.VOGid           = ''  # e.g., VOG0334
        self.VOGannotation   = ''
        self.accnList        = []  # list of accession numbers 
        self.geneCount       = 0
        self.peptideCount    = 0   # ??? are these different?

    def printAll(self):
        print("VOG identifier:",self.VOGid)
        print("VOG annotation:",self.VOGannotation)
        print("accnList:",self.accnList)
        print("geneCount:",self.geneCount)
        print("peptideCount:",self.peptideCount)

class VOGs(object):

    def __init__(self):
        self.databaseName           = "VOGs"
        self.downloadDate           = "July 2020"  # default; master script should set this
        self.version                = "vog99"      # default; master script should set this
        self.VOGmapFile             = ""           # vog.members.tsv
        self.VOGannotationFile      = ""           # vog.annotations.tsv
        self.VOGgeneFastaFile       = ""           # vog.genes.all.fa
        self.VOGproteinFastaFile    = ""           # vog.proteins.all.fa
        self.VOGlist                = []           # list of VOG objects 
        self.VOGobj                 = VOG()
        self.fastaSet               = dbPrep_fastaSequence.multiFasta()  # using this class for headers only (no sequence)
        self.VOGgeneFastaOutFile    = ""
        self.VOGproteinFastaOutFile = ""

    def tagVogFastas(self,kvargs):
        print("dbPrep_vog says, Reading input parameters.")
        self.readParameters(kvargs)
 
        # Create VOG class objects and load with VOG identifiers
        print("dbPrep_vog says, Loading VOG map data and creating VOG objects, from file",self.VOGmapFile)
        vogMap_h = open(self.VOGmapFile,"r")
        self.loadVOGs(vogMap_h)
        vogMap_h.close()

        # Add annotations to VOG objects
        print("dbPrep_vog says, Adding annotations to VOG objects, from file", self.VOGannotationFile)
        vogAnnot_h = open(self.VOGannotationFile)
        self.loadAnnotations(vogAnnot_h)
        vogAnnot_h.close()

        # Insert sequence headers into fastaSet
        print("dbPrep_vog says, Inserting accessions.")
        self.insertAccns() 

        # Link fasta objects to VOG identifiers; Add VOG identifiers to fasta headers
        print("dbPrep_vog says, Linking fastas to VOGs. Please be patient. This may take a while...")
        startTime = datetime.datetime.now()
        self.linkFastas2VOGs()
        endTime = datetime.datetime.now()
        executionTime = endTime - startTime
        print("dbPrep_vog says, Whew! That took a long time: ",executionTime)

        # Identify sequences from files, modify headers, and write to new files.

        # Gene sequences
        if DO_GENE:
            print("dbPrep_vog says, Adding gene sequences to fasta objects, from file", self.VOGgeneFastaFile)
            print("dbPrep_vog says, Again, this may take a while...")
            vogInFile_h  = open(self.VOGgeneFastaFile,"r")
            vogOutFile_h = open(self.VOGgeneFastaOutFile,"w")
            startTime = datetime.datetime.now()
            self.writeVOGtaggedFastaFile(vogInFile_h,vogOutFile_h,"gene") 
            endTime = datetime.datetime.now()
            executionTime = endTime - startTime
            vogInFile_h.close()
            vogOutFile_h.close()
            print("dbPrep_vog says, Finally...after this much time: ",executionTime)

        # Protein sequences
        if DO_PROTEIN:
            print("dbPrep_vog says, Adding protein sequences to fasta objects, from file", self.VOGproteinFastaFile)
            print("dbPrep_vog says, Yeah, this one is really slow too...")
            vogInFile_h  = open(self.VOGproteinFastaFile,"r")
            print("Opening self.VOGproteinFastaOutFile for write:",self.VOGproteinFastaOutFile)
            vogOutFile_h = open(self.VOGproteinFastaOutFile,"w")
            startTime = datetime.datetime.now()
            self.writeVOGtaggedFastaFile(vogInFile_h,vogOutFile_h,"protein")
            endTime = datetime.datetime.now()
            executionTime = endTime - startTime
            vogInFile_h.close()
            vogOutFile_h.close()
            print("dbPrep_vog says, That took a while too: ",executionTime)

        return

    def linkFastas2VOGs(self):
        for vog in self.VOGlist:
            for accn in vog.accnList:
                for fasta in self.fastaSet.fastaList:
                    if fasta.header == accn:
                        vogString = vog.VOGid + '|'
                        fasta.customHeader += vogString
                        if VERBOSE:
                            print("dbPrep_vog says, Tagging ",fasta.customHeader)
        return

    def insertAccns(self):
        # Create a non-redundant list of accessions (also fasta headers)
        accnList = []
        for vog in self.VOGlist:
            for accn in vog.accnList:
                if accn not in accnList:
                    accnList.append(accn)
                    if VERBOSE:
                        print("dbPrep_vog says, Appending accession ",accn)
        # For each accession, create a fasta object and add it to the fasta set
        for accn in accnList:
            newFasta = copy.deepcopy(fastaObj)
            newFasta.header = accn
            self.fastaSet.fastaList.append(newFasta)
        if VERBOSE:
            print("dbPrep_vog says, Fasta list has been created.")
        return

    def readParameters(self,kvargs):
        if "databaseName" in kvargs.keys():
            self.databaseName = kvargs["databaseName"] 
        if "downloadDate" in kvargs.keys():
            self.downloadData = kvargs["downloadDate"] 
        if "version" in kvargs.keys():
            self.version = kvargs["version"] 
        if "VOGmapFile" in kvargs.keys():
            self.VOGmapFile = kvargs["VOGmapFile"] 
        if "VOGannotationFile" in kvargs.keys():
            self.VOGannotationFile = kvargs["VOGannotationFile"] 
        if "VOGgeneFastaFile" in kvargs.keys():
            self.VOGgeneFastaFile = kvargs["VOGgeneFastaFile"] 
        if "VOGproteinFastaFile" in kvargs.keys():
            self.VOGproteinFastaFile = kvargs["VOGproteinFastaFile"] 
        if "VOGgeneFastaOutFile" in kvargs.keys():
            self.VOGgeneFastaOutFile = kvargs["VOGgeneFastaOutFile"] 
        if "VOGproteinFastaOutFile" in kvargs.keys():
            self.VOGproteinFastaOutFile = kvargs["VOGproteinFastaOutFile"] 
        return

    def loadVOGs(self,vogMapFile_h):
        if VERBOSE:
            print("dbPrep_vog says, Loading VOGs.")
        vogCount = 0
        p_comment = re.compile("^#")
        VOGID_COL = 0; ACCN_COL = 4
        vogID = ""; accnListString = ""; accnList = []
        accnCount = 0; totalCount = 0
        fLines = vogMapFile_h.read().splitlines()
        for fLine in fLines:
            match_comment = re.search(p_comment,fLine)
            if not match_comment:
                vogCount += 1
                if VERBOSE:
                    if vogCount >= 1000:
                        print("dbPrep_vog says, Still working...") 
                        vogCount = 0
                fields = fLine.split('\t')
                vogID  = fields[VOGID_COL]
                accnListString = fields[ACCN_COL]
                newVOGobj = copy.deepcopy(self.VOGobj)
                newVOGobj.VOGid = vogID
                accnList = accnListString.split(',')
                for accn in accnList:
                    accnCount += 1
                    totalCount += 1
                    newVOGobj.accnList.append(accn)
                self.VOGlist.append(newVOGobj)
                accnCount = 0
        return

    def loadAnnotations(self,annotFile_h):
        p_comment = re.compile("^#")
        VOGID_COL = 0; ANNOT_COL = 4
        vogID = ""; annot = ""
        annotHash = {}  # key=vogID, value=annotation
        annotCount = 0

        # First, create a hash for easy annotation lookup
        aLines = annotFile_h.read().splitlines()
        for aLine in aLines:
            match_comment = re.search(p_comment,aLine)
            if not match_comment:
                annotCount += 1
                fields = aLine.split('\t')
                vogID = fields[VOGID_COL]
                annot = fields[ANNOT_COL]
                annotHash[vogID] = annot

        # Next, walk through the self.vogList, adding in the annotations
        annotCount = 0
        print("dbPrep_vog says, Walking through self.VOGlist. There are this many VOGs:",len(self.VOGlist))
        for vogObj in self.VOGlist:
            annot = annotHash[vogObj.VOGid]
            vogObj.VOGannotation = annot
            if annot != "":
                annotCount += 1
        print("dbPrep_vog says, In all there are this many VOGs:",len(self.VOGlist))
        print("dbPrep_vog says, This many VOGs were found to have annotation:",annotCount)
        return

    def writeVOGtaggedFastaFile(self,vogInFile_h,vogOutFile_h,seqType):  # seqType is 'gene' or 'protein'
        headerCount    = 0
        header   = ""
        sequence = ""
        FIRST    = True
        p_header = re.compile('^>')
        iLines = vogInFile_h.read().splitlines()
        for iLine in iLines:
            match_header   = re.search(p_header,iLine)
            if match_header:
                (junk,matchHeader) = iLine.split('>')
                if not FIRST: # Process previous fasta sequence
                    for fasta in self.fastaSet.fastaList:
                        if fasta.header == matchHeader:
                            newHeader = '>' + fasta.customHeader + fasta.header
                            headerCount += 1
                            if VERBOSE:
                                if headerCount >= 1000:
                                    print ("working...")
                                    headerCount = 0
                            vogOutFile_h.write("%s\n%s\n" % (newHeader,sequence))
                # Reset
                header = ""
                sequence = ""
                FIRST = False
            else:
                sequence += iLine

        # Process last sequence
        for fasta in self.fastaSet.fastaList:
            if fasta.header == header:
                newHeader = fasta.customHeader + fasta.header
                vogOutFile_h.write("%s\n%s\n" % (newHeader,sequence))
        return

    def getSequence(self,accn,seqFile_h):
        sequence = ""
        return sequence

    def getVogCount(self):
        return len(self.VOGlist)

    def getAccessionCount(self):
        if VERBOSE:
            print("dbPrep_vog says, Counting accessions")
        accessionCount = 0
        for VOG in self.VOGlist:
            accessionCount += len(VOG.accessionList)
        return accessionCount

    def printAll(self):
        print("VOGs object content:")
        print("databaseName:",self.databaseName)
        print("version:",self.version)
        print("VOGmapFile:",self.VOGmapFile)
        print("VOGannotationFile:",self.VOGannotationFile)
        print("VOGgeneFastaFile:",self.VOGgeneFastaFile)
        print("VOGproteinFastaFile:",self.VOGproteinFastaFile)
        for vogObj in self.VOGlist:
            vogObj.printAll()
        for fasta in self.fastaSet:
            fasta.printAll()

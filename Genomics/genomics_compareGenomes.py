#################################################################################
# Module: genomics_compareGenomes.py
#
# Programmer:  Carol L. Ecale Zhou
#
# Most recent update: 01 May 2020
#
# Module comprising classes and data structures for comparing genomes
#
# Classes and Methods:
#    genomes
#    genome
#    paralogSet
#    gene
#    protein
#
#################################################################################

# This code was developed by Carol. L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.PDF FOR DETAILS.

import re, os, copy

CODE_BASE_DIR    = ""
OUTPUT_DIR       = ""

# Verbosity

PHATE_WARNINGS = False
PHATE_MESSAGES = False
PHATE_PROGRESS = False
"""
PHATE_WARNINGS_STRING = os.environ["PHATE_PHATE_WARNINGS"]
PHATE_MESSAGES_STRING = os.environ["PHATE_PHATE_MESSAGES"]
PHATE_PROGRESS_STRING = os.environ["PHATE_PHATE_PROGRESS"]

if PHATE_WARNINGS_STRING.lower() == 'true':
    PHATE_WARNINGS = True

if PHATE_MESSAGES_STRING.lower() == 'true':
    PHATE_MESSAGES = True

if PHATE_PROGRESS_STRING.lower() == 'true':
    PHATE_PROGRESS = True
"""
DEBUG = True
#DEBUG = False

# Templates

# Class genomes represents the set of genomes that have been compared using CGP
class genomes(object):
    
    def __init__(self):
        self.name                 = ""     # Can assign a set to this genome comparison
        self.directories          = []     # List of Results directories to be processed
        self.referenceGenome      = ""     # name of genome assigned as the reference
        self.genomeCount          = 0      # number of genomes in set
        self.genomeList           = []     # set of genome class objects
        self.commonCore           = []     # list of genes that are common among all genomes; how to define exactly???

    def loadData(self):
        if PHATE_PROGRESS:
            print("genomics_compareGenomes says, Loading data from CGP reports.")
        self.readDirectories()
        self.parseDirectories()
        if PHATE_PROGRESS:
            print("genomics_compareGenomes says, Data loading complete.")
        return

    def readDirectories(self):
        dirList = []
        #*** CONTINUE HERE
        if PHATE_MESSAGES:
            print("genomics_compareGenomes says, The following directories have been read:")
        return dirList

    def parseDirectories(self,dirList):
        for directory in dirList:
            if PHATE_PROGRESS:
                print("genomics_compareGenomes says, Loading data from directory,",directory)

            # Query current (home) directory
            command = "pwd"
            homeDir = os.system(command)

            # Change to working directory
            command = "cd " + directory
            os.system(command)

            # Extract data from report file
            #*** CONTINUE HERE

            # Return to starting directory
            command = "cd " + homeDir
            os.system(command)

        return

    # For development purposes
    def checkUnique(self):
        tempList = []
        for genome in genomeList:
            if genome.name in tempList:
                print("genomics_compareGenomes genomes obj says, WARNING: ",genome.name," occurs more than once in genomes object.")
            else:
                tempList.append(genome.name)
        return

    def countGenomes(self):
        self.genomeCount = len(self.genomeList)
        return

    def computeCoreGenome(self):
        return

    def printReport(self):
        return

    def printAll(self):
        print("=========Genome Set=============")
        print("name:",self.name)
        print("referenceGenome:",self.referenceGenome)
        print("genomeCount:",self.genomeCount)
        print("Set of genomes:")
        for genome in self.genomeList:
            genome.printAll()
        print("=========End Genome Set=========")
        return

# Class genome stores metadata about a given genome (no sequences here)
class genome(object):

    def __init__(self):
        self.name                 = ""     # Name of this genome (e.g., Lambda)
        self.species              = ""     # Species name for this genome (e.g., Ecoli_Lambda_phage)
        self.isReference          = False  # True if designated a reference genome
        self.contigList           = []     # Set of contig names (fasta headers)
        self.geneList             = []     # List of gene objects 
        self.proteinList          = []     # List of protein objects 
        self.paralogList          = []     # List of paralogSet objects

    # For development purposes
    def checkUnique(self):
        tempList = []
        for gene in self.geneList:
            if gene.identifier in tempList:
                print("genomics_compareGenomes genome obj says, WARNING: ",gene.identifier," is not unique")
            else:
                tempList.append(gene.identifier)
            
        tempList = []
        for gene in self.geneList:
            if gene.cgpHeader in tempList:
                print("genomics_compareGenomes genome obj says, WARNING: ",gene.cgpHeader," is not unique")
            else:
                tempList.append(gene.cgpHeader)
            
        tempList = []
        for protein in self.proteinList:
            if protein.identifier in tempList:
                print("genomics_compareGenomes genome obj says, WARNING: ",protein.identifier," is not unique")
            else:
                tempList.append(protein.identifier)

        tempList = []
        for protein in self.proteinList:
            if protein.cgpHeader in tempList:
                print("genomics_compareGenomes genome obj says, WARNING: ",protein.cgpHeader," is not unique")
            else:
                tempList.append(protein.cgpHeader)
        return

    def identifyParalogs(self):
        if PHATE_PROGRESS:
            print("genomics_compareGenomes says, Identifying paralogs for genome",self.name)
        return

    def printReport(self):
        return

    def printAll(self):
        print("=========Genome========")
        print("name:",self.name)
        print("species:",self.species)
        print("isReference:",self.isReference)
        print("contigs:",)
        for contig in contigList:
            print(contig)
        for gene in self.geneList:
            gene.printAll()
        for protein in self.proteinList:
            protein.printAll()
        for paralogSet in self.paralogList:
            paralogSet.printAll()
        print("=========End Genome====")
        return

# Class paralog organizes genes and proteins within a given genome that are considered paralogs.
class paralogSet(object):
    
    def __init__(self):
        self.paralogType          = "unknown" # "gene" or "protein"
        self.paralogList          = []        # List of either gene or protein unique identifiers
        self.setSize              = 0         # number of paralogs in this set

    def countParalogs(self):
        self.setSize = len(self.paralogList)
        return

    def printReport(self):
        return

    def printAll(self):
        print("paralogType:",self.paralogType)
        print("setSize:",self.setSize)
        print("paralogList:",self.paralogList)
        return

# Class gene stores meta-data about a gene
class gene(object):

    def __init__(self):
        self.name                 = ""        # Ex: "phanotate_5"
        self.identifier           = ""        # Unique: genome + contig + name
        self.cgpHeader            = ""        # CGP-assigned header; Ex: cgp5/+/72/485/
        self.number               = 0         # Unique: Set to correspond to protein 
        self.parentGenome         = ""        # Ex: "Lambda" - corresponds to a genome object's name
        self.contigName           = ""        # Ex: "Lambda_contig_1"; read from CGP report
        self.annotation           = ""        # The full annotation string, read from CGP report 
        self.isLoner              = False     # True if this gene has no correspondences
        self.mutualBestHitList    = []        # list of gene identifiers that are mutual best hits
        self.singularBestHitList  = []        # list of gene identifiers that are best hits, relative to this gene
        self.correspondenceList   = []        # List of corresponding genes (mutual or best-hit): list of gene identifiers 
        self.paralogList          = []        # List of paralogs within its own parent genome: list of gene identifiers

    def addMutualBestHit(self,hit):
        self.mutualBestHistList.append(hit)
        return

    def addSingularBestHit(self,hit):
        self.singularBestHitList.append(hit)
        return

    # For development purposes
    def verifyLoner(self):
        if self.mutualBestHitList:
            return False
        if self.singularBestHistList:
            return False
        return True

    def printReport(self):
        return

    # For development purposes
    def printAll(self):
        print("===name:",self.name)
        print("identifier:",self.identifier)
        print("cgpHeader:",self.cgpHeader)
        print("number:",self.number)
        print("parentGenome:",self.parentGenome)
        print("contigName:",self.contigName)
        print("annotation:",self.annotation)
        print("isLoner:",self.isLoner)
        print("mutualBestHitList:",self.mutualList)
        print("singularBestHitList:",self.bestHitList)
        print("correspondenceList:",self.correpondenceList)
        print("paralogList:",self.paralogList)
        return

# Class protein stores meta-data about a protein 
class protein(object):

    def __init__(self):
        self.name                 = ""        # Ex: "phanotate_5"
        self.identifier           = ""        # Unique: genome + contig + name
        self.cgpHeader            = ""        # CGP-assigned header; Ex: cgp5/+/72/485/
        self.number               = 0         # Unique: Set to correspond to gene
        self.parentgenome         = ""        # Ex: "Lambda" - corresponds to a genome object's name
        self.contigName           = ""        # Ex: "Lambda_contig_1"; read from CGP report
        self.annotation           = ""        # The full annotation string, read from CGP report
        self.isLoner              = False     #
        self.mutualBestHitList    = []        # list of protein identifiers that are mutual best hits
        self.singularBestHitList  = []        # list of protein identifiers that are best hits, relative to this protein 
        self.correspondenceList   = []        # List of corresponding genes (mutual or best-hit): list of protein identifiers 
        self.paralogList          = []        # List of paralogs within its own parent genome: list of protein identifiers

    def addMutualBestHit(self,hit):
        self.mutualBestHitList.append(hit)
        return

    def addSingularBestHit(self,hit):
        self.singularBestHistList.append(hit)
        return

    # For development purposes
    def verifyLoner(self):
        if self.mutualBestHitList:
            return False
        if self.singularBestHistList:
            return False
        return True

    def printReport(self):
        return

    # For development purposes
    def printAll(self):
        print("===name:",self.name)
        print("identifier:",self.identifier)
        print("cgpHeader:",self.cgpHeader)
        print("number:",self.number)
        print("parentGenome:",self.parentGenome)
        print("contigName:",self.contigName)
        print("annotation:",self.annotation)
        print("isLoner:",self.isLoner)
        print("mutualList:",self.mutualList)
        print("bestHitList:",self.bestHitList)
        print("correspondenceList:",self.correpondenceList)
        print("paralogList:",self.paralogList)
        return


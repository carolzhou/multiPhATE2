#################################################################################
# Module: genomics_compareGenomes.py
#
# Programmer:  Carol L. Ecale Zhou
#
# Most recent update: 06 May 2020
#
# Module comprising classes and data structures for comparing genomes
#
# Classes and Methods:
#    class comparison
#    class genomes
#       performComparison()
#       loadData()
#       readDirectories()
#       parseDirectories()
#       checkUnique()
#       identifyParalogs()
#       computeCoreGenome()
#       computeHomologyGroupes()
#       saveHomologyGroups()
#       countGenomes()
#       checkUnique()
#       printReport()
#       printAll()
#    class genome
#       checkUnique()
#       identifyParalogs()
#       printReport()
#       printAll()
#    class paralogSet
#       countParalogs()
#       printReport()
#       printAll()
#    class gene_protein
#       addMutualBestHit()
#       addSingularBestHit()
#       addGroupMember()
#       constructGroup()
#       verifyLoner()
#       printReport()
#       printAll()
#
#################################################################################

# This code was developed by Carol. L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.PDF FOR DETAILS.

import re, os, copy
import subprocess

CODE_BASE_DIR    = ""
OUTPUT_DIR       = ""

# Verbosity

PHATE_WARNINGS = False
PHATE_MESSAGES = False
PHATE_PROGRESS = False

PHATE_WARNINGS_STRING = os.environ["PHATE_PHATE_WARNINGS"]
PHATE_MESSAGES_STRING = os.environ["PHATE_PHATE_MESSAGES"]
PHATE_PROGRESS_STRING = os.environ["PHATE_PHATE_PROGRESS"]

if PHATE_WARNINGS_STRING.lower() == 'true':
    PHATE_WARNINGS = True

if PHATE_MESSAGES_STRING.lower() == 'true':
    PHATE_MESSAGES = True

if PHATE_PROGRESS_STRING.lower() == 'true':
    PHATE_PROGRESS = True

# Override during development
PHATE_PROGRESS = True
PHATE_WARNINGS = True
PHATE_MESSAGES = True

DEBUG = True
#DEBUG = False

# Locations and files
CGP_RESULTS_DIR      = os.environ["PHATE_CGP_RESULTS_DIR"]
GENOMICS_RESULTS_DIR = os.environ["PHATE_GENOMICS_RESULTS_DIR"]
CGP_LOG_FILE         = "compareGeneProfiles_main.log"           # Lists genome fasta and annotation files
CGP_REPORT_FILE      = "compareGeneProfiles_main.report"        # Tabbed file with gene-gene listing
CGP_OUT_FILE         = "compareGeneProfiles_main.out"           # Same data as report file, but in python list format
CGP_SUMMARY_FILE     = "compareGeneProfiles_main.summary"       # High-level summary of genome-genome comparison

# Class comparison organizes all genomic data and performs comparisons, ultimately yielding
#   homology groups for each gene/protein in a reference genome.
class comparison(object):
    
    def __init__(self):
        self.name                 = "multiPhATE2"     # Can assign a set to this genome comparison
        self.directories          = []                # List of Results directories to be processed
        self.referenceGenome      = ""                # name of genome assigned as the reference
        self.genomeCount          = 0                 # number of genomes in set
        self.genomeList           = []                # set of genome class objects
        self.commonCore           = []                # list of genes that are common among all genomes; how to define exactly???
        # Templates
        self.genomeTemplate       = genome()
        self.geneTemplate         = gene_protein()
        self.geneTemplate.type    = "gene"
        self.proteinTemplate      = gene_protein()
        self.proteinTemplate.type = "protein"

    def performComparison(self):
        if PHATE_PROGRESS:
            print("genomics_compareGenomes says, Loading data.")
        self.loadData()
        if PHATE_PROGRESS:
            print("genomics_compareGenomes says, Identifying paralogs.")
        self.identifyParalogs()
        if PHATE_PROGRESS:
            print("genomics_compareGenomes says, Computing core genome.")
        self.computeCoreGenome()
        if PHATE_PROGRESS:
            print("genomics_compareGenomes says, Computing homology groups.")
        self.computeHomologyGroups()
        if PHATE_PROGRESS:
            print("genomics_compareGenomes says, Writing homology groups to files.")
        self.saveHomologyGroups()
        self.printReport()

    def loadData(self):
        dirList = self.readDirectories()
        self.parseDirectories(dirList)
        if PHATE_PROGRESS:
            print("genomics_compareGenomes says, Data loading complete.")
        return

    def readDirectories(self):
        dirList = []
        # Query Results files in CGP Results directory
        command = "ls " + CGP_RESULTS_DIR
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        (rawresult, err) = proc.communicate()
        result = rawresult.decode('utf-8')
        fileList = str(result).split('\n')   # Python3
        for filename in fileList:
            match_resultdir = re.search('Results_',filename)
            if match_resultdir:
                dirList.append(filename)
        if PHATE_MESSAGES:
            print("genomics_compareGenomes says, The following directories have been read:")
            for directory in dirList:
                print("directory:",directory)
        return dirList

    def parseDirectories(self,dirList):
        genomeList = []; genomeFastaList = []
        for directory in dirList:
            cgp_directory = os.path.join(CGP_RESULTS_DIR,directory)
            if PHATE_PROGRESS:
                print("genomics_compareGenomes says, Loading data from directory,",cgp_directory)

            # Read log file; determine two genomes
            CGP_LOG = os.path.join(cgp_directory,CGP_LOG_FILE)
            cgpLog = open(CGP_LOG,'r')
            cLines = cgpLog.read().splitlines()
            for cLine in cLines:
                match_fileLine = re.search('genome\sfile',cLine)
                if match_fileLine:
                    (preamble,filePath) = cLine.split(': ')
                    print("filePath is",filePath)
                    genomeFasta = os.path.basename(filePath)
                    if genomeFasta not in genomeFastaList:
                        genomeFastaList.append(genomeFasta)
                    (genome,extension) = genomeFasta.split('.')
                    print("genome is",genome)
                    if genome not in genomeList:
                        genomeList.append(genome)
            cgpLog.close()
            print("genomeList:",genomeList)
            # First genome listed is the reference (was listed first in user's config file)
            if genomeList:
                self.referenceGenome = genomeList[0]
            else:
                if PHATE_WARNINGS:
                    print("genomics_compareGenomes says, WARNING: There are no genomes to compare.")

        # Create genome objects, one for each genome found
        #for genome in genomeList:
        for i in range(0,len(genomeList)):
            nextGenome = copy.deepcopy(self.genomeTemplate)
            nextGenome.name = genomeList[i]
            if nextGenome.name == self.referenceGenome:
                nextGenome.isReference = True
            nextGenome.file = genomeFastaList[i]
            self.genomeList.append(nextGenome)
        self.countGenomes()
        print("There are this many genome objects in comparison object:",len(self.genomeList))
        print("List of genome objects in comparison object:")
        self.printAll()

        # Extract data from report files

        # Identify all the genomes that were processed by CGP
        # To be coded... establishes self.genomeList

        # Walk through .report files, add mutual & singular best hits, loners
        # To be coded... establihsed mutualBestHitList, singularBestHit, lonerList for each genome; load data structures 

        if PHATE_PROGRESS:
            print("genomics_compareGenomes says, Data load complete.")

        return

    # Method identifyParalogs inspects the CGP blast output for genome x self to identify paralogs
    def identifyParalogs(self):
        if DEBUG:
            print("genomics_compareGenomes says, DEBUG: Number of genomes in self.genomeList:",len(self.genomeList))
        for genome in self.genomeList:
            if PHATE_PROGRESS:
                print("genomics_compareGenomes says, Identifying paralogs for genome",genome.name)
            genome.identifyParalogs()

        if PHATE_PROGRESS:
            print("genomics_compareGenomes says, Paralog identification complete.")
        return

    def computeCoreGenome(self):

        # Compute core genomes


        if PHATE_PROGRESS:
            print("genomics_compareGenomes says, Core genome computation complete.")
        return

    def computeHomologyGroups(self):

        # Compute homology groups


        if PHATE_PROGRESS:
            print("genomics_compareGenomes says, Homology group computation complete.")
        return

    # For each gene/protein in reference genome, write correspondence set to file
    def saveHomologyGroups(self):
        for genome in self.genomeList:
            if genome.isReference:
                for gene in genome.geneList:
                    pass
                for protein in genome.proteinList:
                    pass
        if PHATE_PROGRESS:
            print("genomics_compareGenomes says, Homology groups have been written to files.")
        return

    def countGenomes(self):
        self.genomeCount = len(self.genomeList)
        if PHATE_MESSAGES:
            print("genomics_compareGenomes says, There are",self.genomeCount,"genomes.")
        if DEBUG:
            print("genomics_compareGenomes says, DEBUG: There are",self.genomeCount,"genomes.")
        return

    # For development purposes
    def checkUnique(self):
        tempList = []
        for genome in genomeList:
            if genome.name in tempList:
                if PHATE_WARNINGS:
                    print("genomics_compareGenomes genomes obj says, WARNING: ",genome.name," occurs more than once in genomes object.")
            else:
                tempList.append(genome.name)
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

# Class genome stores metadat aedDirectories(self):
class genome(object):

    def __init__(self):
        self.name                 = ""     # Name of this genome (e.g., Lambda)
        self.species              = ""     # Species name for this genome (e.g., Ecoli_Lambda_phage)
        self.isReference          = False  # True if designated a reference genome
        self.contigList           = []     # Set of contig names (fasta headers)
        self.geneList             = []     # List of gene_protein objects 
        self.proteinList          = []     # List of gene_protein objects 
        self.paralogList          = []     # List of paralogSet objects

    # For development purposes
    def checkUnique(self):
        tempList = []
        for gene in self.geneList:
            if gene.identifier in tempList:
                if PHATE_WARNINGS:
                    print("genomics_compareGenomes genome obj says, WARNING: ",gene.identifier," is not unique")
            else:
                tempList.append(gene.identifier)
            
        tempList = []
        for gene in self.geneList:
            if gene.cgpHeader in tempList:
                if PHATE_WARNINGS:
                    print("genomics_compareGenomes genome obj says, WARNING: ",gene.cgpHeader," is not unique")
            else:
                tempList.append(gene.cgpHeader)
            
        tempList = []
        for protein in self.proteinList:
            if protein.identifier in tempList:
                if PHATE_WARNINGS:
                    print("genomics_compareGenomes genome obj says, WARNING: ",protein.identifier," is not unique")
            else:
                tempList.append(protein.identifier)

        tempList = []
        for protein in self.proteinList:
            if protein.cgpHeader in tempList:
                if PHATE_WARNINGS:
                    print("genomics_compareGenomes genome obj says, WARNING: ",protein.cgpHeader," is not unique")
            else:
                tempList.append(protein.cgpHeader)
        return

    def identifyParalogs(self):

        # 

        if PHATE_PROGRESS:
            print("genomics_compareGenomes says, Paralogs identified for genome",self.name)
        return

    def printReport(self):
        return

    def printAll(self):
        print("=========Genome========")
        print("name:",self.name)
        print("species:",self.species)
        print("isReference:",self.isReference)
        print("contigs:",)
        for contig in self.contigList:
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
        if PHATE_MESSAGES:
            print("genomics_compareGenomes says, There are",self.setSize,"paralogs")
        return

    def printReport(self):
        return

    def printAll(self):
        print("paralogType:",self.paralogType)
        print("setSize:",self.setSize)
        print("paralogList:",self.paralogList)
        return

# Class gene stores meta-data about a gene
class gene_protein(object):

    def __init__(self):
        self.name                 = ""        # Ex: "phanotate_5"
        self.type                 = "gene"    # "gene" or "protein"
        self.identifier           = ""        # Unique: genome + contig + name
        self.cgpHeader            = ""        # CGP-assigned header; Ex: cgp5/+/72/485/
        self.number               = 0         # Unique: Set to correspond to protein 
        self.parentGenome         = ""        # Ex: "Lambda" - corresponds to a genome object's name
        self.contigName           = ""        # Ex: "Lambda_contig_1"; read from CGP report
        self.annotation           = ""        # The full annotation string, read from CGP report 
        self.isLoner              = False     # True if this gene has no correspondences
        self.mutualBestHitList    = []        # list of gene identifiers that are mutual best hits, across genomes
        self.singularBestHitList  = []        # list of gene identifiers that are best hits, relative to this gene, across genomes
        self.correspondenceList   = []        # List of corresponding genes (mutual or best-hit) across genomes: list of gene identifiers 
        self.groupList            = []        # list of all corresponding genes plus paralogs
        self.paralogList          = []        # List of paralogs within its own parent genome: list of gene identifiers <data is redundant; should pull paralog list as list of lists at genome level.

    def addMutualBestHit(self,hit):
        self.mutualBestHistList.append(hit)
        return

    def addSingularBestHit(self,hit):
        self.singularBestHitList.append(hit)
        return

    def addGroupMember(self,member):
        self.groupList.append(member)

    # Construct a set of genes/proteins that are related by homology
    def constructGroup(self):
        # self.groupList

        if PHATE_MESSAGES:
            print("genomics_compareGenomes says, Group constructed for gene_protein",self.identifier,";",self.type)
        return

    # For development purposes
    def verifyLoner(self):
        if self.mutualBestHitList:
            if PHATE_WARNINGS:
                print("genomics_compareGenomes says, Not a loner due to mutual best hit:",self.identifier,";",self.type)
            return False
        if self.singularBestHistList:
            if PHATE_WARNINGS:
                print("genomics_compareGenomes says, Not a loner due to singular best hit:",self.identifier,";",self.type)
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


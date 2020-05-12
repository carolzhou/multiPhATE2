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

# Report file fields
SORT_POSITION        = 0
GENOME_TYPE          = 1
QUERY_START          = 2
QUERY_END            = 3
SUBJECT_START        = 4
SUBJECT_END          = 5
IDENTITY             = 6
E_VALUE              = 7
G1_Q_S_LONER         = 8
G1_HEADER            = 9
G1_CONTIG            = 10
G1_ANNOTATIONS       = 11
G1_GENE_START        = 12
G1_GENE_END          = 13
G1_STRAND            = 14
G2_Q_S_LONER         = 15
G2_HEADER            = 16
G2_CONTIG            = 17
G2_ANNOTATIONS       = 18
G2_GENE_START        = 19
G2_GENE_END          = 20
G2_STRAND            = 21
G1_COVERAGE          = 22
G2_COVERAGE          = 23
ALIGNMENT_LENGTH     = 24
GAP_OPENS            = 25
G1_SPAN              = 26
G1_LENGTH            = 27
G2_SPAN              = 28
G2_LENGTH            = 29

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
                    genomeFasta = os.path.basename(filePath)
                    if genomeFasta not in genomeFastaList:
                        genomeFastaList.append(genomeFasta)
                    (genome,extension) = genomeFasta.split('.')
                    if genome not in genomeList:
                        genomeList.append(genome)
            cgpLog.close()

            # First genome listed is the reference (was listed first in user's config file)
            if genomeList:
                self.referenceGenome = genomeList[0]
            else:
                if PHATE_WARNINGS:
                    print("genomics_compareGenomes says, WARNING: There are no genomes to compare.")

        # Create genome objects, one for each genome found

        for i in range(0,len(genomeList)):
            nextGenome = copy.deepcopy(self.genomeTemplate)
            nextGenome.name = genomeList[i]
            if nextGenome.name == self.referenceGenome:
                nextGenome.isReference = True
            nextGenome.file = genomeFastaList[i]
            self.genomeList.append(nextGenome)
        #self.countGenomes()

        # Walk through .report files, add mutual & singular best hits, loners
        if PHATE_PROGRESS:
            print("genomics_compareGenomes says, Invoking self.parseReportFiles with dirList,",dirList)
        self.parseReportFiles(dirList)
        
        return

    def parseReportFiles(self,dirList):

        if PHATE_PROGRESS:
            print("genomics_compareGenomes says, Parsing Report files.")
        for directory in dirList:
            if PHATE_PROGRESS:
                print("genomics_compareGenomes says, Identifying Results directories in directory",directory)
            # First, determine which 2 genomes's data are listed in this directory,
            # and which is genome1, genome2.
            genome1 = ""; genome2 = ""
            logFile = os.path.join(CGP_RESULTS_DIR,directory,CGP_LOG_FILE)
            LOG_H = open(logFile,"r")
            lLines = LOG_H.read().splitlines()
            for lLine in lLines:
                match_genome1 = re.search('genome\sfile\s#1:',lLine)
                match_genome2 = re.search('genome\sfile\s#2:',lLine)
                if match_genome1:
                    (preamble,genomePath) = lLine.split(': ')
                    genome1fasta = os.path.basename(genomePath)
                    (genome1,extension) = genome1fasta.split('.') 
                if match_genome2:
                    (preamble,genomePath) = lLine.split(': ')
                    genome2fasta = os.path.basename(genomePath)
                    (genome2,extension) = genome2fasta.split('.') 
            LOG_H.close()

            # Parse report file; extract hits into gene and protein objects
            reportFile = os.path.join(CGP_RESULTS_DIR,directory,CGP_REPORT_FILE)
            if PHATE_PROGRESS:
                print("genomics_compareGenomes says, Parsing reportFile",reportFile)
            DONE = False
            REPORT_H = open(reportFile,"r")
            rLines = REPORT_H.read().splitlines()
            for rLine in rLines:
                fields = []; genomeNum = ""; hitType = ""
                dataArgs = {
                    "genome1"   : genome1,   # the fasta file name
                    "genome2"   : genome2,   # the fasta file name
                    "contig1"   : "",
                    "contig2"   : "",
                    "gene1"     : "",
                    "gene2"     : "",
                    "protein1"  : "",
                    "protein2"  : "",
                    "geneCall1" : "",
                    "geneCall2" : "",
                    "hitType"   : "",
                    }
                match_comment = re.search('^#',rLine)
                match_protein = re.search('^#\sPROTEIN',rLine)
                if match_protein:
                    DONE = True
                    break
                if match_comment:
                    continue
                fields = rLine.split('\t')
                (genomeNum,hitType) = fields[GENOME_TYPE].split('_')
                if genomeNum == "genome1":
                    dataArgs["hitType"] = hitType + '1'
                elif genomeNum == "genome2":
                    dataArgs["hitType"] = hitType + '2'
                else:
                    if PHATE_WARNINGS:
                        print("genomics_compareGenomes says, WARNING: Unrecognized GENOME_TYPE")
                if hitType == "mutual" or hitType == "singular":
                    dataArgs["contig1"]   = fields[G1_CONTIG]
                    dataArgs["contig2"]   = fields[G2_CONTIG]
                    dataArgs["gene1"]     = fields[G1_HEADER]
                    dataArgs["gene2"]     = fields[G2_HEADER]
                    annotations1          = fields[G1_ANNOTATIONS]
                    annotations2          = fields[G2_ANNOTATIONS]
                    dataArgs["geneCall1"] = annotations1.split(';')[0]
                    dataArgs["geneCall2"] = annotations1.split(':')[0]
                elif hitType == "loner":
                    if genomeNum == "genome1":
                        dataArgs["contig1"]   = fields[G1_CONTIG]
                        dataArgs["gene1"]     = fields[G1_HEADER]
                        annotations           = fields[G1_ANNOTATIONS]
                        dataArgs["geneCall1"] = annotations.split(';')[0] 
                    elif genomeNum == "genome2":
                        dataArgs["contig2"]   = fields[G2_CONTIG]
                        dataArgs["gene2"]     = fields[G2_HEADER]
                        annotations           = fields[G2_ANNOTATIONS]
                        dataArgs["geneCall2"] = annotations.split(';')[0] 
                else:
                    if PHATE_WARNINGS:
                        print("genomics_compareGenomes says, WARNING: Unrecognized hitType",hitType,"in report file",reportFile)
                self.addHit2genome(dataArgs)
            REPORT_H.close()

        if PHATE_PROGRESS:
            print("genomics_compareGenomes says, Data load complete.")

        return

    def findGenomeObject(self,genomeName):

        for genome_object in self.genomeList:
            if genome_object.name == genomeName:
                return genome_object
        return

    def isNewGene(self,genomeObj,contig,geneHeader):  # N^2 search; is there a way to do NlogN search here, or better?
        for geneObj in genomeObj.geneList:
            if geneObj.contig == contig and geneObj.name == geneHeader:
                return True 
        return False 

    def addHit2genome(self,dataArgs):  #*** FIX: NOT CAPTURING THE FULL ANNOTATION FIELD FROM CGP REPORT
        #print("CHECKING: dataArgs is",dataArgs)

        GENE_CALL_NUMBER = -3  # Index for accessing gene-call number, as third-to-last item of the gene-call string
        GENOME_ONE = False; GENOME_TWO = False # Should exist data for both genomes if mutual or singular hit; one if loner

        # Data needed for "mutual" or "singular" hit entry
        # Determine which genome's hit to process, or both.
        # Find the two genome objects to which mutual hits will be added
        # Compute unique identifiers for query and subject genes 
        if dataArgs["gene1"] != "":
            GENOME_ONE = True
            genome1_obj = self.findGenomeObject(dataArgs["genome1"])
            gene1id = dataArgs["genome1"] + ':' + dataArgs["contig1"] + ':' + dataArgs["gene1"]
        if dataArgs["gene2"] != "":
            GENOME_TWO = True
            genome2_obj = self.findGenomeObject(dataArgs["genome2"])
            gene2id = dataArgs["genome2"] + ':' + dataArgs["contig2"] + ':' + dataArgs["gene2"]

        # Process mutual-best hit 
        #if dataArgs["hitType"] == "mutual":
        if GENOME_ONE:
            # Enter hit data for query gene (gene1)
            # Search for query gene in genome1 (if exists)
            GENE_FOUND = False
            for geneObj in genome1_obj.geneList:
                if geneObj.contigName == dataArgs["contig1"] and geneObj.name == dataArgs["gene1"]:
                    GENE_FOUND = True
                    # Add hit
                    match_mutual   = re.search("mutual",  dataArgs["hitType"])
                    match_singular = re.search("singular",dataArgs["hitType"])
                    match_loner    = re.search("loner",   dataArgs["hitType"])
                    if match_mutual:
                        geneObj.mutualBestHitList.append(gene2id)         #*** FIX: if we call geneObj, newGene same name, can consolidate code
                        geneObj.isLoner = False
                        genome1_obj.geneList.append(geneObj)
                    elif match_singular:
                        geneObj.singularBestHitList.append(gene2id)
                        geneObj.isLoner = False
                        genome1_obj.geneList.append(geneObj)
                    elif match_loner:
                        geneObj.isLoner = True
                        genome1_obj.geneList.append(geneObj)
                    else:
                        if PHATE_WARNINGS:
                            print("genomics_compareGenomes says, WARNING: Unrecognized hitType,",dataArgs["hitType"]) 
                    continue
            # Or create and populate a new gene object, as needed
            if not GENE_FOUND:  
                newGene = copy.deepcopy(self.geneTemplate)
                newGene.name         = dataArgs["gene1"]
                newGene.type         = "gene"
                newGene.identifier   = gene1id 
                newGene.cgpHeader    = dataArgs["gene1"]
                geneCallFields       = dataArgs["geneCall1"].split('_')
                newGene.number       = geneCallFields[GENE_CALL_NUMBER]
                newGene.parentGenome = dataArgs["genome1"]
                newGene.contigName   = dataArgs["contig1"]
                newGene.annotation   = dataArgs["geneCall1"]  
                # Add hit
                match_mutual   = re.search("mutual",  dataArgs["hitType"])
                match_singular = re.search("singular",dataArgs["hitType"])
                match_loner    = re.search("loner",   dataArgs["hitType"])
                if match_mutual:
                    newGene.mutualBestHitList.append(gene2id) 
                    newGene.isLoner = False
                    genome1_obj.geneList.append(newGene)
                elif match_singular:
                    newGene.singularBestHitList.append(gene2id)
                    newGene.isLoner = False
                    genome1_obj.geneList.append(newGene)
                elif match_loner:
                    newGene.isLoner = True
                    genome1_obj.geneList.append(newGene)
                else:
                    if PHATE_WARNINGS:
                        print("genomics_compareGenomes says, WARNING: Unrecognized hitType,",dataArgs["hitType"]) 

        if GENOME_TWO:
            # Enter hit data for subject gene (gene2)
            # Search for query gene in genome2 (if exists)
            GENE_FOUND = False
            for geneObj in genome2_obj.geneList:
                if geneObj.contigName == dataArgs["contig2"] and geneObj.name == dataArgs["gene2"]:
                    GENE_FOUND = True
                    # Add hit
                    match_mutual   = re.search("mutual",  dataArgs["hitType"])
                    match_singular = re.search("singular",dataArgs["hitType"])
                    match_loner    = re.search("loner",   dataArgs["hitType"])
                    if match_mutual:
                        newGene.mutualBestHitList.append(gene2id) 
                        newGene.isLoner = False
                        genome2_obj.geneList.append(newGene)
                    elif match_singular:  #*** ERROR in CGP: not printing genome2 singular hits in report
                        newGene.singularBestHitList.append(gene2id)
                        newGene.isLoner = False
                        genome2_obj.geneList.append(newGene)
                    elif match_loner:
                        newGene.isLoner = True
                        genome2_obj.geneList.append(newGene)
                    else:
                        if PHATE_WARNINGS:
                            print("genomics_compareGenomes says, WARNING: Unrecognized hitType,",dataArgs["hitType"]) 
                    continue
            # Or create and populate a new gene object, as needed
            if not GENE_FOUND:  
                newGene = copy.deepcopy(self.geneTemplate)
                newGene.name         = dataArgs["gene2"]
                newGene.type         = "gene"
                newGene.identifier   = gene2id 
                newGene.cgpHeader    = dataArgs["gene2"]
                geneCallFields       = dataArgs["geneCall2"].split('_')
                newGene.number       = geneCallFields[GENE_CALL_NUMBER]
                newGene.parentGenome = dataArgs["genome2"]
                newGene.contigName   = dataArgs["contig2"]
                newGene.annotation   = dataArgs["geneCall2"] 
                # Add hit
                match_mutual   = re.search("mutual",  dataArgs["hitType"])
                match_singular = re.search("singular",dataArgs["hitType"])
                match_loner    = re.search("loner",   dataArgs["hitType"])
                if match_mutual:
                    newGene.mutualBestHitList.append(gene2id) 
                    newGene.isLoner = False
                    genome2_obj.geneList.append(newGene)
                elif match_singular:  #*** ERROR in CGP: not printing genome2 singular hits in report
                    newGene.singularBestHitList.append(gene2id)
                    newGene.isLoner = False
                    genome2_obj.geneList.append(newGene)
                elif match_loner:
                    newGene.isLoner = True
                    genome2_obj.geneList.append(newGene)
                else:
                    if PHATE_WARNINGS:
                        print("genomics_compareGenomes says, WARNING: Unrecognized hitType,",dataArgs["hitType"]) 
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


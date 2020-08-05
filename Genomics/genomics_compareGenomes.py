#################################################################################
# Module: genomics_compareGenomes.py
#
# Programmer:  Carol L. Ecale Zhou
#
# Most recent update: 03 June 2020
#
# Module comprising classes and data structures for comparing genomes
#
# Classes and Methods:
#    class comparison
#       performComparison()
#       ---data input
#       loadData()
#       readDirectories()
#       parseDirectories()
#       parseReportFiles()
#       findGenomes()
#       loadBestHits()
#       getGeneCallString()
#       addHit2genome()
#       findGenomeObject()
#       ---comparison
#       identifyParalogs()
#       computeHomologyGroups()
#       saveHomologyGroups()
#       ---verification
#       countGenomes()
#       checkUnique()
#       ---reporting
#       writeCoreGenome()
#       writeCorrespondences()
#       writeMutualBestHitList()
#       writeSingularBestHitList()
#       writeLonerList()
#       printReport()
#       printAll()
#    class genome
#       ---comparison
#       identifyParalogs()
#       ---verification
#       checkUnique()
#       ---reporting
#       printReport()
#       printAll()
#    class paralogSet
#       ---verification
#       countParalogs()
#       ---reporting
#       printReport()
#       printAll()
#    class gene_protein
#       ---data input
#       addMutualBestHit()
#       addSingularBestHit()
#       addGroupMember()
#       ---verification
#       verifyLoner()
#       ---reporting
#       writeMutualBestHitList()
#       writeSingularBestHitList()
#       writeLonerList()
#       printReport()
#       printAll()
#
#################################################################################

# This code was developed by Carol. L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.PDF FOR DETAILS.

import re, os, copy
import subprocess
import ast

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

# Override
PHATE_PROGRESS = True
PHATE_WARNINGS = True
PHATE_MESSAGES = True

DEBUG = False
#DEBUG = True

# Locations and files
CGP_RESULTS_DIR                 = os.environ["PHATE_CGP_RESULTS_DIR"]
CGP_LOG_FILE                    = "compareGeneProfiles_main.log"           # Lists genome fasta and annotation files
CGP_REPORT_FILE                 = "compareGeneProfiles_main.report"        # Tabbed file with gene-gene listing
CGP_PARALOG_FILE                = "compareGeneProfiles_main.paralogs"      # Paralogs detected by CGP; text format 
CGP_OUT_FILE                    = "compareGeneProfiles_main.out"           # Same data as report file, but in python list format
CGP_SUMMARY_FILE                = "compareGeneProfiles_main.summary"       # High-level summary of genome-genome comparison
GENOMICS_RESULTS_DIR            = os.environ["PHATE_GENOMICS_RESULTS_DIR"]
GENOMICS_HOMOLOGY_GROUPS_DIR    = GENOMICS_RESULTS_DIR + "HOMOLOGY_GROUPS"
GENOMICS_MUTUAL_BEST_HIT_FILE   = "genomics_mutualBestHits.out"
GENOMICS_SINGULAR_BEST_HIT_FILE = "genomics_singularBestHits.out"
GENOMICS_CORE_GENOME_FILE       = "genomics_coreGenome.out"
GENOMICS_CORRESPONDENCE_FILE    = "genomics_correspondences.out"
GENOMICS_LONER_FILE             = "genomics_loners.out"
GENOMICS_PARALOG_FILE           = "genomics_paralogs.out"
GENOMICS_HOMOLOGY_GROUP_FILE    = "genomics_homologyGroups.out"
GENOMICS_HOMOLOGY_FASTA_FILE    = "genomics.homGrp"
HOMOLOGY_PREFIX                 = "homologyGroup_"
HOMOLOGY_ANNOT_PREFIX           = "homologyGroup_"
MUTUAL_BEST_HIT_FILE            = os.path.join(GENOMICS_RESULTS_DIR,GENOMICS_MUTUAL_BEST_HIT_FILE)
SINGULAR_BEST_HIT_FILE          = os.path.join(GENOMICS_RESULTS_DIR,GENOMICS_SINGULAR_BEST_HIT_FILE)
CORE_GENOME_FILE                = os.path.join(GENOMICS_RESULTS_DIR,GENOMICS_CORE_GENOME_FILE)
CORRESPONDENCE_FILE             = os.path.join(GENOMICS_RESULTS_DIR,GENOMICS_CORRESPONDENCE_FILE)
LONER_FILE                      = os.path.join(GENOMICS_RESULTS_DIR,GENOMICS_LONER_FILE)
PARALOG_FILE                    = os.path.join(GENOMICS_RESULTS_DIR,GENOMICS_PARALOG_FILE)
HOMOLOGY_GROUP_FILE             = os.path.join(GENOMICS_RESULTS_DIR,GENOMICS_HOMOLOGY_GROUP_FILE)

PHATE_PIPELINE_OUTPUT_DIR       = os.environ["PHATE_PIPELINE_OUTPUT_DIR"] 
CGP_GENES_FILENAME              = "cgp_gene.fnt"
CGP_PROTEINS_FILENAME           = "cgp_protein.fnt"   # This should be ".faa"

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

        self.name                  = "multiPhATE2"  # Can assign a set to this genome comparison
        self.directories           = []             # List of Results directories to be processed
        self.referenceGenome       = ""             # name of genome assigned as the reference
        self.genomeCount           = 0              # number of genomes in set
        self.genomeList            = []             # set of genome class objects
        self.commonCore_gene       = []             # list of genes that are common among all genomes: all around mutual best hits
        self.commonCore_protein    = []             # list of protein that are common among all genomes: all around mutual best hits
        self.geneHomologyGroups    = []             # closely related genes for hmmbuild (list of lists) #*** CHECK THIS - is this used at comparison level?
        self.proteinHomologyGroups = []             # closely related proteins for hmmbuild (list of lists) #*** CHECK THIS
        # Templates
        self.genomeTemplate        = genome()
        self.geneTemplate          = gene_protein()
        self.geneTemplate.type     = "gene"
        self.proteinTemplate       = gene_protein()
        self.proteinTemplate.type  = "protein"

    def performComparison(self):

        if PHATE_PROGRESS:
            print("genomics_compareGenomes says, Loading data.")
        self.loadData()

        # Do data checking
        self.runDataChecks()

        if PHATE_PROGRESS:
            print("genomics_compareGenomes says, Computing homology groups.")
        self.computeHomologyGroups()

        if PHATE_PROGRESS:
            print("genomics_compareGenomes says, Creating homology group fasta files.")
            print("...for genes...")
        try:
            os.stat(GENOMICS_HOMOLOGY_GROUPS_DIR)
        except:
            os.mkdir(GENOMICS_HOMOLOGY_GROUPS_DIR)
        self.createGeneHomologyFastaFiles(GENOMICS_HOMOLOGY_GROUPS_DIR)
        if PHATE_PROGRESS:
            print("...for proteins...")
        self.createProteinHomologyFastaFiles(GENOMICS_HOMOLOGY_GROUPS_DIR)

        if PHATE_PROGRESS:
            print("genomics_compareGenomes says, Writing final reports.")
        self.printReports2files()

    #===== COMPARISON DATA LOAD METHODS

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
        # Walk through the CGP output directory; determine which genomes were processed; add them to a non-redundant list
        genomeList = []; genomeFastaList = []
        for directory in dirList:
            cgp_directory = os.path.join(CGP_RESULTS_DIR,directory)
            if PHATE_PROGRESS:
                print("genomics_compareGenomes says, Preparing to load data from directory,",cgp_directory)

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

        # Create genome objects, one for each genome found in the CGP output directory
        for i in range(0,len(genomeList)):
            nextGenome = copy.deepcopy(self.genomeTemplate)
            nextGenome.name = genomeList[i]
            if nextGenome.name == self.referenceGenome:
                nextGenome.isReference = True
            nextGenome.file = genomeFastaList[i]
            self.genomeList.append(nextGenome)

        # Walk through .report files, add mutual & singular best hits, loners
        if PHATE_PROGRESS:
            print("genomics_compareGenomes says, Invoking parseReportFiles with dirList,",dirList)
        self.parseReportFiles(dirList)
        return

    # Method parseReportFiles pulls data from CGP Results_* directories into memory in order to determine
    # a core genome, gene-gene correspondences among all genomes, and list loner genes per genome.
    def parseReportFiles(self,dirList):
        if PHATE_PROGRESS:
            print("genomics_compareGenomes says, Parsing Report files.")
        # Load data for mutual and singular best hits and loners for all CGP binary genome comparisons.
        for i in range(0,len(dirList)):
            directory = dirList[i]
            (genome1,genome2) = self.findGenomes(directory)
            reportFile  = os.path.join(CGP_RESULTS_DIR,directory,CGP_REPORT_FILE)
            paralogFile = os.path.join(CGP_RESULTS_DIR,directory,CGP_PARALOG_FILE)
            if PHATE_PROGRESS:
                print("genomics_compareGenomes says, Parsing report file",reportFile,"for mutual and singular best hits, and loners.")
            self.loadBestHits(genome1,genome2,reportFile)
            self.loadParalogs(paralogFile)

        return

    def findGenomes(self,directory):
        # Determine which 2 genomes's data are listed in this directory,
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
        return(genome1,genome2)

    # Method loadBestHits reads the CGP report file in a single pass, picking up mutual-best and singular-best hits
    # plus loners. The hits are then forwarded to method addhit2genome for recording in the genome's gene_protein
    # data structure.
    def loadBestHits(self,genome1,genome2,reportFile):
        PROTEIN = False; hitCount = 0
        REPORT_H = open(reportFile,"r")
        rLines = REPORT_H.read().splitlines()

        # Select and process mutual-best-hit data lines
        # Gene hits are loaded first, then protein
        for rLine in rLines:
            fields = []; genomeNum = ""; hitType = ""
            # Skip lines not to be processed in this method
            match_comment  = re.search('^#',rLine)
            match_protein  = re.search('^#PROTEIN\sHITS',rLine)
            match_dataLine = re.search('^\d+',rLine)
            if match_protein: 
                PROTEIN = True  # Prepare for loading protein hits next
                continue
            if match_comment:
                continue
            if match_dataLine:
                fields = rLine.split('\t')
                (genomeNum,hitType) = fields[GENOME_TYPE].split('_')
            else:
                if PHATE_WARNINGS:
                    print("genomics_compareGenomes says, WARNING: expected dataLine: ",rLine)
                    continue
            # Data structure for passing parameters to self.addHit2genome()
            dataArgs = { 
                "genome1"      : genome1,   # the fasta file name without extension
                "genome2"      : genome2,   # the fasta file name without extension
                "contig1"      : "",
                "contig2"      : "",
                "gene1"        : "",
                "gene2"        : "",
                "protein1"     : "",
                "protein2"     : "",
                "geneCall1"    : "",
                "geneCall2"    : "",
                "annotations1" : "",
                "annotations2" : "",
                "hitType"      : "",
                "hitFlavor"    : "",
                }

            # Prepare arguments for passing to addHit2genome method
            if genomeNum == "genome1":
                dataArgs["hitType"] = hitType + '1'
            elif genomeNum == "genome2":
                dataArgs["hitType"] = hitType + '2'
            else:
                if PHATE_WARNINGS:
                    print("genomics_compareGenomes says, WARNING: Unrecognized GENOME_TYPE")

            # Record data from reportFile dataline
            # Fields pertaining to gene or protein
            dataArgs["contig1"]          = fields[G1_CONTIG]
            dataArgs["contig2"]          = fields[G2_CONTIG]
            dataArgs["annotations1"]     = fields[G1_ANNOTATIONS]
            if dataArgs["annotations1"] != "":
                dataArgs["geneCall1"]    = self.getGeneCallString(fields[G1_ANNOTATIONS]) 
            dataArgs["annotations2"]     = fields[G2_ANNOTATIONS]
            if dataArgs["annotations2"] != "":
                dataArgs["geneCall2"]    = self.getGeneCallString(fields[G2_ANNOTATIONS]) 

            # Protein- and gene-specific fields
            if PROTEIN:
                dataArgs["hitFlavor"]    = "protein"
                dataArgs["protein1"]     = fields[G1_HEADER]
                dataArgs["protein2"]     = fields[G2_HEADER]
            else: # gene
                dataArgs["hitFlavor"]    = "gene"
                dataArgs["gene1"]        = fields[G1_HEADER]
                dataArgs["gene2"]        = fields[G2_HEADER]

            # Insert data record into genome object
            self.addHit2genome(dataArgs)

        REPORT_H.close()
        return

    def loadParalogs(self,paralogFile):
        # Select and process mutual-best-hit data lines
        fields = []; genomeNum = ""; paralogType = ""; genomeName = ""
        GENE = False; PROTEIN = False
        query = ""; subject = ""

        # First, reset dataArgs data structure, for passing parameters to self.addHit2genome
        dataArgs = { 
            "genome1"      : "",   # the fasta file name
            "genome2"      : "",   
            "contig1"      : "",   
            "contig2"      : "",   
            "gene1"        : "",
            "gene2"        : "",
            "protein1"     : "",
            "protein2"     : "",
            "geneCall1"    : "",
            "geneCall2"    : "",
            "annotations1" : "",
            "annotations2" : "",
            "hitType"      : "",
            }

        # Parse data from paralogs report file
        PARALOG_H = open(paralogFile,"r")
        pLines = PARALOG_H.read().splitlines()
        for pLine in pLines:

            # Skip lines not to be processed in this method
            match_genome   = re.search('PARALOGS',          pLine)
            match_gene     = re.search('Gene\sParalogs',    pLine)
            match_protein  = re.search('Protein\sParalogs', pLine)
            match_header   = re.search('^header:',           pLine)
            match_hitLine  = re.search('^query:',            pLine)
            match_coverage = re.search('^coverage:',        pLine)

            # Determine which genome's paralogs are being reported
            if match_genome:
                match_path = re.search('PARALOGS\sfor\sgenome\s(.*)',pLine)
                genomePathString = match_path.group(1)
                if genomePathString:
                    genomeFileString = os.path.basename(genomePathString)
                    (genomeName,extension) = genomeFileString.split('.')
                    dataArgs["genome1"] = genomeName
                else:
                    if PHATE_WARNINGS:
                        print("genomics_compareGenomes says, WARNING: Cannot read genome path from paralogs file,",paralogFile)

            # Determine whether it's a gene versus protein paralog
            elif match_gene:
                GENE = True
            elif match_protein:
                PROTEIN = True

            # Parse paralog data; add to paralog list 
            elif match_hitLine:
                fields = pLine.split('\t')
                queryString        = fields[0]
                subjectString      = fields[1]
                hitType            = fields[2]
                identity           = fields[3]
                alignLength        = fields[4]
                mismatches         = fields[5]
                gapopens           = fields[6]
                queryStartEnd      = fields[7]
                subjectStartEnd    = fields[8]
                (preamble,query)   = queryString.split(':')
                (preamble,subject) = subjectString.split(':')
                if GENE:
                    dataArgs["gene1"]    = query
                    dataArgs["gene2"]    = subject
                    dataArgs["hitType"]  = "gene_paralog"
                elif PROTEIN:
                    dataArgs["protein1"] = query
                    dataArgs["protein2"] = subject
                    dataArgs["hitType"]  = "protein_paralog"

            # Read coverage (last data line per paralog) and add this paralog to genome's paralogList; Then, reset
            elif match_coverage:
                (preamble,coverage) = match_coverage.group(0).split(':')

                # Insert data record into genome object
                self.addParalog2genome(dataArgs)

                # Reset data structures
                dataArgs = { 
                    "genome1"      : "",   # the fasta file name
                    "genome2"      : "",   # not used for paralog hit 
                    "contig1"      : "",   
                    "contig2"      : "",   
                    "gene1"        : "",
                    "gene2"        : "",
                    "protein1"     : "",
                    "protein2"     : "",
                    "geneCall1"    : "",
                    "geneCall2"    : "",
                    "annotations1" : "",
                    "annotations2" : "",
                    "hitType"      : "",
                    }
                GENE       = False
                PROTEIN    = False
                genomeName = ""
                query      = ""
                subject    = ""
                coverage   = 0.0

        PARALOG_H.close()
        return

    def addParalog2genome(self,dataArgs):
        if "genome1" in dataArgs.keys():
            genome = dataArgs["genome1"]
        else:
            if PHATE_WARNINGS:
                print("genomics_compareGenomes says, WARNING: Expected genome1 name in addParalog2genome")
                return
        for genome_obj in self.genomeList:
            if genome_obj.name == genome:
                genome_obj.addParalog(dataArgs)
                break
        return

    # Method getGeneCallString extracts the genecall name from the annotation string
    def getGeneCallString(self,inString):  
        inList = ast.literal_eval(inString) # Convert string representation of a list to an actual list
        geneCallString = inList[0][0]          # First element of the list is the genecall string
        return geneCallString

    # Method addHit2genome inserts a new gene hit into a genome's geneList or proteinList, or records a gene as a
    # loner with respect to another genome. Note that every gene is assumed a loner until proven otherwise. Any 
    # mutual or singular hit to a gene/protein of another genome nullifies that gene's status as a loner.
    def addHit2genome(self,dataArgs): 

        # First, create gene or protein object
        hitFlavor = "gene"; GENE = False; PROTEIN = False
        gene1id = ""; gene2id = ""; protein1id = ""; protein2id = ""
        gene_obj = None; protein_obj = None   # a gene_protein object (gene or protein template)

        # Create gene_protein object, according to flavor. hitFlavor is "gene" or "protein".
        if "hitFlavor" in dataArgs.keys():
            hitFlavor = dataArgs["hitFlavor"]
        else:
            if PHATE_WARNINGS:
                print("genomics_compareGenomes says, WARNING: hitFlavor not specified")
                return
        if hitFlavor == "gene":
            gene_obj = self.geneTemplate
            GENE = True
        elif hitFlavor == "protein":
            protein_obj = self.proteinTemplate 
            PROTEIN = True
        else:
            if PHATE_WARNINGS:
                print("genomics_compareGenomes says, WARNING: Unrecognized hitFlavor, ",dataArgs["hitFlavor"])
                return

        MUTUAL = False; SINGULAR_ONE = False; SINGULAR_TWO = False; LONER_ONE = False; LONER_TWO = False
        # Data needed for "mutual" or "singular" hit entry
        # Determine which genome's hit to process, or both.
        # Find the two genome objects to which mutual hits will be added
        # Compute unique identifiers for query and subject genes 

        # Determine type of hit and whether the query is genome1 versus genome2
        match_mutual   = re.search("mutual",  dataArgs["hitType"])
        match_singular = re.search("singular",dataArgs["hitType"])
        match_loner    = re.search("loner",   dataArgs["hitType"])
        match_query    = re.search("(\d)",    dataArgs["hitType"])
        if match_mutual:  # mutual hit is always genome 1 as query
            MUTUAL = True
        if match_singular and match_query.group(1) == '1':
            SINGULAR_ONE = True
        if match_singular and match_query.group(1) == '2':
            SINGULAR_TWO = True
        if match_loner and match_query.group(1) == '1':
            LONER_ONE = True
        if match_loner and match_query.group(1) == '2':
            LONER_TWO = True

        if GENE:
            if dataArgs["gene1"] != "":
                genome1_obj = self.findGenomeObject(dataArgs["genome1"])
                gene1id = dataArgs["genome1"] + ':' + dataArgs["contig1"] + ':' + dataArgs["gene1"]
            if dataArgs["gene2"] != "":
                genome2_obj = self.findGenomeObject(dataArgs["genome2"])
                gene2id = dataArgs["genome2"] + ':' + dataArgs["contig2"] + ':' + dataArgs["gene2"]
        elif PROTEIN:
            if dataArgs["protein1"] != "":
                genome1_obj = self.findGenomeObject(dataArgs["genome1"])
                protein1id = dataArgs["genome1"] + ':' + dataArgs["contig1"] + ':' + dataArgs["protein1"]
            if dataArgs["protein2"] != "":
                genome2_obj = self.findGenomeObject(dataArgs["genome2"])
                protein2id = dataArgs["genome2"] + ':' + dataArgs["contig2"] + ':' + dataArgs["protein2"]

        if MUTUAL or SINGULAR_ONE or LONER_ONE:
            # Search for query gene in genome1 (if exists)
            if GENE:
                GENE_FOUND = False; NEW = False
                for gene_obj in genome1_obj.geneList:
                    if gene_obj.contigName == dataArgs["contig1"] and gene_obj.name == dataArgs["gene1"]:
                        GENE_FOUND = True
                        break
                if not GENE_FOUND: 
                    # Create and populate a new gene object
                    NEW = True
                    gene_obj = copy.deepcopy(self.geneTemplate)
                    gene_obj.name         = dataArgs["gene1"]
                    gene_obj.type         = "gene"
                    gene_obj.identifier   = gene1id 
                    gene_obj.cgpHeader    = dataArgs["gene1"]
                    geneCallFields        = dataArgs["gene1"].split('/')   # format: cds#/strand/start/stop/
                    (cds,geneNumber)      = geneCallFields[0].split('cds')
                    gene_obj.number       = geneNumber
                    gene_obj.parentGenome = dataArgs["genome1"]
                    gene_obj.contigName   = dataArgs["contig1"]
                    gene_obj.annotation   = dataArgs["annotations1"] 

                # Add hit
                if MUTUAL:
                    gene_obj.mutualBestHitList.append(gene2id) 
                    gene_obj.isLoner = False
                if SINGULAR_ONE:
                    gene_obj.singularBestHitList.append(gene2id)
                    gene_obj.isLoner = False
                if LONER_ONE:
                    gene_obj.lonerList.append(dataArgs["genome2"])   # To record a loner, append the name of the genome that this gene is a loner wrt

                # If this gene object needed to be newly created, then add to geneList; else it's already there!
                if NEW:
                    genome1_obj.geneList.append(gene_obj)

            elif PROTEIN:
                PROTEIN_FOUND = False; NEW = False
                for protein_obj in genome1_obj.proteinList:
                    if protein_obj.contigName == dataArgs["contig1"] and protein_obj.name == dataArgs["protein1"]:
                        PROTEIN_FOUND = True
                        break
                if not PROTEIN_FOUND: 
                    # Create and populate a new gene object
                    NEW = True
                    protein_obj = copy.deepcopy(self.proteinTemplate)
                    protein_obj.name         = dataArgs["protein1"]
                    protein_obj.type         = "protein"
                    protein_obj.identifier   = protein1id 
                    protein_obj.cgpHeader    = dataArgs["protein1"]
                    proteinNameFields        = dataArgs["protein1"].split('/')  # format: cds#/strand/start/stop/
                    (cds,proteinNumber)      = proteinNameFields[0].split('cds')
                    protein_obj.number       = proteinNumber
                    protein_obj.parentGenome = dataArgs["genome1"]
                    protein_obj.contigName   = dataArgs["contig1"]
                    protein_obj.annotation   = dataArgs["annotations1"] 

                # Add hit
                if MUTUAL:
                    protein_obj.mutualBestHitList.append(protein2id) 
                    protein_obj.isLoner = False
                if SINGULAR_ONE:
                    protein_obj.singularBestHitList.append(protein2id)
                    protein_obj.isLoner = False
                if LONER_ONE:
                    protein_obj.lonerList.append(dataArgs["genome2"])

                # If this protein object needed to be newly created, then add to proteinList; else it's already there!
                if NEW:
                    genome1_obj.proteinList.append(protein_obj)

        if MUTUAL or SINGULAR_TWO or LONER_TWO:
            # Search for query gene in genome2 (if exists)
            if GENE:
                GENE_FOUND = False; NEW = False
                for gene_obj in genome2_obj.geneList:
                    if gene_obj.contigName == dataArgs["contig2"] and gene_obj.name == dataArgs["gene2"]:
                        GENE_FOUND = True
                        break
                if not GENE_FOUND: 
                    # Create and populate a new gene object
                    NEW = True
                    gene_obj = copy.deepcopy(self.geneTemplate)
                    gene_obj.name         = dataArgs["gene2"]
                    gene_obj.type         = "gene"
                    gene_obj.identifier   = gene2id 
                    gene_obj.cgpHeader    = dataArgs["gene2"]
                    geneCallFields        = dataArgs["gene2"].split('/')
                    (cds,geneNumber)      = geneCallFields[0].split('cds')
                    gene_obj.number       = geneNumber
                    gene_obj.parentGenome = dataArgs["genome2"]
                    gene_obj.contigName   = dataArgs["contig2"]
                    gene_obj.annotation   = dataArgs["annotations2"] 

                # Add hit
                if MUTUAL:
                    gene_obj.mutualBestHitList.append(gene1id) 
                    gene_obj.isLoner = False
                if SINGULAR_TWO:
                    gene_obj.singularBestHitList.append(gene1id)
                    gene_obj.isLoner = False
                if LONER_TWO:
                    gene_obj.lonerList.append(dataArgs["genome1"])
                if NEW:
                    genome2_obj.geneList.append(gene_obj)

            elif PROTEIN:
                PROTEIN_FOUND = False; NEW = False
                for protein_obj in genome2_obj.proteinList:
                    if protein_obj.contigName == dataArgs["contig2"] and protein_obj.name == dataArgs["protein2"]:
                        PROTEIN_FOUND = True
                        break
                if not PROTEIN_FOUND: 
                    # Create and populate a new gene object
                    NEW = True
                    protein_obj = copy.deepcopy(self.proteinTemplate)
                    protein_obj.name         = dataArgs["protein2"]
                    protein_obj.type         = "protein"
                    protein_obj.identifier   = protein2id 
                    protein_obj.cgpHeader    = dataArgs["protein2"]
                    proteinNameFields        = dataArgs["protein2"].split('/')
                    (cds,proteinNumber)      = proteinNameFields[0].split('cds')
                    protein_obj.number       = proteinNumber
                    protein_obj.parentGenome = dataArgs["genome2"]
                    protein_obj.contigName   = dataArgs["contig2"]
                    protein_obj.annotation   = dataArgs["annotations2"] 

                # Add hit
                if MUTUAL:
                    protein_obj.mutualBestHitList.append(protein1id) 
                    protein_obj.isLoner = False
                if SINGULAR_TWO:
                    protein_obj.singularBestHitList.append(protein1id)
                    protein_obj.isLoner = False
                if LONER_TWO:
                    protein_obj.lonerList.append(dataArgs["genome1"])
                if NEW:
                    genome2_obj.proteinList.append(protein_obj)
        return

    def findGenomeObject(self,genomeName):
        for genome_object in self.genomeList:
            if genome_object.name == genomeName:
                return genome_object
        return

    def addMutualBestHit(self,hitList,hit,item1,item2,item3,item4,item5):
        print("genomics_comparGenomes says, Adding hit ",hit," to mutualBestHitList:")
        print("other data: ",item1,' ',item2,' ',item3,' ',item4,' ',item5)
        for hitString in hitList:
            print("   ",hitString)
        hitList.append(hit)
        print(" yielding,")
        for hitString in hitList:
            print("   ",hitString)
        print('\n')
        return

    def addSingularBestHit(self,hitList,hit,item1,item2,item3,item4,item5):
        print("genomics_comparGenomes says, Adding hit ",hit," to singularBestHitList ",hitList)
        print("other data: ",item1,' ',item2,' ',item3,' ',item4,' ',item5)
        for hitString in hitList:
            print("   ",hitString)
        hitList.append(hit)
        print(" yielding,")
        for hitString in hitList:
            print("   ",hitString)
        print('\n')
        return

    def addLoner(self,lonerList,loner):
        print("genomics_compareGenomes says, Adding loner ",loner," to lonerList ",lonerList)
        lonerList.append(loner)
        print("resulting list: ",lonerList,'\n')
        return

    #===== COMPARISON GENOMIC METHODS

    def computeHomologyGroups(self):
        # Compute homology groups  for each reference gene/protein by combining homologous genes 
        # laterally (across genomes) and vertically (paralogs); save to reference genome object.
        refGeneList    = []  # Names of reference genes that have been combined into a group 
        refProteinList = []  # Names of reference proteins that have been combined into a group 
        for genome in self.genomeList:
            if genome.isReference:

                # Compute gene homology lists
                for gene in genome.geneList:
                    if gene.identifier in refGeneList: # Already processed this gene
                        continue 
                    else:
                        # Account for the gene itself
                        refGeneList.append(gene.identifier)
                        # Account for each mutual best hit
                        for homolog in gene.mutualBestHitList:
                            gene.homologyList.append(homolog)
                        # Account for each singular best hit
                        for homolog in gene.singularBestHitList:
                            gene.homologyList.append(homolog)
                        # Account for each paralog of the gene itself
                        for homolog in gene.paralogList:
                            gene.homologyList.append(homolog)
                            refGeneList.append(homolog) # record here so doesn't provoke another homology group 
                            # Account for each mutual and singular best hit of each paralog
                            paralog_obj = self.findGeneParalog(genome.name,homolog)
                            for gHomolog in paralog_obj.mutualBestHitList:
                                gene.homologyList.append(gHomolog)
                            for gHomolog in paralog_obj.singularBestHitList:
                                gene.homologyList.append(gHomolog)

                # Compute protein homology lists
                for protein in genome.proteinList:
                    if protein.identifier in refProteinList: # Already processed this gene
                        continue 
                    else:
                        # Account for the protein itself
                        refProteinList.append(protein.identifier)
                        # Account for each mutual best hit
                        for homolog in protein.mutualBestHitList:
                            protein.homologyList.append(homolog)
                        # Account for each singular best hit
                        for homolog in protein.singularBestHitList:
                            protein.homologyList.append(homolog)
                        # Account for each paralog of the protein itself
                        for homolog in protein.paralogList:
                            protein.homologyList.append(homolog)
                            refProteinList.append(homolog) # record here so doesn't provoke another homology group
                            # Account for each mutual and singular best hit of each paralog
                            paralog_obj = self.findProteinParalog(genome.name,homolog)
                            for pHomolog in paralog_obj.mutualBestHitList:
                                protein.homologyList.append(pHomolog)
                            for pHomolog in paralog_obj.singularBestHitList:
                                protein.homologyList.append(pHomolog)
        if PHATE_PROGRESS:
            print("genomics_compareGenomes says, Homology group computation complete.")
        return

    def findGeneParalog(self,genomeName,geneIdentifier):
        gene_obj = self.geneTemplate
        for genome in self.genomeList:
            if genome.name == genomeName:
                for gene_obj in genome.geneList:
                    if gene_obj.identifier == geneIdentifier:
                        return gene_obj
        return gene_obj

    def findProteinParalog(self,genomeName,proteinIdentifier):
        protein_obj = self.proteinTemplate
        for genome in self.genomeList:
            if genome.name == genomeName:
                for protein_obj in genome.proteinList:
                    if protein_obj.identifier == proteinIdentifier:
                        return protein_obj 
        return protein_obj 

    # Method createGeneHomologyFastaFiles creates homoloyg fasta files (for the reference genome)
    def createGeneHomologyFastaFiles(self,directory='./'):
        geneFasta   = ""
        groupNumber = 0
        for genome in self.genomeList:
            if genome.isReference:
                for gene in genome.geneList:

                    if gene.homologyList:
                        # Construct full path/filename string for next homology group fasta file name
                        groupNumber += 1
                        homFileName = HOMOLOGY_PREFIX + str(groupNumber) + '.fnt' 
                        nextHomologyFastaFile = os.path.join(directory,homFileName)
                        # Get fasta sequence for this reference protein
                        filePathName = self.getGeneFile(genome.name) 
                        geneSeq = self.getSequence(gene.name,filePathName)
                        FILE_H = open(nextHomologyFastaFile,'w')
                        fastaString = '>' + gene.identifier + '\n' + geneSeq + '\n'
                        FILE_H.write("%s" % (fastaString))

                        # Get fasta sequence for each member of this reference genes homology group 
                        for seqID in gene.homologyList:
                            (genomeName,contigName,cgpHeader) = seqID.split(':')
                            filePathName = self.getGeneFile(genomeName)
                            geneFasta = self.getSequence(cgpHeader,filePathName)
                            fastaString = '>' + seqID + '\n' + geneSeq + '\n'
                            FILE_H.write("%s" % (fastaString))
                        FILE_H.close()
        return geneFasta

    # Method createProteinHomologyFastaFiles creates not only homology fasta files, but also annotation files
    def createProteinHomologyFastaFiles(self,directory='./'):
        groupNumber  = 0
        # Identify reference genome and construct homology fasta files and homology annotation files
        for genome in self.genomeList:
            if genome.isReference:  # Reference genome is basis for homology groups
                for protein in genome.proteinList:
                    # Walk through protein list, gather fasta sequences for each homology group
                    if protein.homologyList:  # Remember that homology groups include paralogs and their correspondences

                        # Construct full path/filename string for next homology group fasta file name
                        groupNumber += 1
                        homFileName = HOMOLOGY_PREFIX + str(groupNumber) + '.faa'          
                        nextHomologyFastaFile = os.path.join(directory,homFileName) 

                        # Construct full path/filename string for next homology group annotation file name
                        annotFileName = HOMOLOGY_ANNOT_PREFIX + str(groupNumber) + '.annot'
                        nextHomologyAnnotFile = os.path.join(directory,annotFileName)

                        # Open fasta and annotation files
                        FASTA_FILE_H = open(nextHomologyFastaFile,'w')
                        ANNOT_FILE_H = open(nextHomologyAnnotFile,'w')

                        # Get fasta sequence for this reference protein and write to file
                        filePathName = self.getProteinFile(genome.name) 
                        proteinSeq = self.getSequence(protein.name,filePathName)
                        fastaString = '>' + protein.identifier + '\n' + proteinSeq + '\n'
                        FASTA_FILE_H.write("%s" % (fastaString))

                        # Get fasta sequence for each member of this reference protein's homology group and write to file
                        for seqID in protein.homologyList:
                            (genomeName,contigName,cgpHeader) = seqID.split(':')
                            filePathName = self.getProteinFile(genomeName)
                            proteinSeq = self.getSequence(cgpHeader,filePathName)
                            fastaString = '>' + seqID + '\n' + proteinSeq + '\n'
                            FASTA_FILE_H.write("%s" % (fastaString))

                        # Get annotation for this reference protein and write to file
                        fullAnnotation = ast.literal_eval(protein.annotation)
                        annotation = fullAnnotation[0][1:] 
                        ANNOT_FILE_H.write("%s\t%s\n" % (protein.identifier,annotation))

                        # Get annotation for each member of the reference protein's homology group and write to file
                        for seqID in protein.homologyList:
                            try:
                                fullAnnotation = ast.literal_eval(self.getAnnotation(seqID,"protein"))
                                annotation = fullAnnotation[0][1:]
                            except:
                                if PHATE_WARNINGS:
                                    print("genomics_compareGenomes says, WARNING: No annotation found for ",seqID)
                                annotation = "no annotation found"
                            ANNOT_FILE_H.write("%s\t%s\n" % (seqID,annotation))

                        # Clean up
                        FASTA_FILE_H.close()
                        ANNOT_FILE_H.close()

                        # Get annotations for this reference protein
        return

    def getGeneFile(self,genomeName):
        filePathName = ""
        cgpGenesFilename = genomeName + '_' + CGP_GENES_FILENAME
        filePathName = os.path.join(PHATE_PIPELINE_OUTPUT_DIR,genomeName,cgpGenesFilename)
        return filePathName 

    def getProteinFile(self,genomeName):
        filePathName = ""
        cgpProteinsFilename = genomeName + '_' + CGP_PROTEINS_FILENAME
        filePathName = os.path.join(PHATE_PIPELINE_OUTPUT_DIR,genomeName,cgpProteinsFilename)
        return filePathName 

    def getSequence(self,cgpHeader,filePathName):
        seq = ""
        try:
            file_h = open(filePathName,'r')
        except:
            print("genomics_compareGenomes says, ERROR: Cannot open file ",filePathName)
            return seq 
        fLines = file_h.read().splitlines()
        for i in range(0,len(fLines)):
            fLine = fLines[i]
            # Note: using the cgpHeader as a search string fails if it has '+' (for strand)
            # Thus, need to substitute the '+' in order to find the original string in file
            new_cgpHeader = re.sub('\+','\\+',cgpHeader)
            match_header = re.search(new_cgpHeader,fLine)
            if match_header:
                seq = fLines[i+1]
                break
        if seq == "":
            if PHATE_WARNINGS:
                print("genomics_compareGenomes says, WARNING: Failed to find sequence for cgpHeader ",cgpHeader)
        file_h.close()
        return seq 

    def getAnnotation(self,identifier,seqType):
        annot = ""
        (genomeName,contigName,cgpHeader) = identifier.split(':')
        for genome in self.genomeList:
            # Identify the genome in question
            if genome.name == genomeName:
                # Identify the gene or protein, and capture it's annotation
                if seqType.lower() == "gene" or seqType.lower() == "nucleotide" or seqType.lower() == "nt":
                    for gene in genome.geneList:
                        if gene.cgpHeader == cgpHeader:
                            annot = gene.annotation 
                elif seqType.lower() == "protein" or seqType.lower() == "peptide" or seqType.lower() == "aa" or seqType.lower() == "prot":
                    for protein in genome.proteinList:
                        if protein.cgpHeader == cgpHeader:
                            annot = protein.annotation
        return annot

    #===== COMPARISON DATA CHECK METHODS

    def runDataChecks(self):
        self.countGenomes()
        self.checkUnique()
        self.checkMutualBestHitLists()
        self.checkSingularBestHitLists()
        return

    def countGenomes(self):
        self.genomeCount = len(self.genomeList)
        if PHATE_MESSAGES:
            print("genomics_compareGenomes says, There are",self.genomeCount,"genomes.")
        return

    def checkMutualBestHitLists(self):
        genomeCount = len(self.genomeList)
        for genome in self.genomeList:
            genome.checkMutualBestHitList(genomeCount)
        return

    def checkSingularBestHitLists(self):
        genomeCount = len(self.genomeList)
        for genome in self.genomeList:
            genome.checkSingularBestHitList(genomeCount)
        return

    def checkUnique(self):
        tempList = []
        for genome in self.genomeList:
            if genome.name in tempList:
                if PHATE_WARNINGS:
                    print("genomics_compareGenomes genomes obj says, WARNING: ",genome.name," occurs more than once in genomes object.")
            else:
                tempList.append(genome.name)
        return

    #===== COMPARISON PRINT METHODS

    def writeCorrespondences2file(self,FILE_H):
        for genome in self.genomeList:
            if genome.isReference:
                FILE_H.write("%s\n" % ("GENE CORRESPONDENCES"))
                # For each reference gene, print non-redundant list of mutual or singular best hit wrt each other genome
                for gene in genome.geneList:
                    correspondenceList = []
                    FILE_H.write("%s%s\n" % ("Genes corresponding to",gene.identifier))
                    for hit in gene.mutualBestHitList:
                        if hit not in correspondenceList:
                            correspondenceList.append(hit)
                    for hit in gene.singularBestHitList:
                        if hit not in correspondenceList:
                            correspondenceList.append(hit)
                    for hit in correspondenceList:
                        FILE_H.write("%s%s\n" % ("  ",hit))
                FILE_H.write("%s\n" % ("PROTEIN CORRESPONDENCES"))
                # For each reference protein, print non-redundant list of mutual or singular best hit wrt each other genome
                for protein in genome.proteinList:
                    correspondenceList = []
                    FILE_H.write("%s%s\n" % ("Proteins corresponding to",protein.identifier))
                    for hit in protein.mutualBestHitList:
                        if hit not in correspondenceList:
                            correspondenceList.append(hit)
                    for hit in protein.singularBestHitList:
                        if hit not in correspondenceList:
                            correspondenceList.append(hit)
                    for hit in correspondenceList:
                        FILE_H.write("%s%s\n" % ("  ",hit))
        return

    def writeCorrespondences(self):
        for genome in self.genomeList:
            if genome.isReference:
                print("********** GENE CORRESPONDENCES **********")
                # For each reference gene, print non-redundant list of mutual or singular best hit wrt each other genome
                for gene in genome.geneList:
                    correspondenceList = []
                    print("Genes corresponding to",gene.identifier)
                    for hit in gene.mutualBestHitList:
                        if hit not in correspondenceList:
                            correspondenceList.append(hit)
                    for hit in gene.singularBestHitList:
                        if hit not in correspondenceList:
                            correspondenceList.append(hit)
                    for hit in correspondenceList:
                        print("  ",hit)
                print("********** End Gene Correspondences")
                print("********** PROTEIN CORRESPONDENCES **********")
                # For each reference protein, print non-redundant list of mutual or singular best hit wrt each other genome
                for protein in genome.proteinList:
                    correspondenceList = []
                    print("Genes corresponding to",protein.identifier)
                    for hit in protein.mutualBestHitList:
                        if hit not in correspondenceList:
                            correspondenceList.append(hit)
                    for hit in protein.singularBestHitList:
                        if hit not in correspondenceList:
                            correspondenceList.append(hit)
                    for hit in correspondenceList:
                        print("  ",hit)
                print("********** End Protein Correspondences")
        return

    def writeCoreGenome2file(self,FILE_H):
        genomeCount = len(self.genomeList)
        for genome in self.genomeList:
            if genome.isReference:
                count = 0
                FILE_H.write("%s\n" % ("***** CORE GENOME: GENE"))
                # For each reference gene that matches genes in all other genomes, list the gene set
                for gene in genome.geneList:
                    if len(gene.mutualBestHitList) == genomeCount-1:
                        count += 1
                        FILE_H.write("%s%s%s\n" % ("Gene Set #",count,':'))
                        FILE_H.write("%s\n" % (gene.identifier))
                        for hit in gene.mutualBestHitList:
                            FILE_H.write("%s\n" % (hit))
                count = 0
                FILE_H.write("%s\n" % ("***** CORE GENOME: PROTEIN"))
                # For each reference protein that matches proteins in all other genomes, list the protein set
                for protein in genome.proteinList:
                    if len(protein.mutualBestHitList) == genomeCount-1:
                        count += 1
                        FILE_H.write("%s%s%s\n" % ("Protein Set #",count,':'))
                        FILE_H.write("%s\n" % (protein.identifier))
                        for hit in protein.mutualBestHitList:
                            FILE_H.write("%s\n" % (hit))
        return

    def writeCoreGenome(self):
        genomeCount = len(self.genomeList)
        count = 0
        for genome in self.genomeList:
            if genome.isReference:
                print("********** CORE GENOME: GENE **********")
                # For each reference gene that matches genes in all other genomes, list the gene set
                for gene in genome.geneList:
                    if len(gene.mutualBestHitList) == genomeCount-1:
                        count += 1
                        print("Gene Set #",count,':')
                        print(gene.identifier)
                        for hit in gene.mutualBestHitList:
                            print(hit)
                print("********** End Core Genome: GENE **********")
                print("********** CORE GENOME: PROTEIN **********")
                # For each reference protein that matches protein in all other genomes, list the protein set
                for protein in genome.proteinList:
                    if len(protein.mutualBestHitList) == genomeCount-1:
                        count += 1
                        print("Protein Set #",count,':')
                        print(protein.identifier)
                        for hit in protein.mutualBestHitList:
                            print(hit)
                print("********** End Core Genome: PROTEIN **********")
        return

    def writeMutualBestHitList2file(self,FILE_H):
        for genome in self.genomeList:
            count = 0
            #if genome.isReference:
            FILE_H.write("%s\n" % ("*******************************"))
            FILE_H.write("%s\n" % (genome.name))
            FILE_H.write("%s\n" % ("*** Mutual Best Hits: Gene"))
            for gene in genome.geneList:
                count += 1
                FILE_H.write("%s%s%s\n" % (count,") gene ",gene.identifier))
                gene.writeMutualBestHitList2file(FILE_H)
            FILE_H.write("%s\n" % ("*** Mutual Best Hits: Protein"))
            count = 0
            for protein in genome.geneList:
                count += 1
                FILE_H.write("%s%s%s\n" % (count,") protein ",protein.identifier))
                protein.writeMutualBestHitList2file(FILE_H)
        return

    def writeMutualBestHitList(self):
        for genome in self.genomeList:
            if genome.isReference:
                print("************** ",genome.name," **************")
                print("************** Mutual Best Hits: Gene ")
                for gene in genome.geneList:
                    print("** gene ",gene.identifier)
                    gene.writeMutualBestHitList()
                print("************** End of Mutual Best Hits List: Gene ")
                print("************** Mutual Best Hits: Protein ")
                for protein in genome.proteinList:
                    print("** protein ",protein.identifier)
                    protein.writeMutualBestHitList()
                print("************** End of Mutual Best Hits List: Protein ")
        return

    def writeSingularBestHitList2file(self,FILE_H):
        for genome in self.genomeList:
            FILE_H.write("%s\n" % ("*************************"))
            FILE_H.write("%s%s\n" % ("*** Singular Best Hits for Genome, ",genome.name))
            FILE_H.write("%s\n" % ("*** Singular Best Hits: Gene"))
            count = 0
            for gene in genome.geneList:
                count += 1
                FILE_H.write("%s%s%s\n" % (count,") gene ",gene.identifier))
                gene.writeSingularBestHitList2file(FILE_H)
            FILE_H.write("%s\n" % ("*** Singular Best Hits: Protein"))
            count = 0
            for protein in genome.proteinList:
                count += 1
                FILE_H.write("%s%s%s\n" % (count,") protein ",protein.identifier))
                protein.writeSingularBestHitList2file(FILE_H)
        return

    def writeSingularBestHitList(self):
        for genome in self.genomeList:
            print("************** Singular Best Hits for Genome, ",genome.name," **************")
            print("************** Singular Best Hits ")
            for gene in genome.geneList:
                print("** gene ",gene.identifier)
                gene.writeSingularBestHitList()
            print("************** End of Singular Best Hits List ")
        return

    def writeGeneCorrespondences2file(self,FILE_H):
        for genome in self.genomeList:
            if genome.isReference:
                FILE_H.write("%s%s\n" % ("***** Gene and Protein Correspondences for reference genome, ",genome.name))
                runningList = []
                for gene in self.geneList:
                    for hit in self.mutualBestHitList_gene:
                        if hit not in runningList:
                            runningList.append(hit)
                    for hit in self.singularBestHitList_gene:
                        if hit not in runningList:
                            runningList.append(hit)
                    FILE_H.write("%s%s%s\n" % ("genes corresponding to",gene.identifier,':'))
                    for hit in runningList:
                        FILE_H.write("%s%s\n" % ("     ",hit))
                runningList = []
                for protein in self.proteinList:
                    for hit in self.mutualBestHitList_protein:
                        if hit not in runningList:
                            runningList.append(hit)
                    for hit in self.singularBestHitList_protein:
                        if hit not in runningList:
                            runningList.append(hit)
                    FILE_H.write("%s%s%s\n" % ("proteins corresponding to",protein.identifier,':'))
                    for hit in runningList:
                        FILE_H.write("%s%s\n" % ("     ",hit))
        return

    def writeGeneCorrespondences(self):
        for genome in self.genomeList:
            if genome.isReference:
                runningList = []
                print("********** Gene and Protein Correspondences for genome, ",genome.name," **********")
                for gene in self.geneList:
                    for hit in self.mutualBestHitList_gene:
                        if hit not in runningList:
                            runningList.append(hit)
                    for hit in self.singularBestHitList_gene:
                        if hit not in runningList:
                            runningList.append(hit)
                    print("genes corresponding to",gene.identifier,':')
                    for hit in runningList:
                        print("     ",hit)
                print("********** End of gene correspondences ")
                runningList = []
                for protein in self.proteinList:
                    for hit in self.mutualBestHitList_protein:
                        if hit not in runningList:
                            runningList.append(hit)
                    for hit in self.singularBestHitList_protein:
                        if hit not in runningList:
                            runningList.append(hit)
                    print("proteins corresponding to",protein.identifier,':')
                    for hit in runningList:
                        print("     ",hit)
                print("********** End of protein correspondences ")
        return

    def writeLonerList2file(self,FILE_H):
        for genome in self.genomeList:
            FILE_H.write("%s%s\n" % ("***** ",genome.name))
            FILE_H.write("%s\n" % ("*** Gene Loners"))
            for gene in genome.geneList:
                if gene.isLoner:
                    FILE_H.write("%s\n" % (gene.identifier))
            FILE_H.write("%s\n" % ("*** Protein Loners"))
            for protein in genome.proteinList:
                if protein.isLoner:
                    FILE_H.write("%s\n" % (protein.identifier))
        return

    def writeLonerList(self):
        for genome in self.genomeList:
            print("********** ",genome.name," **********")
            print("********** Gene Loners ")
            for gene in genome.geneList:
                if gene.isLoner:
                    print(gene.identifier)
            print("********** End of gene loner list ")
            print("********** Protein Loners ")
            for protein in genome.proteinList:
                if protein.isLoner:
                    print(protein.identifier)
            print("********** End of protein loner list ")
        return

    def writeParalogs2file(self,FILE_H):
        for genome in self.genomeList:
            FILE_H.write("%s%s\n" % ("***** ",genome.name))
            FILE_H.write("%s\n" % ("*** Paralogs"))
            for gene in genome.geneList:
                if gene.paralogList:
                    FILE_H.write("%s%s%s\n" % ("Paralog(s) for gene ",gene.identifier,":"))
                    for paralog in gene.paralogList:
                        FILE_H.write("%s\n" % (paralog)) 
            for protein in genome.proteinList:
                if protein.paralogList:
                    FILE_H.write("%s%s%s\n" % ("Paralog(s) for protein ",protein.identifier,":"))
                    for paralog in protein.paralogList:
                        FILE_H.write("%s\n" % (paralog))
        return

    def writeParalogs(self):
        for genome in self.genomeList:
            print("********** ",genome.name," **********")
            print("********** Paralogs ")
            for gene in genome.geneList:
                if gene.paralogList:
                    print("Paralog(s) for gene ",gene.identifier,":")
                    for paralog in gene.paralogList:
                        print(paralog) 
                else:
                    print("No paralogs for this gene")
            for protein in genome.proteinList:
                if protein.paralogList:
                    print("Paralog(s) for protein ",protein.identifier,":")
                    for paralog in protein.paralogList:
                        print(paralog)
                else:
                    print("No paralogs for this protein")
            print("********** End of paralogs list ")
        return

    def writeHomologyGroups2file(self,FILE_H):
        count = 0
        for genome in self.genomeList:
            if genome.isReference:
                FILE_H.write("%s%s\n" % ("***** Homology Groups for Reference Genome, ",genome.name))
                FILE_H.write("%s\n" % ("*** Gene Homology Groups"))
                for gene in genome.geneList:
                    if gene.homologyList:
                        count += 1
                        FILE_H.write("%s%s%s%s%s\n" % ("Homology group No. ",count," for gene ",gene.identifier,":"))
                        FILE_H.write("%s%s\n" % ("  ",gene.homologyList))
                FILE_H.write("%s\n" % ("*** Protein Homology Groups"))
                for protein in genome.proteinList:
                    if protein.homologyList:
                        count += 1
                        FILE_H.write("%s%s%s%s%s\n" % ("Homology group No. ",count," for protein ",protein.identifier,":"))
                        FILE_H.write("%s%s\n" % ("  ",protein.homologyList))
        return

    def writeHomologyGroups(self):
        count = 0
        for genome in self.genomeList:
            if genome.isReference:
                print("********** Homology Groups for Reference Genome, ",genome.name," **********")
                print("********** Gene Homology Groups ")
                for gene in genome.geneList:
                    if gene.homologyList:
                        count += 1
                        print("Homology group No. ",count," for gene ",gene.identifier,":")
                        print("  ",gene.homologyList)
                print("********** End of Gene Homology Groups ")
                for protein in genome.proteinList:
                    if protein.homologyList:
                        count += 1
                        print("Homology group No. ",count," for protein ",protein.identifier,":")
                        print("  ",protein.homologyList)
                print("********** End of Protein Homology Groups ")
        return

    def printReports2files(self):
        FILE_H = open(MUTUAL_BEST_HIT_FILE,"w")
        self.writeMutualBestHitList2file(FILE_H)
        FILE_H.close()

        FILE_H = open(SINGULAR_BEST_HIT_FILE,"w")
        self.writeSingularBestHitList2file(FILE_H)
        FILE_H.close()

        FILE_H = open(LONER_FILE,"w")
        self.writeLonerList2file(FILE_H)
        FILE_H.close()

        FILE_H = open(CORE_GENOME_FILE,"w")
        self.writeCoreGenome2file(FILE_H)
        FILE_H.close()

        FILE_H = open(PARALOG_FILE,"w")
        self.writeParalogs2file(FILE_H)
        FILE_H.close()

        FILE_H = open(CORRESPONDENCE_FILE,"w")
        self.writeCorrespondences2file(FILE_H)
        FILE_H.close()

        FILE_H = open(HOMOLOGY_GROUP_FILE,"w")
        self.writeHomologyGroups2file(FILE_H)
        FILE_H.close()

        return

    def printReport2file(self,FILE_H):
        self.writeMutualBestHitList2file(FILE_H)
        self.writeSingularBestHitList2file(FILE_H)
        self.writeLonerList2file(FILE_H)
        self.writeCoreGenome2file(FILE_H)
        self.writeParalogs2file(FILE_H)
        self.writeCorrespondences2file(FILE_H)
        self.writeHomologyGroups2file(FILE_H)
        return

    def printReport(self):
        self.writeMutualBestHitList()
        self.writeSingularBestHitList()
        self.writeLonerList()
        self.writeCoreGenome()
        self.writeParalogs()
        self.writeCorrespondences()
        self.writeHomologyGroups()
        return

    def printAll2file(self,FILE_H):
        FILE_H.write("%s\n" % ("=========Genome Set============="))
        FILE_H.write("%s%s\n" % ("name:",self.name))
        FILE_H.write("%s%s\n" % ("referenceGenome:",self.referenceGenome))
        FILE_H.write("%s%s\n" % ("genomeCount:",self.genomeCount))
        FILE_H.write("%s\n" % ("Set of genomes:"))
        for genome in self.genomeList:
            genome.printAll2file(FILE_H)
        FILE_H.write("%s\n" % ("=========End Genome Set========="))
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

#############################################################################################################
# Class genome stores metadata 
#############################################################################################################
class genome(object):

    def __init__(self):
        self.name                 = ""     # Name of this genome (e.g., Lambda)
        self.species              = ""     # Species name for this genome (e.g., Ecoli_Lambda_phage)
        self.isReference          = False  # True if designated a reference genome
        self.contigList           = []     # Set of contig names (fasta headers)
        self.geneList             = []     # List of gene_protein objects 
        self.proteinList          = []     # List of gene_protein objects 
        self.paralogList          = []     # List of paralogSet objects

    def addParalog(self,dataArgs):
        hitType = ""; gene1 = ""; gene2 = ""; protein1 = ""; protein2 = ""
        geneCall1 = ""; geneCall2 = ""

        # Read input parameters; not all these are being used just yet
        if "hitType" in dataArgs.keys():
            hitType = dataArgs["hitType"]
        else:
            if PHATE_WARNINGS:
                print("genomics_compareGenomes says, WARNING: hitType not defined")
                return
        match_gene    = re.search('gene',   hitType)
        match_protein = re.search('protein',hitType)
        if match_gene:
            if "gene1" in dataArgs.keys():
                gene1 = dataArgs["gene1"]
            if "gene2" in dataArgs.keys():
                gene2 = dataArgs["gene2"]
        elif match_protein:
            if "protein1" in dataArgs.keys():
                protein1 = dataArgs["protein1"]
            if "protein2" in dataArgs.keys():
                protein2 = dataArgs["protein2"]
        if "contig1" in dataArgs.keys():
            contig1 = dataArgs["contig1"]
        if "contig2" in dataArgs.keys():
            contig2 = dataArgs["contig2"]
        if "geneCall1" in dataArgs.keys():
            geneCall1 = dataArgs["geneCall1"]
        if "geneCall2" in dataArgs.keys():
            geneCall2 = dataArgs["geneCall2"]
        if "annotation1" in dataArgs.keys():
            annotation1 = dataArgs["annotation1"]
        if "annotation2" in dataArgs.keys():
            annotation2 = dataArgs["annotation2"]

        # Find gene1 and gene2 in genome; get gene identifiers
        if hitType == "gene_paralog":
            for gene_obj1 in self.geneList:
                if gene_obj1.cgpHeader == gene1:
                    for gene_obj2 in self.geneList:
                        if gene_obj2.cgpHeader == gene2:
                            paralogID = gene_obj2.identifier
                            if paralogID not in gene_obj1.paralogList:
                                gene_obj1.paralogList.append(paralogID)

        # Find protein1 and protein2 in genome; get protein identifiers
        elif hitType == "protein_paralog":
            for protein_obj1 in self.proteinList:
                if protein_obj1.cgpHeader == protein1:
                    for protein_obj2 in self.proteinList:
                        if protein_obj2.cgpHeader == protein2:
                            paralogID = protein_obj2.identifier
                            if paralogID not in protein_obj1.paralogList:
                                protein_obj1.paralogList.append(paralogID)
        return

    #===== GENOME DATA CHECK METHODS

    def checkMutualBestHitList(self,genomeCount):
        for gene in self.geneList:
            if len(gene.mutualBestHitList) > genomeCount-1:
                if PHATE_WARNINGS:
                    print("genomics_compareGenomes says, WARNING: mutualBestHitList for gene",gene.identifier,"incorrect")
                    print("   Length is ",len(gene.mutualBestHitList))
                    print("    ",gene.mutualBestHitList)
            if gene.identifier in gene.mutualBestHitList:
                if PHATE_WARNINGS:
                    print("genomics_compareGenomes says, WARNING: gene,",gene.identifier," is listed in its own mutualBestHitList.")

        for protein in self.proteinList:
            if len(protein.mutualBestHitList) > genomeCount-1:
                if PHATE_WARNINGS:
                    print("genomics_compareGenomes says, WARNING: mutualBestHitList for protein",protein.identifier,"incorrect")
                    print("   Length is ",len(protein.mutualBestHitList))
                    print("    ",protein.mutualBestHitList)
            if protein.identifier in protein.mutualBestHitList:
                if PHATE_WARNINGS:
                    print("genomics_compareGenomes says, WARNING: protein,",protein.identifier," is listed in its own mutualBestHitList.")
        return

    def checkSingularBestHitList(self,genomeCount):
        for gene in self.geneList:
            if gene.identifier in gene.singularBestHitList:
                if PHATE_WARNINGS:
                    print("genomics_compareGenomes says, WARNING: gene,",gene.identifier," is listed in its own singularBestHitList.")
        for protein in self.proteinList:
            if protein.identifier in protein.singularBestHitList:
                if PHATE_WARNINGS:
                    print("genomics_compareGenomes says, WARNING: protein,",protein.identifier," is listed in its own singularBestHitList.")
        return

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

    #===== GENOME PRINT METHODS

    def x_printReport2file(self,FILE_H):
        return

    def x_printReport(self):
        return

    def printAll2file(self,FILE_H):
        FILE_H.write("%s\n" % ("=========Genome========"))
        FILE_H.write("%s%s\n" % ("name:",self.name))
        FILE_H.write("%s%s\n" % ("species:",self.species))
        FILE_H.write("%s%s\n" % ("isReference:",self.isReference))
        FILE_H.write("%s\n" % ("contigs:"))
        if self.contigList:
            for contig in self.contigList:
                FILE_H.write("%s\n" % (contig))
        else:
            FILE_H.write("%s\n" % ("There are no contigs."))
        if self.geneList:
            for gene in self.geneList:
                gene.printAll2file(FILE_H)
        else:
            FILE_H.write("%s\n" % ("There are no genes."))
        if self.proteinList:
            for protein in self.proteinList:
                protein.printAll2file(FILE_H)
        else:
            FILE_H.write("%s\n" % ("There are no proteins."))
        if self.paralogList:
            for paralogSet in self.paralogList:
                paralogSet.printAll2file(FILE_H)
        else:
            FILE_H.write("%s\n" % ("There are no paralogs"))
        FILE_H.write("%s\n" % ("=========End Genome===="))
        return

    def printAll(self):
        print("=========Genome========")
        print("name:",self.name)
        print("species:",self.species)
        print("isReference:",self.isReference)
        print("contigs:",)
        if self.contigList:
            for contig in self.contigList:
                print(contig)
        else:
            print("There are no contigs.")
        if self.geneList:
            for gene in self.geneList:
                gene.printAll()
        else:
            print("There are no genes.")
        if self.proteinList:
            for protein in self.proteinList:
                protein.printAll()
        else:
            print("There are no proteins.")
        if self.paralogList:
            for paralogSet in self.paralogList:
                paralogSet.printAll()
        else:
            print("There are no paralogs")
        print("=========End Genome====")
        return

#####################################################################################################
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

    #===== PARALOG PRINT METHODS

    def x_printReport2file(self,FILE_H):
        return

    def x_printReport(self):
        return

    def printAll2file(self,FILE_H):
        FILE_H.write("%s%s\n" % ("paralogType:",self.paralogType))
        FILE_H.write("%s%s\n" % ("setSize:",self.setSize))
        FILE_H.write("%s%s\n" % ("paralogList:",self.paralogList))
        return

    def printAll(self):
        print("paralogType:",self.paralogType)
        print("setSize:",self.setSize)
        print("paralogList:",self.paralogList)
        return

######################################################################################################
# Class gene_protein stores meta-data about a gene or protein.
class gene_protein(object):

    def __init__(self):
        self.type                 = "gene"    # "gene" or "protein"
        self.name                 = ""        # Ex: "phanotate_5"
        self.identifier           = ""        # Unique: genome + contig + name
        self.cgpHeader            = ""        # CGP-assigned header; Ex: cgp5/+/72/485/
        self.number               = 0         # Unique: Set to correspond to protein 
        self.parentGenome         = ""        # Ex: "Lambda" - corresponds to a genome object's name
        self.contigName           = ""        # Ex: "Lambda_contig_1"; read from CGP report
        self.annotation           = ""        # The full annotation string, read from CGP report 
        self.isLoner              = True      # True if this gene has no correspondences
        self.mutualBestHitList    = []        # list of gene identifiers that are mutual best hits, across genomes
        self.singularBestHitList  = []        # list of gene identifiers that are best hits, relative to this gene, across genomes
        self.correspondenceList   = []        # List of corresponding genes (mutual or best-hit) across genomes: list of gene identifiers 
        self.homologyList         = []        # List of corresponding genes + paralogs and their corresponding genes
        self.lonerList            = []        # List of the genome names that this gene/protein is a loner with respect to
        self.groupList            = []        # list of all corresponding genes plus paralogs
        self.paralogList          = []        # List of paralogs within its own parent genome: list of gene identifiers <data is redundant; should pull paralog list as list of lists at genome level.

    def addMutualBestHit(self,hit):
        self.mutualBestHitList.append(hit)
        return

    def addSingularBestHit(self,hit):
        self.singularBestHitList.append(hit)
        return

    def addGroupMember(self,member):
        self.groupList.append(member)

    # Construct a set of genes/proteins that are related by homology
    def x_constructGroup(self):
        # self.groupList

        if PHATE_MESSAGES:
            print("genomics_compareGenomes says, Group constructed for gene_protein",self.identifier,";",self.type)
        return

    #===== GENE_PROTEIN DATA CHECK METHODS

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

    #===== GENE_PROTEIN PRINT METHODS

    def writeMutualBestHitList2file(self,FILE_H):
        for hit in self.mutualBestHitList:
            FILE_H.write("%s%s\n" % ('   ',hit))
        return

    def writeMutualBestHitList(self):
        for hit in self.mutualBestHitList:
            print(hit)
        return

    def writeSingularBestHitList2file(self,FILE_H):
        for hit in self.singularBestHitList:
            FILE_H.write("%s%s\n" % ('   ',hit))
        return

    def writeSingularBestHitList(self):
        for hit in self.singularBestHitList:
            print(hit)
        return

    def writeLonerList2file(self,FILE_H):
        for genomeString in self.lonerList:
            FILE_H.write("%s\n" % (genomeString))
        return

    def writeLonerList(self):
        for genomeString in self.lonerList:
            print(genomeString)
        return

    def printReport2file(self,FILE_H):
        return

    def printReport(self):
        return

    def printAll2file(self,FILE_H):
        FILE_H.write("%s%s\n" % ("===type:",self.type))
        FILE_H.write("%s%s\n" % ("name:",self.name))
        FILE_H.write("%s%s\n" % ("identifier:",self.identifier))
        FILE_H.write("%s%s\n" % ("cgpHeader:",self.cgpHeader))
        FILE_H.write("%s%s\n" % ("number:",self.number))
        FILE_H.write("%s%s\n" % ("parentGenome:",self.parentGenome))
        FILE_H.write("%s%s\n" % ("contigName:",self.contigName))
        FILE_H.write("%s%s\n" % ("annotation:",self.annotation))
        FILE_H.write("%s%s\n" % ("isLoner:",self.isLoner))
        FILE_H.write("%s%s\n" % ("mutualBestHitList:",self.mutualBestHitList))
        FILE_H.write("%s%s\n" % ("singularBestHitList:",self.singularBestHitList))
        FILE_H.write("%s%s\n" % ("correspondenceList:",self.correspondenceList))
        FILE_H.write("%s%s\n" % ("paralogList:",self.paralogList))
        FILE_H.write("%s\n" % ("==="))
        return

    def printAll(self):
        print("===type:",self.type)
        print("name:",self.name)
        print("identifier:",self.identifier)
        print("cgpHeader:",self.cgpHeader)
        print("number:",self.number)
        print("parentGenome:",self.parentGenome)
        print("contigName:",self.contigName)
        print("annotation:",self.annotation)
        print("isLoner:",self.isLoner)
        print("mutualBestHitList:",self.mutualBestHitList)
        print("singularBestHitList:",self.singularBestHitList)
        print("correspondenceList:",self.correspondenceList)
        print("paralogList:",self.paralogList)
        print("===")
        return

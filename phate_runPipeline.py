#!/usr/bin/env python

################################################################
#
# Program Title:  phate_runPipeline.py ()
#
# Most recent update:  08 August 2020
#
# Description: Runs the phate annotation pipeline.  This code runs under Python 3.7, and requires
#    dependent packages.
#
# Input: a json string, formatted by multiPhATE2.py
#
# Usage:  python phate_runPipeline.py myGenome.json
#
# Setup:
#     CompareCalls/         - code for comparing gene caller results
#     DatabasePrep/         - code for preparing custom databases
#     GeneCalling/          - mini-pipeline runs gene-call programs
#     SequenceAnnotation/   - PhATE sequence annotation codes
#     phate_runPipeline.py  - pipeline driver
#     PipelineInput/        - contains myGenome.fasta (and optionally, myGenome.psat)
#     PipelineOutput/       - output files are written here to a subdirectory specified in config file
#     myPhATE.config        - configuration file (copy/modify sample.config)
#
# Programmer's Notes:
#    This code uses a running log; need to occasionally clean it out
#
# Programmers: 
#    Carol E. Zhou - pipeline programmer: CompareCalls/, DatabasePrep/, SequenceAnnotation/, phate_runPipeline.py
#    Jeff Kimbrel, Carol Zhou  - GeneCalling/
#
################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.PDF FOR DETAILS.

import sys, os, re, string, copy, time, datetime
import subprocess
import json

# CONSTANTS and CONFIGURABLES

seqProgramSet = {'blastp','phmmer','jackhmmer'}  # codes for searches against sequence/blast databases
hmmProgramSet = {'hmmscan'}                      # codes for hmm searches against hmm profile databases

# Output file names 
GENE_FILE              = 'gene.fnt'                              #
PROTEIN_FILE           = 'protein.faa'                           #

##### ENVIRONMENT VARIABLES - These are set by multiPhate.py

BASE_DIR                      = os.environ["PHATE_BASE_DIR"]
DATABASE_DIR                  = os.environ["PHATE_DATABASE_DIR"]
SOFTWARE_DIR                  = os.environ["PHATE_SOFTWARE_DIR"]
PIPELINE_INPUT_DIR            = os.environ["PHATE_PIPELINE_INPUT_DIR"]

EMBOSS_HOME                   = os.environ["PHATE_EMBOSS_PHATE_HOME"] 
PIPELINE_DIR                  = os.environ["PHATE_PIPELINE_DIR"] 

# Data sets
NCBI_VIRUS_BASE_DIR           = os.environ["PHATE_NCBI_VIRUS_BASE_DIR"] 
NCBI_VIRUS_GENOME_BLAST_HOME  = os.environ["PHATE_NCBI_VIRUS_GENOME_BLAST_HOME"]
NCBI_VIRUS_PROTEIN_BLAST_HOME = os.environ["PHATE_NCBI_VIRUS_PROTEIN_BLAST_HOME"]
NCBI_TAXON_DIR                = os.environ["PHATE_NCBI_TAXON_DIR"] 
REFSEQ_PROTEIN_BASE_DIR       = os.environ["PHATE_REFSEQ_PROTEIN_BASE_DIR"]
REFSEQ_PROTEIN_BLAST_HOME     = os.environ["PHATE_REFSEQ_PROTEIN_BLAST_HOME"]
REFSEQ_GENE_BASE_DIR          = os.environ["PHATE_REFSEQ_GENE_BASE_DIR"]
REFSEQ_GENE_BLAST_HOME        = os.environ["PHATE_REFSEQ_GENE_BLAST_HOME"]
PVOGS_BASE_DIR                = os.environ["PHATE_PVOGS_BASE_DIR"]
PVOGS_BLAST_HOME              = os.environ["PHATE_PVOGS_BLAST_HOME"]
VOGS_BASE_DIR                 = os.environ["PHATE_VOGS_BASE_DIR"]
VOGS_BLAST_HOME               = os.environ["PHATE_VOGS_BLAST_HOME"]
PHANTOME_BASE_DIR             = os.environ["PHATE_PHANTOME_BASE_DIR"]
PHANTOME_BLAST_HOME           = os.environ["PHATE_PHANTOME_BLAST_HOME"]
PHAGE_ENZYME_BASE_DIR         = os.environ["PHATE_PHAGE_ENZYME_BASE_DIR"] 
PHAGE_ENZYME_BLAST_HOME       = os.environ["PHATE_PHAGE_ENZYME_BLAST_HOME"]
KEGG_VIRUS_BASE_DIR           = os.environ["PHATE_KEGG_VIRUS_BASE_DIR"] 
KEGG_VIRUS_BLAST_HOME         = os.environ["PHATE_KEGG_VIRUS_BLAST_HOME"]
PFAM_BASE_DIR                 = os.environ["PHATE_PFAM_BASE_DIR"]
PFAM_BLAST_HOME               = os.environ["PHATE_PFAM_BLAST_HOME"]
SMART_BASE_DIR                = os.environ["PHATE_SWISSPROT_BASE_DIR"]
SMART_BLAST_HOME              = os.environ["PHATE_SWISSPROT_BLAST_HOME"]
SWISSPROT_BASE_DIR            = os.environ["PHATE_SWISSPROT_BASE_DIR"]
SWISSPROT_BLAST_HOME          = os.environ["PHATE_SWISSPROT_BLAST_HOME"]
UNIPROT_BASE_DIR              = os.environ["PHATE_UNIPROT_BASE_DIR"]
UNIPROT_BLAST_HOME            = os.environ["PHATE_UNIPROT_BLAST_HOME"]
NR_BLAST_BASE_DIR             = os.environ["PHATE_NR_BLAST_BASE_DIR"]
NR_BLAST_HOME                 = os.environ["PHATE_NR_BLAST_HOME"]

# Gene calling
PRODIGAL_PATH                 = os.environ["PHATE_PRODIGAL_PATH"]
GLIMMER_PATH                  = os.environ["PHATE_GLIMMER_PATH"]
GENEMARKS_PATH                = os.environ["PHATE_GENEMARKS_PATH"]  
PHANOTATE_PATH                = os.environ["PHATE_PHANOTATE_PATH"]
CGC_PATH                      = os.environ["PHATE_CGC_PATH"]

# Blast
BLAST_HOME                    = os.environ["PHATE_BLAST_HOME"]
MIN_BLASTP_IDENTITY           = os.environ["PHATE_MIN_BLASTP_IDENTITY"]
MIN_BLASTN_IDENTITY           = os.environ["PHATE_MIN_BLASTN_IDENTITY"]
MAX_BLASTP_HIT_COUNT          = os.environ["PHATE_MAX_BLASTP_HIT_COUNT"] 
MAX_BLASTN_HIT_COUNT          = os.environ["PHATE_MAX_BLASTN_HIT_COUNT"] 
BLASTP_IDENTITY_DEFAULT       = os.environ["PHATE_BLASTP_IDENTITY_DEFAULT"]
BLASTN_IDENTITY_DEFAULT       = os.environ["PHATE_BLASTN_IDENTITY_DEFAULT"]
BLASTP_HIT_COUNT_DEFAULT      = os.environ["PHATE_BLASTP_HIT_COUNT_DEFAULT"]
BLASTN_HIT_COUNT_DEFAULT      = os.environ["PHATE_BLASTN_HIT_COUNT_DEFAULT"]

# HMM
HMMER_HOME                    = os.environ["PHATE_HMMER_HOME"]

# Global control: verbosity and error capture
CLEAN_RAW_DATA                = os.environ["PHATE_CLEAN_RAW_DATA"]
PHATE_WARNINGS                = os.environ["PHATE_PHATE_WARNINGS"]
PHATE_MESSAGES                = os.environ["PHATE_PHATE_MESSAGES"]
PHATE_PROGRESS                = os.environ["PHATE_PHATE_PROGRESS"]
PHATE_ERR                     = os.environ["PHATE_PHATE_ERR"]
PHATE_OUT                     = os.environ["PHATE_PHATE_OUT"]
CGC_WARNINGS                  = os.environ["PHATE_CGC_WARNINGS"]
CGC_MESSAGES                  = os.environ["PHATE_CGC_MESSAGES"]
CGC_PROGRESS                  = os.environ["PHATE_CGC_PROGRESS"]
CGP_WARNINGS                  = os.environ["PHATE_CGP_WARNINGS"]
CGP_MESSAGES                  = os.environ["PHATE_CGP_MESSAGES"]
CGP_PROGRESS                  = os.environ["PHATE_CGP_PROGRESS"]

# Constants 

CODE_BASE   = "phate_runPipeline"
CODE        = CODE_BASE + ".py"
CONFIG_FILE = "phate.config"  # by default, but user should name their own, ending in ".config"

# Subordinate codes
GENECALL_CODE_DIR           = BASE_DIR + "GeneCalling/"         # Here resides GENECALL_CODE plus CGC (Compare Gene Calls) codes
SEQANNOT_CODE_DIR           = BASE_DIR + "SequenceAnnotation/"  # Performs sequence annotation via blast and incorporates PSAT output
GENECALL_CODE               = GENECALL_CODE_DIR + "phate_genecallPhage.py"
SEQANNOTATION_CODE_BASE     = "phate_sequenceAnnotation_main"
SEQANNOTATION_CODE          = SEQANNOT_CODE_DIR + SEQANNOTATION_CODE_BASE + ".py"

# Configurables set by values in phate.config input file
DEFAULT_PIPELINE_INPUT_DIR  = 'PipelineInput/'                        # Default
DEFAULT_PIPELINE_OUTPUT_DIR = 'PipelineOutput/'                       # Default
PIPELINE_INPUT_DIR          = BASE_DIR + DEFAULT_PIPELINE_INPUT_DIR   # Default
PIPELINE_OUTPUT_DIR         = BASE_DIR + DEFAULT_PIPELINE_OUTPUT_DIR  # Default
PIPELINE_OUTPUT_SUBDIR      = 'unknown'                               # Will be read in from phate.config file 
GENOME_FILE                 = 'unknown'                               # Will be read in from phate.config file
GENOME_DIR                  = PIPELINE_INPUT_DIR                      # Default; can be overridden via phate.config

# In/out files 

logfile = PIPELINE_OUTPUT_DIR + CODE_BASE + ".log" #*** Should be converted to a generic pipeline log
outfile = PIPELINE_OUTPUT_DIR + CODE_BASE + ".out"

# Create log and out files, if not already there
if not os.path.exists(logfile):
    open(logfile,"w").close()
if not os.path.exists(outfile):
    open(outfile,"w").close()

LOGFILE = open(logfile,"a")  # running log; nead to clean it out occasionally
LOGFILE.write("%s%s%s\n" % ("Begin log file ",datetime.datetime.now(), " **************************** "))
OUTFILE = open(outfile,"a")  # eventually this file will contain instructions for where to find the various outputs from each module
OUTFILE.write("%s%s%s\n" % ("Begin out file ",datetime.datetime.now(), " **************************** "))
runLog = "" # Pipeline log file for current pipeline run (written to user's output subdirectory; this is created once we know the subdir name)

# DEBUG messages control (local)
DEBUG = False
#DEBUG = True

##### HELP STRINGS

HELP_STRING = """This code, """ + CODE + """, runs a phage annotation pipeline, comprising 1) gene calling by 4 gene callers (PHANOTATE, GeneMarkS, Glimmer3, and Prodigal), followed by identification of closest phage genome by means of blast against an NCBI-phage database, and sequence-based functional annotation by means of blastp against several peptide databases (NR, NCBI virus protein, KEGG-virus, Phantome, pVOGs, Swissprot, Refseq protein), and HMM search against these same protein databases. If a PSAT output file is provided, then those annotations are merged with the blast results.\nType: python """ + CODE + """ usage - for more information about constructing the command line.\nType: python """ + CODE + """ detail - for more information about how this code can be run.\n"""

INPUT_STRING = """The input files and other parameters for running this code are specified in a configuration file, which is provided as the only input parameter. See sample configuration file (phate.config.sample) for details on how to set up the configuration file.\n"""

USAGE_STRING = """Usage: python """ + CODE + """ phate.config\n"""

DETAIL_STRING = """Currently the PSAT module is run separately as a web service. In order to incorporate PSAT output into your annotations, you should first run this pipeline specifying "translation_only" in the configuration file. Then, use the generated peptide/protein fasta file as input for PSAT processing. Once you have the PSAT output, save it to the pipeline input directory, and re-run this pipeline, specifying that translation_only is false.\n"""

# First, set defaults; note: setting these values in multiPhate.config file is optional

geneFile                = GENE_FILE
proteinFile             = PROTEIN_FILE
PIPELINE_OUTPUT_SUBDIR  = 'unknown'                               # Will be read in from phate.config file 

##### GET INPUT PARAMETERS

# patterns
p_json   = re.compile("json")
p_input  = re.compile("input")
p_usage  = re.compile("usage")
p_detail = re.compile("detail")

jsonFile = ""

if len(sys.argv) != 2:
    print(HELP_STRING)
    dateTime = os.popen('date')
    LOGFILE.write("%s%s%s%s\n" % ("Incorrect number of input parameters: ", len(sys.argv), ". End log ",dateTime))
    LOGFILE.close(); exit(0)
else:
    inParam = sys.argv[1].lower()
    match_json   = re.search(p_json,   inParam)
    match_input  = re.search(p_input,  inParam)
    match_usage  = re.search(p_usage,  inParam)
    match_detail = re.search(p_detail, inParam)
    if match_json:
        jsonFile = sys.argv[1]
        LOGFILE.write("%s%s\n" % ("Json file is ",jsonFile))
    else:
        if match_input:
            print(INPUT_STRING)
        elif match_usage:
            print(USAGE_STRING)
        elif match_detail:
            print(DETAIL_STRING)
        else:
            print(HELP_STRING)
        LOGFILE.write("%s%s\n" % ("A help string was provided to user; End log ",datetime.datetime.now()))
        LOGFILE.close(); exit(0)

# Parameters

genomeNumber             = 0 
genomeFile               = 'unknown' 
genomeType               = 'unknown' 
genomeSpecies            = 'unknown'
genomeName               = 'unknown' 
outputSubdir             = 'unknown' 
geneticCode              = 11 
translateOnly            = True
phanotateCalls           = False
prodigalCalls            = False
glimmerCalls             = False
genemarksCalls           = False
customGeneCalls          = False
customGeneCallerName     = 'unknown' # Set in user's config file (e.g., "custom" or "rast")
#customGeneCallerOutfile  = 'unknown'
primaryCalls             = 'unknown' 
primaryCallsFile         = 'unknown'
blastpSearch             = False
phmmerSearch             = False
jackhmmerSearch          = False
blastnIdentity           = 60
blastpIdentity           = 60
blastnHitCount           = 5
blastpHitCount           = 5
blastThreads             = 1
ncbiVirusGenomeBlast     = False
ncbiVirusProteinBlast    = False
refseqGeneBlast          = False
refseqProteinBlast       = False
pvogsBlast               = False
vogsBlast                = False #*** To be deprecated
vogGeneBlast             = False
vogProteinBlast          = False
phantomeBlast            = False
phageEnzymeBlast         = False
keggVirusBlast           = False
pfamBlast                = False
smartBlast               = False
swissprotBlast           = False
uniprotBlast             = False
nrBlast                  = False
cazyBlast                = False
customGenomeBlast        = False
customGenomeDBname       = 'unknown' 
customGenomeDBpath       = 'unknown' 
customGeneBlast          = False
customGeneDBname         = 'unknown' 
customGeneDBpath         = 'unknown' 
customProteinBlast       = False 
customProteinDBname      = 'unknown'
customProteinDBpath      = 'unknown'
hmmscan                  = False
ncbiVirusGenomeHmm       = False
ncbiVirusProteinHmm      = False
refseqProteinHmm         = False
refseqGeneHmm            = False
pvogsHmm                 = False 
vogsHmm                  = False 
phantomeHmm              = False
phageEnzymeHmm           = False
keggVirusHmm             = False
pfamHmm                  = False
smartHmm                 = False
swissprotHmm             = False
uniprotHmm               = False
nrHmm                    = False
customHmm                = False 
customHmmDBname          = 'unknown' 
customHmmDBpath          = 'unknown'

# Read json configuration file; record parameter values

outputSubdir = ''
with open(jsonFile, 'r') as jsonParameters:
    parameters               = json.load(jsonParameters)
    genomeNumber             = parameters["genomeNumber"]
    genomeFile               = parameters["genomeFile"]
    genomeType               = parameters["genomeType"]
    genomeSpecies            = parameters["genomeSpecies"]
    genomeName               = parameters["genomeName"]
    outputSubdir             = parameters["outputSubdir"]
    geneticCode              = parameters["geneticCode"]
    translateOnly            = parameters["translateOnly"]
    phanotateCalls           = parameters["phanotateCalls"]
    prodigalCalls            = parameters["prodigalCalls"]
    glimmerCalls             = parameters["glimmerCalls"]
    genemarksCalls           = parameters["genemarksCalls"]
    customGeneCalls          = parameters["customGeneCalls"]
    customGeneCallerName     = parameters["customGeneCallerName"]
    #customGeneCallerOutfile  = parameters["customGeneCallerOutfile"]
    primaryCalls             = parameters["primaryCalls"]
    primaryCallsFile         = parameters["primaryCallsFile"]
    blastpSearch             = parameters["blastpSearch"]
    phmmerSearch             = parameters["phmmerSearch"]
    jackhmmerSearch          = parameters["jackhmmerSearch"]
    blastnIdentity           = parameters["blastnIdentity"]
    blastpIdentity           = parameters["blastpIdentity"]
    blastnHitCount           = parameters["blastnHitCount"]
    blastpHitCount           = parameters["blastpHitCount"]
    blastThreads             = parameters["blastThreads"]
    ncbiVirusGenomeBlast     = parameters["ncbiVirusGenomeBlast"]
    ncbiVirusProteinBlast    = parameters["ncbiVirusProteinBlast"]
    refseqProteinBlast       = parameters["refseqProteinBlast"]
    refseqGeneBlast          = parameters["refseqGeneBlast"]
    pvogsBlast               = parameters["pvogsBlast"]
    vogsBlast                = parameters["vogsBlast"]
    vogGeneBlast             = parameters["vogGeneBlast"]
    vogProteinBlast          = parameters["vogProteinBlast"]
    phantomeBlast            = parameters["phantomeBlast"]
    phageEnzymeBlast         = parameters["phageEnzymeBlast"]
    keggVirusBlast           = parameters["keggVirusBlast"]
    pfamBlast                = parameters["pfamBlast"]
    smartBlast               = parameters["smartBlast"]
    swissprotBlast           = parameters["swissprotBlast"]
    uniprotBlast             = parameters["uniprotBlast"]
    nrBlast                  = parameters["nrBlast"]
    cazyBlast                = parameters["cazyBlast"]
    customGenomeBlast        = parameters["customGenomeBlast"]
    customGenomeDBname       = parameters["customGenomeDBname"]
    customGenomeDBpath       = parameters["customGenomeDBpath"]
    customGeneBlast          = parameters["customGeneBlast"]
    customGeneDBname         = parameters["customGeneDBname"]
    customGeneDBpath         = parameters["customGeneDBpath"]
    customProteinBlast       = parameters["customProteinBlast"]
    customProteinDBname      = parameters["customProteinDBname"]
    customProteinDBpath      = parameters["customProteinDBpath"]
    hmmscan                  = parameters["hmmscan"]
    ncbiVirusGenomeHmm       = parameters["ncbiVirusGenomeHmm"]
    ncbiVirusProteinHmm      = parameters["ncbiVirusProteinHmm"]
    refseqProteinHmm         = parameters["refseqProteinHmm"]
    refseqGeneHmm            = parameters["refseqGeneHmm"]
    pvogsHmm                 = parameters["pvogsHmm"]
    vogsHmm                  = parameters["vogsHmm"]
    phantomeHmm              = parameters["phantomeHmm"]
    phageEnzymeHmm           = parameters["phageEnzymeHmm"]
    keggVirusHmm             = parameters["keggVirusHmm"]
    pfamHmm                  = parameters["pfamHmm"]
    smartHmm                 = parameters["smartHmm"]
    swissprotHmm             = parameters["swissprotHmm"]
    uniprotHmm               = parameters["uniprotHmm"]
    nrHmm                    = parameters["nrHmm"]
    customHmm                = parameters["customHmm"]
    customHmmDBname          = parameters["customHmmDBname"]
    customHmmDBpath          = parameters["customHmmDBpath"]

jsonParameters.close()

# Capture user's configured values
#
#PIPELINE_OUTPUT_SUBDIR = PIPELINE_OUTPUT_DIR + outputSubdir 
PIPELINE_OUTPUT_SUBDIR = os.path.join(PIPELINE_OUTPUT_DIR,outputSubdir) 
#runLog = PIPELINE_OUTPUT_SUBDIR + "runPhATE.log"
runLog = os.path.join(PIPELINE_OUTPUT_SUBDIR,"runPhATE.log")
LOGFILE.write("%s%s\n" % ("PIPELINE_OUTPUT_SUBDIR is ",PIPELINE_OUTPUT_SUBDIR))
GENOME_FILE = genomeName + '.fasta' 
LOGFILE.write("%s%s\n" % ("GENOME_FILE is ",GENOME_FILE))

# Create user's output subdirectory, if doesn't already exist
try:
    os.stat(PIPELINE_OUTPUT_SUBDIR)
except:
    os.mkdir(PIPELINE_OUTPUT_SUBDIR)

# Open log file for this PhATE run in User's subdir
if runLog == "":
    LOGFILE.write("%s\n" % ("WARNING: runLog filename is empty!"))
else:
    LOGFILE.write("%s%s%s%s\n" % ("Opening runLog file ",runLog," ",datetime.datetime.now()))
    fileError = False
    try:
        RUNLOG = open(runLog, "w")
        RUNLOG.write("%s%s\n" % ("Opening runLog at ",datetime.datetime.now()))
    except IOError as e:
        fileError = True
        print(e)
        LOGFILE.write("%s%s%s%s\n" % ("ERROR: could not open runLog file: ",runLog, " at ",datetime.datetime.now()))
        LOGFILE.close(); exit(0)

if genomeType.lower() == 'phage' and primaryCallsFile != 'phanotate.cgc':
    if genomeType.lower() != 'superset' and genomeType.lower() != 'consensus' and genomeType.lower() != 'coregenome':
        if PHATE_WARNINGS == 'True':
            print("phate_runPipeline says, WARNING: If genome type is phage, the primary gene-call file is usually phanotate.cgc! Yours is", primaryCallsFile)
        LOGFILE.write("%s%s\n" % ("WARNING:  User has selected genome type as phage, but primary gene-call file as ", primaryCallsFile))
        RUNLOG.write("%s%s\n" % ("WARNING:  User has selected genome type as phage, but primary gene-call file as ", primaryCallsFile))

if PHATE_MESSAGES == 'True':
    print("PIPELINE_INPUT_DIR is", PIPELINE_INPUT_DIR)
    print("PIPELINE_OUTPUT_DIR is", PIPELINE_OUTPUT_DIR)
    print("PIPELINE_OUTPUT_SUBDIR is", PIPELINE_OUTPUT_SUBDIR)
    print("genomeFile is", genomeFile)
    print("geneFile is", geneFile)
    print("proteinFile is", proteinFile)
    print("genomeNumber is", genomeType) 
    print("genomeType is", genomeType) 
    print("genomeName is", genomeName) 
    print("genomeSpecies is", genomeSpecies) 
    print("geneticCode is", geneticCode) 
    print("Status of boolean translateOnly is", translateOnly)
    print("phanotateCalls is", phanotateCalls)
    print("prodigalCalls is", prodigalCalls)
    print("glimmerCalls is", glimmerCalls)
    print("genemarksCalls is", genemarksCalls)
    print("primaryCalls is", primaryCalls)
    print("primaryCallsFile is", primaryCallsFile)
    print("customGeneCalls is", customGeneCalls)
    print("customGeneCallerName is", customGeneCallerName)
    #print("customGeneCallerOutfile is", customGeneCallerOutfile)
    print("primaryCalls is", primaryCalls) 
    print("blastpSearch is", blastpSearch)
    print("phmmerSearch is", phmmerSearch)
    print("jackhmmerSearch is", jackhmmerSearch)
    print("blastpIdentity is", blastpIdentity) 
    print("blastnIdentity is", blastnIdentity) 
    print("blastpHitCount is", blastpHitCount) 
    print("blastnHitCount is", blastnHitCount)
    print("blastThreads is", blastThreads)
    print("ncbiVirusGenomeBlast is", ncbiVirusGenomeBlast)
    print("ncbiVirusProteinBlast is", ncbiVirusProteinBlast)
    print("refseqProteinBlast is", refseqProteinBlast)
    print("refseqGeneBlast is", refseqGeneBlast)
    print("pvogsBlast is", pvogsBlast)
    print("vogsBlast is", vogsBlast)
    print("phantomeBlast is", phantomeBlast)
    print("phageEnzymeBlast is", phageEnzymeBlast)
    print("keggVirusBlast is", keggVirusBlast)
    print("pfamBlast is", pfamBlast)
    print("smartBlast is", smartBlast)
    print("swissprotBlast is", swissprotBlast)
    print("uniprotBlast is", uniprotBlast)
    print("nrBlast is", nrBlast)
    print("cazyBlast is", cazyBlast)
    print("customGenomeBlast is", customGenomeBlast)
    print("customGenomeDBname is", customGenomeDBname)
    print("customGenomeDBpath is", customGenomeDBpath)
    print("customGeneBlast is", customGeneBlast)
    print("customGeneDBname is", customGeneDBname)
    print("customGeneDBpath is", customGeneDBpath)
    print("customProteinBlast is", customProteinBlast)
    print("customProteinDBname is", customProteinDBname)
    print("customProteinDBpath is", customProteinDBpath)
    print("hmmscan is", hmmscan)
    print("ncbiVirusGenomeHmm is", ncbiVirusGenomeHmm)
    print("ncbiVirusProteinHmm is", ncbiVirusProteinHmm)
    print("refseqProteinHmm is", refseqProteinHmm)
    print("refseqGeneHmm is", refseqGeneHmm)
    print("pvogsHmm is", pvogsHmm)
    print("vogsHmm is", vogsHmm)
    print("phantomeHmm is", phantomeHmm)
    print("phageEnzymeHmm is", phageEnzymeHmm)
    print("pfamHmm is", pfamHmm)
    print("smartHmm is", smartHmm)
    print("swissprotHmm is", swissprotHmm)
    print("uniprotHmm is", uniprotHmm)
    print("nrHmm is", nrHmm)
    print("customHmm is", customHmm)
    print("customHmmDBname is", customHmmDBname)
    print("customHmmDBpath is", customHmmDBpath)

RUNLOG.write("%s\n" % ("Input parameters:"))
RUNLOG.write("%s%s\n" % ("   PIPELINE_INPUT_DIR: ", PIPELINE_INPUT_DIR))
RUNLOG.write("%s%s\n" % ("   PIPELINE_OUTPUT_DIR: ", PIPELINE_OUTPUT_DIR))
RUNLOG.write("%s%s\n" % ("   PIPELINE_OUTPUT_SUBDIR: ", PIPELINE_OUTPUT_SUBDIR))
RUNLOG.write("%s%s\n" % ("   genomeFile: ", genomeFile))
RUNLOG.write("%s%s\n" % ("   geneFile: ", geneFile))
RUNLOG.write("%s%s\n" % ("   proteinFile: ", proteinFile))
RUNLOG.write("%s%s\n" % ("   genomeType is ", genomeType))
RUNLOG.write("%s%s\n" % ("   genomeName is ", genomeName))
RUNLOG.write("%s%s\n" % ("   genomeSpecies is ", genomeSpecies))
RUNLOG.write("%s%s\n" % ("   geneticCode: ", geneticCode))
RUNLOG.write("%s%s\n" % ("   Status of boolean translateOnly is ",translateOnly))
RUNLOG.write("%s%s\n" % ("   phanotateCalls is ",phanotateCalls))
RUNLOG.write("%s%s\n" % ("   prodigalCalls is ",prodigalCalls))
RUNLOG.write("%s%s\n" % ("   glimmerCalls is ",glimmerCalls))
RUNLOG.write("%s%s\n" % ("   genemarksCalls is ",genemarksCalls))
RUNLOG.write("%s%s\n" % ("   primaryCalls is ",primaryCalls))
RUNLOG.write("%s%s\n" % ("   primaryCallsFile is ",primaryCallsFile))
RUNLOG.write("%s%s\n" % ("   customGeneCalls is ",customGeneCalls))
RUNLOG.write("%s%s\n" % ("   customGeneCallerName is ",customGeneCallerName))
#RUNLOG.write("%s%s\n" % ("   customGeneCallerOutfile is ",customGeneCallerOutfile))
RUNLOG.write("%s%s\n" % ("   blastpSearch is ",blastpSearch))
RUNLOG.write("%s%s\n" % ("   phmmerSearch is ",phmmerSearch))
RUNLOG.write("%s%s\n" % ("   jackhmmerSearch is ",jackhmmerSearch))
RUNLOG.write("%s%s\n" % ("   blastpIdentity is ",blastpIdentity))
RUNLOG.write("%s%s\n" % ("   blastnIdentity is ",blastnIdentity))
RUNLOG.write("%s%s\n" % ("   blastpHitCount is ",blastpHitCount))
RUNLOG.write("%s%s\n" % ("   blastnHitCount is ",blastnHitCount))
RUNLOG.write("%s%s\n" % ("   blastThreads is ",blastThreads))
RUNLOG.write("%s%s\n" % ("   ncbiVirusGenomeBlast is ",ncbiVirusGenomeBlast))
RUNLOG.write("%s%s\n" % ("   ncbiVirusProteinBlast is ",ncbiVirusProteinBlast))
RUNLOG.write("%s%s\n" % ("   refseqProteinBlast is ",refseqProteinBlast))
RUNLOG.write("%s%s\n" % ("   refseqGeneBlast is ",refseqGeneBlast))
RUNLOG.write("%s%s\n" % ("   pvogsBlast is ",pvogsBlast))   #*** To be deprecated
RUNLOG.write("%s%s\n" % ("   vogsBlast is ",vogsBlast))
RUNLOG.write("%s%s\n" % ("   vogGeneBlast is ",vogGeneBlast))
RUNLOG.write("%s%s\n" % ("   vogProteinBlast is ",vogProteinBlast))
RUNLOG.write("%s%s\n" % ("   phantomeBlast is ",phantomeBlast))
RUNLOG.write("%s%s\n" % ("   phageEnzymeBlast is ",phageEnzymeBlast))
RUNLOG.write("%s%s\n" % ("   keggVirusBlast is ",keggVirusBlast))
RUNLOG.write("%s%s\n" % ("   pfamBlast is ",pfamBlast))
RUNLOG.write("%s%s\n" % ("   smartBlast is ",smartBlast))
RUNLOG.write("%s%s\n" % ("   swissprotBlast is ",swissprotBlast))
RUNLOG.write("%s%s\n" % ("   uniprotBlast is ",uniprotBlast))
RUNLOG.write("%s%s\n" % ("   nrBlast is ",nrBlast))
RUNLOG.write("%s%s\n" % ("   cazyBlast is ",cazyBlast))
RUNLOG.write("%s%s\n" % ("   hmmscan is ",hmmscan))
RUNLOG.write("%s%s\n" % ("   ncbiVirusGenomeHmm is ",ncbiVirusGenomeHmm))
RUNLOG.write("%s%s\n" % ("   ncbiVirusProteinHmm is ",ncbiVirusProteinHmm))
RUNLOG.write("%s%s\n" % ("   refseqProteinHmm is ",refseqProteinHmm))
RUNLOG.write("%s%s\n" % ("   refseqGeneHmm is ",refseqGeneHmm))
RUNLOG.write("%s%s\n" % ("   pvogsHmm is ",pvogsHmm))
RUNLOG.write("%s%s\n" % ("   vogsHmm is ",vogsHmm))
RUNLOG.write("%s%s\n" % ("   phantomeHmm is ",phantomeHmm))
RUNLOG.write("%s%s\n" % ("   phageEnzymeHmm is ",phageEnzymeHmm))
RUNLOG.write("%s%s\n" % ("   keggVirusHmm is ",keggVirusHmm))
RUNLOG.write("%s%s\n" % ("   pfamHmm is ",pfamHmm))
RUNLOG.write("%s%s\n" % ("   smartHmm is ",smartHmm))
RUNLOG.write("%s%s\n" % ("   swissprotHmm is ",swissprotHmm))
RUNLOG.write("%s%s\n" % ("   uniprotHmm is ",uniprotHmm))
RUNLOG.write("%s%s\n" % ("   nrHmm is ",nrHmm))
RUNLOG.write("%s%s\n" % ("   customGenomeBlast is ",customGenomeBlast))
RUNLOG.write("%s%s\n" % ("   customGenomeDBname is ",customGenomeDBname))
RUNLOG.write("%s%s\n" % ("   customGenomeDBpath is ",customGenomeDBpath))
RUNLOG.write("%s%s\n" % ("   customGeneBlast is ",customGeneBlast))
RUNLOG.write("%s%s\n" % ("   customGeneDBname is ",customGeneDBname))
RUNLOG.write("%s%s\n" % ("   customGeneDBpath is ",customGeneDBpath))
RUNLOG.write("%s%s\n" % ("   customProteinBlast is ",customProteinBlast))
RUNLOG.write("%s%s\n" % ("   customProteinDBname is ",customProteinDBname))
RUNLOG.write("%s%s\n" % ("   customProteinDBpath is ",customProteinDBpath))

# Open and check input file(s)

inputDir             = PIPELINE_INPUT_DIR
genomeFile           = PIPELINE_INPUT_DIR + GENOME_FILE
primaryCallsPathFile = PIPELINE_OUTPUT_SUBDIR + primaryCallsFile 
# Reassign geneFile and proteinFile as path/filename
geneFile             = PIPELINE_OUTPUT_SUBDIR + GENE_FILE
proteinFile          = PIPELINE_OUTPUT_SUBDIR + PROTEIN_FILE
outputDir            = PIPELINE_OUTPUT_SUBDIR

RUNLOG.write("%s%s\n" % ("inputDir is ",    inputDir))
RUNLOG.write("%s%s\n" % ("genomeFile is ",  genomeFile))
RUNLOG.write("%s%s\n" % ("primaryCallsPathFile is ",primaryCallsPathFile))
RUNLOG.write("%s%s\n" % ("geneFile is ",    geneFile))
RUNLOG.write("%s%s\n" % ("proteinFile is ", proteinFile))

if PHATE_PROGRESS == 'True':
    print("phate_runPipeline says, Checking files...")
RUNLOG.write("%s\n" % ("Checking files..."))
fileError = False

# Check genome file

fileError = False
try:
    GENOME_H = open(genomeFile,"r")
except IOError as e:
    fileError = True
    print(e) 

if fileError:
    print(USAGE_STRING)
    print("phate_runPipeline says, ERROR: Check your genome file,", genomeFile)
    LOGFILE.write("%s%s%s%s\n" % ("ERROR:  Genome file could not be opened: ", genomeFile, "; End log ", datetime.datetime.now()))
    LOGFILE.close(); exit(0)
GENOME_H.close()

# Copy json configuration file to the user's results directory
jsonSave = outputDir + jsonFile
command = "cp " + jsonFile + ' ' + jsonSave
os.system(command)

if PHATE_PROGRESS == 'True':
    print("phate_runPipeline says, Configuration complete.")

##### BEGIN MAIN ########################################################################################

##### Run Gene-calling Module

if PHATE_PROGRESS == 'True':
    print("phate_runPipeline says, Preparing to run genecall module...")
RUNLOG.write("%s\n" % ("Preparing to run genecall module..."))

param2 = outputDir[:-1]  # remove terminal '/' because subordinate code adds it explicitly
param3 = ''  # can't pass a dict; future re-write genecalling module as class

MULTIPLE_CALLERS = False
geneCallerCount = 0

if genemarksCalls:
    geneCallerCount += 1
    param3 += "genemarks_"

if prodigalCalls:
    geneCallerCount += 1
    param3 += "prodigal_"

if glimmerCalls:
    geneCallerCount += 1
    param3 += "glimmer_"

if phanotateCalls:
    geneCallerCount += 1
    param3 += "phanotate_"

if customGeneCalls: 
    customCallFile  = outputSubdir.rstrip('/') + ".custom.gff"  # This is what the user should have named it
    #customPhateFile = "phate_custom.gff" 
    customPhateFile = customCallFile 
    source      = os.path.join(PIPELINE_INPUT_DIR,    customCallFile)
    destination = os.path.join(PIPELINE_OUTPUT_SUBDIR,customPhateFile)
    try:
        command = "cp " + source + ' ' + destination 
        success = os.system(command)
    except:
        if PHATE_WARNINGS:
            print("phate_runPipeline says, WARNING: Cannot move user's custom genecall file to pipeline output subdirectory. Custom gene calls will not be used.")
        if primaryCalls == 'custom':  # If user selected custom as primary calls, and file is in error, then must abort processing.
            print("phate_runPipeline says, ERROR: Cannot use custom calls as primary")
            exit(0)
    geneCallerCount += 1
    param3 += "custom_"
    # If the user selected 'custom' as the primary calls, then reset to the custom calls pathfile
    #if primaryCalls == 'custom': 
    #    #primaryCallsPathFile = destination 
    #    primaryCallsPathFile = "custom.cgc" 

if geneCallerCount >= 2:
    MULTIPLE_CALLERS = True

# If there are multiple callers and if the user is selecting the primary call set as the superset, consensus,
# or common core, then reassign primaryCallsPathFile (from 'unknown') to the CGC-generated gene set in cgc format.
# The superset is the non-redundant set of gene calls produced by all gene callers.
# The consensus set comprises the non-redundant set of gene calls that agree on stop codon.
# The common_core set comprises the non-redundant set of gene calls that are identical.
if MULTIPLE_CALLERS:
    if primaryCalls == 'superset':
        primaryCallsPathFile = PIPELINE_OUTPUT_SUBDIR + 'superset.cgc' 
    if primaryCalls == 'consensus':
        primaryCallsPathFile = PIPELINE_OUTPUT_SUBDIR + 'consensus.cgc' 
    if primaryCalls == 'common_core':
        primaryCallsPathFile = PIPELINE_OUTPUT_SUBDIR + 'common_core.cgc' 

command = "python " + GENECALL_CODE + ' ' + genomeFile + ' ' + param2 + ' ' + param3 + ' ' + primaryCallsPathFile
if PHATE_PROGRESS == 'True':
    print("phate_runPipeline says, Calling the gene-call module.")
if PHATE_MESSAGES == 'True':
    print("phate_runPipeline says, command is,", command)
RUNLOG.write("%s%s\n" % ("Calling the gene-call module. Command is ", command))

# OS system matters; choose alternate system call if you get error message on this line
result = os.system(command)
#result = subprocess.check_output(command,shell=True)

if PHATE_PROGRESS == 'True':
    print("phate_runPipeline says, Gene-call processing complete.")
RUNLOG.write("%s%s\n" % ("Gene-call processing complete at ", datetime.datetime.now()))

##### Run Sequence Annotation Module

RUNLOG.write("%s\n" % ("Preparing to call sequence annotation module..."))

# Construct command line parameter string

if PHATE_PROGRESS == 'True':
    print("phate_runPipeline says, Preparing command strings for homology searches...")

# First, construct string listing the names of databases to be blasted

blastDatabaseParameterString = ''   # databases
if ncbiVirusGenomeBlast:
    blastDatabaseParameterString += '_ncbiVirusGenome'
if ncbiVirusProteinBlast:
    blastDatabaseParameterString += '_ncbiVirusProtein'
if refseqProteinBlast:
    blastDatabaseParameterString += '_refseqProtein'
if refseqGeneBlast:
    blastDatabaseParameterString += '_refseqGene'
if pvogsBlast:
    blastDatabaseParameterString += '_pvogs'
if vogsBlast:
    blastDatabaseParameterString += '_Vogs'  #*** to be deprecated
if vogGeneBlast:
    blastDatabaseParameterString += '_vogGene'
if vogProteinBlast:
    blastDatabaseParameterString += '_vogProtein'
if phantomeBlast:
    blastDatabaseParameterString += '_phantome'
if phageEnzymeBlast:
    blastDatabaseParameterString += '_phageEnzyme'
if keggVirusBlast:
    blastDatabaseParameterString += '_kegg'
if pfamBlast:
    blastDatabaseParameterString += '_pfam'
if smartBlast:
    blastDatabaseParameterString += '_smart'
if swissprotBlast:
    blastDatabaseParameterString += '_swissprot'
if uniprotBlast:
    blastDatabaseParameterString += '_uniprot'
if nrBlast:
    blastDatabaseParameterString += '_nr'
if cazyBlast:
    blastDatabaseParameterString += '_cazy'
if customGenomeBlast:
    blastDatabaseParameterString += '_customGenome'
if customGeneBlast:
    blastDatabaseParameterString += '_customGene'
if customProteinBlast:
    blastDatabaseParameterString += '_customProtein'

blastProgramParameterString = ''   # programs 
if blastpSearch:
    blastProgramParameterString += '_blastp'

blastThreadsParameterString = ''   # number of blast threads to run
if blastThreads > 0:
    blastThreadsParameterString =  str(blastThreads)

# For now, use same databases for hmm search against sequence databases
seqDatabaseParameterString = blastDatabaseParameterString 

hmmProgramParameterString = ''     # programs
if phmmerSearch:
    hmmProgramParameterString   += '_phmmer'
if jackhmmerSearch:
    hmmProgramParameterString   += '_jackhmmer'

# Construct string listing hmm program and databases to be searched

profileDatabaseParameterString = ''  # hmm profile databases
if ncbiVirusGenomeHmm:
    profileDatabaseParameterString += '_ncbiVirusGenomeHmm'
if ncbiVirusProteinHmm:
    profileDatabaseParameterString += '_ncbiVirusProteinHmm'
if refseqProteinHmm:
    profileDatabaseParameterString += '_refseqProteinHmm'
if refseqGeneHmm:
    profileDatabaseParameterString += '_refseqGeneHmm'
if pvogsHmm:
    profileDatabaseParameterString += '_pvogsHmm'
if vogsHmm:
    profileDatabaseParameterString += '_VogsHmm'
if phantomeHmm:
    profileDatabaseParameterString += '_phantomeHmm'
if phageEnzymeHmm:
    profileDatabaseParameterString += '_phageEnzymeHmm'
if keggVirusHmm:
    profileDatabaseParameterString += '_keggVirusHmm'
if pfamHmm:
    profileDatabaseParameterString += '_pfamHmm'
if smartHmm:
    profileDatabaseParameterString += '_smartHmm'
if uniprotHmm:
    profileDatabaseParameterString += '_uniprotHmm'
if swissprotHmm:
    profileDatabaseParameterString += '_swissprotHmm'
if nrHmm:
    profileDatabaseParameterString += '_nrHmm'
if customHmm:
    profileDatabaseParameterString += '_customHmm'

profileProgramParameterString = ''  # hmm programs to run against sequence databases
if hmmscan:
    profileProgramParameterString  += '_hmmscan'

"""
customDatabasePathParameterString = ''
if customGenomeDBpath != '':
    customDatabasePathParameterString += '^customGenomeDBpath:'  + customGenomeDBpath  + ':' + customGenomeDBname
if customGeneDBpath != '':
    customDatabasePathParameterString += '^customGeneDBpath:'    + customGeneDBpath    + ':' + customGeneDBname
if customProteinDBpath != '':
    customDatabasePathParameterString += '^customProteinDBpath:' + customProteinDBpath + ':' + customProteinDBname
if customHmmDBpath != '':
    customDatabasePathParameterString += '^customHmmDBpath:'     + customHmmDBpath     + ':' + customHmmDBname
"""

# Set empty strings (so parameter is passed regardless, not as '')
if blastProgramParameterString == '':
    blastProgramParameterString = 'none'
if blastThreadsParameterString == '':
    blastThreadsParameterString = 'none'
if blastDatabaseParameterString == '':
    blastDatabaseParameterString = 'none'
if hmmProgramParameterString == '':
    hmmProgramParameterString = 'none'
if seqDatabaseParameterString == '':
    seqDatabaseParameterString = 'none'
if profileProgramParameterString == '':
    profileProgramParameterString = 'none'
if profileDatabaseParameterString == '':
    profileDatabaseParameterString = 'none'
#if customDatabasePathParameterString == '':
#    customDatabasePathParameterString = 'none'

commandRoot1  = "python " + SEQANNOTATION_CODE      + " -o " + outputDir  # code and output direction
commandRoot2  = " -G "    + genomeFile              + " -g " + geneFile             + " -p " + proteinFile   # genome files
commandRoot3  = " -c "    + primaryCalls            + " -f " + primaryCallsPathFile                          # gene-call information
commandRoot4  = " -t "    + genomeType              + " -n " + genomeName           + " -s " + genomeSpecies # genome meta-data
commandRoot5  = " -i "    + blastpIdentity          + " -j " + blastnIdentity                                # blast identity parameters 
commandRoot6  = " -h "    + blastpHitCount          + " -H " + blastnHitCount                                # blast hit count parameters
commandRoot7  = " -B "    + blastProgramParameterString   + " -b " + blastDatabaseParameterString            # blast and hmm search of blast/sequence database(s)
commandRoot8  = " -M "    + hmmProgramParameterString     + " -m " + seqDatabaseParameterString              # hmm search of hmm profile database(s)
commandRoot9  = " -R "    + profileProgramParameterString + " -r " + profileDatabaseParameterString          # program and databases for hmm search
commandRoot10 = " -z "    + blastThreadsParameterString
#commandRoot11 = " -C "    + customDatabasePathParameterString
commandRootA  = commandRoot1 + commandRoot2 + commandRoot3 + commandRoot4 + commandRoot5  + commandRoot6
#commandRoot   = commandRootA + commandRoot7 + commandRoot8 + commandRoot9 + commandRoot10 + commandRoot11
commandRoot   = commandRootA + commandRoot7 + commandRoot8 + commandRoot9 + commandRoot10

# As appropriate, append additional parameters
if translateOnly:
    command = commandRoot + " -x true "        # setting TRANSLATE_ONLY to True
else:
    command = commandRoot 

if PHATE_MESSAGES == 'True':
    print("phate_runPipeline says, command is ",command)
# Communicate and execute
if PHATE_PROGRESS == 'True':
    print("phate_runPipeline says, Calling the sequence annotation module.")
if PHATE_MESSAGES == 'True':
    print("phate_runPipeline says, command is,", command)
RUNLOG.write("%s%s\n" % ("Calling the sequence annotation module. Command is ", command))
result = os.system(command)
if PHATE_PROGRESS == 'True':
    print("phate_runPipeline says, Sequence annotation processing is complete.")
RUNLOG.write("%s%s\n" % ("Sequence annotation processing complete at ", datetime.datetime.now()))

##### CLEAN UP

if PHATE_PROGRESS == 'True':
    print("phate_runPipeline says, Code completed at", datetime.datetime.now())
OUTFILE.write("%s%s\n" %("Pipeline output is in output file created by code ",SEQANNOTATION_CODE))
OUTFILE.close()
LOGFILE.write("%s%s\n" % ("Code completed at ", datetime.datetime.now()))
LOGFILE.close()
RUNLOG.write("%s%s\n" % ("Execution complete at ",datetime.datetime.now()))
RUNLOG.close()

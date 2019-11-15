#!/usr/bin/env python

################################################################
#
# Program Title:  phate_runPipeline.py ()
#
# Description: Runs the phate annotation pipeline.  This code runs under Python 3.7, and requires
#    dependent packages.
#
# Usage:  python phate_runPipeline.py myPhATE.config
#    (see sample.config for how to create your configuration file)
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
#    Jeff Kimbrel  - GeneCalling/
#
################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.PDF FOR DETAILS.

import sys, os, re, string, copy, time, datetime
import subprocess
#import logging

# CONSTANTS and CONFIGURABLES

# Output file names 
CONSENSUS_CALLS_FILE     = 'phanotate.cgc' #*** For now this is PHANOTATE calls, though may be consensus calls in future
GENE_FILE                = 'gene.fnt'                              #
PROTEIN_FILE             = 'protein.faa'                           #

# DEFAULTS

# Organism 
GENETIC_CODE_DEFAULT          = 11               # Bacterial 
GENE_CALLER_DEFAULT           = 'phanotate' 
GENOME_TYPE_DEFAULT           = 'phage' 
NAME_DEFAULT                  = 'unknown'
CONTIG_NAME_DEFAULT           = 'unknown'
SPECIES_DEFAULT               = 'unknown'

# Gene-caller defaults
GENEMARKS_CALLS_DEFAULT  = False     # Requires license
PRODIGAL_CALLS_DEFAULT   = False
GLIMMER_CALLS_DEFAULT    = False
PHANOTATE_CALLS_DEFAULT  = False

# Blast parameter defaults
MAX_BLAST_HIT_COUNT      = 100       # maximum number of hits to capture (user should specify far fewer than max)
MIN_BLASTP_IDENTITY      = '20'      # default; sets a lower limit based on value at which a structure model can provide information
MIN_BLASTN_IDENTITY      = '20'      # default; sets a lower limit based on value at which a structure model can provide information
MAX_BLASTP_HIT_COUNT     = '100'     # default; sets an upper limit; user's value should typically be well below this
MAX_BLASTN_HIT_COUNT     = '10'      # default; sets an upper limit

# Blast database defaults
NCBI_VIRUS_BLAST_DEFAULT         = False
NCBI_VIRUS_PROTEIN_BLAST_DEFAULT = False
KEGG_VIRUS_BLAST_DEFAULT         = False     # Requires license
NR_BLAST_DEFAULT                 = False     # Large data set; blast run takes time
REFSEQ_PROTEIN_BLAST_DEFAULT     = False      # Large data set; blast run takes time
PHANTOME_BLAST_DEFAULT           = False
PVOGS_BLAST_DEFAULT              = False
UNIPARC_BLAST_DEFAULT            = False     # Turned 'off' for now; not yet in service
REFSEQ_GENE_BLAST_DEFAULT        = False
SWISSPROT_BLAST_DEFAULT          = False
UNIPROT_BLAST_DEFAULT            = False     # not yet in service
PFAM_BLAST_DEAFULT               = False     # not yet in service
SMART_BLAST_DEAFULT              = False     # not yet in service

# HMM programs and databases
HMM_PROGRAM_DEFAULT              = 'jackhmmer'  # This is the only program currently supported
NCBI_VIRUS_HMM_DEFAULT           = False     # not yet in service
NCBI_VIRUS_PROTEIN_HMM_DEFAULT   = False     # not yet in service
KEGG_VIRUS_HMM_DEFAULT           = False     # Requires license
NR_HMM_DEFAULT                   = False     # Large data set; hmm run takes time
REFSEQ_PROTEIN_HMM_DEFAULT       = False     # Large data set; hmm run takes time
PHANTOME_HMM_DEFAULT             = False     # not yet in service
PVOGS_HMM_DEFAULT                = False     #
UNIPARC_HMM_DEFAULT              = False     # not yet in service
REFSEQ_GENE_HMM_DEFAULT          = False     # not yet in service
SWISSPROT_HMM_DEFAULT            = False     # not yet in service
UNIPROT_HMM_DEFAULT              = False     # not yet in service
PFAM_HMM_DEFAULT                 = False     # not yet in service
SMART_HMM_DEFAULT                = False     # not yet in service

# Other
PSAT_ANNOTATION_DEFAULT          = False     # Requires LLNL processing

##### ENVIRONMENT VARIABLES

BASE_DIR                      = os.environ["BASE_DIR"]
DATABASE_DIR                  = os.environ["DATABASE_DIR"]
SOFTWARE_DIR                  = os.environ["SOFTWARE_DIR"]

EMBOSS_HOME                   = os.environ["EMBOSS_PHATE_HOME"] 
PIPELINE_DIR                  = os.environ["PIPELINE_DIR"] 
PSAT_OUT_DIR                  = os.environ["PSAT_OUT_DIR"] 

# Data sets
KEGG_VIRUS_BASE_DIR           = os.environ["KEGG_VIRUS_BASE_DIR"] 
KEGG_VIRUS_BLAST_HOME         = os.environ["KEGG_VIRUS_BLAST_HOME"]
NCBI_VIRUS_BASE_DIR           = os.environ["NCBI_VIRUS_BASE_DIR"] 
NCBI_VIRUS_BLAST_HOME         = os.environ["NCBI_VIRUS_BLAST_HOME"]
NCBI_VIRUS_PROTEIN_BLAST_HOME = os.environ["NCBI_VIRUS_PROTEIN_BLAST_HOME"]
NCBI_TAXON_DIR                = os.environ["NCBI_TAXON_DIR"] 
PHANTOME_BASE_DIR             = os.environ["PHANTOME_BASE_DIR"]
PHANTOME_BLAST_HOME           = os.environ["PHANTOME_BLAST_HOME"]
PVOGS_BASE_DIR                = os.environ["PVOGS_BASE_DIR"]
PVOGS_BLAST_HOME              = os.environ["PVOGS_BLAST_HOME"]
UNIPARC_BASE_DIR              = os.environ["UNIPARC_BASE_DIR"]
UNIPARC_VIRUS_BLAST_HOME      = os.environ["UNIPARC_VIRUS_BLAST_HOME"]
NR_BLAST_BASE_DIR             = os.environ["NR_BLAST_BASE_DIR"]
NR_BLAST_HOME                 = os.environ["NR_BLAST_HOME"]
REFSEQ_PROTEIN_BASE_DIR       = os.environ["REFSEQ_PROTEIN_BASE_DIR"]
REFSEQ_PROTEIN_BLAST_HOME     = os.environ["REFSEQ_PROTEIN_BLAST_HOME"]
REFSEQ_GENE_BASE_DIR          = os.environ["REFSEQ_GENE_BASE_DIR"]
REFSEQ_GENE_BLAST_HOME        = os.environ["REFSEQ_GENE_BLAST_HOME"]
SWISSPROT_BASE_DIR            = os.environ["SWISSPROT_BASE_DIR"]
SWISSPROT_BLAST_HOME          = os.environ["SWISSPROT_BLAST_HOME"]
UNIPROT_BASE_DIR              = os.environ["UNIPROT_BASE_DIR"]
UNIPROT_BLAST_HOME            = os.environ["UNIPROT_BLAST_HOME"]
PFAM_BASE_DIR                 = os.environ["PFAM_BASE_DIR"]
PFAM_BLAST_HOME               = os.environ["PFAM_BLAST_HOME"]

# Gene calling
PRODIGAL_PATH                 = os.environ["PRODIGAL_PATH"]
GLIMMER_PATH                  = os.environ["GLIMMER_PATH"]
GENEMARKS_PATH                = os.environ["GENEMARKS_PATH"]  
PHANOTATE_PATH                = os.environ["PHANOTATE_PATH"]
CGC_PATH                      = os.environ["CGC_PATH"]

# Blast
BLAST_HOME                    = os.environ["BLAST_HOME"]
MIN_BLASTP_IDENTITY           = os.environ["MIN_BLASTP_IDENTITY"]
MIN_BLASTN_IDENTITY           = os.environ["MIN_BLASTN_IDENTITY"]
MAX_BLASTP_HIT_COUNT          = os.environ["MAX_BLASTP_HIT_COUNT"] 
MAX_BLASTN_HIT_COUNT          = os.environ["MAX_BLASTN_HIT_COUNT"] 
BLASTP_IDENTITY_DEFAULT       = os.environ["BLASTP_IDENTITY_DEFAULT"]
BLASTN_IDENTITY_DEFAULT       = os.environ["BLASTN_IDENTITY_DEFAULT"]
BLASTP_HIT_COUNT_DEFAULT      = os.environ["BLASTP_HIT_COUNT_DEFAULT"]
BLASTN_HIT_COUNT_DEFAULT      = os.environ["BLASTN_HIT_COUNT_DEFAULT"]

# HMM
HMM_HOME                      = os.environ["HMM_HOME"]

# Global control: verbosity and error capture
CLEAN_RAW_DATA                = os.environ["CLEAN_RAW_DATA"]
PHATE_WARNINGS                = os.environ["PHATE_WARNINGS"]
PHATE_MESSAGES                = os.environ["PHATE_MESSAGES"]
PHATE_PROGRESS                = os.environ["PHATE_PROGRESS"]
PHATE_ERR                     = os.environ["PHATE_ERR"]
PHATE_OUT                     = os.environ["PHATE_OUT"]
CGC_WARNINGS                  = os.environ["CGC_WARNINGS"]
CGC_MESSAGES                  = os.environ["CGC_MESSAGES"]
CGC_PROGRESS                  = os.environ["CGC_PROGRESS"]

# Constants 

CODE_BASE   = "phate_runPipeline"
CODE        = CODE_BASE + ".py"
CONFIG_FILE = "phate.config"  # by default, but user should name their own, ending in ".config"

# Subordinate codes
GENECALL_CODE_DIR       = BASE_DIR + "GeneCalling/"         # Here resides GENECALL_CODE plus CGC (Compare Gene Calls) codes
SEQANNOT_CODE_DIR       = BASE_DIR + "SequenceAnnotation/"  # Performs sequence annotation via blast and incorporates PSAT output
GENECALL_CODE           = GENECALL_CODE_DIR + "phate_genecallPhage.py"
SEQANNOTATION_CODE_BASE = "phate_sequenceAnnotation_main"
SEQANNOTATION_CODE      = SEQANNOT_CODE_DIR + SEQANNOTATION_CODE_BASE + ".py"

# Configurables set by values in phate.config input file
DEFAULT_PIPELINE_INPUT_DIR  = 'PipelineInput/'                        # Default
DEFAULT_PIPELINE_OUTPUT_DIR = 'PipelineOutput/'                       # Default
PIPELINE_INPUT_DIR          = BASE_DIR + DEFAULT_PIPELINE_INPUT_DIR   # Default
PIPELINE_OUTPUT_DIR         = BASE_DIR + DEFAULT_PIPELINE_OUTPUT_DIR  # Default
PIPELINE_OUTPUT_SUBDIR      = 'unknown'                               # Will be read in from phate.config file 
GENOME_FILE                 = 'unknown'                               # Will be read in from phate.config file
GENOME_DIR                  = PIPELINE_INPUT_DIR                      # Default; can be overridden via phate.config
PSAT_FILE                   = 'unknown'                               # Will be read in from phate.config file, if specified
PSAT_DIR                    = PIPELINE_INPUT_DIR                      # Default; can be overridden via phate.config 

# In/out files 

logfile = PIPELINE_OUTPUT_DIR + CODE_BASE + ".log" #*** Should be converted to a generic pipeline log
outfile = PIPELINE_OUTPUT_DIR + CODE_BASE + ".out"

# Create log and out files, if not already there
if not os.path.exists(logfile):
    open(logfile,"w").close()
if not os.path.exists(outfile):
    open(outfile,"w").close()

infile  = CONFIG_FILE   # by default, but this is specified by input parameter
LOGFILE = open(logfile,"a")  # running log; nead to clean it out occasionally
LOGFILE.write("%s%s%s\n" % ("Begin log file ",datetime.datetime.now(), " **************************** "))
OUTFILE = open(outfile,"a")  # eventually this file will contain instructions for where to find the various outputs from each module
OUTFILE.write("%s%s%s\n" % ("Begin out file ",datetime.datetime.now(), " **************************** "))
runLog = "" # Pipeline log file for current pipeline run (written to user's output subdirectory; this is created once we know the subdir name)

##### PATTERNS and CONTROL

# PATTERNS

# General
p_comment               = re.compile("^#")
p_blank                 = re.compile("^$")
p_help                  = re.compile("help")
p_input                 = re.compile("input")
p_usage                 = re.compile("usage")
p_detail                = re.compile("detail")
p_config                = re.compile("config")
p_outputSubdir          = re.compile("output_subdir='(.*)'")
p_genomeFile            = re.compile("genome_file='(.*)'")
p_genomeType            = re.compile("genome_type='(.*)'")
p_name                  = re.compile("genomeName='(.*)'")  
p_contig                = re.compile("contig='(.*)'")  #*** For now, finished genome, single contig only
p_species               = re.compile("species='(.*)'")

# Gene calling
p_geneCaller            = re.compile("primary_calls='(.*)'")
p_genemarksCalls        = re.compile("genemarks_calls='(.*)'")
p_glimmerCalls          = re.compile("glimmer_calls='(.*)'")
p_prodigalCalls         = re.compile("prodigal_calls='(.*)'")
p_phanotateCalls        = re.compile("phanotate_calls='(.*)'")
p_geneticCode           = re.compile("genetic_code='(\d+)'")
p_translateOnly         = re.compile("translate_only='(.*)'")
p_custom_geneCalls      = re.compile("custom_gene_calls='(.*)'")
p_custom_geneCallerName = re.compile("custom_gene_caller_name='(.*)'")
p_custom_geneCallerOutfile = re.compile("custom_gene_caller_outfile='(.*)'")

# Blast
p_blastpIdentity        = re.compile("blastp_identity='(\d+)'")   #*** Complete code to distinguish n,p
p_blastnIdentity        = re.compile("blastn_identity='(\d+)'")   #***
p_blastpHitCount        = re.compile("blastp_hit_count='(\d+)'")
p_blastnHitCount        = re.compile("blastn_hit_count='(\d+)'")
p_ncbiVirusBlast        = re.compile("ncbi_virus_genome_blast='(.*)'")
p_ncbiVirusProteinBlast = re.compile("ncbi_virus_protein_blast='(.*)'")
p_keggVirusBlast        = re.compile("kegg_virus_blast='(.*)'")
p_nrBlast               = re.compile("nr_blast='(.*)'")
p_refseqProteinBlast    = re.compile("refseq_protein_blast='(.*)'")
p_refseqGeneBlast       = re.compile("refseq_gene_blast='(.*)'")
p_swissprotBlast        = re.compile("swissprot_blast='(.*)'")
p_phantomeBlast         = re.compile("phantome_blast='(.*)'")
p_pvogsBlast            = re.compile("pvogs_blast='(.*)'")
p_pfamBlast             = re.compile("pfam_blast='(.*)'")
p_smartBlast            = re.compile("smart_blast='(.*)'")
p_uniparcBlast          = re.compile("uniparc_blast='(.*)'")
p_uniprotBlast          = re.compile("uniprot_blast='(.*)'")

# Custom blast processes and parameters
# Custom blast processes to run (true/false) 
p_customGenomeBlast     = re.compile("custom_genome_blast='(.*)'")                
p_customGeneBlast       = re.compile("custom_gene_blast='(.*)'")                   
p_customProteinBlast    = re.compile("custom_protein_blast='(.*)'")                
# Custom blast database names
p_customGenomeDBname    = re.compile("custom_genome_blast_database_name='(.*)'")   
p_customGeneDBname      = re.compile("custom_gene_blast_database_name='(.*)'")     
p_customProteinDBname   = re.compile("custom_protein_blast_database_name='(.*)'")  
# Paths to custom databases
p_customGenomeDBpath    = re.compile("custom_genome_blast_database_path='(.*)'")  
p_customGeneDBpath      = re.compile("custom_gene_blast_database_path='(.*)'")    
p_customProteinDBpath   = re.compile("custom_protein_blast_database_path='(.*)'")  

# HMM
# programs
p_hmmProgram            = re.compile("hmm_program='(.*)'")
p_jackhmmer             = re.compile("jackhmmer='(.*)'")
p_hmmscan               = re.compile("hmmscan='(.*)'")
p_hmmbuild              = re.compile("hmmbuild='(.*)'")
# databases
p_pvogsHmm              = re.compile("pvogs_hmm_profiles='(.*)'")
p_phantomeHmm           = re.compile("phantome_hmm_profiles='(.*)'")
p_swissprotHmm          = re.compile("swissprot_hmm_profiles='(.*)'")
p_pfamHmm               = re.compile("pfam_hmm_profiles='(.*)'")
p_smartHmm              = re.compile("smart_hmm_profiles='(.*)'")
p_uniprotHmm            = re.compile("uniprot_hmm_profiles='(.*)'")
p_ncbiVirusHmm          = re.compile("ncbi_virus_genome_hmm_profiles='(.*)'")
p_ncbiVirusProteinHmm   = re.compile("ncbi_virus_protein_hmm_profiles='(.*)'")
p_keggVirusHmm          = re.compile("kegg_virus_hmm_profiles='(.*)'")
p_nrHmm                 = re.compile("nr_hmm_profiles='(.*)'")
p_refseqProteinHmm      = re.compile("refseq_protein_hmm_profiles='(.*)'")
p_refseqGeneHmm         = re.compile("refseq_gene_hmm_profiles='(.*)'")
p_uniparcHmm            = re.compile("uniparc_hmm_profiles='(.*)'")

# PSAT
p_psatAnnotation        = re.compile("psat_annotation='(.*)'")
p_psatFile              = re.compile("psat_file='(.*)'")

# BOOLEANS 

# DEBUG messages control (local)
DEBUG = True
#DEBUG = False

# Other boolean control
PSAT = False            # This turns True if a psat output file is specified in the config file AND TRANSLATE_ONLY is False
TRANSLATE_ONLY = False  # User will specify 'True' in config file if only generating gene & protein files

##### HELP STRINGS

HELP_STRING = """This code, """ + CODE + """, runs a phage annotation pipeline, comprising 1) gene calling by 4 gene callers (PHANOTATE, GeneMarkS, Glimmer3, and Prodigal), followed by identification of closest phage genome by means of blast against an NCBI-phage database, and sequence-based functional annotation by means of blastp against several peptide databases (NR, NCBI virus protein, KEGG-virus, Phantome, pVOGs, Swissprot, Refseq protein), and HMM search against these same protein databases. If a PSAT output file is provided, then those annotations are merged with the blast results.\nType: python """ + CODE + """ usage - for more information about constructing the command line.\nType: python """ + CODE + """ detail - for more information about how this code can be run.\n"""

INPUT_STRING = """The input files and other parameters for running this code are specified in a configuration file, which is provided as the only input parameter. See sample configuration file (phate.config.sample) for details on how to set up the configuration file.\n"""

USAGE_STRING = """Usage: python """ + CODE + """ phate.config\n"""

DETAIL_STRING = """Currently the PSAT module is run separately as a web service. In order to incorporate PSAT output into your annotations, you should first run this pipeline specifying "translation_only" in the configuration file. Then, use the generated peptide/protein fasta file as input for PSAT processing. Once you have the PSAT output, save it to the pipeline input directory, and re-run this pipeline, specifying that translation_only is false.\n"""

##### GET INPUT PARAMETERS #####

if len(sys.argv) != 2:
    print(HELP_STRING)
    dateTime = os.popen('date')
    LOGFILE.write("%s%s%s%s\n" % ("Incorrect number of input parameters: ", len(sys.argv), ". End log ",dateTime))
    LOGFILE.close(); exit(0)
else:
    match_config = re.search(p_config,sys.argv[1])
    if match_config:
        configFile = sys.argv[1]
        LOGFILE.write("%s%s\n" % ("Config file is ",configFile))
    else: 
        match_input  = re.search(p_input,  sys.argv[1].lower())
        match_usage  = re.search(p_usage,  sys.argv[1].lower())
        match_detail = re.search(p_detail, sys.argv[1].lower())
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

# Open and check input file

fileError = False
try:
    CONFIG = open(configFile,"r")
except IOError as e:
    fileError = True
    print(e)

if fileError:
    print("Check your config file.")
    print(HELP_STRING)
    LOGFILE.write("%s%s\n" % ("A help string was provided to user; End log ",datetime.datetime.now()))
    LOGFILE.close(); exit(0)

##### Read input parameters from configuration file

# First, set as defaults; note: setting these values in config file is optional

geneticCode           = GENETIC_CODE_DEFAULT
geneCaller            = GENE_CALLER_DEFAULT
genomeType            = GENOME_TYPE_DEFAULT
name                  = NAME_DEFAULT
contigName            = CONTIG_NAME_DEFAULT
species               = SPECIES_DEFAULT

blastpIdentity        = BLASTP_IDENTITY_DEFAULT
blastnIdentity        = BLASTN_IDENTITY_DEFAULT
blastpHitCount        = BLASTP_HIT_COUNT_DEFAULT
blastnHitCount        = BLASTN_HIT_COUNT_DEFAULT
ncbiVirusBlast        = NCBI_VIRUS_BLAST_DEFAULT
ncbiVirusProteinBlast = NCBI_VIRUS_PROTEIN_BLAST_DEFAULT
keggVirusBlast        = KEGG_VIRUS_BLAST_DEFAULT
nrBlast               = NR_BLAST_DEFAULT
refseqProteinBlast    = REFSEQ_PROTEIN_BLAST_DEFAULT
refseqGeneBlast       = REFSEQ_GENE_BLAST_DEFAULT
phantomeBlast         = PHANTOME_BLAST_DEFAULT
pvogsBlast            = PVOGS_BLAST_DEFAULT
uniparcBlast          = UNIPARC_BLAST_DEFAULT
uniprotBlast          = UNIPROT_BLAST_DEFAULT
swissprotBlast        = SWISSPROT_BLAST_DEFAULT

hmmProgram            = HMM_PROGRAM_DEFAULT
ncbiVirusHmm          = NCBI_VIRUS_HMM_DEFAULT
ncbiVirusProteinHmm   = NCBI_VIRUS_PROTEIN_HMM_DEFAULT
keggVirusHmm          = KEGG_VIRUS_HMM_DEFAULT
nrHmm                 = NR_HMM_DEFAULT
refseqGeneHmm         = REFSEQ_GENE_HMM_DEFAULT
refseqProteinHmm      = REFSEQ_PROTEIN_HMM_DEFAULT
phantomeHmm           = PHANTOME_HMM_DEFAULT
pvogsHmm              = PVOGS_HMM_DEFAULT
uniparcHmm            = UNIPARC_HMM_DEFAULT
uniprotHmm            = UNIPROT_HMM_DEFAULT
swissprotHmm          = SWISSPROT_HMM_DEFAULT
pfamHmm               = PFAM_HMM_DEFAULT
smartHmm              = SMART_HMM_DEFAULT

genemarksCalls        = GENEMARKS_CALLS_DEFAULT
prodigalCalls         = PRODIGAL_CALLS_DEFAULT
glimmerCalls          = GLIMMER_CALLS_DEFAULT
phanotateCalls        = PHANOTATE_CALLS_DEFAULT
psatAnnotation        = PSAT_ANNOTATION_DEFAULT

# Capture user's configured values

cLines = CONFIG.read().splitlines()
for cLine in cLines:
    match_comment               = re.search(p_comment,cLine)
    match_blank                 = re.search(p_blank,cLine)
    match_outputSubdir          = re.search(p_outputSubdir,cLine)
    match_genomeFile            = re.search(p_genomeFile,cLine)
    match_psatFile              = re.search(p_psatFile,cLine)
    match_geneticCode           = re.search(p_geneticCode,cLine)
    match_translateOnly         = re.search(p_translateOnly,cLine)
    match_genomeType            = re.search(p_genomeType,cLine)
    match_name                  = re.search(p_name,cLine)
    match_contig                = re.search(p_contig,cLine)
    match_species               = re.search(p_species,cLine)

    match_geneCaller            = re.search(p_geneCaller,cLine)
    match_genemarksCalls        = re.search(p_genemarksCalls,cLine)
    match_prodigalCalls         = re.search(p_prodigalCalls,cLine)
    match_glimmerCalls          = re.search(p_glimmerCalls,cLine)
    match_phanotateCalls        = re.search(p_phanotateCalls,cLine)
    match_custom_geneCalls      = re.search(p_custom_geneCalls,cLine)
    match_custom_geneCallerName = re.search(p_custom_geneCallerName,cLine)
    match_custom_geneCallerOutfile = re.search(p_custom_geneCallerOutfile,cLine)

    match_blastpIdentity        = re.search(p_blastpIdentity,cLine)
    match_blastnIdentity        = re.search(p_blastnIdentity,cLine)
    match_blastpHitCount        = re.search(p_blastpHitCount,cLine)
    match_blastnHitCount        = re.search(p_blastnHitCount,cLine)
    match_ncbiVirusBlast        = re.search(p_ncbiVirusBlast,cLine)
    match_ncbiVirusProteinBlast = re.search(p_ncbiVirusProteinBlast,cLine)
    match_keggVirusBlast        = re.search(p_keggVirusBlast,cLine)
    match_nrBlast               = re.search(p_nrBlast,cLine)
    match_refseqProteinBlast    = re.search(p_refseqProteinBlast,cLine)
    match_refseqGeneBlast       = re.search(p_refseqGeneBlast,cLine)
    match_phantomeBlast         = re.search(p_phantomeBlast,cLine)
    match_pvogsBlast            = re.search(p_pvogsBlast,cLine)
    match_uniparcBlast          = re.search(p_uniparcBlast,cLine)
    match_uniprotBlast          = re.search(p_uniprotBlast,cLine)
    match_swissprotBlast        = re.search(p_swissprotBlast,cLine)
    match_pfamBlast             = re.search(p_pfamBlast,cLine)
    match_smartBlast            = re.search(p_smartBlast,cLine)
    match_refseqGeneBlast       = re.search(p_refseqGeneBlast,cLine)

    match_customGenomeBlast     = re.search(p_customGenomeBlast,cLine)
    match_customGeneBlast       = re.search(p_customGeneBlast,cLine)
    match_customProteinBlast    = re.search(p_customProteinBlast,cLine)
    match_customGenomeDBname    = re.search(p_customGenomeDBname,cLine)
    match_customGeneDBname      = re.search(p_customGeneDBname,cLine)
    match_customProteinDBname   = re.search(p_customProteinDBname,cLine)
    match_customGenomeDBpath    = re.search(p_customGenomeDBpath,cLine)
    match_customGeneDBpath      = re.search(p_customGeneDBpath,cLine)
    match_customProteinDBpath   = re.search(p_customProteinDBpath,cLine)

    match_hmmProgram            = re.search(p_hmmProgram,cLine)
    match_ncbiVirusHmm          = re.search(p_ncbiVirusHmm,cLine)
    match_ncbiVirusProteinHmm   = re.search(p_ncbiVirusProteinHmm,cLine)
    match_keggVirusHmm          = re.search(p_keggVirusHmm,cLine)
    match_nrHmm                 = re.search(p_nrHmm,cLine)
    match_refseqProteinHmm      = re.search(p_refseqProteinHmm,cLine)
    match_phantomeHmm           = re.search(p_phantomeHmm,cLine)
    match_pvogsHmm              = re.search(p_pvogsHmm,cLine)
    match_uniparcHmm            = re.search(p_uniparcHmm,cLine)
    match_uniprotHmm            = re.search(p_uniprotHmm,cLine)
    match_swissprotHmm          = re.search(p_swissprotHmm,cLine)
    match_pfamHmm               = re.search(p_pfamHmm,cLine)
    match_smartHmm              = re.search(p_smartHmm,cLine)
    match_refseqGeneHmm         = re.search(p_refseqGeneHmm,cLine)

    match_psatAnnotation        = re.search(p_psatAnnotation,cLine)
 
    if (match_comment or match_blank):
        continue 

    elif match_outputSubdir: #*** Note that if the output dir is not read before subdir; depends on user not changing order in config - Clean this up!
        value = match_outputSubdir.group(1)
        if value != '':
            value = value.rstrip('/')  # be sure that name of subdir ends in exactly one '/' (user might omit the slash)
            PIPELINE_OUTPUT_SUBDIR = PIPELINE_OUTPUT_DIR + value + '/'
            runLog = PIPELINE_OUTPUT_SUBDIR + "runPhATE.log"
        LOGFILE.write("%s%s\n" % ("PIPELINE_OUTPUT_SUBDIR is ",PIPELINE_OUTPUT_SUBDIR))

    elif match_genomeFile:
        value = match_genomeFile.group(1)
        if value != '':
            GENOME_FILE = value 
        LOGFILE.write("%s%s\n" % ("GENOME_FILE is ",GENOME_FILE))

    elif match_psatFile:
        value = match_psatFile.group(1)
        if value != '':
            PSAT_FILE = value 
            PSAT = True   # Yes, a psat file will be passed to subordinate code
        LOGFILE.write("%s%s\n" % ("PSAT_FILE is ",PSAT_FILE))

    elif match_geneticCode:
        value = match_geneticCode.group(1)
        if value != '':
            geneticCode = value

    elif match_translateOnly:
        value = match_translateOnly.group(1)
        if value.lower() == 'yes' or value.lower() == 'true' or value.lower() == 'on':
            TRANSLATE_ONLY = True
        elif value.lower() == 'no' or value.lower() == 'false' or value.lower() == 'off' or value == '':
            TRANSLATE_ONLY = False
        else:
            if PHATE_WARNINGS == 'True':
                print("WARNING:  Invalid string following translate_only parameter in config file:", value)
            LOGFILE.write("%s%s\n" % ("Invalid string following translate_only parameter in config file: ", value))

    elif match_geneCaller:
        value = match_geneCaller.group(1)
        if value.lower() == 'phanotate':
            geneCaller = 'phanotate'
            CONSENSUS_CALLS_FILE = 'phanotate.cgc'
        elif value.lower() == 'consensus':
            geneCaller = 'consensus'
            CONSENSUS_CALLS_FILE = 'consensus.cgc'
        elif value.lower() == 'genemarks' or value.lower() == 'genemark':
            geneCaller = 'genemarks'
            CONSENSUS_CALLS_FILE = 'genemark.cgc'
        elif value.lower() == 'glimmer2':
            geneCaller = 'glimmer2'
            CONSENSUS_CALLS_FILE = 'glimmer.cgc'
        elif value.lower() == 'glimmer3' or value.lower() == 'glimmer':
            geneCaller = 'glimmer3'
            CONSENSUS_CALLS_FILE = 'glimmer.cgc'
        elif value.lower() == 'prodigal':
            geneCaller = 'prodigal'
            CONSENSUS_CALLS_FILE = 'prodigal.cgc'
        elif value.lower() == 'rast':
            geneCaller = 'rast'
            CONSENSUS_CALLS_FILE = 'rast.cgc'
        if PHATE_MESSAGES == 'True':
            print("Gene caller has been set to:", geneCaller)

    elif match_genomeType:
        value = match_genomeType.group(1)
        if value.lower() == 'phage' or value.lower() == 'bacteriophage':
            genomeType = 'phage' 
        elif value.lower() == 'virus' or value.lower() == 'viral' or value.lower() == 'viridae':
            genomeType = 'virus'
        elif value.lower() == 'bacteria' or value.lower() == 'bacterium' or value.lower() == 'bacterial':
            genomeType = 'bacterium' 

    elif match_name:
        value = match_name.group(1)
        name = value

    elif match_contig:
        value = match_contig.group(1)
        contigName = value

    elif match_species:
        value = match_species.group(1)
        species = value

    # Gene Calls

    elif match_genemarksCalls:
        value = match_genemarksCalls.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            genemarksCalls = True
        else:
            genemarksCalls = False

    elif match_prodigalCalls:
        value = match_prodigalCalls.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            prodigalCalls = True
        else:
            prodigalCalls = False

    elif match_glimmerCalls:
        value = match_glimmerCalls.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            glimmerCalls = True
        else:
            glimmerCalls = False

    elif match_phanotateCalls:
        value = match_phanotateCalls.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            phanotateCalls = True
        else:
            phanotateCalls = False

    # BLAST

    elif match_blastpIdentity:
        value = match_blastpIdentity.group(1)
        if int(value) > int(MIN_BLASTP_IDENTITY) and int(value) <= 100:
            blastpIdentity = value

    elif match_blastnIdentity:
        value = match_blastnIdentity.group(1)
        if int(value) > int(MIN_BLASTN_IDENTITY) and int(value) <= 100:
            blastnIdentity = value

    elif match_blastpHitCount:
        value = match_blastpHitCount.group(1)
        if int(value) > 0 and int(value) <= int(MAX_BLASTP_HIT_COUNT):
            blastpHitCount = value

    elif match_blastnHitCount:
        value = match_blastnHitCount.group(1)
        if int(value) > 0 and int(value) <= int(MAX_BLASTN_HIT_COUNT):
            blastnHitCount = value

    elif match_ncbiVirusBlast:
        value = match_ncbiVirusBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            ncbiVirusBlast = True
        else:
            ncbiVirusBlast = False

    elif match_ncbiVirusProteinBlast:
        value = match_ncbiVirusProteinBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            ncbiVirusProteinBlast = True
        else:
            ncbiVirusProteinBlast = False

    elif match_keggVirusBlast:
        value = match_keggVirusBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
             keggVirusBlast = True
        else:
             keggVirusBlast = False

    elif match_nrBlast:
        value = match_nrBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
             nrBlast = True
        else:
             nrBlast = False 

    elif match_refseqProteinBlast:
        value = match_refseqProteinBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            refseqProteinBlast = True
        else:
            refseqProteinBlast = False

    elif match_refseqGeneBlast:
        value = match_refseqGeneBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            refseqGeneBlast = True
        else:
            refseqGeneBlast = False

    elif match_phantomeBlast:
        value = match_phantomeBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            phantomeBlast = True
        else:
            phantomeBlast = False 

    elif match_pvogsBlast:
        value = match_pvogsBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            pvogsBlast = True
        else:
            pvogsBlast = False

    elif match_uniparcBlast:
        value = match_uniparcBlast.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            uniparcBlast = True
        else:
            uniparcBlast = False

    elif match_uniprotBlast:
        value = match_uniprotBlast.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            uniprotBlast = True
        else:
            uniprotBlast = False

    elif match_swissprotBlast:
        value = match_swissprotBlast.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            swissprotBlast = True
        else:
            swissprotBlast = False

    elif match_pfamBlast:
        value = match_pfamBlast.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            pfamBlast = True
        else:
            pfamBlast = False

    elif match_smartBlast:
        value = match_smartBlast.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            smartBlast = True
        else:
            smartBlast = False

    elif match_customGenomeBlast:
        pass
    elif match_customGeneBlast:
        pass
    elif match_customProteinBlast:
        pass
    elif match_customGenomeDBname:
        pass
    elif match_customGeneDBname:
        pass
    elif match_customProteinDBname:
        pass
    elif match_customGenomeDBpath:
        pass
    elif match_customGeneDBpath:
        pass
    elif match_customProteinDBpath:
        pass

    # HMM

    elif match_hmmProgram:
        value = match_hmmProgram.group(1)
        if value.lower() == 'jackhmmer':
            hmmProgram = 'jackhmmer'
        else:
            if PHATE_WARNINGS == 'True':
                print("WARNING: currenly only jackhmmer hmm search is supported; running jackhmmer")
            hmmProgram = HMM_PROGRAM_DEFAULT 

    elif match_ncbiVirusHmm:
        value = match_ncbiVirusHmm.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            ncbiVirusHmm = True
        else:
            ncbiVirusHmm = False

    elif match_ncbiVirusProteinHmm:
        value = match_ncbiVirusProteinHmm.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            ncbiVirusProteinHmm = True
        else:
            ncbiVirusProteinHmm = False

    elif match_keggVirusHmm:
        value = match_keggVirusHmm.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
             keggVirusHmm = True
        else:
             keggVirusHmm = False

    elif match_nrHmm:
        value = match_nrHmm.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
             nrHmm = True
        else:
             nrHmm = False 

    elif match_refseqProteinHmm:
        value = match_refseqProteinHmm.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            refseqProteinHmm = True
        else:
            refseqProteinHmm = False

    elif match_phantomeHmm:
        value = match_phantomeHmm.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            phantomeHmm = True
        else:
            phantomeHmm = False 

    elif match_pvogsHmm:
        value = match_pvogsHmm.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            pvogsHmm = True
        else:
            pvogsHmm = False

    elif match_uniparcHmm:
        value = match_uniparcHmm.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            uniparcHmm = True
        else:
            uniparcHmm = False

    elif match_uniprotHmm:
        value = match_uniprotHmm.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            uniprotHmm = True
        else:
            uniprotHmm = False

    elif match_swissprotHmm:
        value = match_swissprotHmm.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            swissprotHmm = True
        else:
            swissprotHmm = False

    elif match_pfamHmm:
        value = match_pfamHmm.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            pfamHmm = True
        else:
            pfamHmm = False

    elif match_smartHmm:
        value = match_smartHmm.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            smartHmm = True
        else:
            smartHmm = False

    elif match_refseqGeneHmm:
        pass  # Not yet in service

    # PSAT

    elif match_psatAnnotation:
        value = match_psatAnnotation.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            psatAnnotation = True
        else:
            psatAnnotation = False

    else:
        LOGFILE.write("%s%s\n" % ("ERROR: Unrecognized line in config file: ", cLine))
        print("ERROR: unrecognized line in config file:", cLine)

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

# Create objects for passing genecall, blast, and hmm parameters to subordinate codes 

genecallParameters = {
    "genemarksCalls"        : genemarksCalls,
    "prodigalCalls"         : prodigalCalls,
    "glimmerCalls"          : glimmerCalls,
    "phanotateCalls"        : phanotateCalls,
    }

blastParameters = {
    "ncbiVirusBlast"        : ncbiVirusBlast,
    "ncbiVirusProteinBlast" : ncbiVirusProteinBlast,
    "keggVirusBlast"        : keggVirusBlast,
    "nrBlast"               : nrBlast,
    "refseqProteinBlast"    : refseqProteinBlast,
    "refseqGeneBlast"       : refseqGeneBlast,
    "phantomeBlast"         : phantomeBlast,
    "pvogsBlast"            : pvogsBlast,
    "uniparcBlast"          : uniparcBlast,
    "swissprotBlast"        : swissprotBlast,
    "pfamBlast"             : pfamBlast,
    "smartBlast"            : smartBlast,
    }

hmmParameters = {
    "hmmProgram"            : hmmProgram,
    "ncbiVirusHmm"          : ncbiVirusHmm,
    "ncbiVirusProteinHmm"   : ncbiVirusProteinHmm,
    "keggVirusHmm"          : keggVirusHmm,
    "phantomeHmm"           : phantomeHmm,
    "pvogsHmm"              : pvogsHmm,
    "swissprotHmm"          : swissprotHmm,
    "refseqProteinHmm"      : refseqProteinHmm,
    "refseqGeneHmm"         : refseqGeneHmm,
    "pfamHmm"               : pfamHmm,
    "smartHmm"              : smartHmm,
    "nrHmm"                 : nrHmm,
    "uniparcHmm"            : uniparcHmm,
    }

if DEBUG:
    print("TESTING: hmmParameters['pvogsHmm'] is ", hmmParameters['pvogsHmm'])
# Double check: issue warning if necessary, but continue processing assuming this is what the user intends.
#if GENOME_TYPE == 'PHAGE' and CONSENSUS_CALLS_FILE != 'phanotate.cgc':
if genomeType.lower() == 'phage' and CONSENSUS_CALLS_FILE != 'phanotate.cgc':
    if PHATE_WARNINGS == 'True':
        print("phate_runPipeline says, WARNING: If genome type is phage, the consensus gene-call file should be phanotate.cgc! Yours is", CONSENSUS_CALLS_FILE)
    LOGFILE.write("%s%s\n" % ("WARNING:  User has selected genome type as phage, but consensus gene-call file as ", CONSENSUS_CALLS_FILE))
    RUNLOG.write("%s%s\n" % ("WARNING:  User has selected genome type as phage, but consensus gene-call file as ", CONSENSUS_CALLS_FILE))

CONFIG.close()

# Turn PSAT back to False if the user only needs translations (even if they provided a .psat file)
if PSAT and TRANSLATE_ONLY:  
    PSAT = False

if PHATE_MESSAGES == 'True':
    print("PIPELINE_INPUT_DIR is", PIPELINE_INPUT_DIR)
    print("PIPELINE_OUTPUT_DIR is", PIPELINE_OUTPUT_DIR)
    print("PIPELINE_OUTPUT_SUBDIR is", PIPELINE_OUTPUT_SUBDIR)
    print("GENOME_FILE is", GENOME_FILE)
    print("GENE_FILE is", GENE_FILE)
    print("PROTEIN_FILE is", PROTEIN_FILE)
    print("genomeType is", genomeType) 
    print("name is", name) 
    print("contigName is", contigName)
    print("species is", species) 

    print("geneticCode is", geneticCode) 
    print("Status of boolean TRANSLATE_ONLY is", TRANSLATE_ONLY)
    print("geneCaller is", geneCaller) 
    print("genemarksCalls is", genemarksCalls)
    print("prodigalCalls is", prodigalCalls)
    print("glimmerCalls is", glimmerCalls)
    print("phanotateCalls is", phanotateCalls)
    print("CONSENSUS_CALLS_FILE is", CONSENSUS_CALLS_FILE)

    print("blastpIdentity is", blastpIdentity) 
    print("blastpHitCount is", blastpHitCount) 
    print("blastnHitCount is", blastnHitCount)
    print("ncbiVirusBlast is", ncbiVirusBlast)
    print("ncbiVirusProteinBlast is", ncbiVirusProteinBlast)
    print("keggVirusBlast is", keggVirusBlast)
    print("nrBlast is", nrBlast)
    print("refseqProteinBlast is", refseqProteinBlast)
    print("refseqGeneBlast is", refseqGeneBlast)
    print("phantomeBlast is", phantomeBlast)
    print("pvogsBlast is", pvogsBlast)
    print("swissprotBlast is", swissprotBlast)
    print("uniparcBlast is", uniparcBlast)
    print("uniprotBlast is", uniprotBlast)

    print("ncbiVirusHmm is", ncbiVirusHmm)
    print("ncbiVirusProteinHmm is", ncbiVirusProteinHmm)
    print("keggVirusHmm is", keggVirusHmm)
    print("phantomeHmm is", phantomeHmm)
    print("pvogsHmm is", pvogsHmm)
    print("swissprotHmm is", swissprotHmm)
    print("refseqProteinHmm is", refseqProteinHmm)
    print("refseqGeneHmm is", refseqGeneHmm)
    print("nrHmm is", nrHmm)

    print("psatAnnotation is", psatAnnotation)
    if PSAT:
        print("PSAT_FILE is", PSAT_FILE)
    else:
        print("PSAT_FILE was not provided.") 

RUNLOG.write("%s\n" % ("Input parameters:"))
RUNLOG.write("%s%s\n" % ("   PIPELINE_INPUT_DIR: ", PIPELINE_INPUT_DIR))
RUNLOG.write("%s%s\n" % ("   PIPELINE_OUTPUT_DIR: ", PIPELINE_OUTPUT_DIR))
RUNLOG.write("%s%s\n" % ("   PIPELINE_OUTPUT_SUBDIR: ", PIPELINE_OUTPUT_SUBDIR))
RUNLOG.write("%s%s\n" % ("   GENOME_FILE: ", GENOME_FILE))
RUNLOG.write("%s%s\n" % ("   GENE_FILE: ", GENE_FILE))
RUNLOG.write("%s%s\n" % ("   PROTEIN_FILE: ", PROTEIN_FILE))
RUNLOG.write("%s%s\n" % ("   genomeType is ",genomeType))
RUNLOG.write("%s%s\n" % ("   name is ",name))
RUNLOG.write("%s%s\n" % ("   contigName is ",contigName))
RUNLOG.write("%s%s\n" % ("   species is ",species))

RUNLOG.write("%s%s\n" % ("   geneticCode: ", geneticCode))
RUNLOG.write("%s%s\n" % ("   Status of boolean TRANSLATE_ONLY is ",TRANSLATE_ONLY))
RUNLOG.write("%s%s\n" % ("   geneCaller is ",geneCaller))
RUNLOG.write("%s%s\n" % ("   genemarksCalls is ",genemarksCalls))
RUNLOG.write("%s%s\n" % ("   prodigalCalls is ",prodigalCalls))
RUNLOG.write("%s%s\n" % ("   glimmerCalls is ",glimmerCalls))
RUNLOG.write("%s%s\n" % ("   phanotateCalls is ",phanotateCalls))
RUNLOG.write("%s%s\n" % ("   CONSENSUS_CALLS_FILE is ",CONSENSUS_CALLS_FILE))

RUNLOG.write("%s%s\n" % ("   blastpIdentity is ",blastpIdentity))
RUNLOG.write("%s%s\n" % ("   blastpHitCount is ",blastpHitCount))
RUNLOG.write("%s%s\n" % ("   blastnHitCount is ",blastnHitCount))
RUNLOG.write("%s%s\n" % ("   ncbiVirusBlast is ",ncbiVirusBlast))
RUNLOG.write("%s%s\n" % ("   ncbiVirusProteinBlast is ",ncbiVirusProteinBlast))
RUNLOG.write("%s%s\n" % ("   keggVirusBlast is ",keggVirusBlast))
RUNLOG.write("%s%s\n" % ("   nrBlast is ",nrBlast))
RUNLOG.write("%s%s\n" % ("   refseqProteinBlast is ",refseqProteinBlast))
RUNLOG.write("%s%s\n" % ("   refseqGeneBlast is ",refseqGeneBlast))
RUNLOG.write("%s%s\n" % ("   phantomeBlast is ",phantomeBlast))
RUNLOG.write("%s%s\n" % ("   pvogsBlast is ",pvogsBlast))
RUNLOG.write("%s%s\n" % ("   uniparcBlast is ",uniparcBlast))
RUNLOG.write("%s%s\n" % ("   swissprotBlast is ",swissprotBlast))

RUNLOG.write("%s%s\n" % ("   ncbiVirusHmm is ",ncbiVirusHmm))
RUNLOG.write("%s%s\n" % ("   ncbiVirusProteinHmm is ",ncbiVirusProteinHmm))
RUNLOG.write("%s%s\n" % ("   keggVirusHmm is ",keggVirusHmm))
RUNLOG.write("%s%s\n" % ("   nrHmm is ",nrHmm))
RUNLOG.write("%s%s\n" % ("   refseqProteinHmm is ",refseqProteinHmm))
RUNLOG.write("%s%s\n" % ("   refseqGeneHmm is ",refseqGeneHmm))
RUNLOG.write("%s%s\n" % ("   phantomeHmm is ",phantomeHmm))
RUNLOG.write("%s%s\n" % ("   pvogsHmm is ",pvogsHmm))
RUNLOG.write("%s%s\n" % ("   uniparcHmm is ",uniparcHmm))
RUNLOG.write("%s%s\n" % ("   swissprotHmm is ",swissprotHmm))

RUNLOG.write("%s%s\n" % ("   Status of boolean PSAT is ",PSAT))
RUNLOG.write("%s%s\n" % ("   psatAnnotation is ",psatAnnotation))
RUNLOG.write("%s%s\n" % ("   PSAT_FILE: ", PSAT_FILE))


# Open and check input file(s)

inputDir     = PIPELINE_INPUT_DIR
genomeFile   = PIPELINE_INPUT_DIR + GENOME_FILE
genecallFile = PIPELINE_OUTPUT_SUBDIR + CONSENSUS_CALLS_FILE
geneFile     = PIPELINE_OUTPUT_SUBDIR + GENE_FILE
proteinFile  = PIPELINE_OUTPUT_SUBDIR + PROTEIN_FILE
if PSAT:
    psatFile = PIPELINE_INPUT_DIR + PSAT_FILE
else:
    psatFile = ''
outputDir    = PIPELINE_OUTPUT_SUBDIR

RUNLOG.write("%s%s\n" % ("inputDir is ",    inputDir))
RUNLOG.write("%s%s\n" % ("genomeFile is ",  genomeFile))
RUNLOG.write("%s%s\n" % ("genecallFile is ",genecallFile))
RUNLOG.write("%s%s\n" % ("geneFile is ",    geneFile))
RUNLOG.write("%s%s\n" % ("proteinFile is ", proteinFile))
RUNLOG.write("%s%s\n" % ("psatFile is ",    psatFile))

# Check PSAT file

if PHATE_PROGRESS == 'True':
    print("Checking files...")
RUNLOG.write("%s\n" % ("Checking files..."))
fileError = False
if PSAT:
    try:
        PSAT_H = open(psatFile,"r")
    except IOError as e:
        fileError = True
        print(e)

    if fileError:
        print("Check your PSAT file,", psatFile)
        print(USAGE_STRING)
        LOGFILE.write("%s%s%s%s\n" % ("ERROR:  PSAT file could not be opened: ", psatFile, "; End log ", datetime.datetime.now()))
        LOGFILE.close(); exit(0)
    PSAT_H.close()

# Check genome file

fileError = False
try:
    GENOME_H = open(genomeFile,"r")
except IOError as e:
    fileError = True
    print(e) 

if fileError:
    print(USAGE_STRING)
    print("Check your genome file,", genomeFile)
    LOGFILE.write("%s%s%s%s\n" % ("ERROR:  Genome file could not be opened: ", genomeFile, "; End log ", datetime.datetime.now()))
    LOGFILE.close(); exit(0)
GENOME_H.close()

# Copy config file to the user's results directory
configSave = outputDir + configFile
command = "cp " + configFile + ' ' + configSave
os.system(command)

if PHATE_PROGRESS == 'True':
    print("Configuration complete.")

##### BEGIN MAIN ########################################################################################

##### Run Gene-calling Module

if PHATE_PROGRESS == 'True':
    print("Preparing to run genecall module...")
RUNLOG.write("%s\n" % ("Preparing to run genecall module..."))

param2 = outputDir[:-1]  # remove terminal '/' because subordinate code adds it explicitly

param3 = ''  # can't pass a dict; future re-write genecalling module as class
if genecallParameters["genemarksCalls"]:
    param3 += "genemarks_"
if genecallParameters["prodigalCalls"]:
    param3 += "prodigal_"
if genecallParameters["glimmerCalls"]:
    param3 += "glimmer_"
if genecallParameters["phanotateCalls"]:
    param3 += "phanotate_"

#param4 = "2>&1 > genecall.err"
#command = "python " + GENECALL_CODE + ' ' + genomeFile + ' ' + param2 + ' ' + param3 + ' ' + param4
command = "python " + GENECALL_CODE + ' ' + genomeFile + ' ' + param2 + ' ' + param3

if PHATE_PROGRESS == 'True':
    print("Calling the gene-call module.")
if PHATE_MESSAGES == 'True':
    print("Command is,", command)
RUNLOG.write("%s%s\n" % ("Calling the gene-call module. Command is ", command))

# OS system matters; choose alternate system call if you get error message on this line
result = os.system(command)
#result = subprocess.check_output(command,shell=True)

if PHATE_PROGRESS == 'True':
    print("Done!")
RUNLOG.write("%s%s\n" % ("Gene-call processing complete at ", datetime.datetime.now()))

##### Run Sequence Annotation Module

RUNLOG.write("%s\n" % ("Preparing to call sequence annotation module..."))
if DEBUG:
    print("Before constructing command line to invoke sequence annotation code, contigName is", contigName)

# Construct command line parameter string

if PHATE_PROGRESS == 'True':
    print("Preparing command strings for homology searches...")

# First, construct string listing the names of databases to be blasted

blastParameterString = ''
if blastParameters['ncbiVirusBlast']:
    blastParameterString += '_ncbiVirusGenome'
if blastParameters['ncbiVirusProteinBlast']:
    blastParameterString += '_ncbiVirusProtein'
if blastParameters['nrBlast']:
    blastParameterString += '_nr'
if blastParameters['keggVirusBlast']:
    blastParameterString += '_kegg'
if blastParameters['refseqProteinBlast']:
    blastParameterString += '_refseqProtein'
if blastParameters['refseqGeneBlast']:
    blastParameterString += '_refseqGene'
if blastParameters['pvogsBlast']:
    blastParameterString += '_pvogs'
if blastParameters['phantomeBlast']:
    blastParameterString += '_phantome'
if blastParameters['uniparcBlast']:
    blastParameterString += '_uniparc'
if blastParameters['swissprotBlast']:
    blastParameterString += '_swissprot'
if blastParameters['pfamBlast']:
    blastParameterString += '_pfam'
if blastParameters['smartBlast']:
    blastParameterString += '_smart'

# Construct string listing hmm program and databases to be searched

hmmParameterString = ''
if hmmParameters['hmmProgram']:
    hmmParameterString += '_program=' + hmmProgram 
if hmmParameters['ncbiVirusHmm']:
    hmmParameterString += '_ncbiVirusGenomeHmm'
if hmmParameters['ncbiVirusProteinHmm']:
    hmmParameterString += '_ncbiVirusProteinHmm'
if hmmParameters['nrHmm']:
    hmmParameterString += '_nrHmm'
if hmmParameters['keggVirusHmm']:
    hmmParameterString += '_keggHmm'
if hmmParameters['refseqProteinHmm']:
    hmmParameterString += '_refseqProteinHmm'
if hmmParameters['refseqGeneHmm']:
    hmmParameterString += '_refseqGeneHmm'
if hmmParameters['pvogsHmm']:
    hmmParameterString += '_pvogsHmm'
if hmmParameters['phantomeHmm']:
    hmmParameterString += '_phantomeHmm'
if hmmParameters['uniparcHmm']:
    hmmParameterString += '_uniparcHmm'
if hmmParameters['swissprotHmm']:
    hmmParameterString += '_swissprotHmm'
if hmmParameters['pfamHmm']:
    hmmParameterString += '_pfamHmm'
if hmmParameters['smartHmm']:
    hmmParameterString += '_smartHmm'

commandRoot1 = "python " + SEQANNOTATION_CODE + " -o " + outputDir  # code and output direction
commandRoot2 = " -G " + genomeFile        + " -g " + geneFile       + " -p " + proteinFile             # genome files
commandRoot3 = " -c " + geneCaller        + " -f " + genecallFile   + " -C " + contigName              # gene-call information
commandRoot4 = " -t " + genomeType        + " -n " + name           + " -s " + species                 # genome meta-data
commandRoot5 = " -i " + blastpIdentity    + " -j " + blastnIdentity + " -h " + blastpHitCount + " -H " + blastnHitCount # blast parameters
commandRoot6 = " -d " + blastParameterString                                                           # databases to blast against
commandRoot7 = " -m " + hmmParameterString                                                             # program and databases for hmm search
commandRoot = commandRoot1 + commandRoot2 + commandRoot3 + commandRoot4 + commandRoot5 + commandRoot6 + commandRoot7

# As appropriate, append additional parameters
if TRANSLATE_ONLY:
    command = commandRoot + " -x true "        # setting TRANSLATE_ONLY to True
elif PSAT:
    command = commandRoot + " -P " + psatFile  # including a PSAT results file
else:
    command = commandRoot 

# Communicate and execute
if PHATE_PROGRESS == 'True':
    print("Calling the sequence annotation module.")
if PHATE_MESSAGES == 'True':
    print("Command is,", command)
RUNLOG.write("%s%s\n" % ("Calling the sequence annotation module. Command is ", command))
result = os.system(command)
if PHATE_PROGRESS == 'True':
    print("Done!")
RUNLOG.write("%s%s\n" % ("Sequence annotation processing complete at ", datetime.datetime.now()))

##### CLEAN UP

if PHATE_PROGRESS == 'True':
    print("Code completed at", datetime.datetime.now())
OUTFILE.write("%s%s\n" %("Pipeline output is in output file created by code ",SEQANNOTATION_CODE))
OUTFILE.close()
LOGFILE.write("%s%s\n" % ("Code completed at ", datetime.datetime.now()))
LOGFILE.close()
RUNLOG.write("%s%s\n" % ("Execution complete at ",datetime.datetime.now()))
RUNLOG.close()
#logging.info('PhATE execution complete')

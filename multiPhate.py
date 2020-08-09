#!/usr/bin/env python

################################################################
#
# Program Title:  multiPhate2.py (/multiPhate2/)
#
# Programmer:  Carol L. Ecale Zhou
# Last Update:  08 August 2020
#
# Description: Script multiPhate.py runs an annotation pipeline (phate_runPipeline.py) over any
#    number of genomes specified in the user's input configuration file (multPhate.config). It then
#    runs a code (CGP) to compare the proteomes of the user's genomes.
#    Specifically:
#      1) multiPhate.py runs the annotation pipeline over each input phage genome. The first step
#         involves gene calling using up to 4 gene callers plus an optional custom gene caller's
#         output. The results of gene calling are compared, and several outputs are generated:
#         gene calls for each of the (up to 5) gene callers, the superset computed from all gene
#         callers, a consensus gene-call set, and a common-core gene call set. The second step
#         involves performing searches against specified databases to identify homologies for each
#         predicted peptide from each genome. Codes that may be selected are: blastn (genome and gene),
#         blastp, phmmer, jackhmmer, hmmscan. Databases are appropriate to the code (ie, blast,
#         sequence, or hmm profile). Homology search results are combined into tabbed and gff outfiles.
#      2) multiPhate.py runs the comparative genomics module to compare the predicted genes and proteins
#         of the user's genomes. CGP identifies a reference genome (the first listed in the user
#         config), and lists all corresponding predicted genes and proteins from each of the other genomes.
#         Lastly, it identifies corresponding genes and proteins with respect to the reference genome,
#         and creates homology groups.
#
# Programmer's notes:
#
# Usage:  python multiPhate.py myMultiPhate.config
#    (see sample_multiPhate.config for how to create your configuration file)
#
################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.PDF FOR DETAILS.

# DO NOT MODIFY ANYTHING IN THIS FILE EXCEPT ITEMS LABELED AS "USER CONFIGURATION"

import sys, os, re, string, copy, time, datetime
import subprocess
import json
from multiprocessing import Pool

# Because grep works differently in Python 3.x, and works differently on Mac OSX, need to set this environment
# variable so that grep works correctly in SequenceAnnotation/phate_annotation.py:
# Set MAC_OSX either to True or False by commenting in/out accordingly:
#MAC_OSX = True
MAC_OSX = False

MAC_OSX_STRING = ""
if MAC_OSX:
    MAC_OSX_STRING = 'True'
else:
    MAC_OSX_STRING = 'False'

# Parallelism is now controlled in the user's configuration file, input to multiPhATE. 
# Below are the default values for parameters hpc, threads, and blastThreads.
# 1) If you are running multiPhATE on a high-performance computing (HPC) system (e.g., using SLURM), you will need to mute the multiPhate log.
# Set HPC = True to prevent multiPhATE from writing the multiPhate.log file, as each process attempting to write to this one file will cause
# contention for I/O.
HPC = False  # write to log
# HPC = True  # mute the log
# Set THREADS to 'ALL' to use all available processors, or to an integer to limit number of parallel processes
#THREADS = 'ALL'
THREADS = 1
# Set the number of blast threads that blast+ will invoke
BLAST_THREADS = 1
# 2) If you are running under a linux system, set PHATE_OUT and PHATE_ERR to 'True'. This will capture standard errors to files. Cannot
# guarantee this will work under other operating systems.
PHATE_OUT = 'False'
PHATE_ERR = 'True'
# The following boolean will remove accumulated Results_.. directories from previous multiPhATE processing.
CLEAN_PREVIOUS_CGP_RESULTS = True
#
# Default Verbosity; These are normally set as environment variables based on parameters in the config file,
# but defaults take effect if not specified in config.
# Note: defaults are needed to set environmental parameters in case user config file does not specify.
CLEAN_RAW_DATA_DEFAULT   = 'True'   # if 'False', the raw Blast and Hmm outputs will be saved in the PipelineOutput folder
PHATE_WARNINGS_DEFAULT   = 'False'  # multiPhate and phate codes
PHATE_MESSAGES_DEFAULT   = 'False'
PHATE_PROGRESS_DEFAULT   = 'False'
#*** REVISIT THESE CONTROLS: eliminate module-specific controls; subordinate modules should follow PHATE settings
CGC_WARNINGS_DEFAULT     = 'False'  # compareGeneCalls codes
CGC_MESSAGES_DEFAULT     = 'False'
CGC_PROGRESS_DEFAULT     = 'False'
CGP_WARNINGS_DEFAULT     = 'False'  # compareGeneProfiles codes
CGP_MESSAGES_DEFAULT     = 'False'
CGP_PROGRESS_DEFAULT     = 'False'

CLEAN_RAW_DATA = CLEAN_RAW_DATA_DEFAULT
PHATE_PROGRESS = PHATE_PROGRESS_DEFAULT
PHATE_MESSAGES = PHATE_MESSAGES_DEFAULT
PHATE_WARNINGS = PHATE_WARNINGS_DEFAULT
CGP_PROGRESS   = CGP_PROGRESS_DEFAULT

DEBUG = False
#DEBUG = True     # Controls debug settings in this (local) code only

##### SET DEFAULTS for USER-DEFINED PARAMETERS

# general
genomeType               = 'phage'
genomeName               = 'unknown'
genomeSpecies            = 'phage'
geneticCode              = 11  # for translation

# gene calling
primaryCalls             = 'phanotate'
primaryCallsFile         = 'phanotate.cgc'
phanotateCalls           = False
prodigalCalls            = False
glimmerCalls             = False
genemarksCalls           = False
customGeneCalls          = False
customGeneCallerName     = 'custom'
customGeneCallerOutfile  = 'custom.gff'
translateOnly            = True

# parameters and booleans controlling blast processes
blastpIdentity           = 60
blastnIdentity           = 60
blastpHitCount           = 5
blastnHitCount           = 5
ncbiVirusGenomeBlast     = False
ncbiVirusProteinBlast    = False
keggVirusBlast           = False
refseqProteinBlast       = False
refseqGeneBlast          = False
phantomeBlast            = False
pvogsBlast               = False
vogsBlast                = False
vogGeneBlast             = False
vogProteinBlast          = False
pfamBlast                = False
smartBlast               = False
uniprotBlast             = False
swissprotBlast           = False
phageEnzymeBlast         = False
nrBlast                  = False
cazyBlast                = False

# parameters controlling custom blast processes: booleans, custom names, and paths
customGenomeBlast        = False
customGeneBlast          = False
customProteinBlast       = False
customGenomeDBname       = 'unknown'
customGeneDBname         = 'unknown'
customProteinDBname      = 'unknwon'
customGenomeDBpath       = 'pathNotSet'
customGeneDBpath         = 'pathNotSet'
customProteinDBpath      = 'pathNotSet'

# programs, databases, and booleans controlling hmm processes
blastpSearch             = False  # if True, run blastp against fasta blast database(s)
phmmerSearch             = False  # if True, run phmmer against fasta database(s)
jackhmmerSearch          = False  # if True, run jackhmmer against fasta blast database(s)
hmmscan                  = False
hmmbuild                 = False
hmmsearch                = False
runCGP                   = False
runGenomics              = False
# databases to be used (booleans)
ncbiVirusGenomeHmm       = False
ncbiVirusProteinHmm      = False
keggVirusHmm             = False
nrHmm                    = False
refseqGeneHmm            = False
refseqProteinHmm         = False
phantomeHmm              = False
pvogsHmm                 = False
vogsHmm                  = False
uniparcHmm               = False
uniprotHmm               = False
swissprotHmm             = False
pfamHmm                  = False
smartHmm                 = False
phageEnzymeHmm           = False
# database locations
ncbiVirusGenomeHmmDBpath  = "pathNotSet"
ncbiVirusProteinHmmDBpath = "pathNotSet"
keggVirusHmmDBpath        = "pathNotSet"
nrHmmDBpath               = "pathNotSet"
refseqGeneHmmDBpath       = "pathNotSet"
refseqProteinHmmDBpath    = "pathNotSet"
phantomeHmmDBpath         = "pathNotSet"
pvogsHmmDBpath            = "pathNotSet"
vogsHmmDBpath             = "pathNotSet"
uniprotHmmDBpath          = "pathNotSet"
swissprotHmmDBpath        = "pathNotSet"
pfamHmmDBpath             = "pathNotSet"
smartHmmDBpath            = "pathNotSet"
phageEnzymeHmmDBpath      = "pathNotSet"

# parameters controlling custom hmm processes: booleans, custom names, and paths
customHmm                = False
customHmmName            = 'unknown'
customHmmDBname          = 'unknown'
customHmmDBpath          = 'pathNotSet'

# parameters controlling parallelism
threads                  = THREADS        # Value should be a positive integer or ALL
hpc                      = HPC            # True or False
blastThreads             = BLAST_THREADS  # positive integer

# Constants; defaults will apply if not specified in config file
# Leave all this stuff alone!

# Standard directories
BASE_DIR_DEFAULT            = os.path.join(os.getcwd(),"")    # Ex: /Home/MyName/MyCodeDirectory/multiPhATE2/
DATABASE_DIR_DEFAULT        = BASE_DIR_DEFAULT + "Databases/"
SOFTWARE_DIR_DEFAULT        = BASE_DIR_DEFAULT + "ExternalCodes/"
PIPELINE_INPUT_DIR_DEFAULT  = BASE_DIR_DEFAULT + "PipelineInput/"
PIPELINE_OUTPUT_DIR_DEFAULT = BASE_DIR_DEFAULT + "PipelineOutput/"
CGP_DIR_DEFAULT             = BASE_DIR_DEFAULT + "CompareGeneProfiles/"
CGP_RESULTS_DIR             = PIPELINE_OUTPUT_DIR_DEFAULT + "CGP_RESULTS/"
GENOMICS_DIR_DEFAULT        = BASE_DIR_DEFAULT + "Genomics/"
GENOMICS_RESULTS_DIR        = PIPELINE_OUTPUT_DIR_DEFAULT + "GENOMICS_RESULTS/"
JSON_DIR                    = BASE_DIR_DEFAULT + "JSON/"

PHATE_PIPELINE_CODE             = 'phate_runPipeline.py' # The annotaton engine
GENE_FILE                       = 'gene.fnt'             # default filename where gene sequences are written, per genome's PipelineOutput/
PROTEIN_FILE                    = 'protein.faa'          # default filename where protein sequences are written, per genome's PipelineOutput/
CGP_CODE_NAME                   = 'cgp_driver.py'        # top-level, driver program for running CompareGeneProfiles pipeline
CGP_CODE                        = CGP_DIR_DEFAULT + CGP_CODE_NAME # absolute path of top-level, driver for CompareGeneProfiles pipeline
GENOMICS_CODE_NAME              = 'genomics_driver.py'   # top-level, driver program for running Genomics analysis
GENOMICS_CODE                   = GENOMICS_DIR_DEFAULT + GENOMICS_CODE_NAME # absolute path of top-level, drive for Genomics analysis
VOG_PROTEIN_HEADER_FILENAME     = "vog.protein.headers.lst"  # The headers from the vog.proteins.tagged.all.fa file (computed by multiPhate.py)
VOG_PROTEIN_ANNOTATION_FILENAME = "vog.annotations.tsv"  # The annotations associated with VOG identifiers (downloaded from VOG server)

# naming the custom gene caller
# paths to subordinate codes; '' if installed globally (e.g., via conda)
GENEMARKS_HOME              = '' # Available via license
GLIMMER_HOME                = '' # Can install using Conda
PRODIGAL_HOME               = '' # Can install using Conda
PHANOTATE_HOME              = '' # not currently a conda package; acquire from github
HMMER_HOME                  = '' # Can install using Conda; HMMER includes jackhmmer, hmmbuild, hmmscan

# blast parameters
BLASTP_IDENTITY_DEFAULT     = 60
BLASTN_IDENTITY_DEFAULT     = 60
BLASTP_HIT_COUNT_DEFAULT    = 5
BLASTN_HIT_COUNT_DEFAULT    = 5
MIN_BLASTP_IDENTITY         = 5   # default; sets a lower limit based on value at which a structure model can provide information
MIN_BLASTN_IDENTITY         = 5   # default; sets a lower limit based on value at which a structure model can provide information
MAX_BLASTP_HIT_COUNT        = 100 # default; sets an upper limit; user's value should typically be well below this
MAX_BLASTN_HIT_COUNT        = 100 # default; sets an upper limit
SCORE_EDGE_MAX              = 1.0 # currently fixed setting for blast; modify here if need be
OVERHANG_MAX                = 100 # currently fixed setting for blast; modify here if need be
HIT_COUNT_MAX               = 100 #

# ENVIRONMENT VARIABLES
# It is most convenient to locate the supporting software codes and databases in the above-indicated subdirectories.
# However, if any of your supporting databases or softwares reside elsewhere, then explicit locations will need to
# be filled in in the multiPhate.config file. This will likely be the case for large databases that you may already
# have on your compute cluster (e.g, NR), and for software packages, such as EMBOSS or gene finders that you may
# already have installed on your system. Parameters that differ from defaults will be re-assigned based on information
# provided in the users' multiPhate.config file.

os.environ["PHATE_MAC_OSX"]                 = MAC_OSX_STRING

PIPELINE_INPUT_DIR                          = PIPELINE_INPUT_DIR_DEFAULT   # Default
PIPELINE_OUTPUT_DIR                         = PIPELINE_OUTPUT_DIR_DEFAULT  # Default
PHATE_BASE_DIR                              = BASE_DIR_DEFAULT
CGP_DIR                                     = CGP_DIR_DEFAULT
os.environ["CGP_CODE_BASE_DIR"]             = CGP_DIR
os.environ["PHATE_BASE_DIR"]                = BASE_DIR_DEFAULT
os.environ["PHATE_DATABASE_DIR"]            = DATABASE_DIR_DEFAULT
os.environ["PHATE_SOFTWARE_DIR"]            = SOFTWARE_DIR_DEFAULT
os.environ["PHATE_PIPELINE_INPUT_DIR"]      = PIPELINE_INPUT_DIR
os.environ["PHATE_PIPELINE_OUTPUT_DIR"]     = PIPELINE_OUTPUT_DIR
os.environ["PHATE_PHATE_BASE_DIR"]          = PHATE_BASE_DIR
os.environ["PHATE_PIPELINE_DIR"]            = BASE_DIR_DEFAULT
os.environ["PHATE_CGP_RESULTS_DIR"]         = CGP_RESULTS_DIR
os.environ["PHATE_GENOMICS_RESULTS_DIR"]    = GENOMICS_RESULTS_DIR

os.environ["PHATE_CUSTOM_GENECALLER_NAME"]  = ""  # set from user's config file

# Gene calling and other codes
# if installed globally
#os.environ["PHATE_PHANOTATE_PATH"]          = SOFTWARE_DIR_DEFAULT + "PHANOTATE/PHANOTATE-master/"
os.environ["PHATE_PHANOTATE_PATH"]          = ""   # should be installed globally now
os.environ["PHATE_PRODIGAL_PATH"]           = ""   # global, if installed via conda
os.environ["PHATE_GLIMMER_PATH"]            = ""   # global, if installed via conda
os.environ["PHATE_GENEMARKS_PATH"]          = ""   # global, but not a conda package
#os.environ["PHATE_CGC_PATH"]                = ""    # this code is not installed globally
os.environ["PHATE_tRNAscanSE_HOME"]         = ""    # global, if installed via conda
os.environ["PHATE_EMBOSS_PHATE_HOME"]       = ""    # global, if installed via conda
# if installed locally
#os.environ["PHATE_PHANOTATE_PATH"]          = SOFTWARE_DIR_DEFAULT + "Phanotate/phanotate/"
#os.environ["PHATE_PRODIGAL_PATH"]           = SOFTWARE_DIR_DEFAULT + "prodigal.v2_50/"
#os.environ["PHATE_GLIMMER_PATH"]            = SOFTWARE_DIR_DEFAULT + "glimmer3.02/bin/"
#os.environ["PHATE_GENEMARKS_PATH"]          = SOFTWARE_DIR_DEFAULT + "GeneMarkS/genemark_suite_linux_64/gmsuite/"
os.environ["PHATE_CGC_PATH"]                = BASE_DIR_DEFAULT + "CompareCalls/"
#os.environ["PHATE_tRNAscanSE_HOME"]         = /Users/myName/tRNAscanDir/trnascan-se-2.0/bin/tRNAscan-SE"
#EMBOSS_CODE                                 = "emboss"  # Modify this for the version you have--only if using non-global version
#os.environ["PHATE_EMBOSS_PHATE_HOME"]       = SOFTWARE_DIR_DEFAULT + EMBOSS_CODE

# Data sets (Subject to change based on user input to config file)
# BLAST
os.environ["PHATE_NCBI_VIRUS_BASE_DIR"]             = DATABASE_DIR_DEFAULT + "NCBI/"
os.environ["PHATE_NCBI_VIRUS_GENOME_BLAST_HOME"]    = os.environ["PHATE_NCBI_VIRUS_BASE_DIR"] + "Virus_Genome/"  + "viral.genomic.fna"
os.environ["PHATE_NCBI_VIRUS_PROTEIN_BLAST_HOME"]   = os.environ["PHATE_NCBI_VIRUS_BASE_DIR"] + "Virus_Protein/" + "viral.protein.faa"
os.environ["PHATE_NCBI_TAXON_DIR"]                  = os.environ["PHATE_NCBI_VIRUS_BASE_DIR"] + "Virus_Genome/"
os.environ["PHATE_REFSEQ_PROTEIN_BASE_DIR"]         = DATABASE_DIR_DEFAULT + "Refseq/Protein/"
os.environ["PHATE_REFSEQ_PROTEIN_BLAST_HOME"]       = os.environ["PHATE_REFSEQ_PROTEIN_BASE_DIR"] + "refseq_protein"
os.environ["PHATE_REFSEQ_GENE_BASE_DIR"]            = DATABASE_DIR_DEFAULT + "Refseq/Gene/"
os.environ["PHATE_REFSEQ_GENE_BLAST_HOME"]          = os.environ["PHATE_REFSEQ_GENE_BASE_DIR"] + "refseqgene"
os.environ["PHATE_PVOGS_BASE_DIR"]                  = DATABASE_DIR_DEFAULT + "pVOGs/"
os.environ["PHATE_PVOGS_BLAST_HOME"]                = os.environ["PHATE_PVOGS_BASE_DIR"] + "pVOGs.faa"
os.environ["PHATE_PVOGS_HEADER_FILE"]               = os.environ["PHATE_PVOGS_BASE_DIR"] + "pVOGs.headers.lst"
os.environ["PHATE_VOGS_BASE_DIR"]                   = DATABASE_DIR_DEFAULT + "VOGs/"
os.environ["PHATE_VOGS_BLAST_HOME"]                 = os.environ["PHATE_VOGS_BASE_DIR"] + "VOGs.faa"
os.environ["PHATE_VOGS_ANNOTATION_FILE"]            = os.environ["PHATE_VOGS_BASE_DIR"] + "vog.annotations.tsv"
os.environ["PHATE_VOG_GENE_BLAST_HOME"]             = os.environ["PHATE_VOGS_BASE_DIR"] + "VOG_genes.fnt"      #*** FIX
os.environ["PHATE_VOG_PROTEIN_BLAST_HOME"]          = os.environ["PHATE_VOGS_BASE_DIR"] + "VOG_protein.faa"
os.environ["PHATE_VOG_PROTEIN_BASE_DIR"]            = os.environ["PHATE_VOGS_BASE_DIR"]
os.environ["PHATE_VOG_PROTEIN_HEADERS_FILE"]        = ""
os.environ["PHATE_VOG_PROTEIN_ANNOTATION_FILE"]     = ""
os.environ["PHATE_PHANTOME_BASE_DIR"]               = DATABASE_DIR_DEFAULT + "Phantome/"
os.environ["PHATE_PHANTOME_BLAST_HOME"]             = os.environ["PHATE_PHANTOME_BASE_DIR"] + "Phantome_Phage_genes.faa"
os.environ["PHATE_PHAGE_ENZYME_BASE_DIR"]           = DATABASE_DIR_DEFAULT + "PhageEnzyme/" # not yet in service
os.environ["PHATE_PHAGE_ENZYME_BLAST_HOME"]         = os.environ["PHATE_PHAGE_ENZYME_BASE_DIR"] + "phageEnzyme"
os.environ["PHATE_KEGG_VIRUS_BASE_DIR"]             = DATABASE_DIR_DEFAULT + "KEGG/"
os.environ["PHATE_KEGG_VIRUS_BLAST_HOME"]           = os.environ["PHATE_KEGG_VIRUS_BASE_DIR"] + "T40000.pep"
os.environ["PHATE_PFAM_BASE_DIR"]                   = DATABASE_DIR_DEFAULT + "Pfam/" # not yet in service
os.environ["PHATE_PFAM_BLAST_HOME"]                 = os.environ["PHATE_PFAM_BASE_DIR"] + "pfam"
os.environ["PHATE_SMART_BASE_DIR"]                  = DATABASE_DIR_DEFAULT + "Smart/"
os.environ["PHATE_SMART_BLAST_HOME"]                = os.environ["PHATE_SMART_BASE_DIR"] + "smart"
os.environ["PHATE_SWISSPROT_BASE_DIR"]              = DATABASE_DIR_DEFAULT + "Swissprot/"
os.environ["PHATE_SWISSPROT_BLAST_HOME"]            = os.environ["PHATE_SWISSPROT_BASE_DIR"] + "swissprot"
os.environ["PHATE_UNIPROT_BASE_DIR"]                = DATABASE_DIR_DEFAULT + "Uniprot/" # not yet in service
os.environ["PHATE_UNIPROT_BLAST_HOME"]              = os.environ["PHATE_UNIPROT_BASE_DIR"] + "uniprot"
os.environ["PHATE_NR_BLAST_BASE_DIR"]               = DATABASE_DIR_DEFAULT + "NR/"
os.environ["PHATE_NR_BLAST_HOME"]                   = os.environ["PHATE_NR_BLAST_BASE_DIR"] + "nr"
os.environ["PHATE_CAZY_BLAST_BASE_DIR"]             = DATABASE_DIR_DEFAULT + "CAZY/"
os.environ["PHATE_CAZY_BLAST_HOME"]                 = os.environ["PHATE_CAZY_BLAST_BASE_DIR"] + "cazy.faa"
os.environ["PHATE_CAZY_ANNOTATION_PATH"]            = ""  # read from config file
os.environ["PHATE_CUSTOM_GENOME_BLAST_HOME"]        = ""  # read from config file
os.environ["PHATE_CUSTOM_GENE_BLAST_HOME"]          = ""  # read from config file
os.environ["PHATE_CUSTOM_PROTEIN_BLAST_HOME"]       = ""  # read from config file

# HMM (Subject to change based on user input to config file)
os.environ["PHATE_NCBI_VIRUS_PROTEIN_HMM_BASE_DIR"] = DATABASE_DIR_DEFAULT + "NCBI/Virus_Protein/" # not yet in service
os.environ["PHATE_NCBI_VIRUS_PROTEIN_HMM_HOME"]     = os.environ["PHATE_NCBI_VIRUS_PROTEIN_HMM_BASE_DIR"] + "unknown"
os.environ["PHATE_NCBI_VIRUS_GENOME_HMM_BASE_DIR"]  = DATABASE_DIR_DEFAULT + "NCBI/Virus_Genome/" # not yet in service
os.environ["PHATE_NCBI_VIRUS_GENOME_HMM_HOME"]      = os.environ["PHATE_NCBI_VIRUS_GENOME_HMM_BASE_DIR"] + "unknown"
os.environ["PHATE_REFSEQ_PROTEIN_HMM_BASE_DIR"]     = DATABASE_DIR_DEFAULT + "Refseq/Protein_hmm/" # not yet in service
os.environ["PHATE_REFSEQ_PROTEIN_HMM_HOME"]         = os.environ["PHATE_REFSEQ_PROTEIN_HMM_BASE_DIR"] + "unknown"
os.environ["PHATE_REFSEQ_GENE_HMM_BASE_DIR"]        = DATABASE_DIR_DEFAULT + "Refseq/Gene_hmm/" # not yet in service
os.environ["PHATE_REFSEQ_GENE_HMM_HOME"]            = os.environ["PHATE_REFSEQ_GENE_HMM_BASE_DIR"] + "unknown"
os.environ["PHATE_PVOGS_HMM_BASE_DIR"]              = DATABASE_DIR_DEFAULT + "pVOGhmms/"
os.environ["PHATE_PVOGS_HMM_HOME"]                  = os.environ["PHATE_PVOGS_HMM_BASE_DIR"] + "AllvogHMMprofiles/"
os.environ["PHATE_VOGS_HMM_BASE_DIR"]               = DATABASE_DIR_DEFAULT + "VOGhmms/"
os.environ["PHATE_VOGS_HMM_HOME"]                   = os.environ["PHATE_VOGS_HMM_BASE_DIR"] + "VOGsHMMdb/"
os.environ["PHATE_PHANTOME_HMM_BASE_DIR"]           = DATABASE_DIR_DEFAULT + "PHANTOME_hmm/" # not yet in service
os.environ["PHATE_PHANTOME_HMM_HOME"]               = os.environ["PHATE_PHANTOME_HMM_BASE_DIR"] + "unknown" # not yet in service
os.environ["PHATE_PHAGE_ENZYME_HMM_BASE_DIR"]       = DATABASE_DIR_DEFAULT + "phageEnzyme_hmm/"  # not yet in service
os.environ["PHATE_PHAGE_ENZYME_HMM_HOME"]           = os.environ["PHATE_PHAGE_ENZYME_HMM_BASE_DIR"] + "unknown" # not yet in service
os.environ["PHATE_KEGG_VIRUS_HMM_BASE_DIR"]         = DATABASE_DIR_DEFAULT + "KEGG_hmm/" # not yet in service
os.environ["PHATE_KEGG_VIRUS_HMM_HOME"]             = os.environ["PHATE_KEGG_VIRUS_HMM_BASE_DIR"] + "unknown"
os.environ["PHATE_PFAM_HMM_BASE_DIR"]               = DATABASE_DIR_DEFAULT + "Pfam_hmm/"  # not yet in service
os.environ["PHATE_PFAM_HMM_HOME"]                   = os.environ["PHATE_PFAM_HMM_BASE_DIR"] + "unknown" # not yet in service
os.environ["PHATE_SMART_HMM_BASE_DIR"]              = DATABASE_DIR_DEFAULT + "SMART_hmm/" # not yet in service
os.environ["PHATE_SMART_HMM_HOME"]                  = os.environ["PHATE_SMART_HMM_BASE_DIR"] + "unknown" # not yet in service
os.environ["PHATE_SWISSPROT_HMM_BASE_DIR"]          = DATABASE_DIR_DEFAULT + "Swissprot_hmm/"
os.environ["PHATE_SWISSPROT_HMM_HOME"]              = os.environ["PHATE_SWISSPROT_HMM_BASE_DIR"] + "unknown" # not yet in service
os.environ["PHATE_UNIPROT_HMM_BASE_DIR"]            = DATABASE_DIR_DEFAULT + "UNIPROT_hmm/" # not yet in service
os.environ["PHATE_UNIPROT_HMM_HOME"]                = os.environ["PHATE_UNIPROT_HMM_BASE_DIR"] + "unknown" # not yet in service
os.environ["PHATE_NR_HMM_BASE_DIR"]                 = DATABASE_DIR_DEFAULT + "NR_hmm/"  # not yet in service
os.environ["PHATE_NR_HMM_HOME"]                     = os.environ["PHATE_NR_HMM_BASE_DIR"] + "unknown" # not yet in service
os.environ["PHATE_CUSTOM_HMM_BASE_DIR"]             = DATABASE_DIR_DEFAULT + "custom_hmm/"  # not yet in service
#os.environ["PHATE_CUSTOM_HMM_HOME"]                 = os.environ["PHATE_CUSTOM_HMM_BASE_DIR"] + "unknown" # not yet in service
os.environ["PHATE_CUSTOM_HMM_HOME"]                 = ""   # read from config file

# BLAST
#os.environ["PHATE_BLAST_HOME"]                      = SOFTWARE_DIR_DEFAULT + "/ncbi-blast-2.7.1+/bin/"
os.environ["PHATE_BLAST_HOME"]                      = ""    # global, if installed via conda; use blast+, not legacy blast
os.environ["PHATE_MIN_BLASTP_IDENTITY"]             = str(MIN_BLASTP_IDENTITY)
os.environ["PHATE_MIN_BLASTN_IDENTITY"]             = str(MIN_BLASTN_IDENTITY)
os.environ["PHATE_MAX_BLASTP_HIT_COUNT"]            = str(MAX_BLASTP_HIT_COUNT)
os.environ["PHATE_MAX_BLASTN_HIT_COUNT"]            = str(MAX_BLASTN_HIT_COUNT)
os.environ["PHATE_BLASTP_IDENTITY_DEFAULT"]         = str(BLASTP_IDENTITY_DEFAULT)
os.environ["PHATE_BLASTN_IDENTITY_DEFAULT"]         = str(BLASTN_IDENTITY_DEFAULT)
os.environ["PHATE_BLASTP_HIT_COUNT_DEFAULT"]        = str(BLASTP_HIT_COUNT_DEFAULT)
os.environ["PHATE_BLASTN_HIT_COUNT_DEFAULT"]        = str(BLASTN_HIT_COUNT_DEFAULT)
os.environ["PHATE_SCORE_EDGE_MAX"]                  = str(SCORE_EDGE_MAX)
os.environ["PHATE_OVERHANG_MAX"]                    = str(OVERHANG_MAX)
os.environ["PHATE_HIT_COUNT_MAX"]                   = str(HIT_COUNT_MAX)  #*** This is the hit count actually used in phate_blast

# CODES - These should be installed globally (via Conda), but if not, locations go here
os.environ["PHATE_BLAST_HOME"]                      = ""
os.environ["PHATE_EMBOSS_HOME"]                     = ""
os.environ["PHATE_PHANOTATE_HOME"]                  = ""
os.environ["PHATE_PRODIGAL_HOME"]                   = ""
os.environ["PHATE_GLIMMER_HOME"]                    = ""
os.environ["PHATE_GENEMARKS_HOME"]                  = ""
os.environ["PHATE_CGC_HOME"]                        = ""
os.environ["PHATE_tRNAscanSE_HOME"]                 = ""
os.environ["PHATE_HMMER_HOME"]                      = ""

# Global control: verbosity and error capture
os.environ["PHATE_CLEAN_RAW_DATA"]                  = CLEAN_RAW_DATA_DEFAULT
os.environ["PHATE_PHATE_WARNINGS"]                  = PHATE_WARNINGS_DEFAULT  # Print warnings and errors to standard out
os.environ["PHATE_PHATE_MESSAGES"]                  = PHATE_MESSAGES_DEFAULT  # Print helpful messages (may be verbose)
os.environ["PHATE_PHATE_PROGRESS"]                  = PHATE_PROGRESS_DEFAULT  # Print each step in processing a genome
os.environ["PHATE_PHATE_ERR"]                       = PHATE_ERR               # Capture standard errors to files on linux/mac machine
os.environ["PHATE_PHATE_OUT"]                       = PHATE_OUT               # Capture standard errors to files on linux/mac machine
os.environ["PHATE_CGC_WARNINGS"]                    = CGC_WARNINGS_DEFAULT
os.environ["PHATE_CGC_MESSAGES"]                    = CGC_MESSAGES_DEFAULT
os.environ["PHATE_CGC_PROGRESS"]                    = CGC_PROGRESS_DEFAULT
os.environ["PHATE_CGP_WARNINGS"]                    = CGP_WARNINGS_DEFAULT
os.environ["PHATE_CGP_MESSAGES"]                    = CGP_MESSAGES_DEFAULT
os.environ["PHATE_CGP_PROGRESS"]                    = CGP_PROGRESS_DEFAULT

# Constants

CODE_BASE          = "multiPhate"
CODE               = CODE_BASE + ".py"
CONFIG_FILE        = "multiPhate.config"     # by default, but user should name their own, ending in ".config"
SAMPLE_CONFIG_FILE = "sample_" + CONFIG_FILE # Sample config file; user should copy, then modify.
CGP_CONFIG_FILE    = PIPELINE_OUTPUT_DIR + "cgpNxN.config"

# HELP STRINGS

HELP_STRING = """This code, """ + CODE + """, runs the phage annotation pipeline (phate_runPipeline.py) over multiple genomes. The configuration file input to this code specifies a list of genomes to be processed and the parameters for pipeline execution. The pipeline performs 1) gene calling by 4 gene callers (PHANOTATE, GeneMarkS, Glimmer3, and Prodigal), followed by identification of closest phage genome by means of blast against an NCBI-phage database, and sequence-based functional annotation by means of blastp against several peptide databases (NR, NCBI virus protein, KEGG-virus, Phantome, pVOGs, Swissprot, Refseq protein), and HMM search against the pVOG database. \nType: python """ + CODE + """ usage - for more information about constructing the command line.\nType: python """ + CODE + """ input - for more information about how this code can be run, or consult the README file for extensive information about installation and code usage.\n"""

INPUT_STRING = """The input files and other parameters for running this code are specified in a configuration file, which is provided as the only input parameter. See sample configuration file (""" + SAMPLE_CONFIG_FILE + """) for details on how to customize your configuration file. Copy that file and then modify accordingly.\n"""

USAGE_STRING = """Usage: python """ + CODE + """ """ + CONFIG_FILE + """\n"""

DETAIL_STRING = """For further details please consult the README file.\n"""

##### PATTERNS #####

# LOCATIONS
p_phateDir                    = re.compile("phate_dir='(.*)'")
p_databaseDir                 = re.compile("database_dir='(.*)'")
p_softwareDir                 = re.compile("software_dir='(.*)'")

# GENERAL
p_comment                     = re.compile("^#")
p_blank                       = re.compile("^$")
p_help                        = re.compile("help")
p_input                       = re.compile("input")
p_usage                       = re.compile("usage")
p_detail                      = re.compile("detail")
p_config                      = re.compile("config")
p_outputSubdir                = re.compile("output_subdir='(.*)'")
p_genomeFile                  = re.compile("genome_file='(.*)'")
p_genomeType                  = re.compile("genome_type='(.*)'")
p_genomeName                  = re.compile("genome_name='(.*)'")
p_genomeSpecies               = re.compile("genome_species='(.*)'")

# GENOME INFORMATION
p_genomeList                  = re.compile("Genome\sList")       # non-case-sensitive "Genome List"
p_genomeNumber                = re.compile("Genome\s+(\d+)")     # genome number
p_root                        = re.compile("([\w\d_-]+)\.fasta") # captures the root name of the fasta file (e.g., takes 'P2' from P2.fasta)
p_end                         = re.compile("END")

# GENE CALLING

# Gene call parameters (value|true/false)
p_primaryCalls                = re.compile("primary_calls='(.*)'")
p_geneticCode                 = re.compile("genetic_code='(\d+)'")
p_translateOnly               = re.compile("translate_only='(.*)'")    # true/false
# Gene calling to be performed (true/false)
p_genemarksCalls              = re.compile("genemarks_calls='(.*)'")
p_glimmerCalls                = re.compile("glimmer_calls='(.*)'")
p_prodigalCalls               = re.compile("prodigal_calls='(.*)'")
p_phanotateCalls              = re.compile("phanotate_calls='(.*)'")
# Custom gene calling parameters
p_custom_geneCalls            = re.compile("custom_gene_calls='(.*)'") # true/false 
p_custom_geneCallerName       = re.compile("custom_gene_caller_name='(.*)'")    # the name of the caller (e.g., RAST)
#p_custom_geneCallerOutfile    = re.compile("custom_gene_caller_outfile='(.*)'") 

# BLAST PROCESSING
p_blastpSearch                = re.compile("blastp='(.*)'")             
# Blast Parameters (value)
p_blastpIdentity              = re.compile("blastp_identity='(\d+)'")  
p_blastnIdentity              = re.compile("blastn_identity='(\d+)'")   
p_blastpHitCount              = re.compile("blastp_hit_count='(\d+)'")
p_blastnHitCount              = re.compile("blastn_hit_count='(\d+)'")
# Blast Processes to run (true/false)
p_ncbiVirusGenomeBlast        = re.compile("ncbi_virus_genome_blast='(.*)'")
p_ncbiVirusProteinBlast       = re.compile("ncbi_virus_protein_blast='(.*)'")
p_refseqProteinBlast          = re.compile("refseq_protein_blast='(.*)'")
p_refseqGeneBlast             = re.compile("refseq_gene_blast='(.*)'")
p_pvogsBlast                  = re.compile("pvogs_blast='(.*)'")
p_vogsBlast                   = re.compile("vogs_blast='(.*)'")  #*** To be deprecated
p_vogGeneBlast                = re.compile("vog_gene_blast='(.*)'")
p_vogProteinBlast             = re.compile("vog_protein_blast='(.*)'")
p_phantomeBlast               = re.compile("phantome_blast='(.*)'")
p_phageEnzymeBlast            = re.compile("phage_enzyme_blast='(.*)'")
p_keggVirusBlast              = re.compile("kegg_virus_blast='(.*)'")
p_pfamBlast                   = re.compile("pfam_blast='(.*)'")
p_smartBlast                  = re.compile("smart_blast='(.*)'")
p_swissprotBlast              = re.compile("swissprot_blast='(.*)'")
p_uniprotBlast                = re.compile("uniprot_blast='(.*)'")    # not in service
p_nrBlast                     = re.compile("nr_blast='(.*)'")
p_cazyBlast                   = re.compile("cazy_blast='(.*)'")
# Blast Database Locations; these may be read in from config file, but are set as environment vars
p_ncbiVirusGenomeDBpath       = re.compile("ncbi_virus_genome_database_path='(.*)'")
p_ncbiVirusProteinDBpath      = re.compile("ncbi_virus_protein_database_path='(.*)'")
p_refseqProteinDBpath         = re.compile("refseq_protein_database_path='(.*)'")
p_refseqGeneDBpath            = re.compile("refseq_gene_database_path='(.*)'")
p_pvogsDBpath                 = re.compile("pvogs_database_path='(.*)'")
p_vogsDBpath                  = re.compile("vogs_database_path='(.*)'")
p_vogGeneDBpath               = re.compile("vog_gene_database_path='(.*)'")
p_vogProteinDBpath            = re.compile("vog_protein_database_path='(.*)'")
p_phantomeDBpath              = re.compile("phantome_database_path='(.*)'")
p_phageEnzymeDBpath           = re.compile("phage_enzyme_database_path='(.*)'") # not yet in service
p_keggVirusDBpath             = re.compile("kegg_virus_database_path='(.*)'")
p_pfamDBpath                  = re.compile("pfam_database_path='(.*)'")
p_smartDBpath                 = re.compile("smart_database_path='(.*)'")
p_swissprotDBpath             = re.compile("swissprot_database_path='(.*)'")
p_uniprotDBpath               = re.compile("uniprot_database_path='(.*)'")     # not in service
p_nrDBpath                    = re.compile("nr_database_path='(.*)'")
p_cazyDBpath                  = re.compile("cazy_database_path='(.*)'")
p_cazyAnnotationPath          = re.compile("cazy_annotation_path='(.*)'")

# Custom blast processes and parameters
# Custom blast processes to run (true/false)
p_customGenomeBlast           = re.compile("custom_genome_blast='(.*)'")                 # not yet in service
p_customGeneBlast             = re.compile("custom_gene_blast='(.*)'")                   # not yet in service
p_customProteinBlast          = re.compile("custom_protein_blast='(.*)'")                # not yet in service
# Custom blast database names
p_customGenomeDBname          = re.compile("custom_genome_blast_database_name='(.*)'")   # not yet in service
p_customGeneDBname            = re.compile("custom_gene_blast_database_name='(.*)'")     # not yet in service
p_customProteinDBname         = re.compile("custom_protein_blast_database_name='(.*)'")  # not yet in service
# Paths to custom databases
p_customGenomeDBpath          = re.compile("custom_genome_blast_database_path='(.*)'")   # not yet in service
p_customGeneDBpath            = re.compile("custom_gene_blast_database_path='(.*)'")     # not yet in service
p_customProteinDBpath         = re.compile("custom_protein_blast_database_path='(.*)'")  # not yet in service

# HMM PROCESSING

# HMM Processing to be done (true/false)
p_phmmerSearch                = re.compile("phmmer='(.*)'")           # for hmm search of fasta database(s)
p_jackhmmerSearch             = re.compile("jackhmmer='(.*)'")        # for hmm search of fasta database(s)
p_hmmscan                     = re.compile("hmmscan='(.*)'")          # for hmm search of hmm profile database(s) 
# HMM Databases to be used (true/false)
p_ncbiVirusGenomeHmm          = re.compile("ncbi_virus_genome_hmm_profiles='(.*)'")  # not in service
p_ncbiVirusProteinHmm         = re.compile("ncbi_virus_protein_hmm_profiles='(.*)'") # not in service
p_refseqProteinHmm            = re.compile("refseq_protein_hmm_profiles='(.*)'")     # not in service
p_refseqGeneHmm               = re.compile("refseq_gene_hmm_profiles='(.*)'")        # not in service
p_pvogsHmm                    = re.compile("pvogs_hmm_profiles='(.*)'")
p_vogsHmm                     = re.compile("vogs_hmm_profiles='(.*)'")
p_phantomeHmm                 = re.compile("phantome_hmm_profiles='(.*)'")           # not in service
p_phageEnzymeHmm              = re.compile("phage_enzyme_hmm_profiles='(.*)'")       # not in service
p_keggVirusHmm                = re.compile("kegg_virus_hmm_profiles='(.*)'")         # not in service
p_pfamHmm                     = re.compile("pfam_hmm_profiles='(.*)'")               # not in service
p_smartHmm                    = re.compile("smart_hmm_profiles='(.*)'")              # not in service
p_swissprotHmm                = re.compile("swissprot_hmm_profiles='(.*)'")          # not in service
p_uniprotHmm                  = re.compile("uniprot_hmm_profiles='(.*)'")            # not in service (will this be combined with swissprot?)
p_nrHmm                       = re.compile("nr_hmm_profiles='(.*)'")
# HMM profile databases (locations; string)
p_ncbiVirusGenomeHmmDBpath    = re.compile("ncbi_virus_genome_hmm_profiles_database_path='(.*)'")
p_ncbiVirusProteinHmmDBpath   = re.compile("ncbi_virus_protein_hmm_profiles_database_path='(.*)'")
p_refseqProteinHmmDBpath      = re.compile("refseq_protein_hmm_profiles_database_path='(.*)'")
p_refseqGeneHmmDBpath         = re.compile("refseq_gene_hmm_profiles_database_path='(.*)'")
p_pvogsHmmDBpath              = re.compile("pvogs_hmm_profiles_database_path='(.*)'")
p_vogsHmmDBpath               = re.compile("vogs_hmm_profiles_database_path='(.*)'")
p_phantomeHmmDBpath           = re.compile("phantome_hmm_profiles_database_path='(.*)'")
p_phageEnzymeHmmDBpath        = re.compile("phage_enzyme_hmm_profiles_database_path='(.*)'")
p_keggVirusHmmDBpath          = re.compile("kegg_virus_hmm_profiles_database_path='(.*)'")
p_pfamHmmDBpath               = re.compile("pfam_hmm_profiles_database_path='(.*)'")
p_smartHmmDBpath              = re.compile("smart_hmm_profiles_database_path='(.*)'")
p_swissprotHmmDBpath          = re.compile("swissprot_hmm_profiles_database_path='(.*)'")
p_uniprotHmmDBpath            = re.compile("uniprot_hmm_profiles_database_path='(.*)'")
p_nrHmmDBpath                 = re.compile("nr_hmm_profiles_database_path='(.*)'")
# HMM custom processing ### not yet in service
p_customHmm                   = re.compile("custom_hmm_profiles='(.*)'")    # not yet in service
p_customHmmDBname             = re.compile("custom_hmm_profiles_database_name='(.*)'")  # value/string
p_customHmmDBpath             = re.compile("custom_hmm_profiles_database_path='(.*)'")  # value/string

# COMPARATIVE GENOMICS
p_runCGP                      = re.compile("CGP='(.*)'")              # run CGP followed by Genomics module
p_hmmbuild                    = re.compile("hmmbuild='(.*)'")         # not yet in service; hmmbuild creates profiles from homology groups
p_hmmsearch                   = re.compile("hmmsearch='(.*)'")        # not yet in service; hmmscan is the program, hmmsearch uses hmmscan w/profile query 

# DEPENDENT CODE LOCATIONS: These codes should be installed globally, except for phanotate.
p_blastPlusHome               = re.compile("blast_plus_home='(.*)'")
p_embossHome                  = re.compile("emboss_home='(.*)'")
p_tRNAscanSEhome              = re.compile("tRNAscanSE_home='(.*)'")
p_glimmerHome                 = re.compile("glimmer_home='(.*)'")
p_prodigalHome                = re.compile("prodigal_home='(.*)'")
p_phanotateHome               = re.compile("phanotate_home='(.*)'")
p_genemarkHome                = re.compile("genemark_home='(.*)'")
p_hmmerHome                   = re.compile("hmmer_home='(.*)'")

# PARALLELISM
p_threads                     = re.compile("threads='(.*)'")
p_hpc                         = re.compile("HPC='(.*)'")
p_blastThreads                = re.compile("blast_threads='(.*)'")

# VERBOSTIY
p_phateWarnings               = re.compile("phate_warnings='(.*)'")
p_phateMessages               = re.compile("phate_messages='(.*)'")
p_phateProgress               = re.compile("phate_progress='(.*)'")
p_cgcWarnings                 = re.compile("cgc_warnings='(.*)'")    # CGC warnings, messages, and progress are now linked w/ those of phate
p_cgcMessages                 = re.compile("cgc_messages='(.*)'")
p_cgcProgress                 = re.compile("cgc_progress='(.*)'")
p_cgpWarnings                 = re.compile("cgp_warnings='(.*)'")
p_cgpMessages                 = re.compile("cgp_messages='(.*)'")
p_cgpProgress                 = re.compile("cgp_progress='(.*)'")
p_cleanRawData                = re.compile("clean_raw_data='(.*)'")

##### GET INPUT PARAMETERS #####

# Open log file
logfile = BASE_DIR_DEFAULT + CODE_BASE + ".log"
if not HPC:
    LOG = open(logfile,'w')
    LOG.write("%s%s\n" % ("Begin log file ",datetime.datetime.now()))

if len(sys.argv) != 2:
    print(USAGE_STRING)
    dateTime = os.popen('date')
    if not HPC:
        LOG.write("%s%s%s%s\n" % ("Incorrect number of input parameters: ", len(sys.argv), ". End log ",dateTime))
        LOG.close()
    exit(0)
else:
    match_config = re.search(p_config,sys.argv[1])
    if match_config:
        configFile = sys.argv[1]
        if not HPC:
            LOG.write("%s%s\n" % ("Config file is ",configFile))
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
        if not HPC:
            LOG.write("%s%s\n" % ("A help string was provided to user; End log ",datetime.datetime.now()))
            LOG.close()
        exit(0)

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
    if not HPC:
        LOG.write("%s%s\n" % ("A help string was provided to user; End log ",datetime.datetime.now()))
        LOG.close()
    exit(0)

##### Read input parameters from configuration file; capture user's parameter choices

FIRST_GENOME = True
DATA_ITEMS_NUM = 6
genomeDataDict = {
    "genomeNumber"  : "",
    "genomeFile"    : "",
    "genomeType"    : "",
    "genomeSpecies" : "",
    "genomeName"    : "",
    "outputSubdir"  : "",
    }
genomeList      = []  # List of genomeData objects
genomeNumber    = ""  # Number of current genome; assigned by user; could be a string
genomeDataItems = 0   # Number of data items collected for current genome; should be DATA_ITEMS_NUM
BEGIN_GENOME_LIST = False
nextGenomeData = genomeDataDict

cLines = CONFIG.read().splitlines()
for cLine in cLines:
    match_comment                   = re.search(p_comment,cLine)
    match_blank                     = re.search(p_blank,cLine)

    # genome info
    match_genomeList                = re.search(p_genomeList,cLine)
    match_genomeNumber              = re.search(p_genomeNumber,cLine)
    match_genomeFile                = re.search(p_genomeFile,cLine)
    match_genomeType                = re.search(p_genomeType,cLine)
    match_genomeSpecies             = re.search(p_genomeSpecies,cLine)
    match_genomeName                = re.search(p_genomeName,cLine)
    match_outputSubdir              = re.search(p_outputSubdir,cLine)
    match_end                       = re.search(p_end,cLine)

    # gene calling
    match_genemarksCalls            = re.search(p_genemarksCalls,cLine)
    match_prodigalCalls             = re.search(p_prodigalCalls,cLine)
    match_glimmerCalls              = re.search(p_glimmerCalls,cLine)
    match_phanotateCalls            = re.search(p_phanotateCalls,cLine)
    match_customGeneCalls           = re.search(p_custom_geneCalls,cLine)
    match_customGeneCallerName      = re.search(p_custom_geneCallerName,cLine)   # Name of program (eg, RAST)
    #match_customGeneCallerOutfile   = re.search(p_custom_geneCallerOutfile,cLine)

    # translation info
    match_primaryCalls              = re.search(p_primaryCalls,cLine)
    match_geneticCode               = re.search(p_geneticCode,cLine)
    match_translateOnly             = re.search(p_translateOnly,cLine)

    # codes to run against sequence databases
    match_blastpSearch              = re.search(p_blastpSearch,cLine)    # for blast search
    match_phmmerSearch              = re.search(p_phmmerSearch,cLine)    # for hmm search
    match_jackhmmerSearch           = re.search(p_jackhmmerSearch,cLine) # for iterative hmm search

    # blast parameters
    match_blastpIdentity            = re.search(p_blastpIdentity,cLine)
    match_blastnIdentity            = re.search(p_blastnIdentity,cLine)
    match_blastpHitCount            = re.search(p_blastpHitCount,cLine)
    match_blastnHitCount            = re.search(p_blastnHitCount,cLine)

    # sequence/blast databases to be used
    match_ncbiVirusGenomeBlast      = re.search(p_ncbiVirusGenomeBlast,cLine)
    match_ncbiVirusProteinBlast     = re.search(p_ncbiVirusProteinBlast,cLine)
    match_refseqProteinBlast        = re.search(p_refseqProteinBlast,cLine)
    match_refseqGeneBlast           = re.search(p_refseqGeneBlast,cLine)
    match_pvogsBlast                = re.search(p_pvogsBlast,cLine)
    match_vogsBlast                 = re.search(p_vogsBlast,cLine)     #*** To be deprecated
    match_vogGeneBlast              = re.search(p_vogGeneBlast,cLine)
    match_vogProteinBlast           = re.search(p_vogProteinBlast,cLine)
    match_phantomeBlast             = re.search(p_phantomeBlast,cLine)
    match_phageEnzymeBlast          = re.search(p_phageEnzymeBlast,cLine)
    match_keggVirusBlast            = re.search(p_keggVirusBlast,cLine)
    match_pfamBlast                 = re.search(p_pfamBlast,cLine)
    match_smartBlast                = re.search(p_smartBlast,cLine)
    match_swissprotBlast            = re.search(p_swissprotBlast,cLine)
    match_uniprotBlast              = re.search(p_uniprotBlast,cLine)
    match_nrBlast                   = re.search(p_nrBlast,cLine)
    match_cazyBlast                 = re.search(p_cazyBlast,cLine)

    match_customGenomeBlast         = re.search(p_customGenomeBlast,cLine)
    match_customGenomeDBname        = re.search(p_customGenomeDBname,cLine)
    match_customGeneBlast           = re.search(p_customGeneBlast,cLine)
    match_customGeneDBname          = re.search(p_customGeneDBname,cLine)
    match_customProteinBlast        = re.search(p_customProteinBlast,cLine)
    match_customProteinDBname       = re.search(p_customProteinDBname,cLine)

    # hmm codes to run against profile databases
    #match_hmmProgram                = re.search(p_hmmProgram,cLine)      # for blast search
    match_hmmscan                   = re.search(p_hmmscan,cLine)         # for hmm search

    # True/False running a custom hmm profile database
    match_customHmm                 = re.search(p_customHmm,cLine)       # for hmm search

    # hmm profiles
    #*** some of the hmms listed may not be implemented, but are included just in case
    match_ncbiVirusGenomeHmm        = re.search(p_ncbiVirusGenomeHmm,cLine)    # ? big; is this used?
    match_ncbiVirusProteinHmm       = re.search(p_ncbiVirusProteinHmm,cLine)   # ? big; is this used?
    match_refseqProteinHmm          = re.search(p_refseqProteinHmm,cLine)      # ? big; is this used?
    match_refseqGeneHmm             = re.search(p_refseqGeneHmm,cLine)         # ? big; is this used?
    match_pvogsHmm                  = re.search(p_pvogsHmm,cLine)
    match_vogsHmm                   = re.search(p_vogsHmm,cLine)
    match_phantomeHmm               = re.search(p_phantomeHmm,cLine)
    match_phageEnzymeHmm            = re.search(p_phageEnzymeHmm,cLine)
    match_keggVirusHmm              = re.search(p_keggVirusHmm,cLine)          # ? is there an hmm profiles db for kegg? build it?
    match_pfamHmm                   = re.search(p_pfamHmm,cLine)
    match_smartHmm                  = re.search(p_smartHmm,cLine)
    match_swissprotHmm              = re.search(p_swissprotHmm,cLine)
    match_uniprotHmm                = re.search(p_uniprotHmm,cLine)            # will this be combined w/Swissprot?
    match_nrHmm                     = re.search(p_nrHmm,cLine)                 # ? big; is this used?
    match_customHmmDBname           = re.search(p_customHmmDBname,cLine)

    # gene calls
    match_genemarksCalls            = re.search(p_genemarksCalls,cLine)
    match_prodigalCalls             = re.search(p_prodigalCalls,cLine)
    match_glimmerCalls              = re.search(p_glimmerCalls,cLine)
    match_phanotateCalls            = re.search(p_phanotateCalls,cLine)
    match_customGeneCalls           = re.search(p_custom_geneCalls,cLine)
    match_customGeneCallerName      = re.search(p_custom_geneCallerName,cLine)
    #match_customGeneCallerOutfile   = re.search(p_custom_geneCallerOutfile,cLine)

    # directories; should be unnecessary to read from config; using standard predefined directory locations
    match_phateDir                  = re.search(p_phateDir,cLine)
    match_databaseDir               = re.search(p_databaseDir,cLine)
    match_softwareDir               = re.search(p_softwareDir,cLine)

    # 3rd party codes
    #*** The below parameters are no longer in use; assuming global installs of 3rd party software
    # phanotate should be located to the external softwares directory
    match_blastPlusHome             = re.search(p_blastPlusHome,cLine)
    match_embossHome                = re.search(p_embossHome,cLine)
    match_tRNAscanSEhome            = re.search(p_tRNAscanSEhome,cLine)
    match_glimmerHome               = re.search(p_glimmerHome,cLine)
    match_prodigalHome              = re.search(p_prodigalHome,cLine)
    match_phanotateHome             = re.search(p_phanotateHome,cLine)
    match_genemarkHome              = re.search(p_genemarkHome,cLine)
    match_hmmerHome                 = re.search(p_hmmerHome,cLine)

    # paths to databases
    # blast
    match_ncbiVirusGenomeDBpath     = re.search(p_ncbiVirusGenomeDBpath,cLine)
    match_ncbiVirusProteinDBpath    = re.search(p_ncbiVirusProteinDBpath,cLine)
    match_refseqGeneDBpath          = re.search(p_refseqGeneDBpath,cLine)
    match_refseqProteinDBpath       = re.search(p_refseqProteinDBpath,cLine)
    match_pvogsDBpath               = re.search(p_pvogsDBpath,cLine)
    match_vogsDBpath                = re.search(p_vogsDBpath,cLine)
    match_vogGeneDBpath             = re.search(p_vogGeneDBpath,cLine)
    match_vogProteinDBpath          = re.search(p_vogProteinDBpath,cLine)
    match_phantomeDBpath            = re.search(p_phantomeDBpath,cLine)
    match_phageEnzymeDBpath         = re.search(p_phageEnzymeDBpath,cLine)
    match_keggVirusDBpath           = re.search(p_keggVirusDBpath,cLine)
    match_pfamDBpath                = re.search(p_pfamDBpath,cLine)
    match_smartDBpath               = re.search(p_smartDBpath,cLine)
    match_swissprotDBpath           = re.search(p_swissprotDBpath,cLine)
    match_uniprotDBpath             = re.search(p_uniprotDBpath,cLine)
    match_nrDBpath                  = re.search(p_nrDBpath,cLine)
    match_cazyDBpath                = re.search(p_cazyDBpath,cLine)
    match_cazyAnnotationPath        = re.search(p_cazyAnnotationPath,cLine)
    match_customGenomeDBpath        = re.search(p_customGenomeDBpath,cLine)
    match_customGeneDBpath          = re.search(p_customGeneDBpath,cLine)
    match_customProteinDBpath       = re.search(p_customProteinDBpath,cLine)
    # hmm
    match_ncbiVirusGenomeHmmDBpath  = re.search(p_ncbiVirusGenomeHmmDBpath,cLine)
    match_ncbiVirusProteinHmmDBpath = re.search(p_ncbiVirusProteinHmmDBpath,cLine)
    match_refseqGeneHmmDBpath       = re.search(p_refseqGeneHmmDBpath,cLine)
    match_refseqProteinHmmDBpath    = re.search(p_refseqProteinHmmDBpath,cLine)
    match_pvogsHmmDBpath            = re.search(p_pvogsHmmDBpath,cLine)
    match_vogsHmmDBpath             = re.search(p_vogsHmmDBpath,cLine)
    match_phantomeHmmDBpath         = re.search(p_phantomeHmmDBpath,cLine)
    match_phageEnzymeHmmDBpath      = re.search(p_phageEnzymeHmmDBpath,cLine)
    match_keggVirusHmmDBpath        = re.search(p_keggVirusHmmDBpath,cLine)
    match_pfamHmmDBpath             = re.search(p_pfamHmmDBpath,cLine)
    match_smartHmmDBpath            = re.search(p_smartHmmDBpath,cLine)
    match_swissprotHmmDBpath        = re.search(p_swissprotHmmDBpath,cLine)
    match_uniprotHmmDBpath          = re.search(p_uniprotHmmDBpath,cLine)
    match_nrHmmDBpath               = re.search(p_nrHmmDBpath,cLine)
    match_customHmmDBpath           = re.search(p_customHmmDBpath,cLine)

    # Comparative genomics
    match_runCGP                    = re.search(p_runCGP,cLine)
    match_hmmbuild                  = re.search(p_hmmbuild,cLine)
    match_hmmsearch                 = re.search(p_hmmsearch,cLine)

    # parallelism
    match_threads                   = re.search(p_threads,cLine)
    match_hpc                       = re.search(p_hpc,cLine)
    match_blastThreads              = re.search(p_blastThreads,cLine)

    # verbosity
    match_phateWarnings             = re.search(p_phateWarnings,cLine)
    match_phateMessages             = re.search(p_phateMessages,cLine)
    match_phateProgress             = re.search(p_phateProgress,cLine)
    match_cgcWarnings               = re.search(p_cgcWarnings,cLine)
    match_cgcMessages               = re.search(p_cgcMessages,cLine)
    match_cgcProgress               = re.search(p_cgcProgress,cLine)
    match_cgpWarnings               = re.search(p_cgpWarnings,cLine)
    match_cgpMessages               = re.search(p_cgpMessages,cLine)
    match_cgpProgress               = re.search(p_cgpProgress,cLine)
    match_cleanRawData              = re.search(p_cleanRawData,cLine)

    ##### Capture list of genomes and associated data #####

    if (match_comment or match_blank):
        pass

    elif match_genomeList:  # Capture all genomes listed; for each, gather genome file, genome type, species, name, output subdir
        pass  # line in config not needed; present for user's benefit only

    elif match_genomeNumber:  # The next genome's data
        # First, record previous genome data, if this is not the first genome
        if not FIRST_GENOME:
            if not HPC:
                LOG.write("%s\n" % ("Appending a genome data set"))
            if genomeDataItems != DATA_ITEMS_NUM:  # If record appears incomplete, flag a problem
                if not HPC:
                    LOG.write("%s%s\n" % ("WARNING: check config file for possible incorrect data items: ", genomeDataItems))
            genomeList.append(nextGenomeData)

        # Next, begin collecting next genome's data
        if not HPC:
            LOG.write("%s%s\n" % ("Creating a new genome data set for ", match_genomeNumber.group(0)))
        genomeDataItems = 0
        nextGenomeData = copy.deepcopy(genomeDataDict)  # make new genome data object
        genomeNumber = match_genomeNumber.group(1)
        nextGenomeData["genomeNumber"] = genomeNumber
        genomeDataItems += 1
        FIRST_GENOME = False

    elif match_genomeFile:
        value = match_genomeFile.group(1)
        if value != '':
            GENOME_FILE = value
        else:
            GENOME_FILE = "unknown"
            if PHATE_WARNINGS:
                print("multiPhate says, WARNING:  GENOME_FILE is", GENOME_FILE)
        nextGenomeData["genomeFile"] = GENOME_FILE
        if not HPC:
            LOG.write("%s%s\n" % ("GENOME_FILE is ",GENOME_FILE))
        genomeDataItems += 1

    elif match_genomeType:
        value = match_genomeType.group(1)
        if value.lower() == 'phage' or value.lower() == 'bacteriophage':
            nextGenomeData["genomeType"] = 'phage'
        elif value.lower() == 'virus' or value.lower() == 'viral' or value.lower() == 'viridae':
            nextGenomeData["genomeType"] = 'virus'
        elif value.lower() == 'bacteria' or value.lower() == 'bacterium' or value.lower() == 'bacterial':
            nextGenomeData["genomeType"] = 'bacterium'
        else:
            nextGenomeData["genomeType"] = 'other'
        if not HPC:
            LOG.write("%s%s\n" % ("genome type is ",nextGenomeData["genomeType"]))
        genomeDataItems += 1

    elif match_genomeSpecies:
        genomeSpecies = match_genomeSpecies.group(1)
        nextGenomeData["genomeSpecies"] = genomeSpecies
        if not HPC:
            LOG.write("%s%s\n" % ("genomeSpecies is ",genomeSpecies))
        genomeDataItems += 1

    elif match_genomeName:
        genomeName = match_genomeName.group(1)
        nextGenomeData["genomeName"] = genomeName
        if not HPC:
            LOG.write("%s%s\n" % ("genome name is ",genomeName))
        genomeDataItems += 1

    elif match_phateDir:
        if match_phateDir.group(1) != '':
            os.environ["PHATE_BASE_DIR"] = match_phateDir.group(1)

    elif match_databaseDir:
        if match_databaseDir.group(1) != '':
            os.environ["PHATE_DATABASE_DIR"] = match_databaseDir.group(1)

    elif match_softwareDir:
        if match_softwareDir.group(1) != '':
            os.environ["PHATE_SOFTWARE_DIR"] = match_softwareDir.group(1)

    elif match_outputSubdir: #*** Note that if the output dir is not read before subdir; depends on user not changing order in config - Clean this up!
        value = match_outputSubdir.group(1)
        if value != '':
            value = value.rstrip('/')  # be sure that name of subdir ends in exactly one '/' (user might omit the slash)
            value = value + '/'
            nextGenomeData["outputSubdir"] = value
        else:
            nextGenomeData["outputSubdir"] = "unknown"
            if PHATE_WARNINGS:
                print("multiPhate says, WARNING: pipeline output subdir is ", "unknown")
        genomeDataItems += 1

    elif match_end:  # List of genomes complete; record last genome's data
        if genomeDataItems != DATA_ITEMS_NUM:
            if not HPC:
                LOG.write("%s%s%s%s%s%s\n" % ("WARNING: check config file for possible incorrect data items: ", genomeDataItems, " for genome ",nextGenomeData["genomeName"],' ',nextGenomeData["genomeNumber"]))
        genomeList.append(nextGenomeData)
        if not HPC:
            LOG.write("%s%s\n" % ("END: Length of genomeList is ",len(genomeList)))

    ##### Other processing #####

    elif match_geneticCode:
        value = match_geneticCode.group(1)
        if value != '':
            geneticCode = value

    elif match_translateOnly:
        value = match_translateOnly.group(1)
        if value.lower() == 'yes' or value.lower() == 'true' or value.lower() == 'on':
            translateOnly = True
        elif value.lower() == 'no' or value.lower() == 'false' or value.lower() == 'off' or value == '':
            translateOnly = False
        else:
            if PHATE_WARNINGS == 'True':
                print("multiPhate says, WARNING: Invalid string following translate_only parameter in config file:", value)
            if not HPC:
                LOGFILE.write("%s%s\n" % ("Invalid string following translate_only parameter in config file: ", value))

    ##### Gene Calls #####

    elif match_primaryCalls:
        value = match_primaryCalls.group(1)
        if value.lower() == 'phanotate':
            primaryCalls = 'phanotate'
        elif value.lower() == 'consensus':
            primaryCalls = 'consensus'
        elif value.lower() == 'genemarks' or value.lower() == 'genemark':
            primaryCalls = 'genemarks'
        elif value.lower() == 'glimmer2':
            primaryCalls = 'glimmer2'
        elif value.lower() == 'glimmer3' or value.lower() == 'glimmer':
            primaryCalls = 'glimmer3'
        elif value.lower() == 'prodigal':
            primaryCalls = 'prodigal'
        elif value.lower() == 'rast':
            primaryCalls = 'rast'
        elif value.lower() == 'consensus':
            primaryCalls = 'consensus'
        elif value.lower() == 'superset':
            primaryCalls = 'superset'
        elif value.lower() == 'commoncore' or value.lower() == 'common-core' or value.lower() == 'common_core':
            primaryCalls = 'commoncore'
        elif value.lower() == 'custom':
            primaryCalls = 'custom'
        else:
            if PHATE_WARNINGS == 'True':
                print("multiPhate says, WARNING: Invalid string for primary calls in config file:", value)
        if primaryCalls == 'glimmer2' or primaryCalls == 'glimmer3':
            primaryCallsFile = 'glimmer' + '.cgc'
        else:
            primaryCallsFile = primaryCalls + '.cgc'

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

    elif match_customGeneCalls:
        value = match_customGeneCalls.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            customGeneCalls = True
        else:
            customGeneCalls = False

    elif match_customGeneCallerName:   # ie, the name of the program that generated the calls (e.g., RAST)
        value = match_customGeneCallerName.group(1)
        customGeneCallerName = value
        os.environ["PHATE_CUSTOM_GENECALLER_NAME"] = customGeneCallerName

    # To be deprecated: user must use standard file naming for custom calls (ie, genomeName.custom)
    #elif match_customGeneCallerOutfile:
    #    value = match_customGeneCallerOutfile.group(1)
    #    customGeneCallerOutfile = value

    ##### BLAST #####

    elif match_blastpSearch:     # blastp search against fasta blast database(s)
        value = match_blastpSearch.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            blastpSearch = True

    elif match_phmmerSearch:     # phmmer search against fasta blast database(s)
        value = match_phmmerSearch.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            phmmerSearch = True

    elif match_jackhmmerSearch:     # jackhmmer search against fasta blast database(s)
        value = match_jackhmmerSearch.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            jackhmmerSearch = True

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

    elif match_ncbiVirusGenomeBlast:
        value = match_ncbiVirusGenomeBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            ncbiVirusGenomeBlast = True

    elif match_ncbiVirusProteinBlast:
        value = match_ncbiVirusProteinBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            ncbiVirusProteinBlast = True

    elif match_keggVirusBlast:
        value = match_keggVirusBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
             keggVirusBlast = True

    elif match_refseqProteinBlast:
        value = match_refseqProteinBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            refseqProteinBlast = True

    elif match_refseqGeneBlast:
        value = match_refseqGeneBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            refseqGeneBlast = True

    elif match_pvogsBlast:
        value = match_pvogsBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            pvogsBlast = True

    elif match_vogsBlast:   #*** To be deprecated
        value = match_vogsBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            vogsBlast = True

    elif match_vogGeneBlast:
        value = match_vogGeneBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            vogGeneBlast = True

    elif match_vogProteinBlast:
        value = match_vogProteinBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            vogProteinBlast = True

    elif match_phantomeBlast:
        value = match_phantomeBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            phantomeBlast = True

    elif match_phageEnzymeBlast:
        value = match_phageEnzymeBlast.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            phageEnzymeBlast = True

    elif match_pfamBlast:
        value = match_pfamBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            pfamBlast = True

    elif match_smartBlast:
        value = match_smartBlast.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            smartBlast = True

    elif match_swissprotBlast:
        value = match_swissprotBlast.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            swissprotBlast = True

    elif match_uniprotBlast:
        value = match_uniprotBlast.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            uniprotBlast = True

    elif match_nrBlast:
        value = match_nrBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
             nrBlast = True

    elif match_cazyBlast:
        value = match_cazyBlast.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
             cazyBlast = True

    elif match_customGenomeBlast:
        value = match_customGenomeBlast.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            customGenomeBlast = True

    elif match_customGenomeDBname:
        value = match_customGenomeDBname.group(1)
        if value != '':
            customGenomeDBname = value
    
    #elif match_customGenomeDBpath:
    #    value = match_customGenomeDBpath.group(1)
    #    if value != '':
    #        customGenomeDBpath = value
    
    elif match_customGeneBlast:
        value = match_customGeneBlast.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            customGeneBlast = True

    elif match_customGeneDBname:
        value = match_customGeneDBname.group(1)
        if value != '':
            customGeneDBname = value
    
    #elif match_customGeneDBpath:
    #    value = match_customGeneDBpath.group(1)
    #    if value != '':
    #        customGeneDBpath = value
    
    elif match_customProteinBlast:
        value = match_customProteinBlast.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            customProteinBlast = True

    elif match_customProteinDBname:
        value = match_customProteinDBname.group(1)
        if value != '':
            customProteinDBname = value
    
    #elif match_customProteinDBpath:
    #    value = match_customProteinDBpath.group(1)
    #    if value != '':
    #        customProteinDBpath = value
    
    ##### HMM #####

    # HMM programs

    elif match_hmmscan:          # hmmscan search against hmm profile database(s)
        value = match_hmmscan.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            hmmscan = True

    # HMM databases

    elif match_ncbiVirusGenomeHmm:
        value = match_ncbiVirusGenomeHmm.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            ncbiVirusGenomeHmm = True

    elif match_ncbiVirusProteinHmm:
        value = match_ncbiVirusProteinHmm.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            ncbiVirusProteinHmm = True

    elif match_refseqGeneHmm:
        value = match_refseqGeneHmm.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            refseqGeneHmm = True

    elif match_refseqProteinHmm:
        value = match_refseqProteinHmm.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            refseqProteinHmm = True

    elif match_pvogsHmm:
        value = match_pvogsHmm.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            pvogsHmm = True

    elif match_vogsHmm:
        value = match_vogsHmm.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            vogsHmm = True

    elif match_phantomeHmm:
        value = match_phantomeHmm.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            phantomeHmm = True

    elif match_phageEnzymeHmm:
        value = match_phageEnzymeHmm.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            phageEnzymeHmm = True

    elif match_keggVirusHmm:
        value = match_keggVirusHmm.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
             keggVirusHmm = True

    elif match_pfamHmm:
        value = match_pfamHmm.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            pfamHmm = True

    elif match_smartHmm:
        value = match_smartHmm.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            smartHmm = True

    elif match_swissprotHmm:
        value = match_swissprotHmm.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            swissprotHmm = True

    elif match_uniprotHmm:
        value = match_uniprotHmm.group(1).lower()
        if value == 'true' or value == 'yes' or value == 'on':
            uniprotHmm = True

    elif match_nrHmm:
        value = match_nrHmm.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
             nrHmm = True

    elif match_customHmm:
        value = match_customHmm.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            customHmm = True

    elif match_customHmmDBname:
        value = match_customHmmDBname.group(1)
        if value != '':
            customHmmDBname = value
    
    #elif match_customHmmDBpath:
    #    value = match_customHmmDBpath.group(1)
    #    if value != '':
    #        customHmmDBpath = value
    
    # Comparative Genomics

    elif match_runCGP:
        value = match_runCGP.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            runCGP      = True
            runGenomics = True

    elif match_hmmbuild:         # create hmm profiles
        value = match_hmmbuild.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            hmmbuild = True

    elif match_hmmsearch:         # search hmm profiles against fasta database(s)
        value = match_hmmsearch.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            hmmsearch = True

    ##### DEPENDENT CODE LOCATIONS #####

    elif match_blastPlusHome:
        if match_blastPlusHome.group(1) != '':
            os.environ["PHATE_BLAST_HOME"] = match_blastPlusHome.group(1)

    elif match_embossHome:
        if match_embossHome.group(1) != '':
            os.environ["PHATE_EMBOSS_PHATE_HOME"] = match_embossHome.group(1)

    elif match_tRNAscanSEhome:
        if match_tRNAscanSEhome.group(1) != '':
            os.environ["PHATE_tRNAscanSE_HOME"] = match_tRNAscanSEhome.group(1)

    elif match_glimmerHome:
        if match_glimmerHome.group(1) != '':
            os.environ["PHATE_GLIMMER_PATH"] = match_glimmerHome.group(1)

    elif match_prodigalHome:
        if match_prodigalHome.group(1) != '':
            os.environ["PHATE_PRODIGAL_PATH"] = match_prodigalHome.group(1)

    elif match_phanotateHome:
        if match_phanotateHome.group(1) != '':
            os.environ["PHATE_PHANOTATE_PATH"] = match_phanotateHome.group(1)

    elif match_genemarkHome:
        if match_genemarkHome.group(1) != '':
            os.environ["PHATE_GENEMARKS_PATH"] = match_genemarkHome.group(1)

    ##### DATABASE LOCATIONS #####

    # Blast database locations

    elif match_ncbiVirusGenomeDBpath:
        if match_ncbiVirusGenomeDBpath.group(1) != '':
            os.environ["PHATE_NCBI_VIRUS_GENOME_BLAST_HOME"] = match_ncbiVirusGenomeDBpath.group(1)

    elif match_ncbiVirusProteinDBpath:
        if match_ncbiVirusProteinDBpath.group(1) != '':
            os.environ["PHATE_NCBI_VIRUS_PROTEIN_BLAST_HOME"] = match_ncbiVirusProteinDBpath.group(1)

    elif match_refseqGeneDBpath:
        if match_refseqGeneDBpath.group(1) != '':
            os.environ["PHATE_REFSEQ_GENE_BLAST_HOME"] = match_refseqGeneDBpath.group(1)

    elif match_refseqProteinDBpath:
        if match_refseqProteinDBpath.group(1) != '':
            os.environ["PHATE_REFSEQ_PROTEIN_BLAST_HOME"] = match_refseqProteinDBpath.group(1)

    elif match_pvogsDBpath:
        if match_pvogsDBpath.group(1) != '':
            os.environ["PHATE_PVOGS_BLAST_HOME"] = match_pvogsDBpath.group(1)

    elif match_vogsDBpath:   #*** To be deprecated
        if match_vogsDBpath.group(1) != '':
            os.environ["PHATE_VOGS_BLAST_HOME"] = match_vogsDBpath.group(1)

    elif match_vogGeneDBpath:   
        if match_vogGeneDBpath.group(1) != '':
            os.environ["PHATE_VOG_GENE_BLAST_HOME"] = match_vogGeneDBpath.group(1)

    elif match_vogProteinDBpath:   
        if match_vogProteinDBpath.group(1) != '':
            os.environ["PHATE_VOG_PROTEIN_BLAST_HOME"] = match_vogProteinDBpath.group(1)
            PHATE_VOG_PROTEIN_BASE_DIR = os.path.dirname(match_vogProteinDBpath.group(1)) + '/'
            os.environ["PHATE_VOG_PROTEIN_BASE_DIR"]        = PHATE_VOG_PROTEIN_BASE_DIR 
            os.environ["PHATE_VOG_PROTEIN_HEADER_FILE"]     = os.path.join(PHATE_VOG_PROTEIN_BASE_DIR,VOG_PROTEIN_HEADER_FILENAME)
            os.environ["PHATE_VOG_PROTEIN_ANNOTATION_FILE"] = os.path.join(PHATE_VOG_PROTEIN_BASE_DIR,VOG_PROTEIN_ANNOTATION_FILENAME)
            # Create Vog Protein headers file
            try:
                command = "grep '>' " + os.environ["PHATE_VOG_PROTEIN_BLAST_HOME"] + ' > ' + os.environ["PHATE_VOG_PROTEIN_HEADER_FILE"]
                success = os.system(command)
            except:
                print ("multiPhate says, ERROR: Could not create PHATE_VOG_PROTEIN_HEADER_FILE")

    elif match_phantomeDBpath:
        if match_phantomeDBpath.group(1) != '':
            os.environ["PHATE_PHANTOME_BLAST_HOME"] = match_phantomeDBpath.group(1)
            os.environ["PHATE_PHANTOME_BASE_DIR"] = os.path.dirname(match_phantomeDBpath.group(1)) + '/'

    elif match_phageEnzymeDBpath:
        if match_phageEnzymeDBpath.group(1) != '':
            os.environ["PHATE_PHAGE_ENZYME_BLAST_HOME"] = match_phageEnzymeDBpath.group(1)

    elif match_keggVirusDBpath:
        if match_keggVirusDBpath.group(1) != '':
            os.environ["PHATE_KEGG_VIRUS_BLAST_HOME"] = match_keggVirusDBpath.group(1)
            os.environ["PHATE_KEGG_VIRUS_BASE_DIR"] = os.path.dirname(match_keggVirusDBpath.group(1)) + '/'

    elif match_pfamDBpath:
        if match_pfamDBpath.group(1) != '':
            os.environ["PHATE_PFAM_BLAST_HOME"] = match_pfamDBpath.group(1)

    elif match_smartDBpath:
        if match_smartDBpath.group(1) != '':
            os.environ["PHATE_SMART_BLAST_HOME"] = match_smartDBpath.group(1)

    elif match_swissprotDBpath:
        if match_swissprotDBpath.group(1) != '':
            os.environ["PHATE_SWISSPROT_BLAST_HOME"] = match_swissprotDBpath.group(1)

    elif match_uniprotDBpath:
        if match_uniprotDBpath.group(1) != '':
            os.environ["PHATE_UNIPROT_BLAST_HOME"] = match_uniprotDBpath.group(1)

    elif match_nrDBpath:
        if match_nrDBpath.group(1) != '':
            os.environ["PHATE_NR_BLAST_HOME"] = match_nrDBpath.group(1)
     
    elif match_cazyDBpath:
        if match_cazyDBpath.group(1) != '':
            phateCAZyBlastHome = match_cazyDBpath.group(1)
            os.environ["PHATE_CAZY_BLAST_HOME"] = phateCAZyBlastHome 
            os.environ["PHATE_CAZY_BASE_DIR"]   = os.path.dirname(phateCAZyBlastHome) + '/'

    elif match_cazyAnnotationPath:
        if match_cazyAnnotationPath.group(1) != '':
            os.environ["PHATE_CAZY_ANNOTATION_PATH"] = match_cazyAnnotationPath.group(1)
     
    elif match_customGenomeDBpath:
        value = match_customGenomeDBpath.group(1)
        if value != '':
            customGenomeDBpath = value
        os.environ["PHATE_CUSTOM_GENOME_BLAST_HOME"] = customGenomeDBpath

    elif match_customGeneDBpath:
        value = match_customGeneDBpath.group(1)
        if value != '':
            customGeneDBpath = value
        os.environ["PHATE_CUSTOM_GENE_BLAST_HOME"] = customGeneDBpath

    elif match_customProteinDBpath:
        value = match_customProteinDBpath.group(1)
        if value != '':
            customProteinDBpath = value
        os.environ["PHATE_CUSTOM_PROTEIN_BLAST_HOME"] = customProteinDBpath
    
    # Hmm database locations

    elif match_ncbiVirusGenomeHmmDBpath:
        if match_ncbiVirusGenomeHmmDBpath.group(1) != '':
            os.environ["PHATE_NCBI_VIRUS_GENOME_HMM_HOME"] = match_ncbiVirusGenomeHmmDBpath.group(1)

    elif match_ncbiVirusProteinHmmDBpath:
        if match_ncbiVirusProteinHmmDBpath.group(1) != '':
            os.environ["PHATE_NCBI_VIRUS_PROTEIN_HMM_HOME"] = match_ncbiVirusProteinHmmDBpath.group(1)

    elif match_refseqGeneHmmDBpath:
        if match_refseqGeneHmmDBpath.group(1) != '':
            os.environ["PHATE_REFSEQ_GENE_HMM_HOME"] = match_refseqGeneHmmDBpath.group(1)

    elif match_refseqProteinHmmDBpath:
        if match_refseqProteinHmmDBpath.group(1) != '':
            os.environ["PHATE_REFSEQ_PROTEIN_HMM_HOME"] = match_refseqProteinHmmDBpath.group(1)

    elif match_pvogsHmmDBpath:
        value = match_pvogsHmmDBpath.group(1)
        if value != '':
            pvogsHmmDBpath = value
        os.environ["PHATE_PVOGS_HMM_HOME"] = pvogsHmmDBpath

    elif match_vogsHmmDBpath:
        value = match_vogsHmmDBpath.group(1)
        if value != '':
            vogsHmmDBpath = value
        os.environ["PHATE_VOGS_HMM_HOME"] = vogsHmmDBpath

    elif match_phantomeHmmDBpath:
        value = match_phantomeHmmDBpath.group(1)
        if value != '':
            phantomeHmmDBpath = value
        os.environ["PHATE_PHANTOME_HMM_HOME"] = phantomeHmmDBpath

    elif match_phageEnzymeHmmDBpath:
        value = match_phageEnzymeHmmDBpath.group(1)
        if value != '':
            phageEnzymeHmmDBpath = value
        os.environ["PHATE_PHAGE_ENZYME_HMM_HOME"] = phageEnzymeHmmDBpath

    elif match_keggVirusHmmDBpath:
        value = match_keggVirusHmmDBpath.group(1)
        if value != '':
            keggVirusHmmDBpath = value
        os.environ["PHATE_KEGG_HMM_HOME"] = keggVirusHmmDBpath

    elif match_pfamHmmDBpath:
        value = match_pfamHmmDBpath.group(1)
        if value != '':
            pfamHmmDBpath = value
        os.environ["PHATE_PFAM_HMM_HOME"] = pfamHmmDBpath

    elif match_smartHmmDBpath:
        value = match_smartHmmDBpath.group(1)
        if value != '':
            smartHmmDBpath = value
        os.environ["PHATE_SMART_HMM_HOME"] = smartHmmDBpath

    elif match_swissprotHmmDBpath:
        value = match_swissprotHmmDBpath.group(1)
        if value != '':
            swissprotHmmDBpath = value
        os.environ["PHATE_SWISSPROT_HMM_HOME"] = swissprotHmmDBpath

    elif match_uniprotHmmDBpath:
        value = match_uniprotHmmDBpath.group(1)
        if value != '':
            uniprotHmmDBpath = value
        os.environ["PHATE_UNIPROT_HMM_HOME"] = uniprotHmmDBpath

    elif match_nrHmmDBpath:
        value = match_nrHmmDBpath.group(1)
        if value != '':
            nrHmmDBpath = value
        os.environ["PHATE_NR_HMM_HOME"] = nrHmmDBpath

    elif match_customHmmDBpath:
        value = match_customHmmDBpath.group(1)
        if value != '':
            customHmmDBpath = value
        os.environ["PHATE_CUSTOM_HMM_HOME"] = customHmmDBpath

    ##### PARALLELISM #####

    elif match_threads:
        value = match_threads.group(1)
        if value == 'ALL' or value == 'All' or value == 'all':
            threads = value
        elif int(value) >= 1:
            threads = int(value) 
        else:
            threads = 1

    elif match_hpc:
        value = match_hpc.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            hpc = True
            HPC = hpc   #*** need some cleanup wrt this variable, now that it's been moved to config file

    elif match_blastThreads:
        value = match_blastThreads.group(1)
        if int(value) >= 1:
            blastThreads = int(value)

    ##### VERBOSITY #####

    elif match_phateWarnings:
        value = match_phateWarnings.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            os.environ["PHATE_PHATE_WARNINGS"] = 'True'
            PHATE_WARNINGS = True
        else:
            os.environ["PHATE_PHATE_WARNINGS"] = 'False'
            PHATE_WARNINGS = False

    elif match_phateMessages:
        value = match_phateMessages.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            os.environ["PHATE_PHATE_MESSAGES"] = 'True'
            PHATE_MESSAGES = True
        else:
            os.environ["PHATE_PHATE_MESSAGES"] = 'False'
            PHATE_MESSAGES = False

    elif match_phateProgress:
        value = match_phateProgress.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            os.environ["PHATE_PHATE_PROGRESS"] = 'True'
            PHATE_PROGRESS = True
        else:
            os.environ["PHATE_PHATE_PROGRESS"] = 'False'
            PHATE_PROGRESS = False

    elif match_cgcWarnings:
        value = match_cgcWarnings.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            os.environ["PHATE_CGC_WARNINGS"] = 'True'
        else:
            os.environ["PHATE_CGC_WARNINGS"] = 'False'

    elif match_cgcMessages:
        value = match_cgcMessages.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            os.environ["PHATE_CGC_MESSAGES"] = 'True'
        else:
            os.environ["PHATE_CGC_MESSAGES"] = 'False'

    elif match_cgcProgress:
        value = match_cgcProgress.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            os.environ["PHATE_CGC_PROGRESS"] = 'True'
        else:
            os.environ["PHATE_CGC_PROGRESS"] = 'False'

    elif match_cgpWarnings:
        value = match_cgpWarnings.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            os.environ["PHATE_CGP_WARNINGS"] = 'True'
        else:
            os.environ["PHATE_CGP_WARNINGS"] = 'False'

    elif match_cgpMessages:
        value = match_cgpMessages.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            os.environ["PHATE_CGP_MESSAGES"] = 'True'
        else:
            os.environ["PHATE_CGP_MESSAGES"] = 'False'

    elif match_cgpProgress:
        value = match_cgpProgress.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            os.environ["PHATE_CGP_PROGRESS"] = 'True'
        else:
            os.environ["PHATE_CGP_PROGRESS"] = 'False'

    elif match_cleanRawData:
        value = match_cleanRawData.group(1)
        if value.lower() == 'true' or value.lower() == 'yes' or value.lower() == 'on':
            os.environ["PHATE_CLEAN_RAW_DATA"] = 'True'
        else:
            os.environ["PHATE_CLEAN_RAW_DATA"] = 'False'

    else:
        if not HPC:
            LOG.write("%s%s\n" % ("ERROR: Unrecognized line in config file: ", cLine))
        print("multiPhate says, ERROR: unrecognized line in config file:", cLine)

if not HPC: # Skip logging if running in high-throughput, else multiPhate.log files will clash
    LOG.write("%s\n" % ("Input parameters and configurables:"))
    LOG.write("%s%s\n" % ("   BASE_DIR is ",os.environ["PHATE_BASE_DIR"]))
    LOG.write("%s%s\n" % ("   PHATE_BASE_DIR is ",os.environ["PHATE_PHATE_BASE_DIR"]))
    LOG.write("%s%s\n" % ("   DATABASE_DIR is ",os.environ["PHATE_DATABASE_DIR"]))
    LOG.write("%s%s\n" % ("   SOFTWARE_DIR is ",os.environ["PHATE_SOFTWARE_DIR"]))
    LOG.write("%s%s\n" % ("   GENE_FILE: ", GENE_FILE))
    LOG.write("%s%s\n" % ("   PROTEIN_FILE: ", PROTEIN_FILE))
    LOG.write("%s%s\n" % ("   geneticCode: ",geneticCode))
    LOG.write("%s%s\n" % ("   Status of boolean translateOnly is ",translateOnly))
    LOG.write("%s%s\n" % ("   primaryCalls is ",primaryCalls))
    LOG.write("%s%s\n" % ("   primaryCallsFile is ",primaryCallsFile))
    LOG.write("%s%s\n" % ("   genemarksCalls is ",genemarksCalls))
    LOG.write("%s%s\n" % ("   prodigalCalls is ",prodigalCalls))
    LOG.write("%s%s\n" % ("   glimmerCalls is ",glimmerCalls))
    LOG.write("%s%s\n" % ("   phanotateCalls is ",phanotateCalls))
    if customGeneCalls:
        LOG.write("%s%s\n" % ("   customGeneCalls is ",customGeneCalls))
        LOG.write("%s%s\n" % ("   customGeneCallerName is ",customGeneCallerName))
        LOG.write("%s%s\n" % ("   customGeneCallerOutfile is ",customGeneCallerOutfile))
    LOG.write("%s%s\n" % ("   blastpSearch is ",blastpSearch))
    LOG.write("%s%s\n" % ("   phmmerSearch is ",phmmerSearch))
    LOG.write("%s%s\n" % ("   jackhmmerSearch is ",jackhmmerSearch))
    LOG.write("%s%s\n" % ("   blastpIdentity is ",blastpIdentity))
    LOG.write("%s%s\n" % ("   blastnIdentity is ",blastnIdentity))
    LOG.write("%s%s\n" % ("   blastpHitCount is ",blastpHitCount))
    LOG.write("%s%s\n" % ("   blastnHitCount is ",blastnHitCount))
    LOG.write("%s%s\n" % ("   ncbiVirusGenomeBlast is ",ncbiVirusGenomeBlast))
    LOG.write("%s%s\n" % ("   ncbiVirusProteinBlast is ",ncbiVirusProteinBlast))
    LOG.write("%s%s\n" % ("   refseqProteinBlast is ",refseqProteinBlast))
    LOG.write("%s%s\n" % ("   refseqGeneBlast is ",refseqGeneBlast))
    LOG.write("%s%s\n" % ("   pvogsBlast is ",pvogsBlast))
    LOG.write("%s%s\n" % ("   vogsBlast is ",vogsBlast)) #*** To be deprecated
    LOG.write("%s%s\n" % ("   vogGeneBlast is ",vogGeneBlast))
    LOG.write("%s%s\n" % ("   vogProteinBlast is ",vogProteinBlast))
    LOG.write("%s%s\n" % ("   phantomeBlast is ",phantomeBlast))
    LOG.write("%s%s\n" % ("   phageEnzymeBlast is ",phageEnzymeBlast))
    LOG.write("%s%s\n" % ("   keggVirusBlast is ",keggVirusBlast))
    LOG.write("%s%s\n" % ("   pfamBlast is ",pfamBlast))
    LOG.write("%s%s\n" % ("   smartBlast is ",smartBlast))
    LOG.write("%s%s\n" % ("   swissprotBlast is ",swissprotBlast))
    LOG.write("%s%s\n" % ("   uniprotBlast is ",uniprotBlast))
    LOG.write("%s%s\n" % ("   nrBlast is ",nrBlast))
    LOG.write("%s%s\n" % ("   cazyBlast is ",cazyBlast))
    if customGenomeBlast:
        LOG.write("%s%s\n" % ("   customGenomeBlast is ",customGenomeBlast))
        LOG.write("%s%s\n" % ("   customGenomeDBname is ",customGenomeDBname))
        LOG.write("%s%s\n" % ("   customGenomeDBpath is ",customGenomeDBpath))
    if customGeneBlast:
        LOG.write("%s%s\n" % ("   customGeneBlast is ",customGeneBlast))
        LOG.write("%s%s\n" % ("   customGeneDBname is ",customGeneDBname))
        LOG.write("%s%s\n" % ("   customGeneDBpath is ",customGeneDBpath))
    if customProteinBlast:
        LOG.write("%s%s\n" % ("   customProteinBlast is ",customProteinBlast))
        LOG.write("%s%s\n" % ("   customProteinDBname is ",customProteinDBname))
        LOG.write("%s%s\n" % ("   customProteinDBpath is ",customProteinDBpath))
    LOG.write("%s%s\n" % ("   hmmscan is ",hmmscan))
    LOG.write("%s%s\n" % ("   hmmbuild is ",hmmbuild))
    LOG.write("%s%s\n" % ("   ncbiVirusGenomeHmm is ",ncbiVirusGenomeHmm))
    LOG.write("%s%s\n" % ("   ncbiVirusProteinHmm is ",ncbiVirusProteinHmm))
    LOG.write("%s%s\n" % ("   refseqProteinHmm is ",refseqProteinHmm))
    LOG.write("%s%s\n" % ("   refseqGeneHmm is ",refseqGeneHmm))
    LOG.write("%s%s\n" % ("   pvogsHmm is ",pvogsHmm))
    LOG.write("%s%s\n" % ("   vogsHmm is ",vogsHmm))
    LOG.write("%s%s\n" % ("   phantomeHmm is ",phantomeHmm))
    LOG.write("%s%s\n" % ("   phageEnzymeHmm is ",phageEnzymeHmm))
    LOG.write("%s%s\n" % ("   keggVirusHmm is ",keggVirusHmm))
    LOG.write("%s%s\n" % ("   pfamHmm is ",pfamHmm))
    LOG.write("%s%s\n" % ("   smartHmm is ",smartHmm))
    LOG.write("%s%s\n" % ("   swissprotHmm is ",swissprotHmm))
    LOG.write("%s%s\n" % ("   uniprotHmm is ",uniprotHmm))
    LOG.write("%s%s\n" % ("   nrHmm is ",nrHmm))
    if customHmm:
        LOG.write("%s%s\n" % ("   customHmm is ",customHmm))
        LOG.write("%s%s\n" % ("   customHmmName is ",customHmmName))
        LOG.write("%s%s\n" % ("   customHmmDBname is ",customHmmDBname))
        LOG.write("%s%s\n" % ("   customHmmDBpath is ",customHmmDBpath))
    LOG.write("%s%s\n" % ("   blast+ home is ",os.environ["PHATE_BLAST_HOME"]))
    LOG.write("%s%s\n" % ("   emboss home is ",os.environ["PHATE_EMBOSS_PHATE_HOME"]))
    LOG.write("%s%s\n" % ("   tRNAscanSE home is ",os.environ["PHATE_tRNAscanSE_HOME"]))
    LOG.write("%s%s\n" % ("   glimmer home is ",os.environ["PHATE_GLIMMER_PATH"]))
    LOG.write("%s%s\n" % ("   prodigal home is ",os.environ["PHATE_PRODIGAL_PATH"]))
    LOG.write("%s%s\n" % ("   phanotate home is ",os.environ["PHATE_PHANOTATE_PATH"]))
    LOG.write("%s%s\n" % ("   genemark home is ",os.environ["PHATE_GENEMARKS_PATH"]))
    LOG.write("%s%s\n" % ("   ncbi virus genome database is located in ",os.environ["PHATE_NCBI_VIRUS_GENOME_BLAST_HOME"]))
    LOG.write("%s%s\n" % ("   ncbi virus protein databases is located in ",os.environ["PHATE_NCBI_VIRUS_PROTEIN_BLAST_HOME"]))
    LOG.write("%s%s\n" % ("   refseq gene database is located in ",os.environ["PHATE_REFSEQ_GENE_BLAST_HOME"]))
    LOG.write("%s%s\n" % ("   refseq protein database is located in ",os.environ["PHATE_REFSEQ_PROTEIN_BLAST_HOME"]))
    LOG.write("%s%s\n" % ("   pVOGs database is located in ",os.environ["PHATE_PVOGS_BLAST_HOME"]))
    LOG.write("%s%s\n" % ("   VOGs database is located in ",os.environ["PHATE_VOGS_BLAST_HOME"]))
    LOG.write("%s%s\n" % ("   VOG Genes database is located in ",os.environ["PHATE_VOG_GENE_BLAST_HOME"]))
    LOG.write("%s%s\n" % ("   VOG Proteins database is located in ",os.environ["PHATE_VOG_PROTEIN_BLAST_HOME"]))
    LOG.write("%s%s\n" % ("   phantome database is located in ",os.environ["PHATE_PHANTOME_BLAST_HOME"]))
    LOG.write("%s%s\n" % ("   phage enzyme database is located in ",os.environ["PHATE_PHAGE_ENZYME_BLAST_HOME"]))
    LOG.write("%s%s\n" % ("   kegg virus database is located in ",os.environ["PHATE_KEGG_VIRUS_BLAST_HOME"]))
    LOG.write("%s%s\n" % ("   pfam database is located in ",os.environ["PHATE_PFAM_BLAST_HOME"]))
    LOG.write("%s%s\n" % ("   smart database is located in ",os.environ["PHATE_SMART_BLAST_HOME"]))
    LOG.write("%s%s\n" % ("   swissprot database is located in ",os.environ["PHATE_SWISSPROT_BLAST_HOME"]))
    LOG.write("%s%s\n" % ("   uniprot database is located in ",os.environ["PHATE_UNIPROT_BLAST_HOME"]))
    LOG.write("%s%s\n" % ("   NR database is located in ",os.environ["PHATE_NR_BLAST_HOME"]))
    LOG.write("%s%s\n" % ("   CAZy database is located in ",os.environ["PHATE_CAZY_BLAST_HOME"]))
    if customGenomeBlast:
        LOG.write("%s%s\n" % ("   Custom genome blast database is located in ",os.environ["PHATE_CUSTOM_GENOME_BLAST_HOME"]))
    if customGeneBlast:
        LOG.write("%s%s\n" % ("   Custom gene blast database is located in ",os.environ["PHATE_CUSTOM_GENE_BLAST_HOME"]))
    if customProteinBlast:
        LOG.write("%s%s\n" % ("   Custom protein blast database is located in ",os.environ["PHATE_CUSTOM_PROTEIN_BLAST_HOME"]))
    LOG.write("%s%s\n" % ("   ncbi virus genome hmm profiles are located in ",ncbiVirusGenomeHmmDBpath))
    LOG.write("%s%s\n" % ("   ncbi virus protein hmm profiles are located in ",ncbiVirusProteinHmmDBpath))
    LOG.write("%s%s\n" % ("   refseq gene hmm profiles are located in ",refseqGeneHmmDBpath))
    LOG.write("%s%s\n" % ("   refseq protein hmm profiles are located in ",refseqProteinHmmDBpath))
    LOG.write("%s%s\n" % ("   pvogs hmm profiles are located in ",pvogsHmmDBpath))
    LOG.write("%s%s\n" % ("   vogs hmm profiles are located in ",vogsHmmDBpath))
    LOG.write("%s%s\n" % ("   phantome hmm profiles are located in ",phantomeHmmDBpath))
    LOG.write("%s%s\n" % ("   phage enzyme hmm profiles are located in ",phageEnzymeHmmDBpath))
    LOG.write("%s%s\n" % ("   kegg virus hmm profiles are located in ",keggVirusHmmDBpath))
    LOG.write("%s%s\n" % ("   pfam hmm profiles are located in ",pfamHmmDBpath))
    LOG.write("%s%s\n" % ("   smart hmm profiles are located in ",smartHmmDBpath))
    LOG.write("%s%s\n" % ("   swissprot hmm profiles are located in ",swissprotHmmDBpath))
    LOG.write("%s%s\n" % ("   uniprot hmm profiles are located in ",uniprotHmmDBpath))
    LOG.write("%s%s\n" % ("   nr hmm profiles are located in ",nrHmmDBpath))
    if customHmm:
        LOG.write("%s%s\n" % ("   Custom hmm database is located in ",os.environ["PHATE_CUSTOM_HMM_HOME"]))
    LOG.write("%s%s\n" % ("   runCGP is ",runCGP))
    LOG.write("%s%s\n" % ("   runGenomics is ",runGenomics))
    LOG.write("%s%s\n" % ("   hmmbuild is ",hmmbuild))
    LOG.write("%s%s\n" % ("   hmmsearch is ",hmmsearch))
    LOG.write("%s%s\n" % ("   threads is ",threads))
    LOG.write("%s%s\n" % ("   HPC is ",HPC))
    LOG.write("%s%s\n" % ("   hpc is ",hpc))
    LOG.write("%s%s\n" % ("   blastThreads is ",blastThreads))
    LOG.write("%s%s\n" % ("   phate warnings is set to ",os.environ["PHATE_PHATE_WARNINGS"]))
    LOG.write("%s%s\n" % ("   phate messages is set to ",os.environ["PHATE_PHATE_MESSAGES"]))
    LOG.write("%s%s\n" % ("   phate progress is set to ",os.environ["PHATE_PHATE_PROGRESS"]))
    LOG.write("%s%s\n" % ("   cgc warnings is set to ",os.environ["PHATE_CGC_WARNINGS"]))
    LOG.write("%s%s\n" % ("   cgc messages is set to ",os.environ["PHATE_CGC_MESSAGES"]))
    LOG.write("%s%s\n" % ("   cgc progress is set to ",os.environ["PHATE_CGC_PROGRESS"]))
    LOG.write("%s%s\n" % ("   cgp warnings is set to ",os.environ["PHATE_CGP_WARNINGS"]))
    LOG.write("%s%s\n" % ("   cgp messages is set to ",os.environ["PHATE_CGP_MESSAGES"]))
    LOG.write("%s%s\n" % ("   cgp progress is set to ",os.environ["PHATE_CGP_PROGRESS"]))
    LOG.write("%s%s\n" % ("   clean raw data is set to ",os.environ["PHATE_CLEAN_RAW_DATA"]))
    LOG.write("%s%s\n" % ("Number of genomes to be processed: ",len(genomeList)))
    LOG.write("%s\n" % ("List of genomes to be processed:"))
for genome in genomeList:
    if not HPC:
        LOG.write("%s%c%s%c%s%c%s%c%s%c%s\n" % (genome["genomeNumber"],' ',genome["genomeName"],' ',genome["genomeType"],' ',genome["genomeSpecies"],' ',genome["genomeFile"],' ',genome["outputSubdir"]))

##### BEGIN MAIN ########################################################################################

# For each genome, create a phate.json file for running phate_runPipeline.py

nextJsonFile = ""
jsonList = []  # List of json filenames
for genome in genomeList:
    match_root = re.search(p_root,genome["genomeFile"])
    if match_root:
        genomeRoot = match_root.group(1)
        nextJsonFile = genomeRoot + '.json'
        jsonList.append(nextJsonFile)
        genomeFile = PIPELINE_INPUT_DIR + genome["genomeFile"]
    else:
        if PHATE_WARNINGS:
            print("multiPhate says, WARNING: Fasta filename not recognized for genome", genome["genomeName"])
            print("  Expected fasta filename extension: .fasta")
        if not HPC:
            LOG.write("%s%s%s\n" % ("WARNING: problem with fasta file: ",genome["genomeFile"]," This genome not processed!"))
        continue

    NEXT_JSON = open(nextJsonFile,"w")

    parameters = {
            "genomeNumber":genome["genomeNumber"],
            "genomeFile":genome["genomeFile"],
            "genomeType":genome["genomeType"],
            "genomeSpecies":genome["genomeSpecies"],
            "genomeName":genome["genomeName"],
            "outputSubdir":genome["outputSubdir"],
            "geneticCode":geneticCode,
            "translateOnly":translateOnly,
            "phanotateCalls":phanotateCalls,
            "prodigalCalls":prodigalCalls,
            "glimmerCalls":glimmerCalls,
            "genemarksCalls":genemarksCalls,
            "customGeneCalls":customGeneCalls,
            "customGeneCallerName":customGeneCallerName,
            "customGeneCallerOutfile":customGeneCallerOutfile,
            "primaryCalls":primaryCalls,
            "primaryCallsFile":primaryCallsFile,
            "blastnIdentity":blastnIdentity,
            "blastpIdentity":blastpIdentity,
            "blastnHitCount":blastnHitCount,
            "blastpHitCount":blastpHitCount,
            "ncbiVirusGenomeBlast":ncbiVirusGenomeBlast,
            "ncbiVirusProteinBlast":ncbiVirusProteinBlast,
            "refseqGeneBlast":refseqGeneBlast,
            "refseqProteinBlast":refseqProteinBlast,
            "pvogsBlast":pvogsBlast,
            "vogsBlast":vogsBlast,
            "vogGeneBlast":vogGeneBlast,
            "vogProteinBlast":vogProteinBlast,
            "phantomeBlast":phantomeBlast,
            "phageEnzymeBlast":phageEnzymeBlast,
            "keggVirusBlast":keggVirusBlast,
            "pfamBlast":pfamBlast,
            "smartBlast":smartBlast,
            "swissprotBlast":swissprotBlast,
            "uniprotBlast":uniprotBlast,
            "nrBlast":nrBlast,
            "cazyBlast":cazyBlast,
            "customGenomeBlast":customGenomeBlast,
            "customGenomeDBname":customGenomeDBname,
            "customGenomeDBpath":customGenomeDBpath,
            "customGeneBlast":customGeneBlast,
            "customGeneDBname":customGeneDBname,
            "customGeneDBpath":customGeneDBpath,
            "customProteinBlast":customProteinBlast,
            "customProteinDBname":customProteinDBname,
            "customProteinDBpath":customProteinDBpath,
            "ncbiVirusGenomeHmm":ncbiVirusGenomeHmm,
            "ncbiVirusProteinHmm":ncbiVirusProteinHmm,
            "refseqGeneHmm":refseqGeneHmm,
            "refseqProteinHmm":refseqProteinHmm,
            "pvogsHmm":pvogsHmm,
            "vogsHmm":vogsHmm,
            "phantomeHmm":phantomeHmm,
            "phageEnzymeHmm":phageEnzymeHmm,
            "keggVirusHmm":keggVirusHmm,
            "pfamHmm":pfamHmm,
            "smartHmm":smartHmm,
            "swissprotHmm":swissprotHmm,
            "uniprotHmm":uniprotHmm,
            "nrHmm":nrHmm,
            "blastpSearch":blastpSearch,
            "jackhmmerSearch":jackhmmerSearch,
            "phmmerSearch":phmmerSearch,
            "hmmscan":hmmscan,
            "customHmm":customHmm,
            "customHmmDBname":customHmmDBname,
            "customHmmDBpath":customHmmDBpath,
            "blastThreads":blastThreads,
            }

    # Write next json file
    json.dump(parameters, NEXT_JSON)

    NEXT_JSON.close()

##### Run the pipeline (phate_runPipeline.py) over each genome;
##### The code below runs in serial;
##### Modify this section to implement in parallel on your cluster

if not HPC:
    LOG.write("%s%s\n" % ("Processing genomes through PhATE. Begin processing at ",datetime.datetime.now()))

# Run the phate pipeline over each genome, as specified in the multiPhate.config file
def phate_threaded(jsonFile):
    print(f'multiPhate says, Running {jsonFile} on PID {os.getpid()}')
    if not HPC:
        LOG.write("%s%s\n" % ("Running PhATE using genome json file ",jsonFile))
    command = "python " + PHATE_BASE_DIR + PHATE_PIPELINE_CODE + " " + jsonFile
    if not HPC:
        LOG.write("%s%s\n" % ("Command is ",command))
        LOG.write("%s%s\n" % ("Begin PhATE processing at ",datetime.datetime.now()))
    result = os.system(command)
    if not HPC:
        LOG.write("%s%s\n" % ("End PhATE processing at ",datetime.datetime.now()))

# Set up to run processing in parallel on multiple threads
if threads is 'ALL':
    threads = os.cpu_count()
elif threads > os.cpu_count():
    threads = os.cpu_count()

#print(f'Using {THREADS} threads')
print(f'Using {threads} threads')

#pool = Pool(int(THREADS))
pool = Pool(int(threads))
pool.map(phate_threaded, jsonList)

pool.close()

# Run the comparative genomics module to compare proteoms among the user's specified genomes
GFF_OUTFILE = "phate_sequenceAnnotation_main.gff"
if runCGP and not translateOnly and (len(genomeList) > 1):
    if not HPC:
        LOG.write("%s%s\n" % ("Starting CompareGeneProfiles at ",datetime.datetime.now()))
    if PHATE_PROGRESS:
        print("multiPhATE says, Starting CompareGeneProfiles...")

    # Prepare the configuration file for input to CGP driver code
    cgpOutputDir = PIPELINE_OUTPUT_DIR
    try:
        if PHATE_PROGRESS:
            print("multiPhate says, Opening CGP_CONFIG_FILE:", CGP_CONFIG_FILE)
        CGP_CONFIG_H = open(CGP_CONFIG_FILE,"w")
        # Print output directory to config file
        CGP_CONFIG_H.write("%s\n" % (cgpOutputDir))
        for genome in genomeList:
            genomeFile = PIPELINE_INPUT_DIR  + genome["genomeFile"]
            gffFile    = PIPELINE_OUTPUT_DIR + genome["outputSubdir"] + GFF_OUTFILE
            # Print data line for each genome
            dataLine = genomeFile + ' ' + gffFile
            CGP_CONFIG_H.write("%s\n" % (dataLine))
        CGP_CONFIG_H.close()
    except:
        if not HPC:
            LOG.write("%s%s\n" % ("Error preparing CompareGeneProfiles configuration file,", CGP_CONFIG_FILE))
        print ("multiPhate says, ERROR preparing CompareGeneProfiles configuration file")

    # Run CompareGeneProfiles
    try:
        command = "python " + CGP_CODE + ' ' + CGP_CONFIG_FILE
        if PHATE_MESSAGES:
            print("multiPhate says, command is:", command)
        result = os.system(command)
    except:
        if not HPC:
            LOG.write("%s\n" % ("ERROR: Problem running CompareGeneProfiles"))
        print ("multiPhate says, ERROR running CompareGeneProfiles")

    if PHATE_PROGRESS:
        print("multiPhate says, Completed CompareGeneProfiles. Moving results to CGP_RESULTS/ directory")

    if runCGP:

        # Create output directory for loading CompareGeneProfiles results
        try:
            os.stat(CGP_RESULTS_DIR)
        except:
            os.mkdir(CGP_RESULTS_DIR)

        # Move CompareGeneProfiles results directories and files to the CGP_RESULTS_DIR directory
        if CLEAN_PREVIOUS_CGP_RESULTS:
            if PHATE_PROGRESS:
                print ("multiPhate says, Removing previous CGP Results_ directories")
            try:
                command = "rm -r " + CGP_RESULTS_DIR + "Results_* "
                result = os.system(command)
            except:
                if PHATE_WARNINGS:
                    print("multiPhate says, Failure in removing previous CGP Results directories")
        else:
            if PHATE_WARNINGS:
                print("multiPhate says, CAUTION: Not removing previous CGP Results_ directories--data files will accumulate!")

        try:
            # move directories and files from main working directory to CGP results directory
            command = "mv " + PIPELINE_OUTPUT_DIR + "Results_* "                 + CGP_RESULTS_DIR + '.'
            result = os.system(command)
            command = "mv " + PIPELINE_OUTPUT_DIR + "cgpNxN.config "             + CGP_RESULTS_DIR + '.'
            result = os.system(command)
            command = "mv " + PIPELINE_OUTPUT_DIR + "cgp_wrapper.config "        + CGP_RESULTS_DIR + '.'
            result = os.system(command)

            # remove excess copies of CGP out files from final comparison (these are duplicates)
            command = "rm " + PIPELINE_OUTPUT_DIR + "compareGeneProfiles_main.log.copy"
            result = os.system(command)
            command = "rm " + PIPELINE_OUTPUT_DIR + "compareGeneProfiles_main.out"
            result = os.system(command)
            command = "rm " + PIPELINE_OUTPUT_DIR + "compareGeneProfiles_main.report"
            result = os.system(command)
            command = "rm " + PIPELINE_OUTPUT_DIR + "compareGeneProfiles_main.summary"
            result = os.system(command)
        except:
            print("multiPhate says, ERROR: cannot move CGP results to output directory")

else:
    if PHATE_PROGRESS:
        print("multiPhate says, Skipping CompareGeneProfiles.")
    LOG.write("%s\n" % ("Skipping CompareGeneProfiles."))

# Turn off Genomics processing if user selected translate only
if translateOnly == True:
    runGenomics = False

# Override
#runGenomics = False
if runGenomics:
    if PHATE_PROGRESS:
        print("multiPhate says, Performing gene-correspondence analysis.")

    # Invoke the genomics module
    command = "python " + GENOMICS_CODE
    result = os.system(command)
    if PHATE_PROGRESS:
        print("multiPhate says, Result of genomics processing is:",result)

    # Clean up
    if PHATE_PROGRESS:
        print("multiPhate says, Genomics completed.")

else:
    if PHATE_PROGRESS:
        print("multiPhate says, Skipping Genomics.")
    LOG.write("%s\n" % ("Skipping Genomics."))

##### CLEAN UP

# Move json files to a subdir to keep main dir tidy
try:
    os.stat(JSON_DIR)
except:
    os.mkdir(JSON_DIR)
command = "mv *.json " + JSON_DIR + '.'
result = os.system(command)

if not HPC:
    LOG.write("%s%s\n" % ("End log file ",datetime.datetime.now()))
    LOG.close()

if PHATE_MESSAGES:
    print("multiPhate says, Processing complete. Thank you for using multiPhATE2.")

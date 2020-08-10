#!/usr/bin/env python

################################################################
#
# phate_sequenceAnnotation_main.py
#
# Description:  Performs blast and/or hmm searches on a given input gene or protein databases.
#    Databases may be blast, sequence, or hmm profiles. Current search programs supported are:
#    blastn, blastp, jackhmmer, phmmer, hmmscan.
#
# Programmer: CEZhou
#
# Latest Update: 07 August 2020
# Version 1.5
#
################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.PDF FOR DETAILS.

import sys, os, re, string, copy
import time, datetime
from subprocess import call

DEBUG = False 
#DEBUG = True

# Defaults/Parameters
PRIMARY_CALLS          = 'phanotate'   # Default; can be configured by user
GENOME_TYPE            = 'phage'       # Default; can be configured by user
EVALUE_MIN             = 10            # This and the following blast parameters might be parameterized, eventually
EVALUE_SELECT          = 10
XML_OUT_FORMAT         = 5
LIST_OUT_FORMAT        = 7
SCORE_EDGE             = 0.1
OVERHANG               = 0.1
GENOME_IDENTITY_MIN    = 20
GENOME_IDENTITY_SELECT = 20
GENETIC_CODE           = 11
TOP_HIT_COUNT          = 5             # Limit on number of hits to report

# Constants
CODE_BASE = "phate_sequenceAnnotation_main"
CODE = CODE_BASE + ".py"
CODE_OUT_FILE = CODE_BASE + ".out"
LOGFILE = CODE_BASE + ".log"
OUTFILE = CODE_BASE + ".out"
GFFFILE = CODE_BASE + ".gff"

# Get environment variables (set by phate_runPipeline.py)
# These are required by this code and modules that are imported (below)
# General
PIPELINE_DIR                  = os.environ["PHATE_PIPELINE_DIR"]
CODE_BASE_DIR                 = os.environ["PHATE_PIPELINE_DIR"]
BASE_DIR                      = os.environ["PHATE_BASE_DIR"]
DATABASE_DIR                  = os.environ["PHATE_DATABASE_DIR"]
SOFTWARE_DIR                  = os.environ["PHATE_SOFTWARE_DIR"]

# Codes
BLAST_HOME                    = os.environ["PHATE_BLAST_HOME"]
EMBOSS_HOME                   = os.environ["PHATE_EMBOSS_PHATE_HOME"]
PHANOTATE_HOME                = os.environ["PHATE_PHANOTATE_PATH"]
PRODIGAL_HOME                 = os.environ["PHATE_PRODIGAL_PATH"]
GLIMMER_HOME                  = os.environ["PHATE_GLIMMER_PATH"]
GENEMARKS_HOME                = os.environ["PHATE_GENEMARKS_PATH"]
CGC_HOME                      = os.environ["PHATE_CGC_PATH"]
tRNAscanSE_HOME               = os.environ["PHATE_tRNAscanSE_HOME"]
HMMER_HOME                    = os.environ["PHATE_HMMER_HOME"]

# BLAST parameters - needed to formulate parameters input to blast module
MIN_BLASTP_IDENTITY           = os.environ["PHATE_MIN_BLASTP_IDENTITY"]   # Sets a lower limit
MIN_BLASTN_IDENTITY           = os.environ["PHATE_MIN_BLASTN_IDENTITY"]   # Sets a lower limit
MAX_BLASTP_HIT_COUNT          = os.environ["PHATE_MAX_BLASTP_HIT_COUNT"]  # Sets an upper limit
MAX_BLASTN_HIT_COUNT          = os.environ["PHATE_MAX_BLASTN_HIT_COUNT"]  # Sets an upper limit
BLASTP_IDENTITY_DEFAULT       = os.environ["PHATE_BLASTP_IDENTITY_DEFAULT"]
BLASTN_IDENTITY_DEFAULT       = os.environ["PHATE_BLASTN_IDENTITY_DEFAULT"]
BLASTP_HIT_COUNT_DEFAULT      = os.environ["PHATE_BLASTP_HIT_COUNT_DEFAULT"]
BLASTN_HIT_COUNT_DEFAULT      = os.environ["PHATE_BLASTN_HIT_COUNT_DEFAULT"]

# BLAST database locations 
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
#VOG_GENE_BASE_DIR             = os.environ["PHATE_VOG_GENE_BASE_DIR"]
VOG_GENE_BLAST_HOME           = os.environ["PHATE_VOG_GENE_BLAST_HOME"]
#VOG_PROTEIN_BASE_DIR          = os.environ["PHATE_VOG_PROTEIN_BASE_DIR"]
VOG_PROTEIN_BLAST_HOME        = os.environ["PHATE_VOG_PROTEIN_BLAST_HOME"]
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
CAZY_BLAST_BASE_DIR           = os.environ["PHATE_CAZY_BASE_DIR"]
CAZY_BLAST_HOME               = os.environ["PHATE_CAZY_BLAST_HOME"]
CUSTOM_GENOME_BLAST_HOME      = os.environ["PHATE_CUSTOM_GENOME_BLAST_HOME"]
CUSTOM_GENE_BLAST_HOME        = os.environ["PHATE_CUSTOM_GENE_BLAST_HOME"]
CUSTOM_PROTEIN_BLAST_HOME     = os.environ["PHATE_CUSTOM_PROTEIN_BLAST_HOME"]

# HMM database locations
NCBI_VIRUS_GENOME_HMM_HOME    = os.environ["PHATE_NCBI_VIRUS_GENOME_HMM_HOME"]
NCBI_VIRUS_PROTEIN_HMM_HOME   = os.environ["PHATE_NCBI_VIRUS_PROTEIN_HMM_HOME"]
REFSEQ_PROTEIN_HMM_HOME       = os.environ["PHATE_REFSEQ_PROTEIN_HMM_HOME"]
REFSEQ_GENE_HMM_HOME          = os.environ["PHATE_REFSEQ_GENE_HMM_HOME"]
PVOGS_HMM_HOME                = os.environ["PHATE_PVOGS_HMM_HOME"]
VOGS_HMM_HOME                 = os.environ["PHATE_VOGS_HMM_HOME"]
PHANTOME_HMM_HOME             = os.environ["PHATE_PHANTOME_HMM_HOME"]
PHAGE_ENZYME_HMM_HOME         = os.environ["PHATE_PHAGE_ENZYME_BLAST_HOME"]
PFAM_HMM_HOME                 = os.environ["PHATE_PFAM_BLAST_HOME"]
SMART_HMM_HOME                = os.environ["PHATE_SMART_BLAST_HOME"]
SWISSPROT_HMM_HOME            = os.environ["PHATE_SWISSPROT_BLAST_HOME"]
UNIPROT_HMM_HOME              = os.environ["PHATE_UNIPROT_BLAST_HOME"]
NR_HMM_HOME                   = os.environ["PHATE_NR_HMM_HOME"]
CUSTOM_HMM_HOME               = os.environ["PHATE_CUSTOM_HMM_HOME"]

# Verbosity
CLEAN_RAW_DATA                = os.environ["PHATE_CLEAN_RAW_DATA"]
PHATE_WARNINGS                = os.environ["PHATE_PHATE_WARNINGS"]
PHATE_MESSAGES                = os.environ["PHATE_PHATE_MESSAGES"]
PHATE_PROGRESS                = os.environ["PHATE_PHATE_PROGRESS"]

import phate_fastaSequence  # generic fasta sequence module
import phate_genomeSequence # manages genomes to be annotated
import phate_annotation     # records annotations, including gene-call info and secondary annotations
import phate_blast          # runs blast search against blast database(s) 
import phate_hmm            # runs hmm search against specified sequence database(s) 
import phate_profile        # runs hmm search against specified hmm profile database(s)

##### FILES

logfile = ""
outfile = ""  # will be constructed using user's specified output subdir   

outputDir           = ""         # user-specified output subdirectory
infile_genome       = ""         # user-provided file containing the genome sequence that was gene-called
infile_geneCall     = ""         # gene-call file (output from a gene caller; PHANOTATE for now)
infile_gene         = ""         # genes to be blasted (not yet in service)
infile_protein      = ""         # user-provided file containing protein fasta sequences
infile_primaryCalls = ""         # will be assigned later; depends on user's choice

##### USER-SPECIFIED META-DATA and PARAMTERS

genomeType      = GENOME_TYPE              # user-provided, typically 'phage'
genomeName      = "unknown"                # user-provided, name of current genome, e.g., 'LYP215'
genomeSpecies   = "unknown"                # user-provided, e.g., 'YpPhage_LYP215'
configFile      = "unknown"                # must be passed by calling code (phate_runPipeline.py)
primaryCalls    = PRIMARY_CALLS            # default, unless changed
blastpIdentity  = BLASTP_IDENTITY_DEFAULT  # integer percent identity cutoff
blastnIdentity  = BLASTN_IDENTITY_DEFAULT  # integer percent identity cutoff
blastpHitCount  = BLASTP_HIT_COUNT_DEFAULT # number of top hits to capture
blastnHitCount  = BLASTN_HIT_COUNT_DEFAULT # number of top blastn hits to capture 
geneticCode     = GENETIC_CODE             # default, unless changed

##### BOOLEANS:  These may be changed by input parameters

# If True, annotate only through translation of gene calls, then stop 
TRANSLATE_ONLY           = False

# Annotation processes to perform
RUN_BLAST                = False         # any flavor of blast
RUN_HMM_SEARCH           = False         # any flavor of hmm search
RUN_PROFILE_SEARCH       = False         # any flavor of hmm profile search
RUN_GENOME_BLAST         = False         # blasting against genome database
RUN_GENE_BLAST           = False         # blasting against gene database
RUN_PROTEIN_BLAST        = False         # blasting against any protein database
BLASTP_SEARCH            = False         # If True, run blastp against fasta blast DB(s)
BLAST_THREADS            = 1             # Positive integer; assume serial execution unless indicated otherwise
PHMMER_SEARCH            = False         # if True, run phmmer against fasta blast DB(s)
JACKHMMER_SEARCH         = False         # if True, run jackhmmer against fasta blast DB(s)
HMMSCAN                  = False         # Default; can be configured by user ###*** to be deprecated

# Databases: assume turned 'off', unless input string indicates otherwise
# BLAST Fasta Databases
NCBI_VIRUS_GENOME_BLAST  = False
NCBI_VIRUS_PROTEIN_BLAST = False
REFSEQ_PROTEIN_BLAST     = False
REFSEQ_GENE_BLAST        = False
PVOGS_BLAST              = False
VOGS_BLAST               = False  #*** To be deprecated
VOG_GENE_BLAST           = False
VOG_PROTEIN_BLAST        = False
PHANTOME_BLAST           = False
PHAGE_ENZYME_BLAST       = False
KEGG_VIRUS_BLAST         = False
PFAM_BLAST               = False
SMART_BLAST              = False
SWISSPROT_BLAST          = False
UNIPROT_BLAST            = False
NR_BLAST                 = False
CAZY_BLAST               = False
CUSTOM_GENOME_BLAST      = False
CUSTOM_GENE_BLAST        = False
CUSTOM_PROTEIN_BLAST     = False
# HMM Profiles Databases
NCBI_VIRUS_GENOME_HMM    = False
NCBI_VIRUS_PROTEIN_HMM   = False
REFSEQ_PROTEIN_HMM       = False
REFSEQ_GENE_HMM          = False
KEGG_VIRUS_HMM           = False
PVOGS_HMM                = False
VOGS_HMM                 = False
PHANTOME_HMM             = False
PHAGE_ENZYME_HMM         = False
SWISSPROT_HMM            = False
PFAM_HMM                 = False
SMART_HMM                = False
SWISSPROT_HMM            = False
UNIPROT_HMM              = False
NR_HMM                   = False
CUSTOM_HMM               = False

##### PATTERNS and CONTROL

# Input parameters 

# parameter tags
p_outputDirParam             = re.compile('^-o')   # outout directory (e.g., 'LYP215')
p_genomeFileParam            = re.compile('^-G')   # Genome with a capital 'G'
p_geneFileParam              = re.compile('^-g')   # gene with a lower-case 'g'
p_proteinFileParam           = re.compile('^-p')   # protein or peptide
p_geneticCodeParam           = re.compile('^-e')   # genetic code
p_primaryCallsParam          = re.compile('^-c')   # gene caller = primary calls
p_primaryCallsPathFileParam  = re.compile('^-f')   # primary gene calls file
p_genomeTypeParam            = re.compile('^-t')   # genome type (e.g., 'phage', 'bacterium')
p_genomeNameParam            = re.compile('^-n')   # genome name (e.g., 'my_fave_Yp_genome')
p_genomeSpeciesParam         = re.compile('^-s')   # species (e.g., 'Y_pestis') 
p_blastpIdentityParam        = re.compile('^-i')   # blastp identity cutoff
p_blastnIdentityParam        = re.compile('^-j')   # blastn identity cutoff
p_blastpHitCountParam        = re.compile('^-h')   # blastp top hit count
p_blastnHitCountParam        = re.compile('^-H')   # blastn top hit count
p_blastThreadsParam          = re.compile('^-z')   # number of blast threads to execute
p_translateOnlyParam         = re.compile('^-x')   # if user passes 'true' => get genes, translate, compare, then stop before annotation
p_blastDatabaseStringParam   = re.compile('^-b')   # string listing database(s) to blast against
p_blastProgramStringParam    = re.compile('^-B')   # string listing hmm program and database(s) to search against
p_seqDatabaseStringParam     = re.compile('^-m')   # string listing hmm database(s) to search against
p_hmmProgramStringParam      = re.compile('^-M')   # string listing hmm program(s) to search sequence databases 
p_profileDatabaseStringParam = re.compile('^-r')   # string listing hmm profile database(s) to search
p_profileProgramStringParam  = re.compile('^-R')   # string listing hmm program(s) for searching hmm profile databases
p_customDatabasePathParam    = re.compile('^-C')   # string listing path/filenames to custom database(s)

# Parts of input strings naming databases to blast or hmm-search against

# blast DBs
p_ncbiVirusGenome      = re.compile('ncbiVirusGenome')     
p_ncbiVirusProtein     = re.compile('ncbiVirusProtein')    
p_refseqProtein        = re.compile('refseqP')             
p_refseqGene           = re.compile('refseqG')             
p_pvogs                = re.compile('pvogs')               
p_vogs                 = re.compile('Vogs')     # *** To be deprecated
p_vogGene              = re.compile('vogGene')
p_vogProtein           = re.compile('vogProtein')
p_phantome             = re.compile('phantome')            
p_phageEnzyme          = re.compile('phageEnzyme')               
p_keggVirus            = re.compile('kegg')                
p_pfam                 = re.compile('pfam')                
p_smart                = re.compile('smart')               
p_swissprot            = re.compile('swissprot')           
p_uniprot              = re.compile('uniprot')             
p_nr                   = re.compile('nr') 
p_cazy                 = re.compile('cazy')
p_customGenome         = re.compile('customGenome')
p_customGene           = re.compile('customGene')
p_customProtein        = re.compile('customProtein')

# hmm profiles DBs
p_ncbiVirusGenomeHmm   = re.compile('ncbiVirusGenomeHmm')  
p_ncbiVirusProteinHmm  = re.compile('ncbiVirusProteinHmm') 
p_refseqProteinHmm     = re.compile('refseqPhmm')          
p_refseqGeneHmm        = re.compile('refseqGhmm')          
p_pvogsHmm             = re.compile('pvogsHmm')            
p_vogsHmm              = re.compile('VogsHmm')            
p_phantomeHmm          = re.compile('phantomeHmm')         
p_phageEnzymeHmm       = re.compile('phageEnzymeHmm')            
p_keggVirusHmm         = re.compile('keggHmm')             
p_pfamHmm              = re.compile('pfamHmm')             
p_smartHmm             = re.compile('smartHmm')            
p_swissprotHmm         = re.compile('swissprotHmm')        
p_uniprotHmm           = re.compile('uniprotHmm')          
p_nrHmm                = re.compile('nrHmm')  
p_customHmm            = re.compile('customHmm')

# custom DB Names 
p_customGenomeDBname   = re.compile('customGenomeDBname:')
p_customGeneDBname     = re.compile('customGeneDBname:')
p_customProteinDBname  = re.compile('customProteinDBname:')
p_customHmmDBname      = re.compile('customHmmDBname:')

# custom DB Paths
p_customGenomeDBpath   = re.compile('customGenomeDBpath:')
p_customGeneDBpath     = re.compile('customGeneDBpath:')
p_customProteinDBpath  = re.compile('customProteinDBpath:')
p_customHmmDBpath      = re.compile('customHmmDBpath:')

# programs 
p_blastp               = re.compile('blastp')
p_phmmer               = re.compile('phmmer')
p_jackhmmer            = re.compile('jackhmmer')
p_hmmscan              = re.compile('hmmscan')

# Other patterns
p_comment  = re.compile('^#')
p_blank    = re.compile("^\s*$")

##### GET INPUT PARAMETERS #####

argList  = sys.argv
argCount = len(argList)

for i in range(0,argCount):

    # Look for parameter tags
    match_outputDirParam             = re.search(p_outputDirParam,             argList[i])
    match_genomeFileParam            = re.search(p_genomeFileParam,            argList[i])
    match_geneFileParam              = re.search(p_geneFileParam,              argList[i])
    match_proteinFileParam           = re.search(p_proteinFileParam,           argList[i])
    match_primaryCallsParam          = re.search(p_primaryCallsParam,          argList[i])
    match_primaryCallsPathFileParam  = re.search(p_primaryCallsPathFileParam,  argList[i])
    match_geneticCodeParam           = re.search(p_geneticCodeParam,           argList[i])
    match_genomeTypeParam            = re.search(p_genomeTypeParam,            argList[i])
    match_genomeNameParam            = re.search(p_genomeNameParam,            argList[i])
    match_genomeSpeciesParam         = re.search(p_genomeSpeciesParam,         argList[i])
    match_translateOnlyParam         = re.search(p_translateOnlyParam,         argList[i]) # if True, stop after gene translation

    match_blastpIdentityParam        = re.search(p_blastpIdentityParam,        argList[i])
    match_blastnIdentityParam        = re.search(p_blastnIdentityParam,        argList[i])
    match_blastpHitCountParam        = re.search(p_blastpHitCountParam,        argList[i])
    match_blastnHitCountParam        = re.search(p_blastnHitCountParam,        argList[i])
    match_blastThreadsParam          = re.search(p_blastThreadsParam,          argList[i])

    match_blastDatabaseStringParam   = re.search(p_blastDatabaseStringParam,   argList[i]) # blast databases
    match_blastProgramStringParam    = re.search(p_blastProgramStringParam,    argList[i]) # blast programs for blast database search
    match_seqDatabaseStringParam     = re.search(p_seqDatabaseStringParam,     argList[i]) # sequence databases (same as for blast)
    match_hmmProgramStringParam      = re.search(p_hmmProgramStringParam,      argList[i]) # hmm programs for sequence database search
    match_profileDatabaseStringParam = re.search(p_profileDatabaseStringParam, argList[i]) # hmm profile databases
    match_profileProgramStringParam  = re.search(p_profileProgramStringParam,  argList[i]) # hmm programs for hmm profile database search
    match_customDatabasePathParam    = re.search(p_customDatabasePathParam,    argList[i]) # path/filenames to custom database(s)

    ### Capture parameters 

    # Filenames that are tagged 

    if match_outputDirParam:
        if i < argCount:
            outputDir = argList[i+1]
            logfile = outputDir + LOGFILE
            outfile = outputDir + OUTFILE
            gfffile = outputDir + GFFFILE

    if match_genomeFileParam:
        if i < argCount:
            infile_genome = argList[i+1]
        
    if match_geneFileParam:
        if i < argCount:
            infile_gene = argList[i+1] 

    if match_proteinFileParam:
        if i < argCount:
            infile_protein = argList[i+1] 

    # Other parameterized arguments

    if match_primaryCallsParam:
        if i < argCount:
            primaryCalls = argList[i+1]

    if match_primaryCallsPathFileParam:
        if i < argCount:
            infile_primaryCalls = argList[i+1]

    if match_genomeTypeParam:
        if i < argCount:
            genomeType = argList[i+1]

    if match_geneticCodeParam:
        if i < argCount:
            geneticCode = argList[i+1]

    if match_genomeNameParam:
        if i < argCount:
            genomeName = argList[i+1]
 
    if match_genomeSpeciesParam:
        if i < argCount:
            genomeSpecies = argList[i+1]
  
    if match_translateOnlyParam:
        if i < argCount:
            value = argList[i+1]
            if value.lower() == 'yes' or value.lower() == 'true' or value.lower() == 'on' or value.lower() == 'y':
                TRANSLATE_ONLY = True
                translateOnly = True
            else:
                TRANSLATE_ONLY = False 
                translateOnly = False 

    # Blast parameters

    if match_blastpIdentityParam:
        if i < argCount:
            value = argList[i+1]
            if int(value) > int(MIN_BLASTP_IDENTITY) and int(value) <= 100:
                blastpIdentity = int(value)

    if match_blastnIdentityParam:
        if i < argCount:
            value = argList[i+1]
            if int(value) > int(MIN_BLASTN_IDENTITY) and int(value) <= 100:
                blastnIdentity = int(value)

    if match_blastpHitCountParam:
        if i < argCount:
            value = argList[i+1]
            if int(value) > 0 and int(value) <= int(MAX_BLASTP_HIT_COUNT):
                blastpHitCount = int(value)

    if match_blastnHitCountParam:
        if i < argCount:
            value = argList[i+1]
            if int(value) > 0 and int(value) <= int(MAX_BLASTN_HIT_COUNT):
                blastnHitCount = int(value)

    if match_blastThreadsParam:
        if i < argCount:
            blastThreads = int(argList[i+1])

    # Blast, Hmm, and Profile database processing

    if match_blastDatabaseStringParam: # blast databases to use (blast DBs are seq DBs formatted with makeblastdb)
        if i < argCount:
            value = argList[i+1]
            match_ncbiVirusGenome  = re.search(p_ncbiVirusGenome,value)
            match_ncbiVirusProtein = re.search(p_ncbiVirusProtein,value)
            match_refseqProtein    = re.search(p_refseqProtein,value)
            match_refseqGene       = re.search(p_refseqGene,value)
            match_pvogs            = re.search(p_pvogs,value)
            match_vogs             = re.search(p_vogs,value)
            match_vogGene          = re.search(p_vogGene,value)
            match_vogProtein       = re.search(p_vogProtein,value)
            match_phantome         = re.search(p_phantome,value)
            match_phageEnzyme      = re.search(p_phageEnzyme,value)
            match_keggVirus        = re.search(p_keggVirus,value)
            match_pfam             = re.search(p_pfam,value)
            match_smart            = re.search(p_smart,value)
            match_swissprot        = re.search(p_swissprot,value)
            match_uniprot          = re.search(p_uniprot,value)
            match_nr               = re.search(p_nr,value)
            match_cazy             = re.search(p_cazy,value)
            match_customGenome     = re.search(p_customGenome,value)
            match_customGene       = re.search(p_customGene,value)
            match_customProtein    = re.search(p_customProtein,value)
            if match_ncbiVirusGenome:
                NCBI_VIRUS_GENOME_BLAST = True
            if match_ncbiVirusProtein:
                NCBI_VIRUS_PROTEIN_BLAST = True
            if match_refseqProtein:
                REFSEQ_PROTEIN_BLAST = True
            if match_refseqGene:
                REFSEQ_GENE_BLAST = True
            if match_pvogs:
                PVOGS_BLAST = True 
            if match_vogs:
                VOGS_BLAST = True 
            if match_vogGene:
                VOG_GENE_BLAST = True
            if match_vogProtein:
                VOG_PROTEIN_BLAST = True
            if match_phantome:
                PHANTOME_BLAST = True
            if match_phageEnzyme:
                PHAGE_ENZYME_BLAST = True
            if match_keggVirus:
                KEGG_VIRUS_BLAST = True
            if match_pfam:
                PFAM_BLAST = True
            if match_smart:
                SMART_BLAST = True
            if match_swissprot:
                SWISSPROT_BLAST = True
            if match_uniprot:
                UNIPROT_BLAST = True
            if match_nr:
                NR_BLAST = True
            if match_cazy:
                CAZY_BLAST = True
            if match_customGenome:
                CUSTOM_GENOME_BLAST = True
            if match_customGene:
                CUSTOM_GENE_BLAST = True
            if match_customProtein:
                CUSTOM_PROTEIN_BLAST = True

    if match_blastProgramStringParam: # blast programs to run against (protein) blast databases
        if i < argCount:
            value = argList[i+1]
            match_blastp = re.search(p_blastp, value)
            if match_blastp:
                BLASTP_SEARCH = True

    if match_seqDatabaseStringParam: # sequence databases to search using hmm programs (seq DBs are in blast DB folders)
        if i < argCount:
            value = argList[i+1]
            match_ncbiVirusGenome  = re.search(p_ncbiVirusGenome,value)
            match_ncbiVirusProtein = re.search(p_ncbiVirusProtein,value)
            match_refseqProtein    = re.search(p_refseqProtein,value)
            match_refseqGene       = re.search(p_refseqGene,value)
            match_pvogs            = re.search(p_pvogs,value)
            match_vogs             = re.search(p_vogs,value)
            match_vogGene          = re.search(p_vogGene,value)
            match_vogProtein       = re.search(p_vogProtein,value)
            match_phantome         = re.search(p_phantome,value)
            match_phageEnzyme      = re.search(p_phageEnzyme,value)
            match_keggVirus        = re.search(p_keggVirus,value)
            match_pfam             = re.search(p_pfam,value)
            match_smart            = re.search(p_smart,value)
            match_swissprot        = re.search(p_swissprot,value)
            match_uniprot          = re.search(p_uniprot,value)
            match_nr               = re.search(p_nr,value)
            match_cazy             = re.search(p_cazy,value)
            match_customGenome     = re.search(p_customGenome,value)
            match_customGene       = re.search(p_customGene,value)
            match_customProtein    = re.search(p_customProtein,value)
            if match_ncbiVirusGenome:
                NCBI_VIRUS_GENOME_BLAST = True
            if match_ncbiVirusProtein:
                NCBI_VIRUS_PROTEIN_BLAST = True
            if match_refseqProtein:
                REFSEQ_PROTEIN_BLAST = True
            if match_refseqGene:
                REFSEQ_GENE_BLAST = True
            if match_pvogs:
                PVOGS_BLAST = True 
            if match_vogs:
                VOGS_BLAST = True 
            if match_vogGene:
                VOG_GENE_BLAST = True 
            if match_vogProtein:
                VOG_PROTEIN_BLAST = True 
            if match_phantome:
                PHANTOME_BLAST = True
            if match_phageEnzyme:
                PHAGE_ENZYME_BLAST = True
            if match_keggVirus:
                KEGG_VIRUS_BLAST = True
            if match_pfam:
                PFAM_BLAST = True
            if match_smart:
                SMART_BLAST = True
            if match_swissprot:
                SWISSPROT_BLAST = True
            if match_uniprot:
                UNIPROT_BLAST = True
            if match_nr:
                NR_BLAST = True
            if match_cazy:
                CAZY_BLAST = True
            if match_customGenome:
                CUSTOM_GENOME_BLAST = True
            if match_customGene:
                CUSTOM_GENE_BLAST = True
            if match_customProtein:
                CUSTOM_PROTEIN_BLAST = True

    if match_hmmProgramStringParam: # hmm programs for searching sequence databases
        if i < argCount:
            value = argList[i+1]
            match_phmmer    = re.search(p_phmmer,value)
            match_jackhmmer = re.search(p_jackhmmer,value)
            if match_phmmer:
                PHMMER_SEARCH = True
            if match_jackhmmer:
                JACKHMMER_SEARCH = True

    if match_profileDatabaseStringParam: # hmm profile databases to be searched
        if i < argCount:
            value = argList[i+1]
            match_ncbiVirusGenomeHmm  = re.search(p_ncbiVirusGenomeHmm,value)
            match_ncbiVirusProteinHmm = re.search(p_ncbiVirusProteinHmm,value)
            match_refseqProteinHmm    = re.search(p_refseqProteinHmm,value)
            match_refseqGeneHmm       = re.search(p_refseqGeneHmm,value)
            match_pvogsHmm            = re.search(p_pvogsHmm,value)
            match_vogsHmm             = re.search(p_vogsHmm,value)
            match_phantomeHmm         = re.search(p_phantomeHmm,value)
            match_phageEnzymeHmm      = re.search(p_phageEnzymeHmm,value)
            match_keggVirusHmm        = re.search(p_keggVirusHmm,value)
            match_pfamHmm             = re.search(p_pfamHmm,value)
            match_smartHmm            = re.search(p_smartHmm,value)
            match_swissprotHmm        = re.search(p_swissprotHmm,value)
            match_uniprotHmm          = re.search(p_uniprotHmm,value)
            match_nrHmm               = re.search(p_nrHmm,value)
            match_customHmm           = re.search(p_customHmm,value)
            if match_ncbiVirusGenomeHmm:
                NCBI_VIRUS_GENOME_HMM = True
            if match_ncbiVirusProteinHmm:
                NCBI_VIRUS_PROTEIN_HMM = True
            if match_refseqProteinHmm:
                REFSEQ_PROTEIN_HMM = True
            if match_refseqGeneHmm:
                REFSEQ_GENE_HMM = True
            if match_pvogsHmm:
                PVOGS_HMM = True 
            if match_vogsHmm:
                VOGS_HMM = True 
            if match_phantomeHmm:
                PHANTOME_HMM = True
            if match_phageEnzymeHmm:
                PHAGE_ENZYME_HMM = True
            if match_keggVirusHmm:
                KEGG_VIRUS_HMM = True
            if match_pfamHmm:
                PFAM_HMM = True
            if match_smartHmm:
                SMART_HMM = True
            if match_swissprotHmm:
                SWISSPROT_HMM = True
            if match_uniprotHmm:
                UNIPROT_HMM = True
            if match_nrHmm:
                NR_HMM = True
            if match_customHmm:
                CUSTOM_HMM = True

    if match_profileProgramStringParam: # hmm programs for searching hmm profile databases
        if i < argCount:
            value = argList[i+1]
            match_hmmscan = re.search(p_hmmscan,value)
            if match_hmmscan:
                HMMSCAN = True       # search hmm profile dB(s) with hmm program

    #if match_customDatabasePathParam:  # custom database paths and names
    #    if i < argCount:
    #        value = argList[i+1]
    #        paramSet = value.split('^')
    #        for parameter in paramSet:
    #            match_customGenomeDBpath  = re.search(p_customGenomeDBpath,parameter)
    #            match_customGeneDBpath    = re.search(p_customGeneDBpath,parameter)
    #            match_customProteinDBpath = re.search(p_customProteinDBpath,parameter)
    #            match_customHmmDBpath     = re.search(p_customHmmDBpath,parameter)
    #            if match_customGenomeDBpath:
    #                (tag, customGenomeDBpath, customGenomeDBname) = parameter.split(':') 
    #            elif match_customGeneDBpath:
    #                (tag, customGeneDBpath, customGeneDBname) = parameter.split(':')
    #            elif match_customProteinDBpath:
    #                (tag, customProteinDBpath,customProteinDBname) = parameter.split(':') 
    #            elif match_customHmmDBpath:
    #                (tag, customHmmDBpath, customHmmDBname) = parameter.split(':') 

# Blast or Hmm search of blast or sequence databases
# Set local booleans: if any one database has been selected, then user intends to process
blastOutputDir = ''
if NCBI_VIRUS_GENOME_BLAST or CUSTOM_GENOME_BLAST:
    RUN_BLAST         = True
    RUN_HMM_SEARCH    = True
    RUN_GENOME_BLAST  = True
if NCBI_VIRUS_PROTEIN_BLAST:
    RUN_BLAST         = True
    RUN_HMM_SEARCH    = True
    RUN_PROTEIN_BLAST = True
if REFSEQ_GENE_BLAST or VOG_GENE_BLAST or CUSTOM_GENE_BLAST:
    RUN_BLAST         = True
    RUN_HMM_SEARCH    = True
    RUN_GENE_BLAST    = True
if NCBI_VIRUS_PROTEIN_BLAST or REFSEQ_PROTEIN_BLAST or NR_BLAST or CAZY_BLAST:
    RUN_BLAST         = True
    RUN_HMM_SEARCH    = True
    RUN_PROTEIN_BLAST = True
if PVOGS_BLAST or VOGS_BLAST or VOG_PROTEIN_BLAST or PHANTOME_BLAST or PHAGE_ENZYME_BLAST or KEGG_VIRUS_BLAST:
    RUN_BLAST         = True
    RUN_HMM_SEARCH    = True
    RUN_PROTEIN_BLAST = True
if PFAM_BLAST or SMART_BLAST or SWISSPROT_BLAST or UNIPROT_BLAST or CUSTOM_PROTEIN_BLAST:
    RUN_BLAST         = True
    RUN_HMM_SEARCH    = True
    RUN_PROTEIN_BLAST = True

# If user did not also select a code for searching a blast or sequence DB, then turn off.
if not BLASTP_SEARCH:
    RUN_PROTEIN_BLAST = False
if not PHMMER_SEARCH and not JACKHMMER_SEARCH:
    RUN_HMM_SEARCH    = False

# HMM profile search
# Set local booleans: if any one database has been selected, then user intends to process
if NCBI_VIRUS_GENOME_HMM: 
    RUN_PROFILE_SEARCH = True 
if REFSEQ_GENE_HMM:
    RUN_PROFILE_SEARCH = True 
if REFSEQ_PROTEIN_HMM or NR_HMM:
    RUN_PROFILE_SEARCH = True
if PVOGS_HMM or VOGS_HMM or PHANTOME_HMM or PHAGE_ENZYME_HMM or KEGG_VIRUS_HMM:
    RUN_PROFILE_SEARCH = True
if PFAM_HMM or SMART_HMM or SWISSPROT_HMM or UNIPROT_HMM or CUSTOM_HMM:
    RUN_PROFILE_SEARCH = True
# If user did not select hmmscan, then turn off, regardless of whether profile DBs were selected.
if not HMMSCAN:  # Only hmmscan, at least for now.
    RUN_PROFILE_SEARCH = False   

# Open and Check files

fileError = False

try:
    LOGFILE_H = open(logfile,"w")
except IOError as e:
    fileError = True
    if PHATE_WARNINGS == 'True':
        print(e, "logfile,", logfile) 
        print("phate_sequenceAnnotation_main says, ERROR: Cannot write log file; Exit computation at phate_sequenceAnnotation_main.py, attpempting to open log file")
    exit(0)
LOGFILE_H.write("%s%s\n" % ("Processing begun ",datetime.datetime.now()))
LOGFILE_H.write("%s%s\n" % ("sys.argv is ",sys.argv))

try:
    GENOME_FILE = open(infile_genome,"r")
except IOError as e:
    fileError = True
    if PHATE_WARNINGS == 'True':
        print(e, "genome file,", infile_genome)
    LOGFILE_H.write("%s%s%s%s\n" % ("fileError ",e," genome file, ", infile_genome))

try:
    LOGFILE_H.write("%s%s\n" % ("Opening primary calls file, ",infile_primaryCalls))
    PRIMARY_CALLS_FILE_H = open(infile_primaryCalls,"r")
except IOError as e:
    fileError = True
    if PHATE_WARNINGS == 'True':
        print(e, "primary calls file,", infile_primaryCalls)
    LOGFILE_H.write("%s%s%s%s\n" % ("fileError ",e," primary calls file, ", infile_primaryCalls))

try:
    OUTFILE = open(outfile,"w")
except IOError as e:
    fileError = True
    if PHATE_WARNINGS == 'True':
        print(e, "outfile,", outfile)
    LOGFILE_H.write("%s%s%s%s\n" % ("fileError ",e," outfile, ", outfile))

try:
    GFFFILE = open(gfffile,"w")
except IOError as e:
    fileError = True
    if PHATE_WARNINGS == 'True':
        print(e, "gfffile,", gfffile)
    LOGFILE_H.write("%s%s%s%s\n" % ("fileError ",e," gfffile, ", gfffile))

if fileError:
    print("phate_sequenceAnnotation_main says, ERROR: Check the formats of your input file(s):")
    print("  genome file is", infile_genome)
    print("  primary gene call file is", infile_primaryCalls)
    print("  outfile is", outfile)
    print("  gfffile is", gfffile)
    LOGFILE_H.write("%s%s\n" % ("Terminating due to file error at ",datetime.datetime.now()))
    LOGFILE_H.close(); exit(0)

# Communicate to log
LOGFILE_H.write("%s%s\n" % ("Parameters recorded at ",datetime.datetime.now()))
LOGFILE_H.write("%s%s\n" % ("outputDir is", outputDir))
LOGFILE_H.write("%s%s\n" % ("outfile is ",outfile))
LOGFILE_H.write("%s%s\n" % ("gfffile is ",gfffile))
LOGFILE_H.write("%s%s\n" % ("infile_genome is ", infile_genome))
LOGFILE_H.write("%s%s\n" % ("geneticCode is ",geneticCode))
LOGFILE_H.write("%s%s\n" % ("primaryCalls is ",primaryCalls))
LOGFILE_H.write("%s%s\n" % ("infile_primaryCalls is ", infile_primaryCalls))
LOGFILE_H.write("%s%s\n" % ("infile_gene is ", infile_gene))
LOGFILE_H.write("%s%s\n" % ("infile_protein is ", infile_protein))
LOGFILE_H.write("%s%s\n" % ("genomeType is ",genomeType))
LOGFILE_H.write("%s%s\n" % ("genomeName is ",genomeName))
LOGFILE_H.write("%s%s\n" % ("genomeSpecies is ",genomeSpecies))
LOGFILE_H.write("%s%s\n" % ("blastpIdentity is ",blastpIdentity))
LOGFILE_H.write("%s%s\n" % ("blastnIdentity is ",blastnIdentity))
LOGFILE_H.write("%s%s\n" % ("blastpHitCount is ",blastpHitCount))
LOGFILE_H.write("%s%s\n" % ("blastnHitCount is ",blastnHitCount))
LOGFILE_H.write("%s%s\n" % ("blastThreads is ",str(blastThreads)))
if TRANSLATE_ONLY:
    LOGFILE_H.write("%s\n" % ("Translating only; no annotation."))
else:
    LOGFILE_H.write("%s\n" % ("Annotating"))
# Blast search against blast databases
LOGFILE_H.write("%s%s\n" % ("RUN_BLAST is ",RUN_BLAST))
LOGFILE_H.write("%s%s\n" % ("RUN_GENOME_BLAST is ",RUN_GENOME_BLAST))
LOGFILE_H.write("%s%s\n" % ("RUN_GENE_BLAST is ",RUN_GENE_BLAST))
LOGFILE_H.write("%s%s\n" % ("RUN_PROTEIN_BLAST is ",RUN_PROTEIN_BLAST))
LOGFILE_H.write("%s%s\n" % ("BLASTP_SEARCH is ",BLASTP_SEARCH))
# Hmm search against sequence databases
LOGFILE_H.write("%s%s\n" % ("RUN_HMM_SEARCH is ",RUN_HMM_SEARCH))
LOGFILE_H.write("%s%s\n" % ("PHMMER_SEARCH is ",PHMMER_SEARCH))
LOGFILE_H.write("%s%s\n" % ("JACKHMMER_SEARCH is ",JACKHMMER_SEARCH))
# Sequence and blast databases
LOGFILE_H.write("%s%s\n" % ("NCBI_VIRUS_GENOME_BLAST is ",NCBI_VIRUS_GENOME_BLAST))
LOGFILE_H.write("%s%s\n" % ("NCBI_VIRUS_PROTEIN_BLAST is ",NCBI_VIRUS_PROTEIN_BLAST))
LOGFILE_H.write("%s%s\n" % ("REFSEQ_PROTEIN_BLAST is ",REFSEQ_PROTEIN_BLAST))
LOGFILE_H.write("%s%s\n" % ("REFSEQ_GENE_BLAST is ",REFSEQ_GENE_BLAST))
LOGFILE_H.write("%s%s\n" % ("PVOGS_BLAST is ",PVOGS_BLAST))
LOGFILE_H.write("%s%s\n" % ("VOGS_BLAST is ",VOGS_BLAST))    #*** To be deprecated
LOGFILE_H.write("%s%s\n" % ("VOG_GENE_BLAST is ",VOG_GENE_BLAST))
LOGFILE_H.write("%s%s\n" % ("VOG_PROTEIN_BLAST is ",VOG_PROTEIN_BLAST))
LOGFILE_H.write("%s%s\n" % ("PHANTOME_BLAST is ",PHANTOME_BLAST))
LOGFILE_H.write("%s%s\n" % ("PHAGE_ENZYME_BLAST is ",PHAGE_ENZYME_BLAST))
LOGFILE_H.write("%s%s\n" % ("KEGG_VIRUS_BLAST is ",KEGG_VIRUS_BLAST))
LOGFILE_H.write("%s%s\n" % ("PFAM_BLAST is ",PFAM_BLAST))
LOGFILE_H.write("%s%s\n" % ("SMART_BLAST is ",SMART_BLAST))
LOGFILE_H.write("%s%s\n" % ("SWISSPROT_BLAST is ",SWISSPROT_BLAST))
LOGFILE_H.write("%s%s\n" % ("UNIPROT_BLAST is ",UNIPROT_BLAST))
LOGFILE_H.write("%s%s\n" % ("NR_BLAST is ",NR_BLAST))
LOGFILE_H.write("%s%s\n" % ("CAZY_BLAST is ",CAZY_BLAST))
LOGFILE_H.write("%s%s\n" % ("CUSTOM_GENOME_BLAST is ",CUSTOM_GENOME_BLAST))
LOGFILE_H.write("%s%s\n" % ("CUSTOM_GENE_BLAST is ",CUSTOM_GENE_BLAST))
LOGFILE_H.write("%s%s\n" % ("CUSTOM_PROTEIN_BLAST is ",CUSTOM_PROTEIN_BLAST))
# Hmm search against hmm profile databases
LOGFILE_H.write("%s%s\n" % ("RUN_PROFILE_SEARCH is ",RUN_PROFILE_SEARCH))
LOGFILE_H.write("%s%s\n" % ("HMMSCAN is ",HMMSCAN))
# Hmm profile databases
LOGFILE_H.write("%s%s\n" % ("NCBI_VIRUS_GENOME_HMM is ",NCBI_VIRUS_GENOME_HMM))
LOGFILE_H.write("%s%s\n" % ("NCBI_VIRUS_PROTEIN_HMM is ",NCBI_VIRUS_PROTEIN_HMM))
LOGFILE_H.write("%s%s\n" % ("REFSEQ_PROTEIN_HMM is ",REFSEQ_PROTEIN_HMM))
LOGFILE_H.write("%s%s\n" % ("REFSEQ_GENE_HMM is ",REFSEQ_GENE_HMM))
LOGFILE_H.write("%s%s\n" % ("PVOGS_HMM is ",PVOGS_HMM))
LOGFILE_H.write("%s%s\n" % ("VOGS_HMM is ",VOGS_HMM))
LOGFILE_H.write("%s%s\n" % ("PHANTOME_HMM is ",PHANTOME_HMM))
LOGFILE_H.write("%s%s\n" % ("PHAGE_ENZYME_HMM is ",PHAGE_ENZYME_HMM))
LOGFILE_H.write("%s%s\n" % ("KEGG_VIRUS_HMM is ",KEGG_VIRUS_HMM))
LOGFILE_H.write("%s%s\n" % ("PFAM_HMM is ",PFAM_HMM))
LOGFILE_H.write("%s%s\n" % ("SMART_HMM is ",SMART_HMM))
LOGFILE_H.write("%s%s\n" % ("SWISSPROT_HMM is ",SWISSPROT_HMM))
LOGFILE_H.write("%s%s\n" % ("UNIPROT_HMM is ",UNIPROT_HMM))
LOGFILE_H.write("%s%s\n" % ("NR_HMM is ",NR_HMM))
#LOGFILE_H.write("%s%s\n" % ("CUSTOM_GENE_HMM is ",CUSTOM_GENE_HMM))
#LOGFILE_H.write("%s%s\n" % ("CUSTOM_PROTEIN_HMM is ",CUSTOM_PROTEIN_HMM))

# Communicate to user
if PHATE_MESSAGES == 'True':
    print("phate_sequenceAnnotation_main says, Parameters are: ")
    print("  outputDir is", outputDir)
    print("  outfile is", outfile)
    print("  gfffile is", gfffile)
    print("  genome file is", infile_genome)
    print("  genetic code is", geneticCode)
    print("  primary calls is", primaryCalls)
    print("  primary gene call file is", infile_primaryCalls)
    print("  gene file is", infile_gene)
    print("  protein file is", infile_protein)
    print("  genome type is", genomeType)
    print("  genome name is", genomeName)
    #print("  contigName is", contigName)
    print("  genomeSpecies is", genomeSpecies)
    print("  blastp identity is", blastpIdentity)
    print("  blastn identity is", blastnIdentity)
    print("  blastp hit count is", blastpHitCount)
    print("  blastn hit count is", blastnHitCount)
    print("  blastThreads is", blastThreads)
    if TRANSLATE_ONLY:
        print("  Translating only; no annotation.")
    else:
        print("  Annotating")
    print("  RUN_BLAST is", RUN_BLAST)
    print("  RUN_GENOME_BLAST is", RUN_GENOME_BLAST)
    print("  RUN_PROTEIN_BLAST is", RUN_PROTEIN_BLAST)
    print("  BLASTP_SEARCH is", BLASTP_SEARCH)
    print("  RUN_HMM_SEARCH is", RUN_HMM_SEARCH)
    print("  PHMMER_SEARCH is", PHMMER_SEARCH)
    print("  JACKHMMER_SEARCH is", JACKHMMER_SEARCH)
    print("  NCBI_VIRUS_GENOME_BLAST is", NCBI_VIRUS_GENOME_BLAST)
    print("  NCBI_VIRUS_PROTEIN_BLAST is", NCBI_VIRUS_PROTEIN_BLAST)
    print("  REFSEQ_PROTEIN_BLAST is", REFSEQ_PROTEIN_BLAST)
    print("  REFSEQ_GENE_BLAST is", REFSEQ_GENE_BLAST)
    print("  PVOGS_BLAST is", PVOGS_BLAST)
    print("  VOGS_BLAST is", VOGS_BLAST)   #*** To be deprecated
    print("  VOG_GENE_BLAST is", VOG_GENE_BLAST)
    print("  VOG_PROTEIN_BLAST is", VOG_PROTEIN_BLAST)
    print("  PHANTOME_BLAST is", PHANTOME_BLAST)
    print("  PHAGE_ENZYME_BLAST is", PHAGE_ENZYME_BLAST)
    print("  KEGG_VIRUS_BLAST is", KEGG_VIRUS_BLAST)
    print("  PFAM_BLAST is", PFAM_BLAST)
    print("  SMART_BLAST is", SMART_BLAST)
    print("  SWISSRPOT_BLAST is", SWISSPROT_BLAST)
    print("  UNIRPOT_BLAST is", UNIPROT_BLAST)
    print("  NR_BLAST is", NR_BLAST)
    print("  CAZY_BLAST is", CAZY_BLAST)
    print("  CUSTOM_GENOME_BLAST is", CUSTOM_GENOME_BLAST)
    print("  CUSTOM_GENE_BLAST is", CUSTOM_GENE_BLAST)
    print("  CUSTOM_PROTEIN_BLAST is", CUSTOM_PROTEIN_BLAST)
    print("  RUN_PROFILE_SEARCH is", RUN_PROFILE_SEARCH)
    print("  HMMSCAN is", HMMSCAN)
    print("  NCBI_VIRUS_GENOME_HMM is", NCBI_VIRUS_GENOME_HMM)
    print("  NCBI_VIRUS_PROTEIN_HMM is", NCBI_VIRUS_PROTEIN_HMM)
    print("  REFSEQ_PROTEIN_HMM is", REFSEQ_PROTEIN_HMM)
    print("  REFSEQ_GENE_HMM is", REFSEQ_GENE_HMM)
    print("  PVOGS_HMM is", PVOGS_HMM)
    print("  VOGS_HMM is", VOGS_HMM)
    print("  PHANTOME_HMM is", PHANTOME_HMM)
    print("  PHAGE_ENZYME_HMM is", PHAGE_ENZYME_HMM)
    print("  KEGG_VIRUS_HMM is", KEGG_VIRUS_HMM)
    print("  PFAM_HMM is", PFAM_HMM)
    print("  SMART_HMM is", SMART_HMM)
    print("  SWISSPROT_HMM is", SWISSPROT_HMM)
    print("  UNIPROT_HMM is", UNIPROT_HMM)
    print("  NR_HMM is", NR_HMM)
    #print("  CUSTOM_GENE_HMM is", CUSTOM_GENE_HMM)
    #print("  CUSTOM_PROTEIN_HMM is", CUSTOM_PROTEIN_HMM)

if PHATE_PROGRESS == 'True':
    print("phate_sequenceAnnotation_main says, Configuration complete for sequence annotation module.")

##### BEGIN MAIN ################################################################################3

geneCallInfo = {      # For passing info to genomeSequence module  #*** ???
    'primaryCalls'         : primaryCalls,
    'primaryCallsPathFile' : infile_primaryCalls,
    'genomeName'           : genomeName,
}

# Create a genome object and set parameters 

if PHATE_PROGRESS == 'True':
    print("phate_sequenceAnnotation_main says, Preparing for sequence annotation: setting parameters...")
LOGFILE_H.write("%s%s\n" % ("Setting parameters for genome at ",datetime.datetime.now()))
myGenome = phate_genomeSequence.genome()
myGenome.genomeType    = genomeType
myGenome.genomeName    = genomeName
myGenome.genomeSpecies = genomeSpecies 
myGenome.setCodeBaseDir(CODE_BASE_DIR)
myGenome.setOutputDir(outputDir)
LOGFILE_H.write("%s\n" % ("Reading sequence into genome object"))
gLines = GENOME_FILE.read().splitlines()
myGenome.contigSet.addFastas(gLines,'nt')
if PHATE_MESSAGES == 'True':
    print("phate_sequenceAnnotation_main says, contigName is", myGenome.contigSet.contig)

# Create a hash to record contig names and their sequence lengths, for GFF output
contigSeqLen_hash = {}
for i in range(0,len(myGenome.contigSet.fastaList)):
    contigName = myGenome.contigSet.fastaList[i].cleanHeader
    seqLen     = len(myGenome.contigSet.fastaList[i].sequence)
    contigSeqLen_hash[contigName] = seqLen

# Extract gene calls
if PHATE_PROGRESS == 'True':
    print("phate_sequenceAnnotation_main says, Processing gene calls...")
LOGFILE_H.write("%s%s\n" % ("Processing gene calls at ",datetime.datetime.now()))
# Record the gene calls, and pass the hash so contig lengths can also be recorded...
# ...for access when ultimately writing the GFF output file. (oh what we do for want of a pointer.)
if PHATE_PROGRESS == 'True':
    print("phate_sequenceAnnotation_main says, geneCallInfo",geneCallInfo,"infile_primaryCalls",infile_primaryCalls,"contigSeqLen_hash",contigSeqLen_hash)
myGenome.processGeneCalls(geneCallInfo,PRIMARY_CALLS_FILE_H)
myGenome.cleanUpAfterEMBOSS()
PRIMARY_CALLS_FILE_H.close()

# Output the gene and protein sets, if newly created

fastaOut = {
    "mtype"      : "",
    "headerType" : "",
    "filename"   : "",
}

if PHATE_PROGRESS == 'True':
    print("phate_sequenceAnnotation_main says, Writing genes file")
LOGFILE_H.write("%s\n" % ("Writing genes file"))
# Print out newly created gene list
fastaOut["mtype"] = "gene"
fastaOut["headerType"] = "full"  #*** Should this be "compound" ???
fastaOut["filename"] = infile_gene 
myGenome.printFastas2file(fastaOut)

if PHATE_PROGRESS == 'True':
    print("phate_sequenceAnnotation_main says, Writing peptides file")
LOGFILE_H.write("%s\n" % ("Writing peptides file"))
# Print out newly created protein list 
fastaOut["mtype"] = "protein"   #*** Should this be "compound" ???
fastaOut["headerType"] = "full"
fastaOut["filename"] = infile_protein 
myGenome.printFastas2file(fastaOut)

if PHATE_PROGRESS == 'True':
    print("phate_sequenceAnnotation_main says, Gene and protein files created.")
LOGFILE_H.write("%s%s\n" % ("Gene and protein files created at ",datetime.datetime.now()))

# If user specified to translate only, then skip this segment of the pipeline.
if TRANSLATE_ONLY:
    if PHATE_PROGRESS == 'True':
        print("phate_sequenceAnnotation_main says, Translate only: computations completed.")
    LOGFILE_H.write("%s%s\n" % ("Translating only: computations completed at ",datetime.datetime.now()))
else:
    ######################################################## ANNOTATE ####################################################
    # Create a blast object and set parameters
    if RUN_BLAST:
        if PHATE_PROGRESS == 'True':
            print("phate_sequenceAnnotation_main says, Preparing for blast...")
        LOGFILE_H.write("%s\n" % ("Creating a blast object"))
        blast = phate_blast.multiBlast()

        if PHATE_PROGRESS == 'True':
            print("phate_sequenceAnnotation_main says, Preparing to run", blast.blastFlavor)
        if PHATE_MESSAGES == 'True':
            print("phate_sequenceAnnotation_main says, Running at the following settings:")
            blast.printParameters()
        LOGFILE_H.write("%s%s%s%s\n" % (datetime.datetime.now(), " Preparing to run ", blast.blastFlavor, " at the following settings:"))
        blast.printParameters2file(LOGFILE_H)

        ### Create directory for blast output
        blastOutputDir = outputDir + 'BLAST/'
        try:
            os.stat(blastOutputDir)
        except:
            os.mkdir(blastOutputDir)

    if RUN_HMM_SEARCH:

        ### Create an hmm object and set parameters
        if PHATE_PROGRESS == 'True':
            print("phate_sequenceAnnotation_main says, Preparing for hmm...")
        LOGFILE_H.write("%s\n" % ("Creating an hmm object"))
        hmm = phate_hmm.multiHMM()

        # Create directory for hmm output
        hmmOutputDir = outputDir + 'HMM/'
        try:
            os.stat(hmmOutputDir)
        except:
            os.mkdir(hmmOutputDir)

    if RUN_PROFILE_SEARCH:

        ### Create profile search object
        if PHATE_PROGRESS == 'True':
            print("phate_sequenceAnnotation_main says, Preparing for profile search...")
        LOGFILE_H.write("%s\n" % ("Creating a profile object"))
        profile = phate_profile.multiProfile()

        # Create directories for hmm/profile search output
        profileOutputDir = outputDir + 'PROFILE/'
        try:
            os.stat(profileOutputDir)
        except:
            os.mkdir(profileOutputDir)
        proteinProfileOutputDir = profileOutputDir + 'Protein/'
        try:
            os.stat(proteinProfileOutputDir)
        except:
            os.mkdir(proteinProfileOutputDir)

    ##### Run blast against sequences/databases we have:  genome, gene, protein, virus, phage databases

    ##### GENOME BLAST

    LOGFILE_H.write("%s%s\n" % ("Checking RUN_GENOME_BLAST ",RUN_GENOME_BLAST))
    if RUN_GENOME_BLAST:

        # Create blast output directory for genome blast
        genomeBlastOutputDir = blastOutputDir + 'Genome/'
        try:
            os.stat(genomeBlastOutputDir)
        except:
            os.mkdir(genomeBlastOutputDir)

        # Prepare for genome blast
        myParamSet = {
            'identityMin'          : int(blastnIdentity),  
            'identitySelect'       : int(blastnIdentity),   
            'evalueMin'            : EVALUE_MIN,
            'evalueSelect'         : EVALUE_SELECT,
            'topHitCount'          : int(blastnHitCount),
            'outputFormat'         : XML_OUT_FORMAT,  # blast output format (use 5=XML or 7=LIST only) 
            'scoreEdge'            : SCORE_EDGE,
            'overhang'             : OVERHANG,
            'geneCallDir'          : outputDir, 
            'blastOutDir'          : genomeBlastOutputDir,
            'ncbiVirusGenomeBlast' : NCBI_VIRUS_GENOME_BLAST,
            'customGenomeBlast'    : CUSTOM_GENOME_BLAST,
        }
        blast.setBlastParameters(myParamSet)
        blast.setBlastFlavor('blastn') 

        if PHATE_PROGRESS == 'True':
            print("phate_sequenceAnnotation_main says, Preparing to run", blast.blastFlavor)
        if PHATE_MESSAGES == 'True':
            print("phate_sequenceAnnotation_main says, Running", blast.blastFlavor, "at the following settings:")
            blast.printParameters()
        if PHATE_PROGRESS == 'True':
            print("phate_sequenceAnnotation_main says, Running Blast against phage genome database(s)...")

        LOGFILE_H.write("%s%s%s%s\n" % (datetime.datetime.now(), " Preparing to run genome ", blast.blastFlavor, " at the following settings:"))
        blast.printParameters2file(LOGFILE_H)
        LOGFILE_H.write("%s\n" % ("Running Blast against phage genome database(s)..."))
        # Run Genome blast 
        blast.runBlast(myGenome.contigSet,'genome')
        if PHATE_PROGRESS == 'True':
            print("phate_sequenceAnnotation_main says, Genome blast complete")
        LOGFILE_H.write("%s%s\n" % ("Genome blast complete at ", datetime.datetime.now()))

    else:
        if PHATE_PROGRESS == 'True':
            print("phate_sequenceAnnotation_main says, Skipping genome blast")

    ##### GENE BLAST

    if RUN_GENE_BLAST:

        LOGFILE_H.write("%s\n" % ("Preparing to run gene blast"))
        if PHATE_PROGRESS == 'True':
            print("phate_sequenceAnnotation_main says, Preparing to run gene blast...")

        # Create blast output directory for gene blast
        geneBlastOutputDir = blastOutputDir + 'Gene/'
        try:
            os.stat(geneBlastOutputDir)
        except:
            os.mkdir(geneBlastOutputDir)

        # Prepare for gene blast
        myParamSet = {
            'identityMin'         : int(blastnIdentity),   
            'identitySelect'      : int(blastnIdentity),  
            'evalueMin'           : EVALUE_MIN,
            'evalueSelect'        : EVALUE_SELECT,
            'topHitCount'         : int(blastnHitCount),  #*** maybe should parameterize this
            'outputFormat'        : XML_OUT_FORMAT,  # XML=5; LIST=7
            'scoreEdge'           : 0.1,
            'overhang'            : 0.1,
            'geneCallDir'         : outputDir, 
            'blastOutDir'         : geneBlastOutputDir,
            'refseqGeneBlast'     : REFSEQ_GENE_BLAST,   #*** To ge deprecated
            'vogGeneBlast'        : VOG_GENE_BLAST,
            'pvogsOutDir'         : geneBlastOutputDir,
            'vogsOutDir'          : geneBlastOutputDir,
            'customGeneBlast'     : CUSTOM_GENE_BLAST,
            'blastThreads'        : blastThreads,
        }
        blast.setBlastParameters(myParamSet)
        blast.setBlastFlavor('blastn') 

        if PHATE_PROGRESS == 'True':
            print("phate_sequenceAnnotation_main says, Running Blast against gene database(s)...")
        LOGFILE_H.write("%s\n" % ("Running Blast against gene databases at the following settings:"))
        blast.printParameters2file(LOGFILE_H)

        # Run Gene blast
        blast.runBlast(myGenome.geneSet,'gene')

        if PHATE_PROGRESS == 'True':
            print("phate_sequenceAnnotation_main says, Gene blast complete.")
        LOGFILE_H.write("%s%s\n" % ("Gene blast complete at ",datetime.datetime.now()))

    else:
        if PHATE_PROGRESS == 'True':
            print("phate_sequenceAnnotation_main says, Skipping gene blast.")

    ##### PROTEIN BLAST

    if RUN_PROTEIN_BLAST:
   
        # Create blast output directory for protein blast
        proteinBlastOutputDir = blastOutputDir + 'Protein/'
        try:
            os.stat(proteinBlastOutputDir)
        except:
            os.mkdir(proteinBlastOutputDir)

        if PHATE_PROGRESS == 'True':
            print("phate_sequenceAnnotation_main says, Preparing for protein blast...")
        LOGFILE_H.write("%s%s\n" % ("Preparing for protein blast at ",datetime.datetime.now()))

        # Prepare for protein blast
        myParamSet = {
            'blastpSearch'          : BLASTP_SEARCH,
            'identityMin'           : int(blastpIdentity),  
            'identitySelect'        : int(blastpIdentity),  
            'evalueMin'             : EVALUE_MIN,
            'evalueSelect'          : EVALUE_SELECT,
            'topHitCount'           : int(blastpHitCount),  #*** maybe should parameterize this
            'outputFormat'          : XML_OUT_FORMAT,  # XML=5, LIST=7  
            'scoreEdge'             : SCORE_EDGE,
            'overhang'              : OVERHANG,
            'geneCallDir'           : outputDir, 
            'blastOutDir'           : proteinBlastOutputDir,
            'pvogsOutDir'           : proteinBlastOutputDir,
            'vogsOutDir'            : proteinBlastOutputDir,
            'ncbiVirusProteinBlast' : NCBI_VIRUS_PROTEIN_BLAST,
            'refseqProteinBlast'    : REFSEQ_PROTEIN_BLAST,
            'pvogsBlast'            : PVOGS_BLAST,
            'vogsBlast'             : VOGS_BLAST,   #*** To be deprecated
            'vogProteinBlast'       : VOG_PROTEIN_BLAST,
            'phantomeBlast'         : PHANTOME_BLAST,
            'phageEnzymeBlast'      : PHAGE_ENZYME_BLAST,
            'keggVirusBlast'        : KEGG_VIRUS_BLAST,
            'pfamBlast'             : PFAM_BLAST,
            'smartBlast'            : SMART_BLAST,
            'swissprotBlast'        : SWISSPROT_BLAST,
            'uniprotBlast'          : UNIPROT_BLAST,
            'nrBlast'               : NR_BLAST,
            'cazyBlast'             : CAZY_BLAST,
            'customProteinBlast'    : CUSTOM_PROTEIN_BLAST,
            'blastThreads'          : blastThreads,
        }
        blast.setBlastParameters(myParamSet)
        blast.setBlastFlavor('blastp')

        if PHATE_PROGRESS == 'True': 
            print("phate_sequenceAnnotation_main says, Running blastp against protein database(s)...")
        LOGFILE_H.write("%s\n" % ("Running Blast against protein database(s) at the following settings:"))
        blast.printParameters2file(LOGFILE_H)

        # Run protein blast
        blast.runBlast(myGenome.proteinSet,'protein')

        if PHATE_PROGRESS == 'True':
            print("phate_sequenceAnnotation_main says, Protein blast complete.")
        LOGFILE_H.write("%s%s\n" % ("Protein blast complete at ", datetime.datetime.now()))

    else:
        if PHATE_PROGRESS == 'True':
            print("phate_sequenceAnnotation_main says, Skipping protein blast")

    #######################################################################################
    # HMM SEARCH:  Use hmm program to search fasta blast-formatted database(s)

    if RUN_HMM_SEARCH:

        if PHATE_PROGRESS == 'True':
            print("phate_sequenceAnnotation_main says, Preparing for hmm search...")
        LOGFILE_H.write("%s\n" % ("Preparing for hmm search"))

        # Create hmm output directory for genome hmm search
        genomeHmmOutputDir = hmmOutputDir + 'Genome/'
        try:
            os.stat(genomeHmmOutputDir)
        except:
            os.mkdir(genomeHmmOutputDir)
        # Done for now; genome hmm search not yet in service

        # Create hmm output directory for gene hmm search
        geneHmmOutputDir = hmmOutputDir + 'Gene/'
        try:
            os.stat(geneHmmOutputDir)
        except:
            os.mkdir(geneHmmOutputDir)
        # Done for now; gene hmm search not yet in service

        # Create hmm output directory for protein hmm search
        proteinHmmOutputDir = hmmOutputDir + 'Protein/'
        try:
            os.stat(proteinHmmOutputDir)
        except:
            os.mkdir(proteinHmmOutputDir)

        # Prepare for protein hmm search
        # HMM search uses jackhmmer or phmmer against fasta blast databases
        myParamSet = {
            'hmmProgram'            : '', 
            'phmmerSearch'          : PHMMER_SEARCH,
            'jackhmmerSearch'       : JACKHMMER_SEARCH,
            'geneCallDir'           : outputDir,
            'genomeHmmOutDir'       : genomeHmmOutputDir,
            'geneHmmOutDir'         : geneHmmOutputDir,
            'proteinHmmOutDir'      : proteinHmmOutputDir,
            'pVOGsOutDir'           : proteinHmmOutputDir,
            'VOGsOutDir'            : proteinHmmOutputDir,
            'ncbiVirusGenomeBlast'  : NCBI_VIRUS_GENOME_BLAST,
            'ncbiVirusProteinBlast' : NCBI_VIRUS_PROTEIN_BLAST,
            'refseqProteinBlast'    : REFSEQ_PROTEIN_BLAST,
            'refseqGeneBlast'       : REFSEQ_GENE_BLAST,
            'pvogsBlast'            : PVOGS_BLAST,
            'vogsBlast'             : VOGS_BLAST,
            'vogGeneBlast'          : VOG_GENE_BLAST,
            'vogProteinBlast'       : VOG_PROTEIN_BLAST,
            'phantomeBlast'         : PHANTOME_BLAST,
            'phageEnzymeBlast'      : PHAGE_ENZYME_BLAST,
            'keggVirusBlast'        : KEGG_VIRUS_BLAST,
            'pfamBlast'             : PFAM_BLAST,
            'smartBlast'            : SMART_BLAST,
            'swissprotBlast'        : SWISSPROT_BLAST,
            'uniprotBlast'          : UNIPROT_BLAST,
            'nrBlast'               : NR_BLAST,
            'cazyBlast'             : CAZY_BLAST,
            'customGeneBlast'       : CUSTOM_GENE_BLAST,
            'customProteinBlast'    : CUSTOM_PROTEIN_BLAST,
        }
        hmm.setHmmParameters(myParamSet)

        if PHATE_PROGRESS == 'True':
            print("phate_sequenceAnnotation_main says, Checking whether to run Hmm searches against protein database(s)...")

        # Run jackhmmer 
        if JACKHMMER_SEARCH:
            if PHATE_PROGRESS == 'True':
                print("phate_sequenceAnnotation_main says, Running jackhmmer against protein database(s)")
            LOGFILE_H.write("%s%s\n" % ("Running jackhmmer search against protein database(s) at ", datetime.datetime.now()))
            hmm.setHmmProgram('jackhmmer')
            hmm.runHmm(myGenome.proteinSet,'protein')
            LOGFILE_H.write("%s%s\n" % ("HMM jackhmmer search complete at ",datetime.datetime.now()))

        # Run phmmer
        if PHMMER_SEARCH:
            if PHATE_PROGRESS == 'True':
                print("phate_sequenceAnnotation_main says, Running hmm search with phmmer...")
            LOGFILE_H.write("%s%s\n" % ("Running phmmer search against protein database(s) at ", datetime.datetime.now()))
            hmm.setHmmProgram('phmmer')
            hmm.runHmm(myGenome.proteinSet,'protein')
            LOGFILE_H.write("%s%s\n" % ("HMM phmmer search complete at ",datetime.datetime.now()))

    else:
        if PHATE_PROGRESS == 'True':
            print("phate_sequenceAnnotation_main says, Skipping HMM sequence search")

    #################################################################################
    ##### HMM PROFILE SEARCH

    if RUN_PROFILE_SEARCH:

        if PHATE_PROGRESS == 'True':
            print("phate_sequenceAnnotation_main says, Preparing for profile search...")
        LOGFILE_H.write("%s\n" % ("Preparing for profile search"))

        # Create profile output directory for genome profile search
        genomeProfileOutputDir = profileOutputDir + 'Genome/'
        try:
            os.stat(genomeProfileOutputDir)
        except:
            os.mkdir(genomeProfileOutputDir)
        # Done for now; genome profile search not yet in service

        # Create profile output directory for gene profile search
        geneProfileOutputDir = profileOutputDir + 'Gene/'
        try:
            os.stat(geneProfileOutputDir)
        except:
            os.mkdir(geneProfileOutputDir)
        # Done for now; gene profile search not yet in service

        # Create profile output directory for protein profile search and pVOG groups
        proteinProfileOutputDir = profileOutputDir + 'Protein/'
        pVOGsOutputDir          = profileOutputDir + 'Protein/'
        VOGsOutputDir           = profileOutputDir + 'Protein/'
        try:
            os.stat(proteinProfileOutputDir)
        except:
            os.mkdir(proteinProfileOutputDir)

        # profile object already created (see above) 

        if HMMSCAN:
            profileProgram = 'hmmscan'
            hmmscan = True

        myParamSet = {
            'profileProgram'       : profileProgram,
            'geneCallDir'          : outputDir,
            'profileOutDir'        : profileOutputDir,
            'genomeProfileOutDir'  : genomeProfileOutputDir,
            'geneProfileOutDir'    : geneProfileOutputDir,
            'proteinProfileOutDir' : proteinProfileOutputDir,
            'pVOGsOutDir'          : pVOGsOutputDir,
            'VOGsOutDir'           : VOGsOutputDir,
            'hmmscan'              : hmmscan,
            'ncbiVirusProteinHmm'  : NCBI_VIRUS_PROTEIN_HMM,
            'refseqProteinHmm'     : REFSEQ_PROTEIN_HMM,
            'pvogsHmm'             : PVOGS_HMM,
            'vogsHmm'              : VOGS_HMM,  # Vog Hmms are protein
            'phantomeHmm'          : PHANTOME_HMM,
            'phageEnzymeHmm'       : PHAGE_ENZYME_HMM,
            'keggVirusHmm'         : KEGG_VIRUS_HMM,
            'pfamHmm'              : PFAM_HMM,
            'smartHmm'             : SMART_HMM,
            'swissprotHmm'         : SWISSPROT_HMM,
            'uniprotHmm'           : UNIPROT_HMM,
            'nrHmm'                : NR_HMM,
        }

        profile.setProfileParameters(myParamSet)

        # Run protein hmm search
        if PHATE_PROGRESS == 'True':
            print("phate_sequenceAnnotation_main says, Running hmmscan search against profile database(s)...")
        LOGFILE_H.write("%s%s\n" % ("Running hmmscan search against profile database(s) at ", datetime.datetime.now()))
        profile.runProfile(myGenome.proteinSet,'protein')
        LOGFILE_H.write("%s%s\n" % ("Hmmscan  search complete at ",datetime.datetime.now()))

    else:
        if PHATE_PROGRESS == 'True':
            print("phate_sequenceAnnotation_main says, Skipping HMM profile search.")

    ##### REPORT OUT 

    if PHATE_PROGRESS == 'True':
        print("phate_sequenceAnnotation_main says, Reporting final annotations")
    LOGFILE_H.write("%s%s\n" % ("Reporting final annotations at ",datetime.datetime.now()))
    myGenome.printGenomeData2file_tab(OUTFILE)
    myGenome.printGenomeData2file_GFF(GFFFILE)

##### CLEAN UP

if PHATE_PROGRESS == 'True':
    print("phate_sequenceAnnotation_main says, Sequence annotation complete.")
GENOME_FILE.close()
OUTFILE.close()
GFFFILE.close()

LOGFILE_H.write("%s%s\n" % ("Sequence annotation processing complete at ",datetime.datetime.now()))
LOGFILE_H.close()

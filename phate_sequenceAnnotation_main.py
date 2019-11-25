#!/usr/bin/env python

################################################################
#
# phate_sequenceAnnotation_main.py
#
# Description:  Performs blast and/or hmm search on a given input gene or protein fasta set
#    and returns the results. The databases are pre-determined, but can
#    be added to by inserting new specifications in the blast object.
#
# Programmer: CEZhou
#
# Date: 21 November 2019
# Version 1.5
#
################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.PDF FOR DETAILS.

import sys, os, re, string, copy
import time, datetime
from subprocess import call

# DEBUG control
DEBUG = True
#DEBUG = False 

# Defaults/Parameters
GENE_CALLER            = 'phanotate'   # Default; can be configured by user
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
HMM_PROGRAM            = 'jackhmmer'        # Default; can be configured by user
HMM_DATABASE           = 'pvogsHmm'         # Default; can be configured by user

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
REFSEQ_PROTEIN_BASE_DIR       = os.environ["PHATE_REFSEQ_PROTEIN_BASE_DIR"]
REFSEQ_PROTEIN_BLAST_HOME     = os.environ["PHATE_REFSEQ_PROTEIN_BLAST_HOME"]
REFSEQ_GENE_BASE_DIR          = os.environ["PHATE_REFSEQ_GENE_BASE_DIR"]
REFSEQ_GENE_BLAST_HOME        = os.environ["PHATE_REFSEQ_GENE_BLAST_HOME"]
PVOGS_BASE_DIR                = os.environ["PHATE_PVOGS_BASE_DIR"]
PVOGS_BLAST_HOME              = os.environ["PHATE_PVOGS_BLAST_HOME"]
PHANTOME_BASE_DIR             = os.environ["PHATE_PHANTOME_BASE_DIR"]
PHANTOME_BLAST_HOME           = os.environ["PHATE_PHANTOME_BLAST_HOME"]
KEGG_VIRUS_BASE_DIR           = os.environ["PHATE_KEGG_VIRUS_BASE_DIR"]
KEGG_VIRUS_BLAST_HOME         = os.environ["PHATE_KEGG_VIRUS_BLAST_HOME"]
SWISSPROT_BASE_DIR            = os.environ["PHATE_SWISSPROT_BASE_DIR"]
SWISSPROT_BLAST_HOME          = os.environ["PHATE_SWISSPROT_BLAST_HOME"]
PFAM_BASE_DIR                 = os.environ["PHATE_PFAM_BASE_DIR"]
PFAM_BLAST_HOME               = os.environ["PHATE_PFAM_BLAST_HOME"]
SMART_BASE_DIR                = os.environ["PHATE_SWISSPROT_BASE_DIR"]
SMART_BLAST_HOME              = os.environ["PHATE_SWISSPROT_BLAST_HOME"]
PHAGE_ENZYME_BASE_DIR         = os.environ["PHATE_PHAGE_ENZYME_BASE_DIR"]
PHAGE_ENZYME_BLAST_HOME       = os.environ["PHATE_PHAGE_ENZYME_BLAST_HOME"]
UNIPROT_BASE_DIR              = os.environ["PHATE_UNIPROT_BASE_DIR"]
UNIPROT_BLAST_HOME            = os.environ["PHATE_UNIPROT_BLAST_HOME"]
NR_BLAST_BASE_DIR             = os.environ["PHATE_NR_BLAST_BASE_DIR"]
NR_BLAST_HOME                 = os.environ["PHATE_NR_BLAST_HOME"]
NCBI_TAXON_DIR                = os.environ["PHATE_NCBI_TAXON_DIR"]
CUSTOM_GENOME_BLAST_HOME      = os.environ["PHATE_CUSTOM_GENOME_BLAST_HOME"]
CUSTOM_GENE_BLAST_HOME        = os.environ["PHATE_CUSTOM_GENE_BLAST_HOME"]
CUSTOM_PROTEIN_BLAST_HOME     = os.environ["PHATE_CUSTOM_PROTEIN_BLAST_HOME"]

# HMM database locations
HMM_HOME                      = os.environ["PHATE_HMM_HOME"]
REFSEQ_PROTEIN_HMM_HOME       = os.environ["PHATE_REFSEQ_PROTEIN_BLAST_HOME"]
REFSEQ_GENE_HMM_HOME          = os.environ["PHATE_REFSEQ_GENE_BLAST_HOME"]
PVOGS_HMM_HOME                = os.environ["PHATE_PVOGS_BLAST_HOME"]
PHANTOME_HMM_HOME             = os.environ["PHATE_PHANTOME_HMM_HOME"]
SWISSPROT_HMM_HOME            = os.environ["PHATE_SWISSPROT_BLAST_HOME"]
PFAM_HMM_HOME                 = os.environ["PHATE_PFAM_BLAST_HOME"]
SMART_HMM_HOME                = os.environ["PHATE_SMART_BLAST_HOME"]
PHAGE_ENZYME_HMM_HOME         = os.environ["PHATE_PHAGE_ENZYME_BLAST_HOME"]
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
import phate_blast          # uses phate_blastAnalysis to handle blast requests 
import phate_hmm            # runs hmm search against specified database(s) 

##### FILES

logfile = ""
outfile = ""  # will be constructed using user's specified output subdir   

outputDir       = ""         # user-specified output subdirectory
infile_genome   = ""         # user-provided file containing the genome sequence that was gene-called
infile_geneCall = ""         # gene-call file (output from a gene caller; PHANOTATE for now)
infile_gene     = ""         # genes to be blasted (not yet in service)
infile_protein  = ""         # user-provided file containing protein fasta sequences

##### USER-SPECIFIED META-DATA and PARAMTERS

genomeType      = GENOME_TYPE              # user-provided, typically 'phage'
genomeName      = "unknown"                # user-provided, name of current genome, e.g., 'LYP215'
genomeSpecies   = "unknown"                # user-provided, e.g., 'YpPhage_LYP215'
configFile      = "unknown"                # must be passed by calling code (phate_runPipeline.py)
geneCaller      = GENE_CALLER              # default, unless changed
blastpIdentity  = BLASTP_IDENTITY_DEFAULT  # integer percent identity cutoff
blastnIdentity  = BLASTN_IDENTITY_DEFAULT  # integer percent identity cutoff
blastpHitCount  = BLASTP_HIT_COUNT_DEFAULT # number of top hits to capture
blastnHitCount  = BLASTN_HIT_COUNT_DEFAULT # number of top blastn hits to capture 
geneticCode     = GENETIC_CODE             # default, unless changed
hmmProgram      = HMM_PROGRAM              # default, unless changed by input parameter
hmmDatabase     = HMM_DATABASE             # default, unless changed by input parameter

##### BOOLEANS

# Databases: assume turned 'off', unless input string indicates otherwise
# BLAST
NCBI_VIRUS_GENOME_BLAST  = False
NCBI_VIRUS_PROTEIN_BLAST = False
NR_BLAST                 = False
KEGG_VIRUS_BLAST         = False
REFSEQ_PROTEIN_BLAST     = False
REFSEQ_GENE_BLAST        = False
PHANTOME_BLAST           = False
PVOGS_BLAST              = False
SWISSPROT_BLAST          = False
PFAM_BLAST               = False
SMART_BLAST              = False
PHAGE_ENZYME_BLAST       = False
# HMM
NCBI_VIRUS_GENOME_HMM    = False
NCBI_VIRUS_PROTEIN_HMM   = False
NR_HMM                   = False
KEGG_VIRUS_HMM           = False
REFSEQ_PROTEIN_HMM       = False
REFSEQ_GENE_HMM          = False
PHANTOME_HMM             = False
PVOGS_HMM                = False
SWISSPROT_HMM            = False
PFAM_HMM                 = False
SMART_HMM                = False
PHAGE_ENZYME_HMM         = False

##### PATTERNS and CONTROL

# Input parameters 

# parameter tags
p_outputDirParam            = re.compile('^-o')   # outout directory (e.g., 'LYP215')
p_genomeFileParam           = re.compile('^-G')   # Genome with a capital 'G'
p_geneFileParam             = re.compile('^-g')   # gene with a lower-case 'g'
p_proteinFileParam          = re.compile('^-p')   # protein or peptide
p_geneCallerParam           = re.compile('^-c')   # gene caller
p_geneticCodeParam          = re.compile('^-e')   # genetic code
p_geneCallFileParam         = re.compile('^-f')   # gene calls file
p_contigNameParam           = re.compile('^-C')   # Contig name #*** TEMPORARY until code handles draft genomes 
p_genomeTypeParam           = re.compile('^-t')   # genome type (e.g., 'phage', 'bacterium')
p_genomeNameParam           = re.compile('^-n')   # genome name (e.g., 'my_fave_Yp_genome')
p_speciesParam              = re.compile('^-s')   # species (e.g., 'Y_pestis') 
p_blastpIdentityParam       = re.compile('^-i')   # blastp identity cutoff
p_blastnIdentityParam       = re.compile('^-j')   # blastn identity cutoff
p_blastpHitCountParam       = re.compile('^-h')   # blastp top hit count
p_blastnHitCountParam       = re.compile('^-H')   # blastn top hit count
p_translateOnlyParam        = re.compile('^-x')   # if user passes 'true' => get genes, translate, compare, then stop before annotation
p_blastDBsStringParam       = re.compile('^-d')   # string listing database(s) to blast against
p_hmmProgramDBsStringParam  = re.compile('^-m')   # string listing hmm program and database(s) to search against

# Parts of input strings naming databases to blast or hmm-search against
p_hmmProgram           = re.compile('program=([a-zA-Z\d]+)')   # name of hmm program to use in hmm search
# blast DBs
p_ncbiVirus            = re.compile('ncbiVirusGenome')     
p_ncbiVirusProtein     = re.compile('ncbiVirusProtein')    
p_nr                   = re.compile('nr')                  
p_keggVirus            = re.compile('kegg')                
p_refseqProtein        = re.compile('refseqP')             
p_refseqGene           = re.compile('refseqG')             
p_phantome             = re.compile('phantome')            
p_pvogs                = re.compile('pvogs')               
p_uniparc              = re.compile('uniparc')             
p_swissprot            = re.compile('swissprot')           
p_pfam                 = re.compile('pfam')                
p_smart                = re.compile('smart')               
p_phageEnzyme          = re.compile('phageEnzyme')               
# hmm profiles DBs
p_ncbiVirusHmm         = re.compile('ncbiVirusGenomeHmm')  
p_ncbiVirusProteinHmm  = re.compile('ncbiVirusProteinHmm') 
p_nrHmm                = re.compile('nrHmm')               
p_keggVirusHmm         = re.compile('keggHmm')             
p_refseqProteinHmm     = re.compile('refseqPhmm')          
p_refseqGeneHmm        = re.compile('refseqGhmm')          
p_phantomeHmm          = re.compile('phantomeHmm')         
p_pvogsHmm             = re.compile('pvogsHmm')            
p_uniparcHmm           = re.compile('uniparcHmm')          
p_swissprotHmm         = re.compile('swissprotHmm')        
p_pfamHmm              = re.compile('pfamHmm')             
p_smartHmm             = re.compile('smartHmm')            
p_phageEnzymeHmm       = re.compile('phageEnzymeHmm')            

# Initialize
TRANSLATE_ONLY = False
HMM_SEARCH     = False

# Other patterns
p_comment  = re.compile('^#')
p_blank    = re.compile("^\s*$")

##### GET INPUT PARAMETERS #####

argList = sys.argv
argCount = len(argList)

for i in range(0,argCount):

    # Look for parameter tags
    match_outputDirParam           = re.search(p_outputDirParam,           argList[i])
    match_genomeFileParam          = re.search(p_genomeFileParam,          argList[i])
    match_geneFileParam            = re.search(p_geneFileParam,            argList[i])
    match_proteinFileParam         = re.search(p_proteinFileParam,         argList[i])
    match_geneCallerParam          = re.search(p_geneCallerParam,          argList[i])
    match_geneticCodeParam         = re.search(p_geneticCodeParam,         argList[i])
    match_geneCallFileParam        = re.search(p_geneCallFileParam,        argList[i])
    match_contigNameParam          = re.search(p_contigNameParam,          argList[i])
    match_genomeTypeParam          = re.search(p_genomeTypeParam,          argList[i])
    match_genomeNameParam          = re.search(p_genomeNameParam,          argList[i])
    match_speciesParam             = re.search(p_speciesParam,             argList[i])
    match_blastpIdentityParam      = re.search(p_blastpIdentityParam,      argList[i])
    match_blastnIdentityParam      = re.search(p_blastnIdentityParam,      argList[i])
    match_blastpHitCountParam      = re.search(p_blastpHitCountParam,      argList[i])
    match_blastnHitCountParam      = re.search(p_blastnHitCountParam,      argList[i])
    match_translateOnlyParam       = re.search(p_translateOnlyParam,       argList[i])
    match_blastDBsStringParam      = re.search(p_blastDBsStringParam,      argList[i])
    match_hmmProgramDBsStringParam = re.search(p_hmmProgramDBsStringParam, argList[i])

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

    if match_geneCallFileParam:
        if i < argCount:
            infile_geneCall = argList[i+1] 

    # Other parameterized arguments

    if match_geneCallerParam:
        if i < argCount:
            geneCaller = argList[i+1]

    if match_contigNameParam:  #*** TEMPORARY until PHANOTATE reports contig name
        CONTIG = True
        if i < argCount:
            contigName = argList[i+1]

    if match_genomeTypeParam:
        if i < argCount:
            genomeType = argList[i+1]

    if match_geneticCodeParam:
        if i < argCount:
            geneticCode = argList[i+1]

    if match_genomeNameParam:
        if i < argCount:
            genomeName = argList[i+1]
 
    if match_speciesParam:
        if i < argCount:
            species = argList[i+1]
  
    if match_translateOnlyParam:
        if i < argCount:
            value = argList[i+1]
            if value.lower() == 'yes' or value.lower() == 'true' or value.lower() == 'on' or value.lower() == 'y':
                TRANSLATE_ONLY = True
                translateOnly = True
            else:
                TRANSLATE_ONLY = False 
                translateOnly = False 

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

    if match_blastDBsStringParam:
        if i < argCount:
            value = argList[i+1]
            match_ncbiVirus        = re.search(p_ncbiVirus,value)
            match_ncbiVirusProtein = re.search(p_ncbiVirusProtein,value)
            match_nr               = re.search(p_nr,value)
            match_kegg             = re.search(p_keggVirus,value)
            match_refseqProtein    = re.search(p_refseqProtein,value)
            match_refseqGene       = re.search(p_refseqGene,value)
            match_phantome         = re.search(p_phantome,value)
            match_pvogs            = re.search(p_pvogs,value)
            match_uniparc          = re.search(p_uniparc,value)
            match_swissprot        = re.search(p_swissprot,value)
            match_pfam             = re.search(p_pfam,value)
            match_smart            = re.search(p_smart,value)
            match_phageEnzyme      = re.search(p_phageEnzyme,value)
            if match_ncbiVirus:
                NCBI_VIRUS_GENOME_BLAST = True
            if match_ncbiVirusProtein:
                NCBI_VIRUS_PROTEIN_BLAST = True
            if match_nr:
                NR_BLAST = True
            if match_kegg:
                KEGG_VIRUS_BLAST = True
            if match_refseqProtein:
                REFSEQ_PROTEIN_BLAST = True
            if match_refseqGene:
                REFSEQ_GENE_BLAST = True
            if match_phantome:
                PHANTOME_BLAST = True
            if match_pvogs:
                PVOGS_BLAST = True 
            if match_swissprot:
                SWISSPROT_BLAST = True
            if match_pfam:
                PFAM_BLAST = True
            if match_smart:
                SMART_BLAST = True
            if match_phageEnzyme:
                PHAGE_ENZYME_BLAST = True

    if match_hmmProgramDBsStringParam:
        if i < argCount:
            value = argList[i+1]
            match_hmmProgram          = re.search(p_hmmProgram,value)
            match_ncbiVirusHmm        = re.search(p_ncbiVirusHmm,value)
            match_ncbiVirusProteinHmm = re.search(p_ncbiVirusProteinHmm,value)
            match_nrHmm               = re.search(p_nrHmm,value)
            match_keggHmm             = re.search(p_keggVirusHmm,value)
            match_refseqProteinHmm    = re.search(p_refseqProteinHmm,value)
            match_refseqGeneHmm       = re.search(p_refseqGeneHmm,value)
            match_phantomeHmm         = re.search(p_phantomeHmm,value)
            match_pvogsHmm            = re.search(p_pvogsHmm,value)
            match_uniparcHmm          = re.search(p_uniparcHmm,value)
            match_swissprotHmm        = re.search(p_swissprotHmm,value)
            match_pfamHmm             = re.search(p_pfamHmm,value)
            match_smartHmm            = re.search(p_smart,value)
            match_phageEnzymeHmm      = re.search(p_phageEnzyme,value)
            if match_hmmProgram:
                hmmProgram = match_hmmProgram.group(1)
                HMM_SEARCH = True
            if match_ncbiVirusHmm:
                NCBI_VIRUS_HMM = True
            if match_ncbiVirusProteinHmm:
                NCBI_VIRUS_PROTEIN_HMM = True
            if match_nrHmm:
                NR_HMM = True
            if match_keggHmm:
                KEGG_VIRUS_HMM = True
            if match_refseqProteinHmm:
                REFSEQ_PROTEIN_HMM = True
            if match_refseqGeneHmm:
                REFSEQ_GENE_HMM = True
            if match_phantomeHmm:
                PHANTOME_HMM = True
            if match_pvogsHmm:
                PVOGS_HMM = True 
            if match_swissprotHmm:
                SWISSPROT_HMM = True
            if match_pfamHmm:
                PFAM_HMM = True
            if match_smartHmm:
                SMART_HMM = True
            if match_phageEnzymeHmm:
                PHAGE_ENZYME_HMM = True

# Open and Check files

fileError = False

try:
    LOGFILE_H = open(logfile,"w")
except IOError as e:
    fileError = True
    if PHATE_WARNINGS == 'True':
        print(e, "logfile,", logfile) 
        print("Cannot write log file; Exit computation at phate_sequenceAnnotation_main.py, attpempting to open log file")
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
    GENE_CALL_FILE = open(infile_geneCall,"r")
except IOError as e:
    fileError = True
    if PHATE_WARNINGS == 'True':
        print(e, "gene call file,", infile_geneCall)
    LOGFILE_H.write("%s%s%s%s\n" % ("fileError ",e," gene call file, ", infile_geneCall))

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
    print("Check the formats of your input file(s):")
    print("    genome file is", infile_genome)
    print("    gene call file is", infile_geneCall)
    print("    outfile is", outfile)
    print("    gfffile is", gfffile)
    LOGFILE_H.write("%s%s\n" % ("Terminating due to file error at ",datetime.datetime.now()))
    LOGFILE_H.close(); exit(0)

# Communicate to log
LOGFILE_H.write("%s%s\n" % ("Parameters recorded at ",datetime.datetime.now()))
LOGFILE_H.write("%s%s\n" % ("outputDir is", outputDir))
LOGFILE_H.write("%s%s\n" % ("outfile is ",outfile))
LOGFILE_H.write("%s%s\n" % ("gfffile is ",gfffile))
LOGFILE_H.write("%s%s\n" % ("infile_genome is ", infile_genome))
LOGFILE_H.write("%s%s\n" % ("geneCaller is ",geneCaller))
LOGFILE_H.write("%s%s\n" % ("geneticCode is ",geneticCode))
LOGFILE_H.write("%s%s\n" % ("infile_geneCall is ", infile_geneCall))
LOGFILE_H.write("%s%s\n" % ("infile_gene is ", infile_gene))
LOGFILE_H.write("%s%s\n" % ("infile_protein is ", infile_protein))
LOGFILE_H.write("%s%s\n" % ("genomeType is ",genomeType))
LOGFILE_H.write("%s%s\n" % ("genomeName is ",genomeName))
LOGFILE_H.write("%s%s\n" % ("genomeSpecies is ",genomeSpecies))
LOGFILE_H.write("%s%s\n" % ("blastpIdentity is ",blastpIdentity))
LOGFILE_H.write("%s%s\n" % ("blastpHitCount is ",blastpHitCount))
LOGFILE_H.write("%s%s\n" % ("blastnHitCount is ",blastnHitCount))
if TRANSLATE_ONLY:
    LOGFILE_H.write("%s\n" % ("Translating only; no annotation."))
else:
    LOGFILE_H.write("%s\n" % ("Annotating"))
LOGFILE_H.write("%s%s\n" % ("NCBI_VIRUS_GENOME_BLAST is ",NCBI_VIRUS_GENOME_BLAST))
LOGFILE_H.write("%s%s\n" % ("NCBI_VIRUS_PROTEIN_BLAST is ",NCBI_VIRUS_PROTEIN_BLAST))
LOGFILE_H.write("%s%s\n" % ("KEGG_VIRUS_BLAST is ",KEGG_VIRUS_BLAST))
LOGFILE_H.write("%s%s\n" % ("NR_BLAST is ",NR_BLAST))
LOGFILE_H.write("%s%s\n" % ("REFSEQ_PROTEIN_BLAST is ",REFSEQ_PROTEIN_BLAST))
LOGFILE_H.write("%s%s\n" % ("REFSEQ_GENE_BLAST is ",REFSEQ_GENE_BLAST))
LOGFILE_H.write("%s%s\n" % ("PHANTOME_BLAST is ",PHANTOME_BLAST))
LOGFILE_H.write("%s%s\n" % ("PVOGS_BLAST is ",PVOGS_BLAST))
LOGFILE_H.write("%s%s\n" % ("SWISSPROT_BLAST is ",SWISSPROT_BLAST))
LOGFILE_H.write("%s%s\n" % ("PHAGE_ENZYME_BLAST is ",PHAGE_ENZYME_BLAST))
LOGFILE_H.write("%s%s\n" % ("hmmProgram is ",hmmProgram))
LOGFILE_H.write("%s%s\n" % ("NCBI_VIRUS_GENOME_HMM is ",NCBI_VIRUS_GENOME_HMM))
LOGFILE_H.write("%s%s\n" % ("NCBI_VIRUS_PROTEIN_HMM is ",NCBI_VIRUS_PROTEIN_HMM))
LOGFILE_H.write("%s%s\n" % ("KEGG_VIRUS_HMM is ",KEGG_VIRUS_HMM))
LOGFILE_H.write("%s%s\n" % ("NR_HMM is ",NR_HMM))
LOGFILE_H.write("%s%s\n" % ("REFSEQ_PROTEIN_HMM is ",REFSEQ_PROTEIN_HMM))
LOGFILE_H.write("%s%s\n" % ("REFSEQ_GENE_HMM is ",REFSEQ_GENE_HMM))
LOGFILE_H.write("%s%s\n" % ("PHANTOME_HMM is ",PHANTOME_HMM))
LOGFILE_H.write("%s%s\n" % ("PVOGS_HMM is ",PVOGS_HMM))
LOGFILE_H.write("%s%s\n" % ("SWISSPROT_HMM is ",SWISSPROT_HMM))
LOGFILE_H.write("%s%s\n" % ("PHAGE_ENZYME_HMM is ",PHAGE_ENZYME_HMM))

# Communicate to user
if PHATE_MESSAGES == 'True':
    print("outputDir is", outputDir)
    print("outfile is", outfile)
    print("gfffile is", gfffile)
    print("genome file is", infile_genome)
    print("gene call file is", infile_geneCall)
    print("gene caller is", geneCaller)
    print("genetic code is", geneticCode)
    print("gene call file is", infile_geneCall)
    print("gene file is", infile_gene)
    print("protein file is", infile_protein)
    print("genome type is", genomeType)
    print("name is", name)
    print("contigName is", contigName)
    print("species is", species)
    print("blastp identity is", blastpIdentity)
    print("blastp hit count is", blastpHitCount)
    print("blastn hit count is", blastnHitCount)
    if TRANSLATE_ONLY:
        print("Translating only; no annotation.")
    else:
        print("Annotating")
    print("NCBI_VIRUS_GENOME_BLAST is", NCBI_VIRUS_GENOME_BLAST)
    print("NCBI_VIRUS_PROTEIN_BLAST is", NCBI_VIRUS_PROTEIN_BLAST)
    print("NR_BLAST is", NR_BLAST)
    print("KEGG_VIRUS_BLAST is", KEGG_VIRUS_BLAST)
    print("REFSEQ_PROTEIN_BLAST is", REFSEQ_PROTEIN_BLAST)
    print("REFSEQ_GENE_BLAST is", REFSEQ_GENE_BLAST)
    print("PHANTOME_BLAST is", PHANTOME_BLAST)
    print("PVOGS_BLAST is", PVOGS_BLAST)
    print("SWISSRPOT_BLAST is", SWISSPROT_BLAST)
    print("PHAGE_ENZYME_BLAST is", PHAGE_ENZYME_BLAST)
    print("hmmProgram is", hmmProgram)
    print("NCBI_VIRUS_HMM is", NCBI_VIRUS_HMM)
    print("NCBI_VIRUS_PROTEIN_HMM is", NCBI_VIRUS_PROTEIN_HMM)
    print("KEGG_VIRUS_HMM is", KEGG_VIRUS_HMM)
    print("NR_HMM is", NR_HMM)
    print("REFSEQ_PROTEIN_HMM is", REFSEQ_PROTEIN_HMM)
    print("REFSEQ_GENE_HMM is", REFSEQ_GENE_HMM)
    print("PHANTOME_HMM is", PHANTOME_HMM)
    print("PVOGS_HMM is", PVOGS_HMM)
    print("SWISSPROT_HMM is", SWISSPROT_HMM)
    print("PHAGE_ENZYME_HMM is", PHAGE_ENZYME_HMM)

if PHATE_PROGRESS == 'True':
    print("Sequence annotation main says: Configuration complete for sequence annotation module.")

##### BEGIN MAIN 

geneCallInfo = {      # For passing info to genomeSequence module
    'geneCaller'           : geneCaller, 
    'geneCallFile'         : infile_geneCall, 
    'genomeName'           : genomeName,
}

# Create a genome object and set parameters 

if PHATE_PROGRESS == 'True':
    print("Preparing for sequence annotation: setting parameters...")
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
#myGenome.contigSet.assignContig(contigName) #***  TEMPORARY handles only finished genome / single contig for now
if PHATE_MESSAGES == 'True':
    print("Sequence annotation main says: contigName is", myGenome.contigSet.contig)

# Create a hash to record contig names and their sequence lengths, for GFF output
contigSeqLen_hash = {}
for i in range(0,len(myGenome.contigSet.fastaList)):
    contigName = myGenome.contigSet.fastaList[i].cleanHeader
    #seqLen     = myGenome.contigSet.fastaList[i].sequenceLength
    seqLen     = len(myGenome.contigSet.fastaList[i].sequence)
    contigSeqLen_hash[contigName] = seqLen

# Extract gene calls
if PHATE_PROGRESS == 'True':
    print("Sequence annotation main says: Processing gene calls...")
LOGFILE_H.write("%s%s\n" % ("Processing gene calls at ",datetime.datetime.now()))
# Record the gene calls, and pass the hash so contig lengths can also be recorded...
# ...for access when ultimately writing the GFF output file. (oh what we do for want of a pointer.)
print("geneCallInfo",geneCallInfo,"GENE_CALL_FILE",GENE_CALL_FILE,"contigSeqLen_hash",contigSeqLen_hash)
#myGenome.processGeneCalls(geneCallInfo,GENE_CALL_FILE,contigSeqLen_hash)
myGenome.processGeneCalls(geneCallInfo,GENE_CALL_FILE)
myGenome.cleanUpAfterEMBOSS()
#myGenome.contigSet.assignContig2all(contigName) #*** TEMPORARY handles only finished genome /single contig for now
GENE_CALL_FILE.close()

# Output the gene and protein sets, if newly created

fastaOut = {
    "mtype"      : "",
    "headerType" : "",
    "filename"   : "",
}

if PHATE_PROGRESS == 'True':
    print("Writing genes file")
LOGFILE_H.write("%s\n" % ("Writing genes file"))
# Print out newly created gene list
fastaOut["mtype"] = "gene"
fastaOut["headerType"] = "full"  #*** Should this be "compound" ???
fastaOut["filename"] = infile_gene 
myGenome.printFastas2file(fastaOut)

if PHATE_PROGRESS == 'True':
    print("Writing peptides file")
LOGFILE_H.write("%s\n" % ("Writing peptides file"))
# Print out newly created protein list 
fastaOut["mtype"] = "protein"   #*** Should this be "compound" ???
fastaOut["headerType"] = "full"
fastaOut["filename"] = infile_protein 
myGenome.printFastas2file(fastaOut)

if PHATE_PROGRESS == 'True':
    print("Gene and protein files created.")
LOGFILE_H.write("%s%s\n" % ("Gene and protein files created at ",datetime.datetime.now()))

# If user specified to translate only, then skip this segment of the pipeline.
if TRANSLATE_ONLY:
    if PHATE_PROGRESS == 'True':
        print("Translate only: computations completed.")
    LOGFILE_H.write("%s%s\n" % ("Translating only: computations completed at ",datetime.datetime.now()))
else:
    # Create a blast object and set parameters
    if PHATE_PROGRESS == 'True':
        print("Preparing for blast...")
    LOGFILE_H.write("%s\n" % ("Creating a blast object"))
    blast = phate_blast.multiBlast()

    if PHATE_PROGRESS == 'True':
        print("Sequence annotation module says: Preparing to run", blast.blastFlavor)
    if PHATE_MESSAGES == 'True':
        print("Sequence annotation module says: Running at the following settings:")
        blast.printParameters()
    LOGFILE_H.write("%s%s%s%s\n" % (datetime.datetime.now(), " Preparing to run ", blast.blastFlavor, " at the following settings:"))
    blast.printParameters2file(LOGFILE_H)

    # Create directory for blast output
    blastOutputDir = outputDir + '/BLAST/'
    try:
        os.stat(blastOutputDir)
    except:
        os.mkdir(blastOutputDir)

    # Create an hmm object and set parameters
    if PHATE_PROGRESS == 'True':
        print("Preparing for hmm...")
    LOGFILE_H.write("%s\n" % ("Creating an hmm object"))
    hmm = phate_hmm.multiHMM()

    # Create directory for hmm output
    hmmOutputDir = outputDir + 'HMM/'
    try:
        os.stat(hmmOutputDir)
    except:
        os.mkdir(hmmOutputDir)

    ##### Run blast against whatever sequences/databases we have:  genome, gene, protein, phage databases

    # GENOME BLAST

    # Create blast output directory for genome blast
    genomeBlastOutputDir = blastOutputDir + 'Genome/'
    try:
        os.stat(genomeBlastOutputDir)
    except:
        os.mkdir(genomeBlastOutputDir)

    # Prepare for genome blast
    myParamSet = {
        'identityMin'          : GENOME_IDENTITY_MIN,  #*** Should not be hard coded, but should be a low-stringency setting
        'identitySelect'       : GENOME_IDENTITY_SELECT,  #*** this one too 
        'evalueMin'            : EVALUE_MIN,
        'evalueSelect'         : EVALUE_SELECT,
        'topHitCount'          : int(blastnHitCount),
        'outputFormat'         : XML_OUT_FORMAT,  # blast output format (use 5=XML or 7=LIST only) 
        'scoreEdge'            : SCORE_EDGE,
        'overhang'             : OVERHANG,
        'geneCallDir'          : outputDir, 
        'blastOutDir'          : genomeBlastOutputDir,
        'ncbiVirusGenomeBlast' : NCBI_VIRUS_GENOME_BLAST,
    }
    blast.setBlastParameters(myParamSet)
    blast.setBlastFlavor('blastn') 

    LOGFILE_H.write("%s%s%s%s\n" % (datetime.datetime.now(), " Preparing to run ", blast.blastFlavor, " at the following settings:"))
    blast.printParameters2file(LOGFILE_H)

    if PHATE_PROGRESS == 'True':
        print("Sequence annotation main says: Preparing to run", blast.blastFlavor)
    if PHATE_MESSAGES == 'True':
        print("Sequence annotation main says: Running", blast.blastFlavor, "at the following settings:")
        blast.printParameters()

    if PHATE_PROGRESS == 'True':
        print("Sequence annotation main says: Running Blast against phage genome database(s)...")

    LOGFILE_H.write("%s%s%s%s\n" % (datetime.datetime.now(), " Preparing to run ", blast.blastFlavor, " at the following settings:"))
    blast.printParameters2file(LOGFILE_H)
    LOGFILE_H.write("%s\n" % ("Running Blast against phage genome database(s)..."))

    # Run Genome blast 
    blast.runBlast(myGenome.contigSet,'genome')

    if PHATE_PROGRESS == 'True':
        print("Sequence annotation main says: Genome blast complete")
    LOGFILE_H.write("%s%s\n" % ("Genome blast complete at ", datetime.datetime.now()))

    # GENE BLAST

    LOGFILE_H.write("%s\n" % ("Preparing to run gene blast"))
    if PHATE_PROGRESS == 'True':
        print("Sequence annotation main says: Preparing to run gene blast...")

    # Create blast output directory for gene blast
    geneBlastOutputDir = blastOutputDir + 'Gene/'
    try:
        os.stat(geneBlastOutputDir)
    except:
        os.mkdir(geneBlastOutputDir)

    # Prepare for gene blast
    myParamSet = {
        'identityMin'     : blastpIdentity,  #*** Setting as per blastp, for now 
        'identitySelect'  : blastpIdentity,  #*** this one too 
        'evalueMin'       : 10,
        'evalueSelect'    : 10,
        'topHitCount'     : int(blastpHitCount),  #*** maybe should parameterize this
        'outputFormat'    : XML_OUT_FORMAT,  # XML=5; LIST=7
        'scoreEdge'       : 0.1,
        'overhang'        : 0.1,
        'geneCallDir'     : outputDir, 
        'blastOutDir'     : geneBlastOutputDir,
        'refseqGeneBlast' : REFSEQ_GENE_BLAST,
    }
    blast.setBlastParameters(myParamSet)
    blast.setBlastFlavor('blastn') 

    if PHATE_PROGRESS == 'True':
        print("Sequence annotation main says: Running Blast against gene database(s)...")
    LOGFILE_H.write("%s\n" % ("Running Blast against gene databases..."))

    # Run Gene blast
    blast.runBlast(myGenome.geneSet,'gene')

    if PHATE_PROGRESS == 'True':
        print("Sequence annotation main says: Gene blast complete.")
    LOGFILE_H.write("%s%s\n" % ("Gene blast complete at ",datetime.datetime.now()))

    # PROTEIN BLAST

    # Create blast output directory for protein blast
    proteinBlastOutputDir = blastOutputDir + 'Protein/'
    try:
        os.stat(proteinBlastOutputDir)
    except:
        os.mkdir(proteinBlastOutputDir)

    if PHATE_PROGRESS == 'True':
        print("Sequence annotation main says: Preparing for protein blast...")
    LOGFILE_H.write("%s%s\n" % ("Preparing for protein blast at ",datetime.datetime.now()))

    # Prepare for protein blast
    myParamSet = {
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
        'nrBlast'               : NR_BLAST,
        'ncbiVirusProteinBlast' : NCBI_VIRUS_PROTEIN_BLAST,
        'keggVirusBlast'        : KEGG_VIRUS_BLAST,
        'refseqProteinBlast'    : REFSEQ_PROTEIN_BLAST,
        'phantomeBlast'         : PHANTOME_BLAST,
        'pvogsBlast'            : PVOGS_BLAST,
        'swissprotBlast'        : SWISSPROT_BLAST,
        'phageEnzymeBlast'      : PHAGE_ENZYME_BLAST,
    }
    blast.setBlastParameters(myParamSet)
    blast.setBlastFlavor('blastp')

    if PHATE_PROGRESS == 'True': 
        print("Sequence annotation main says: Running Blast against protein database(s)...")
    LOGFILE_H.write("%s\n" % ("Running Blast against protein database(s)..."))

    # Run protein blast
    blast.runBlast(myGenome.proteinSet,'protein')

    if PHATE_PROGRESS == 'True':
        print("Sequence annotation main says: Protein blast complete.")
    LOGFILE_H.write("%s%s\n" % ("Protein blast complete at ", datetime.datetime.now()))

    # Write out pVOGs sequences to enable alignments

    # HMM SEARCH

    if PHATE_PROGRESS == 'True':
        print("Sequence annotation main says: Preparing for hmm search...")
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
    myParamSet = {
        'hmmProgram'            : hmmProgram,
        'geneCallDir'           : outputDir,
        'hmmOutDir'             : proteinHmmOutputDir,
        'pvogsOutDir'           : proteinHmmOutputDir,
        'nrHmm'                 : NR_HMM,
        'keggVirusHmm'          : KEGG_VIRUS_HMM,
        'refseqProteinHmm'      : REFSEQ_PROTEIN_HMM,
        'phantomeHmm'           : PHANTOME_HMM,
        'pvogsHmm'              : PVOGS_HMM,
        'swissprotHmm'          : SWISSPROT_HMM,
        'pfamHmm'               : PFAM_HMM,
        'smartHmm'              : SMART_HMM,
        'phageEnzymeHmm'        : PHAGE_ENZYME_HMM,
        'ncbiVirusGenomeHmm'    : NCBI_VIRUS_GENOME_HMM,
        'ncbiVirusProteinHmm'   : NCBI_VIRUS_PROTEIN_HMM,
    }
    if DEBUG:
        print("DEBUGGING Hmm parameters:", myParamSet)

    hmm.setHmmParameters(myParamSet)

    if PHATE_PROGRESS == 'True':
        print("Sequence annotation main says: Running Hmm search against protein database(s)...")
    LOGFILE_H.write("%s%s\n" % ("Running Hmm search against protein database(s) at ", datetime.datetime.now()))

    # Run protein hmm search
    hmm.runHmm(myGenome.proteinSet,'protein')

    LOGFILE_H.write("%s%s\n" % ("HMM search complete at ",datetime.datetime.now()))

    # REPORT OUT 

    if PHATE_PROGRESS == 'True':
        print("Reporting final annotations")
    LOGFILE_H.write("%s%s\n" % ("Reporting final annotations at ",datetime.datetime.now()))
    myGenome.printGenomeData2file_tab(OUTFILE)
    myGenome.printGenomeData2file_GFF(GFFFILE)

##### CLEAN UP

if PHATE_PROGRESS == 'True':
    print("Sequence annotation main says: Done!")
GENOME_FILE.close()
OUTFILE.close()
GFFFILE.close()

LOGFILE_H.write("%s%s\n" % ("Sequence annotation processing complete at ",datetime.datetime.now()))
LOGFILE_H.close()

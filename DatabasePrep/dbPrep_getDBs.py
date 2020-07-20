#!/usr/bin/env python

#############################################################################
#
# program: dbPrep_getDBs.py
#
# programmer:  C. E. Zhou
#
# Programmer's Notes:
# 1) Move this script to the /DatabasePrep/ folder. Adjust accordingly.
# 2) update_blastdb.pl --showall <to see which databases are available>
#
# Summary:  This script facilitates the downloading of databases to be used with multiPhATE.
#
# Most recent update:  17 July 2020
#
##############################################################################

import os, sys, re, time, datetime
from ftplib import FTP
from pathlib import Path

##############################################################################
# CONSTANTS, BOOLEANS

VOG_VERSION      = "99"
VOG_DOWNLOAD_URL = "http://fileshare.csb.univie.ac.at/vog/vog" + VOG_VERSION + "/"
ftpAddr  = "ftp.ncbi.nlm.nih.gov"
httpAddr = "https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/"
virGenomeFile       = "all.fna"
virProteinFile      = "all.faa"
virGenomeFile_gz    = virGenomeFile  + ".tar.gz"
virProteinFile_gz   = virProteinFile + ".tar.gz"
virGenomeFile_tar   = virGenomeFile  + ".tar"
virProteinFile_tar  = virProteinFile + ".tar"
virGenome_httpAddr  = os.path.join(httpAddr,virGenomeFile_gz) 
virProtein_httpAddr = os.path.join(httpAddr,virProteinFile_gz) 
accn2taxid_httpAddr = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid"
accn2taxid_file     = "nucl_gb.accession2taxid"
accn2taxid_file_gz  = accn2taxid_file + ".gz"
accn2taxid_fileAddr = os.path.join(accn2taxid_httpAddr,accn2taxid_file_gz)
#pvogHmm_httpAddr    = "http://dmk-brain.ecn.uiowa.edu/pVOGs/downloads/All/AllvogHMMprofiles.tar.gz"
pvogHmm_httpAddr    = "https://ftp.ncbi.nlm.nih.gov/pub/kristensen/pVOGs/downloads/All/AllvogHMMprofiles.tar.gz"

BLAST               = False
NCBI_VIRUS_GENOME   = False
NCBI_VIRUS_PROTEIN  = False
REFSEQ_GENE         = False
REFSEQ_PROTEIN      = False
SWISSPROT           = False
NR                  = False
PHANTOME            = False 
PVOGS               = False 
PVOG_HMMS           = False
VOGS                = False
VOG_HMMS            = False

# VARIABLES
blast               = ''
ncbi_virus_genome   = ''
ncbi_virus_protein  = ''
refseq_gene         = ''
refseq_protein      = ''
swissprot           = ''
nr                  = ''
pvogs               = ''
pvog_hmms           = ''
vogs                = ''
vog_hmms            = ''
decision            = ''
blastPath           = '' # path to user's blast installation
cwd                 = '' # current working direcgtory
emailAddr           = '' # user's email address for ftp login

# Set up database directories
cwd                    = os.getcwd()
# Get directory one level up; this is where the Databases directory will go
baseDir                = Path(__file__).resolve().parents[1]
dbDir                  = os.path.join(baseDir,        "Databases/") 
ncbiDir                = os.path.join(dbDir,          "NCBI/")
ncbiGenomeDir          = os.path.join(ncbiDir,        "Virus_Genome/")
ncbiGenomeDownloadDir  = os.path.join(ncbiGenomeDir,  "Download/")
ncbiProteinDir         = os.path.join(ncbiDir,        "Virus_Protein/")
ncbiProteinDownloadDir = os.path.join(ncbiProteinDir, "Download/")
nrDir                  = os.path.join(dbDir,          "NR/")
refseqDir              = os.path.join(dbDir,          "Refseq/")
refseqGeneDir          = os.path.join(refseqDir,      "Gene/")
refseqProteinDir       = os.path.join(refseqDir,      "Protein/")
swissprotDir           = os.path.join(dbDir,          "Swissprot/")
phantomeDir            = os.path.join(dbDir,          "Phantome/")
pVOGsDir               = os.path.join(dbDir,          "pVOGs/")
pVOGhmmsDir            = os.path.join(dbDir,          "pVOGhmms/")
VOGsDir                = os.path.join(dbDir,          "VOGs/")
VOGsBlastDBdir         = os.path.join(VOGsDir,        "BlastDBs/")
VOGhmmsDir             = os.path.join(dbDir,          "VOGhmms/")

# Filenames from source
vogProteins_filename             = "vog.proteins.all.fa.gz"
vogGenes_filename                = "vog.genes.all.fa.gz"
vogAnnotations_filename          = "vog.annotations.tsv.gz"
vogHmms_filename                 = "vog.hmm.tar.gz"
vogHmms_filename_tar             = "vog.hmm.tar"
vogMembers_filename              = "vog.members.tsv.gz"
vogFunctionalCategories_filename = "vog_functional_categories.txt"
pvogHmms_filename                = "pvog.hmm.tar.gz"
pvogHmms_filename_tar            = "pvog.hmm.tar"

# Filenames for processed data
vogGenes                         = "vog.genes.tagged.all.fa"
vogProteins                      = "vog.proteins.tagged.all.fa"
ncbiVirusGenomes                 = "ncbiVirusGenomes.fasta"
ncbiVirusProteins                = "ncbiVirusProteins.faa"

if not os.path.exists(dbDir):
    os.mkdir(dbDir)
if not os.path.exists(ncbiDir):
    os.mkdir(ncbiDir)
if not os.path.exists(ncbiGenomeDir):
    os.mkdir(ncbiGenomeDir)
if not os.path.exists(ncbiGenomeDownloadDir):
    os.mkdir(ncbiGenomeDownloadDir)
if not os.path.exists(ncbiProteinDir):
    os.mkdir(ncbiProteinDir)
if not os.path.exists(ncbiProteinDownloadDir):
    os.mkdir(ncbiProteinDownloadDir)
if not os.path.exists(nrDir):
    os.mkdir(nrDir)
if not os.path.exists(refseqDir):
    os.mkdir(refseqDir)
#if not os.path.exists(refseqGeneDir):  # Skipping Refseq Gene (at least for now)
#    os.mkdir(refseqGeneDir)
if not os.path.exists(refseqProteinDir):
    os.mkdir(refseqProteinDir)
if not os.path.exists(swissprotDir):
    os.mkdir(swissprotDir)
if not os.path.exists(phantomeDir):
    os.mkdir(phantomeDir)
if not os.path.exists(pVOGsDir):
    os.mkdir(pVOGsDir)
if not os.path.exists(pVOGhmmsDir):
    os.mkdir(pVOGhmmsDir)
if not os.path.exists(VOGsDir):
    os.mkdir(VOGsDir)
if not os.path.exists(VOGhmmsDir):
    os.mkdir(VOGhmmsDir)

##############################################################################
# First, determine if user needs to download BLAST+.

print ("Let's download the databases you will need to run multiPhATE")
print ("First, you need blast+ in order to install several of the databases")
print ("Please confirm that you have downloaded and installed blast+: type 'y' (yes) or 'n' (no)")
blast = input()
if re.search('Y|y|Yes|yes|YES', blast):
    BLAST = True 
    print ("That's great! Now tell me where your blast executables are located.") 
    print ("If you installed them within your Conda environment, then you should find")
    print ("them in a path something like, \"/Users/myName/miniconda3/envs/multiPhate/bin/\"")
    print ("That directory should contain executables: blastn, blastp, blastx, and makeblastdb")
    print ("Please input the fully qualified path to blast+ : ")
    #blastPath = input()
elif re.search('N|n|No|no|NO', blast):
    BLAST = False 
    print ("Please consult the README file for how to acquire and install BLAST+.") 
    print ("Note: The easiest way to install BLAST+ is within a Conda environment.")
    print ("Without blast+, we can still download some of the databases.")
    print ("Shall we continue? please respond 'y' or 'n': ")
    toContinue = input()
    if re.search('Y','y','yes','Yes','YES', toContinue):
        pass 
    else:
        print ("Bye")
        exit()
else:
    print ("That was not a correct response; please try again")
    exit()

##############################################################################
# Next, determine which databases to download

print ("Please indicate which databases you would like to download...")
print ("For each, type 'y' or 'n'")

##### NCBI Virus Genome database is downloaded via ftp 
print ("NCBI Virus Genome database: ('y'/'n')")
ncbi_virus_genome = input()
if re.search('Y|y|yes|Yes|YES',ncbi_virus_genome):
    print ("Great, let's download NCBI Virus Genome")
    NCBI_VIRUS_GENOME = True
elif re.search('N|n|no|No|NO',ncbi_virus_genome):
    print ("Ok, we'll skip that one")
else:
    print ("That was not a correct response; please run this script again to download the database.")

##### NCBI Virus Protein database is downloaded via ftp 
print ("NCBI Virus Protein database: ('y'/'n')")
ncbi_virus_protein = input()
if re.search('Y|y|yes|Yes|YES',ncbi_virus_protein):
    print ("Great, let's download NCBI Virus Protein")
    NCBI_VIRUS_PROTEIN = True
elif re.search('N|n|no|No|NO',ncbi_virus_protein):
    print ("Ok, we'll skip that one")
else:
    print ("That was not a correct response; please run this script again to download the database.")

#### If downloading via FTP, get user's email address for login
#if NCBI_VIRUS_GENOME or NCBI_VIRUS_PROTEIN:
#    print ("Please enter your email address for ftp anonymous login: ")
#    emailAddr = input()

##### REFSEQ PROTEIN
print ("Refseq Protein database: ('y'/'n')")
refseq_protein = input()
if re.search('Y|y|yes|Yes|YES',refseq_protein):
    print ("Great, let's download the Refseq Protein database")
    REFSEQ_PROTEIN = True
elif re.search('N|n|no|No|NO',refseq_protein):
    print ("Ok, we'll skip that one")
else:
    print ("That was not a correct response; please run this script again to download the database.")

##### REFSEQ GENE  # Refseq Gene is no longer downloadable via update_blastdb.pl
# Refseq Gene is human centric. I am replacing it with VOG Gene.
"""
print ("Refseq Gene database: ('y'/'n')")
refseq_gene = input()
if re.search('Y|y|yes|Yes|YES',refseq_gene):
    print ("Great, let's download the Refseq Gene database")
    REFSEQ_GENE = True
elif re.search('N|n|no|No|NO',refseq_gene):
    print ("Ok, we'll skip that one")
else:
    print ("That was not a correct response; please run this script again to download the database.")
"""
##### SWISSPROT
print ("Swissprot database: ('y'/'n')")
swissprot = input()
if re.search('Y|y|yes|Yes|YES',swissprot):
    print ("Great, let's download the Swissprot database")
    SWISSPROT = True
elif re.search('N|n|no|No|NO',swissprot):
    print ("Ok, we'll skip that one")
else:
    print ("That was not a correct response; please run this script again to download the database.")

##### VOGS
print ("VOGs database: ('y'/'n')")
vogs = input()
if re.search('Y|y|yes|Yes|YES',vogs):
    print ("Great, let's download the VOGs database. This includes gene and protein sequences.")
    VOGS = True
elif re.search('N|n|no|No|NO',vogs):
    print ("Ok, we'll skip that one")
else:
    print ("That was not a correct response; please run this script again to download the database.")

##### VOG HMMS
print ("VOGhmms database: ('y'/'n')")
voghmms = input()
if re.search('Y|y|yes|Yes|YES',voghmms):
    print ("Great, let's download the VOG HMMs database")
    VOG_HMMS = True
elif re.search('N|n|no|No|NO',voghmms):
    print ("Ok, we'll skip that one")
else:
    print ("That was not a correct response; please run this script again to download the database.")

##### PVOG HMMS
print ("PVOGhmms database: ('y'/'n')")
pvoghmms = input()
if re.search('Y|y|yes|Yes|YES',pvoghmms):
    print ("Great, let's download the PVOG HMMs database")
    PVOG_HMMS = True
elif re.search('N|n|no|No|NO',pvoghmms):
    print ("Ok, we'll skip that one")
else:
    print ("That was not a correct response; please run this script again to download the database.")

##### NR
print ("Caution: the NR database is extremely large. ")
print ("If you already have it on disk, you are advised not to download here. ")
print ("Downloading NR will take a long time. ")
print ("NR database: ('y'/'n')")
nr = input()
if re.search('Y|y|yes|Yes|YES',nr):
    print ("Great, let's download NCBI Virus Genome")
    NR = True
elif re.search('N|n|no|No|NO',nr):
    print ("Ok, we'll skip that one")
else:
    print ("That was not a correct response; please run this script again to download the database.")

##### PVOGS
print ("The PVOGS sequence database is included in the multiPhATE distribution.")
print ("Shall we format the PVOGS database for blast? 'y'/'n'")
pvogs = input()
if re.search('Y|y|yes|Yes|YES',pvogs):
    print ("Great, we will do the PVOGS formatting")
    PVOGS = True
elif re.search('N|n|no|No|NO',pvogs):
    print ("Ok, we'll skip that")
else:
    print ("That was not a correct response; please run this script again to format the PVOGS database.")

##### PHANTOME
print ("The PHANTOME sequence database is included in the multiPhATE distribution.")
print ("Shall we format the PHANTOME database for blast? 'y'/'n'")
phantome = input()
if re.search('Y|y|yes|Yes|YES',phantome):
    print ("Great, we will do the PHANTOME formatting")
    PHANTOME = True
elif re.search('N|n|no|No|NO',phantome):
    print ("Ok, we'll skip that")
else:
    print ("That was not a correct response; please run this script again to format the PHANTOME database.")

###################################################################################################
if (NCBI_VIRUS_GENOME or NCBI_VIRUS_PROTEIN or REFSEQ_PROTEIN or SWISSPROT or VOGS or VOG_HMMS or PVOG_HMMS or NR or PVOGS or PHANTOME):

    if (NCBI_VIRUS_GENOME or NCBI_VIRUS_PROTEIN or REFSEQ_PROTEIN or SWISSPROT or VOGS or VOG_HMMS or PVOG_HMMS or NR):
        print ("Ok, this is what we are going to download: ")
        #if not BLAST:
        #    print ("BLAST plus")
        if NCBI_VIRUS_GENOME:
            print ("NCBI Virus Genome database")
        if NCBI_VIRUS_PROTEIN:
            print ("NCBI Virus Protein database")
        if REFSEQ_PROTEIN:
            print ("Refseq Protein database")
        if SWISSPROT:
            print ("Swissprot database")
        if VOGS:
            print ("VOGS database")
        if VOG_HMMS:
            print ("VOG hmms database")
        if PVOG_HMMS:
            print ("PVOG hmms database")
        if NR:
            print ("NR database")

    if (PVOGS or PHANTOME):
        print ("This is what we are going to format:")
        if PVOGS:
            print ("PVOGS")
        if PHANTOME:
            print ("PHANTOME")

    print ("\nType 'go' to proceed, or 'stop' to reconsider. ")
    print ("(You can always run this script again.) ")
    decision = input()
    if (re.search('go|Go|GO',decision)):
        print ("Ok, let's get started with downloading.")
        print ("Databases will be installed into the Databases/ folder within multiPhATE")
    else:
        print ("Ok, maybe some other time. Bye!")
        exit()
else:
    print ("You have selected no downloads or databases to format. Have a happy day! :-)")
    exit()

##############################################################################
# Install BLAST+

# Install blast+ for user; NOT YET IN SERVICE
#if not BLAST:
#    command1 = "wget \"ftp//ftp.ncbi.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi- blast-${BLAST_VERSION}+-x64-linux.tar.gzi\""
#    command2 = "tar -zxf ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz"
#    command3 = "cd ncbi-blast-${BLAST_VERSION}+/bin"
#    command4 = "pwd"

##############################################################################
# Install NCBI_VIRUS_GENOME

if NCBI_VIRUS_GENOME:
    os.chdir(ncbiGenomeDownloadDir)

    # Download directories containing virus genome fasta files
    try:
        print ("Downloading NCBI Genome fasta files.")
        print ("This may take a while...")
        command = "wget " + virGenome_httpAddr
        success = os.system(command)
    except:
        print ("Command " + command + " failed.")

    # Untar the file
    try:
        print ("Uncompressing downloaded file, ",virGenomeFile_gz)
        command = "gunzip " + virGenomeFile_gz
        success = os.system(command)
        print ("gunzip success: ",success)
        command = "tar -xvf " + virGenomeFile_tar
        success = os.system(command)
        print ("tar success: ",success)
    except:
        print ("Uncompression of genome fasta file was unsuccessful.")

    # Read in a list of the directories
    try:
        print ("Listing directories and concatenating fasta files.")
        dirList = os.listdir(r'./')
        for directory in dirList:
            if os.path.isdir(directory):
                command = "ls " + directory + " > ./fileList.out"
                success = os.system(command)
                if os.path.isfile("./fileList.out"):
                    FILELIST_H = open("./fileList.out",'r')
                    fLines = FILELIST_H.read().splitlines()
                    for fLine in fLines:
                        nextFile = os.path.join(directory,fLine)
                        command = "cat " + nextFile + " >> " + ncbiVirusGenomes
                        success = os.system(command)
                        command = "rm " + nextFile  # Remove file once we are done with it
                        success = os.system(command)
                else:
                    print ("Problem with fileList.out")
                command = "rmdir " + directory  # Remove directory once we are done with it
                success = os.system(command)
        print ("Concatenation is complete!")
    except:
        print ("Concatenation failed.")

    # Move consolidated fasta file up one directory
    try:
        print ("Moving ",ncbiVirusGenomes," to ",ncbiGenomeDir," directory")
        command = "mv ./ncbiVirusGenomes.fasta ../."
        success = os.system(command)
    except:
        print ("Cannot move ncbiVirusGenomes.fastaile to directory ",ncbiGenomeDir)

    # Before leaving the Download directory, remove the ncbi downloaded file
    try:
        print ("Removing the NCBI Virus Genome compressed file that we downloaded.")
        command = "rm " + virGenomeFile_tar
        os.system(command)
        command = "rm fileList.out"
        os.system(command)
    except:
        print ("Problem removing ",ncbiVirusGenomes," or fileList.out file.")

    os.chdir(ncbiGenomeDir)  # Change to NCBI/Virus_Genome directory

    # Download the accession2taxid file
    try:
        print ("Downloading the accession2taxid file from ",accn2taxid_fileAddr)
        command = "wget " + accn2taxid_fileAddr
        success = os.system(command)
    except:
        print ("Download of the accession2taxid file failed.")

    try:
        print ("Unpacking accession2taxid.")
        command = "gunzip " + accn2taxid_file_gz
        success = os.system(command)
        print ("accn2taxid file has been unzipped.")
    except:
        print ("Unpacking of accession2taxid failed.")

    # Finally, format the Virus Genome database for blast
    try:
        print ("Formatting Virus Genome database for blast.")
        command = "makeblastdb -dbtype nucl -in " + ncbiVirusGenomes
        os.system(command)
        print ("Formatting complete.")
    except:
        print ("Formatting of Virus Genome database for blast failed.")

    os.chdir(cwd)  # Return to home directory

##############################################################################
# Install NCBI VIRUS PROTEIN database

if NCBI_VIRUS_PROTEIN:
    os.chdir(ncbiProteinDownloadDir)

    # Download directories containing virus peptide fasta files
    try:
        print ("Downloading NCBI Virus peptide fasta files.")
        print ("This may take a while...")
        command = "wget " + virProtein_httpAddr
        success = os.system(command)
    except:
        print ("Command " + command + " failed.")

    # Uncompress the downloaded file
    try:
        print ("Uncompressing NCBI Virus Protein gz file.")
        command = "gunzip " + virProteinFile_gz
        os.system(command)
        command = "tar -xvf " + virProteinFile_tar
        os.system(command)
        print ("Virus protein data is now uncompressed.")
    except:
        print ("Uncompression of file failed.")

    # Remove the DBV/ directory from download
    # This appears to hold peptides from raw Prokka gene-call predictions
    try:
        print ("Removing the DBV/ subdirectory, which we do not need.")
        command = "rm -r DBV/"
        os.system(command)
        print ("DBV/ has been removed.")
    except:
        print ("DBV directory removal failed.")

    # Read in a list of the directories
    try:
        print ("Listing directories and concatenating fasta files.")
        dirList = os.listdir(r'./')
        for directory in dirList:
            if os.path.isdir(directory):
                command = "ls " + directory + " > ./fileList.out"
                success = os.system(command)
                if os.path.isfile("./fileList.out"):
                    FILELIST_H = open("./fileList.out",'r')
                    fLines = FILELIST_H.read().splitlines()
                    for fLine in fLines:
                        nextFile = os.path.join(directory,fLine)
                        command = "cat " + nextFile + " >> " + ncbiVirusProteins
                        success = os.system(command)
                        command = "rm " + nextFile  # Remove file once done with it
                        success = os.system(command)
                else:
                    print ("Problem with fileLine.out")
                command = "rmdir " + directory  # Remove directory once done with it
                success = os.system(command)  
        print ("Concatenation is complete!")
    except:
        print ("Listing and concatenating files failed.")

    # Move consolidated fasta file up one directory
    try:
        print ("Moving ",ncbiVirusProteins," to ",ncbiProteinDir," directory")
        command = "mv ./ncbiVirusProteins.faa ../."
        success = os.system(command)
    except:
        print ("Cannot move ncbiVirusProteins.faa file to directory ",ncbiProteinDir)
   
    # Before leaving the Download directory, remove the ncbi downloaded file
    try:
        print ("Removing the NCBI Virus Genome compressed file that we downloaded.")
        command = "rm " + virProteinFile_tar
        os.system(command)
        command = "rm fileList.out"
        os.system(command)
        command = "rm ls.out"
        os.system(command)
    except:
        print ("Problem removing ",ncbiVirusProteins," or fileList.out file.")

    os.chdir(ncbiProteinDir)  # Change to NCBI/Virus_Protein directory

    # Finally, format the Virus Protein database for blast
    try:
        print ("Formatting Virus Protein database for blast.")
        command = "makeblastdb -dbtype prot -in " + ncbiVirusProteins
        os.system(command)
        print ("Formatting complete.")
    except:
        print ("Formatting of Virus Protein database for blast failed.")

    os.chdir(cwd)  # Return to home directory

##############################################################################
# Install REFSEQ PROTEIN database

if REFSEQ_PROTEIN:
    filename_root = "refseq_protein"
    os.chdir(refseqProteinDir)
    try:
        print ("Downloading NCBI Refseq Protein database.")
        print ("This may take a while...")
        command = blastPath + "update_blastdb.pl" + ' ' + "refseq_protein"
        success = os.system(command)
        print ("NCBI Refseq Protein database download complete.")
    except BlastError:  
        print ("Command " + command + " failed; please check the location of your blast executables")

    print ("Unpacking files...")
    try:
        command = "ls > ls.out"
        success = os.system(command)
        ls_h = open("ls.out",'r')
        files = ls_h.read().splitlines()
        for filename in files:
            if not re.search('md5', filename):
                command = "gunzip -c " + filename + " | tar xopf -"
                success = os.system(command)
        ls_h.close()
    except Exception:
        print ("Error encountered in unpacking files")

    # Clean up direcgtory
    try:
        print ("Cleaning up directory")
        command = "rm ls.out"
        os.system(command)
    except:
        print ("WARNING: Could not remove unneeded file.")

    os.chdir(cwd)

##############################################################################
# Install REFSEQ GENE database
#*** REFSEQ GENE is no longer downloadable using update_blastdb.pl
#*** NEED TO FIND A SUBSTITUTE COMPRISING VIRUS OR PHAGE GENES
"""
if REFSEQ_GENE:
    os.chdir(refseqGeneDir)
    try:
        print ("Downloading NCBI Refseq Gene database.")
        print ("This may take a while...")
        command = blastPath + "update_blastdb.pl" + ' ' + "refseqgene"
        success = os.system(command)
        print ("NCBI Refseq Gene database download complete.")
    except BlastError:
        print ("Command " + command + " failed; please check the location of your blast executables")

    print ("Unpacking files...")
    try:
        command = "ls > ls.out"
        success = os.system(command)
        ls_h = open("ls.out",'r')
        files = ls_h.read().splitlines()
        for filename in files:
            if not re.search('md5', filename):
                command = "gunzip -c " + filename + " | tar xopf -"
                success = os.system(command)
        ls_h.close()
    except Exception:
        print ("Error encountered in unpacking files")
    os.chdir(cwd)
"""
##############################################################################
# Install SWISSPROT database

if SWISSPROT:
    os.chdir(swissprotDir)
    try:
        print ("Downloading Swissprot database.")
        print ("This may take a while...")
        command = blastPath + "update_blastdb.pl" + ' ' + "swissprot"
        success = os.system(command)
        print ("Swissprot database download complete.")
    except BlastError:
        print ("Command " + command + " failed; please check the location of your blast executables")

    print ("Unpacking files...")
    try:
        command = "ls > ls.out"
        success = os.system(command)
        ls_h = open("ls.out",'r')
        files = ls_h.read().splitlines()
        for filename in files:
            if not re.search('md5', filename):
                command = "gunzip -c " + filename + " | tar xopf -"
                success = os.system(command)
        ls_h.close()
    except Exception:
        print ("Error encountered in unpacking files")

    os.chdir(cwd)

##############################################################################
# Install NR database

if NR:
    os.chdir(nrDir)
    try:
        print ("Downloading NR database.")
        print ("This may take a long time...")
        command = blastPath + "update_blastdb.pl" + ' ' + "nr"
        success = os.system(command)
        print ("NR database download complete.")
    except BlastError:
        print ("Command " + command + " failed; please check the location of your blast executables")

    print ("Unpacking files...")
    try:
        command = "ls > ls.out"
        success = os.system(command)
        ls_h = open("ls.out",'r')
        files = ls_h.read().splitlines()
        for filename in files:
            if not re.search('md5', filename):
                command = "gunzip -c " + filename + " | tar xopf -"
                success = os.system(command)
        ls_h.close()
    except Exception:
        print ("Error encountered in unpacking files")

    os.chdir(cwd)

##############################################################################
# Install VOG databases
# 
if VOGS:
    os.chdir(VOGsDir)

    # First, download the sequences
    OK2FORMAT4VOG = False
    MEMBERS = False; ANNOTATIONS = False
    GENES   = False; PROTEINS    = False
    CATEGORIES  = False
    try:
        print ("Downloading and unzipping VOGS database files.")
        try:
            command = 'wget -O ' + vogMembers_filename + ' "' + VOG_DOWNLOAD_URL + vogMembers_filename + '"' 
            #success = os.system(command)
            command = 'gunzip ' + vogMembers_filename
            #success = os.system(command)
            MEMBERS = True
        except:
            print ("WARNING: Download of vog members file unsuccessful")

        try:
            command = 'wget -O ' + vogAnnotations_filename + ' "' + VOG_DOWNLOAD_URL + vogAnnotations_filename + '"'
            #success = os.system(command)
            command = 'gunzip ' + vogAnnotations_filename
            #success = os.system(command)
            ANNOTATIONS = True
        except:
            print ("WARNING: Download of vog annotations file unsuccessful")

        try:
            command = 'wget -O ' + vogGenes_filename + ' "' + VOG_DOWNLOAD_URL + vogGenes_filename + '"'
            #success = os.system(command)
            command = 'gunzip ' + vogGenes_filename
            #success = os.system(command)
            GENES = True
        except:
            print ("WARNING: Download of vog genes fasta file unsuccessful")

        try:
            command = 'wget -O ' + vogProteins_filename + ' "' + VOG_DOWNLOAD_URL + vogProteins_filename + '"'
            #success = os.system(command)
            command = 'gunzip ' + vogProteins_filename
            #success = os.system(command)
            PROTEINS = True
        except:
            print ("WARNING: Download of vog proteins fasta file unsuccessful")

        try:
            command = 'wget -O ' + vogFunctionalCategories_filename + ' "' + VOG_DOWNLOAD_URL + vogFunctionalCategories_filename + '"'
            #success = os.system(command)
            CATEGORIES = True
        except:
            print ("WARNING: Download of vog members file unsuccessful")

        if MEMBERS and ANNOTATIONS and GENES and PROTEINS and CATEGORIES:
            print ("VOGS database files download complete.")
            OK2FORMAT4VOG = True
        else:
            print ("WARNING: Not all data files were successfully downloaded.")
    except:
        print("Downloading of VOG file(s) failed") 

    # Next, reformat sequence headers with VOG identifiers
    OK2FORMAT4BLAST = False
    os.chdir(cwd)
    if OK2FORMAT4VOG:
        try:
            # The dbPrep_vogTagFastas.py code will tag gene and protein sequence headers
            print ("Reformatting sequence headers with VOG identifiers.")
            print ("This may take a long time (perhap 8 hours).")
            try:
                command = "python dbPrep_vogTagFastas.py " + VOGsDir 
                #result = os.system(command)
                OK2FORMAT4BLAST = True
            except:
                command = "python3 dbPrep_vogTagFastas.py " + VOGsDir
                #result = os.system(command)
                OK2FORMAT4BLAST = True
        except:
            print ("WARNING: Reformatting of sequence headers unsuccessful.")
    os.chdir(VOGsDir)

    # Last, format the sequences for blast
    if OK2FORMAT4BLAST:
        try:
            print ("Formatting VOGs protein database for blast.")
            command = blastPath + "makeblastdb -dbtype prot -in vog.proteins.tagged.all.fa"
            success = os.system(command)
            print ("done")
        except BlastError:
            print ("Command " + command + " failed; please check the location of your blast executables")

        try:
            print ("Formatting VOGs gene database for blast.")
            command = blastPath + "makeblastdb -dbtype nucl -in vog.genes.tagged.all.fa"
            success = os.system(command)
            print ("done")
        except BlastError:
            print ("Command " + command + " failed; please check the location of your blast executables")

    os.chdir(cwd)

##### Install VOG HMMS

if VOG_HMMS:
    os.chdir(VOGhmmsDir)
    try:
        command = 'wget -O ' + vogHmms_filename + ' "' + VOG_DOWNLOAD_URL + vogHmms_filename + '"'
        success = os.system(command)
        command = 'gunzip ' + vogHmms_filename
        success = os.system(command)
        command = 'tar -xvf ' + vogHmms_filename_tar
        success = os.system(command)
        HMMS = True
    except:
        print ("WARNING: Download of vog hmms file unsuccessful")

    if HMMS:
        try:
            print ("Concatenating individual hmm files into one.")
            command = "ls > ./fileList.out"
            success = os.system(command)
            FILELIST_H = open("fileList.out",'r')
            fLines = FILELIST_H.read().splitlines()
            for fLine in fLines:
                if os.path.isfile(fLine):
                    print ("Processing fLine ",fLine)
                    match = re.search("VOG",fLine)
                    if match:
                        command = "cat " + fLine + " >> VOGsHmmProfilesDB.hmm"
                        success = os.system(command)
                        command = "rm " + fLine
                        success = os.system(command)
                else:
                    print ("Next item is not a file:",fLine)
            print ("Concatenation is complete.")
        except:
            print ("Concatenation failed.")

        try:
            print ("Formatting VOG HMM file.")
            command = "hmmpress VOGsHmmProfilesDB.hmm"
            success = os.system(command)
            print ("Formatting complete.")
        except:
            command = "hmmconvert VOGsHmmProfilesDB.hmm VOGsHmmProfilesDB.hmm"
            os.system(command)
            command = "hmmpress VOGsHmmProfilesDB.hmm"
            success = os.system(command)
            print ("Formatting complete.")

        # Clean up the directory
        try:
            print ("Cleaning up unneeded files from the directory.")
            command = "rm fileList.out"
            success = os.system(command)
            command = "rm vog.hmm.tar"
            success = os.system(command)
            command = "rm VOGsHmmProfilesDB.hmm"
            success = os.system(command)
            print ("Directory cleanup completed.")
        except:
            print ("WARNING: Something went wrong cleaning up the VOG HMM directory.")

    os.chdir(cwd)

##############################################################################
# Prepare PVOG HMMs

if PVOG_HMMS:
    os.chdir(pVOGhmmsDir)

    HMMS = False
    try:
        print ("Downloading PVOG HMM profiles file.")
        command = 'wget "' + pvogHmm_httpAddr + '"'
        success = os.system(command)
        print ("gunzip")
        command = 'gunzip AllvogHMMprofiles.tar.gz'
        success = os.system(command)
        print ("tar")
        command = 'tar -xvf AllvogHmmprofiles.tar'
        success = os.system(command)
        HMMS = True
    except:
        print ("Download of PVOG HMM Profiles database unsuccessful.")
  
    if HMMS:
        pVOGhmmsSubDir = "AllvogHMMprofiles/"
        os.chdir(pVOGhmmsSubDir)
        try:
            print ("Concatenating individual hmm files into one.")
            command = "ls > ./fileList.out"
            success = os.system(command)
            FILELIST_H = open("fileList.out",'r')
            fLines = FILELIST_H.read().splitlines()
            for fLine in fLines:
                if os.path.isfile(fLine):
                    print ("Processing fLine ",fLine)
                    match = re.search("VOG",fLine)
                    if match:
                        command = "cat " + fLine + " >> PVOGsHmmProfilesDB.hmm"
                        success = os.system(command)
                        command = "rm " + fLine
                        success = os.system(command)
                else:
                    print ("Next item is not a file:",fLine)
            print ("Concatenation is complete.")
        except:
            print ("Concatenation failed.")

        try:
            print ("Formatting PVOG HMM file.")
            command = "hmmpress PVOGsHmmProfilesDB.hmm"
            success = os.system(command)
            print ("Formatting complete.")
        except:  # Could be older format
            command = "hmmconvert PVOGsHmmProfilesDB.hmm PVOGsHmmProfilesDB.hmm"
            os.system(command)
            command = "hmmpress PVOGsHmmProfilesDB.hmm"
            success = os.system(command)
            print ("Formatting complete.")

        # Clean up the directory
        try:
            print ("Cleaning up unneeded files from the directory.")
            command = "rm fileList.out"
            success = os.system(command)
            command = "rm PVOGsHmmProfilesDB.hmm"
            success = os.system(command)
            print ("Directory cleanup completed.")
        except:
            print ("WARNING: Something went wrong cleaning up the VOG HMM directory.")

        os.chdir(pVOGhmmsDir)

        # Clean up directory
        try:
            print ("Removing file(s) no longer needed.")
            command = "rm AllvogHMMprofiles.tar"
            os.system(command)
            print ("Cleanup successful.")
        except:
            print ("WARNING: Could not clean up pVOGhmms dir")

    os.chdir(cwd)

##############################################################################
# Run makeblastdb on the Phantome and pVOGs databases

if PHANTOME:
    os.chdir(phantomeDir)
    try:
        print ("Formatting Phantome database for blast.")
        command = blastPath + "makeblastdb -dbtype prot -in Phantome_Phage_genes.faa"
        success = os.system(command)
        print ("done")
    except BlastError:
        print ("Command " + command + " failed; please check the location of your blast executables")
    os.chdir(cwd)

if PVOGS:
    os.chdir(pVOGsDir)
    try:
        print ("Formatting pVOGs database for blast.")
        command = blastPath + "makeblastdb -dbtype prot -in pVOGs.faa"
        success = os.system(command)
        print ("done")
    except BlastError:
        print ("Command " + command + " failed; please check the location of your blast executables")
    os.chdir(cwd)

#############################################################################

print ("Done!")

##############################################################################
##############################################################################

#!/usr/bin/env python

#############################################################################
#
# program: getDBs.py
#
# programmer:  C. E. Zhou
#
# Summary:  This script facilitates the downloading of databases to be used with multiPhATE.
#
# Most recent update:  25 March 2019
#
##############################################################################

import os, sys, re, time, datetime
from ftplib import FTP

##############################################################################
# CONSTANTS, BOOLEANS

BLAST               = False
NCBI_VIRUS_GENOME   = False
NCBI_VIRUS_PROTEIN  = False
NCBI_REFSEQ_PROTEIN = False
NCBI_REFSEQ_GENE    = False
NCBI_SWISSPROT      = False
NR                  = False
PHANTOME            = True
PVOGS               = True

# VARIABLES
blast               = ''
ncbi_virus_genome   = ''
ncbi_virus_protein  = ''
ncbi_refseq_protein = ''
ncbi_refseq_gene    = ''
ncbi_swissprot      = ''
nr                  = ''
decision            = ''
blastPath           = '' # path to user's blast installation
cwd                 = '' # current working direcgtory
emailAddr           = '' # user's email address for ftp login

# Set up database directories
cwd              = os.getcwd()
dbDir            = cwd       + "/Databases/"
ncbiDir          = dbDir     + "NCBI/"
ncbiGenomeDir    = ncbiDir   + "Virus_Genome/"
ncbiProteinDir   = ncbiDir   + "Virus_Protein/"
nrDir            = dbDir     + "NR/"
refseqDir        = dbDir     + "Refseq/"
refseqProteinDir = refseqDir + "Protein/"
refseqGeneDir    = refseqDir + "Gene/"
swissprotDir     = dbDir     + "Swissprot/"
phantomeDir      = dbDir     + "Phantome/"
pVOGsDir         = dbDir     + "pVOGs/"

if not os.path.exists(dbDir):
    os.mkdir(dbDir)
if not os.path.exists(ncbiDir):
    os.mkdir(ncbiDir)
if not os.path.exists(ncbiGenomeDir):
    os.mkdir(ncbiGenomeDir)
if not os.path.exists(ncbiProteinDir):
    os.mkdir(ncbiProteinDir)
if not os.path.exists(refseqDir):
    os.mkdir(refseqDir)
if not os.path.exists(refseqProteinDir):
    os.mkdir(refseqProteinDir)
if not os.path.exists(refseqGeneDir):
    os.mkdir(refseqGeneDir)
if not os.path.exists(swissprotDir):
    os.mkdir(swissprotDir)
if not os.path.exists(nrDir):
    os.mkdir(nrDir)

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
    blastPath = input()
    # CHECK THAT THIS IS CORRECT !!!"
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
    print ("That was not a correct response; please try again: type 'y' or 'n'")

##### NCBI Virus Protein database is downloaded via ftp 
print ("NCBI Virus Protein database: ('y'/'n')")
ncbi_virus_protein = input()
if re.search('Y|y|yes|Yes|YES',ncbi_virus_protein):
    print ("Great, let's download NCBI Virus Protein")
    NCBI_VIRUS_PROTEIN = True
elif re.search('N|n|no|No|NO',ncbi_virus_genome):
    print ("Ok, we'll skip that one")
else:
    print ("That was not a correct response; please try again: type 'y' or 'n'")

#### If downloading via FTP, get user's email address for login
if NCBI_VIRUS_GENOME or NCBI_VIRUS_PROTEIN:
    print ("Please enter your email address for ftp anonymous login: ")
    emailAddr = input()

#####
print ("NCBI Refseq Protein database: ('y'/'n')")
ncbi_refseq_protein = input()
if re.search('Y|y|yes|Yes|YES',ncbi_refseq_protein):
    print ("Great, let's download NCBI Virus Genome")
    NCBI_REFSEQ_PROTEIN = True
elif re.search('N|n|no|No|NO',ncbi_refseq_protein):
    print ("Ok, we'll skip that one")
else:
    print ("That was not a correct response; please try again: type 'y' or 'n'")

#####
print ("NCBI Refseq Gene database: ('y'/'n')")
ncbi_refseq_gene = input()
if re.search('Y|y|yes|Yes|YES',ncbi_refseq_gene):
    print ("Great, let's download NCBI Virus Genome")
    NCBI_REFSEQ_GENE = True
elif re.search('N|n|no|No|NO',ncbi_refseq_gene):
    print ("Ok, we'll skip that one")
else:
    print ("That was not a correct response; please try again: type 'y' or 'n'")

#####
print ("NCBI Swissprot database: ('y'/'n')")
ncbi_swissprot = input()
if re.search('Y|y|yes|Yes|YES',ncbi_swissprot):
    print ("Great, let's download NCBI Virus Genome")
    NCBI_SWISSPROT = True
elif re.search('N|n|no|No|NO',ncbi_swissprot):
    print ("Ok, we'll skip that one")
else:
    print ("That was not a correct response; please try again: type 'y' or 'n'")

#####
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
    print ("That was not a correct response; please try again: type 'y' or 'n'")

#####
if (BLAST or NCBI_VIRUS_GENOME or NCBI_REFSEQ_PROTEIN or NCBI_REFSEQ_GENE or NCBI_SWISSPROT or NR):
    print ("Ok, this is what we are going to download. ")
    if not BLAST:
        print ("BLAST plus")
    if NCBI_VIRUS_GENOME:
        print ("NCBI Virus Genome database")
    if NCBI_REFSEQ_PROTEIN:
        print ("NCBI Refseq Protein database")
    if NCBI_REFSEQ_GENE:
        print ("NCBI Refseq Gene database")
    if NCBI_SWISSPROT:
        print ("NCBI Swissprot database")
    if NR:
        print ("NR database")
    print ("Type 'go' to proceed, or 'stop' to reconsider. ")
    print ("(You can always run this script again.) ")
    decision = input()
    if (re.search('go|Go|GO',decision)):
        print ("Ok, let's get started with downloading.")
        print ("Databases will be installed into the Databases/ folder within multiPhATE")
    else:
        print ("Ok, maybe some other time. Bye!")
        exit()
else:
    print ("You have selected no downloads. Have a happy day! :-)")
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

ftpAddr  = "ftp.ncbi.nlm.nih.gov"
file1_1u = "viral.1.1.genomic.fna"
file1_2u = "viral.2.1.genomic.fna"
file1_1  = file1_1u + ".gz"
file1_2  = file1_2u + ".gz"
file2    = "nucl_gb.accession2taxid.gz"

if NCBI_VIRUS_GENOME:
    os.chdir(ncbiGenomeDir)

    ### Download accessory file
    print ("Downloading NCBI accession2taxid accessory file for Virus Genome database.")
    print ("This may take a while...")

    try:
        ftp = FTP(ftpAddr)
        ftp.login(user='anonymous', passwd=emailAddr)
        ftp.cwd('pub/')
        ftp.cwd('taxonomy/')
        ftp.cwd('accession2taxid/')
        localfile = open(file2, 'wb')
        ftp.retrbinary('RETR ' + file2, localfile.write, 1024)
        localfile.close()
        print ("done")
        ftp.quit()
    except Exception:
        print ("FTP download for accessory file failed")

    ### Download two volumes of NCBI Virus Genome database
    print ("Downloading two NCBI Virus Genome database volumes.")
    print ("This may take a while...")

    # Volume 1
    try:
        print ("Volume 1...")
        ftp = FTP(ftpAddr)
        ftp.login(user='anonymous', passwd=emailAddr)
        ftp.cwd('refseq/')
        ftp.cwd('release/')
        ftp.cwd('viral/')
        localfile = open(file1_1, 'wb')
        ftp.retrbinary('RETR ' + file1_1, localfile.write, 1024)
        localfile.close()
        ftp.quit()
        print ("done")
    except Exception:
        print ("FTP download for Volume 1 file failed")
   
    # Volume 2
    try: 
        print ("Volume 2...")
        ftp = FTP(ftpAddr)
        ftp.login(user='anonymous', passwd=emailAddr)
        ftp.cwd('refseq/')
        ftp.cwd('release/')
        ftp.cwd('viral/')
        localfile = open(file1_2, 'wb')
        ftp.retrbinary('RETR ' + file1_2, localfile.write, 1024)
        localfile.close()
        ftp.quit()
        print ("done")
    except Exception:
        print ("FTP download for Volume 2 file failed")

    ### Unpack files

    print ("Unpacking files...")
    try:
        command = "gunzip " + file2
        success = os.system(command)
    except Exception:
        print ("Unpacking of accessory file failed")

    try:
        command = "gunzip " + file1_1
        success = os.system(command)
    except Exception:
        print ("Unpacking of Volume 1 failed")

    try:
        command = "gunzip " + file1_2
        success = os.system(command)
    except Exception:
        print ("Unpacking of Volume 2 failed")

    print ("done")

    ### Format files for blasting

    print ("Performing blast formatting...")

    try:
        command = blastPath + "makeblastdb -dbtype nucl -in " + file1_1u
        success = os.system(command)
    except Exception:
        print ("Formatting of Volume 1 failed")

    try:
        command = blastPath + "makeblastdb -dbtype nucl -in " + file1_2u
        success = os.system(command)
    except Exception:
        print ("Formatting of Volume 2 failed")

    print ("done")
    os.chdir(cwd)

##############################################################################

if NCBI_VIRUS_PROTEIN:

    ftpAddr  = "ftp.ncbi.nlm.nih.gov"
    file3_1u = "viral.1.protein.faa"
    file3_2u = "viral.2.protein.faa"
    file3_1  = file3_1u + ".gz"
    file3_2  = file3_2u + ".gz"

    os.chdir(ncbiProteinDir)

    ### Download two volumes of NCBI Virus Protein database
    print ("Downloading two NCBI Virus Protein database volumes.")
    print ("This may take a while...")

    # Volume 1
    try:
        print ("Volume 1...")
        ftp = FTP(ftpAddr)
        ftp.login(user='anonymous', passwd=emailAddr)
        ftp.cwd('refseq/')
        ftp.cwd('release/')
        ftp.cwd('viral/')
        localfile = open(file3_1, 'wb')
        ftp.retrbinary('RETR ' + file3_1, localfile.write, 1024)
        localfile.close()
        ftp.quit()
        print ("done")
    except Exception:
        print ("FTP download for Volume 1 file failed")

    # Volume 2
    try:
        print ("Volume 2...")
        ftp = FTP(ftpAddr)
        ftp.login(user='anonymous', passwd=emailAddr)
        ftp.cwd('refseq/')
        ftp.cwd('release/')
        ftp.cwd('viral/')
        localfile = open(file3_2, 'wb')
        ftp.retrbinary('RETR ' + file3_2, localfile.write, 1024)
        localfile.close()
        ftp.quit()
        print ("done")
    except Exception:
        print ("FTP download for Volume 2 file failed")

    ### Unpack files

    print ("Unpacking files...")

    try:
        command = "gunzip " + file3_1
        success = os.system(command)
    except Exception:
        print ("Unpacking of Volume 1 failed")

    try:
        command = "gunzip " + file3_2
        success = os.system(command)
    except Exception:
        print ("Unpacking of Volume 2 failed")

    print ("done")

    ### Format files for blasting

    print ("Performing blast formatting...")

    try:
        command = blastPath + "makeblastdb -dbtype prot -in " + file3_1u
        success = os.system(command)
    except Exception:
        print ("Formatting of Volume 1 failed")

    try:
        command = blastPath + "makeblastdb -dbtype nucl -in " + file3_2u
        success = os.system(command)
    except Exception:
        print ("Formatting of Volume 2 failed")

##############################################################################
# Install NCBI_REFSEQ_PROTEIN

if NCBI_REFSEQ_PROTEIN:
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

    os.chdir(cwd)

##############################################################################
# Install NCBI_REFSEQ_GENE

if NCBI_REFSEQ_GENE:
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

##############################################################################
# Install SWISSPROT

if NCBI_SWISSPROT:
    os.chdir(swissprotDir)
    try:
        print ("Downloading Swissprot database.")
        print ("This may take a while...")
        command = blastPath + "update_blastdb.pl" + ' ' + "swissprot"
        success = os.system(command)
        print ("NCBI Swissprot database download complete.")
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
# Install NR

if NR:
    os.chdir(nrDir)
    try:
        print ("Downloading NCBI NR database.")
        print ("This may take a while...")
        command = blastPath + "update_blastdb.pl" + ' ' + "nr"
        success = os.system(command)
        print ("NCBI NR database download complete.")
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

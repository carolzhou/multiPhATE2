#!/usr/bin/env python

#############################################################################
#
# program: dbPrep_getDBs.py
#
# programmer:  C. E. Zhou
#
# Programmer's Notes:
# 1) update_blastdb.pl --showall <to see which databases are available>
#
# Summary:  This script facilitates the downloading of databases to be used with multiPhATE.
#
# Most recent update:  08 August 2020
#
##############################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.pdf FOR DETAILS.

import os, sys, re, time, datetime
#from ftplib import FTP
from pathlib import Path

##############################################################################
# CONSTANTS, BOOLEANS

# Run interactive for user to input which data sets they want to install
INTERACTIVE = False
# Run remote to pre-set download instructions and skip user input
REMOTE = True
# Set verbost to true for remote processing if server is killing idle processes.
# Verbose will write voluminous progress to console, keeping user process non-idle during long computations.
VERBOSE = False 

CODE_BASE = "dbPrep_getDBs"
CODE_NAME = CODE_BASE + ".py"

# OPEN LOG FILES
# Open Log file
logFile = './' + CODE_BASE + '.log'
LOG_H = open(logFile,'a')

# Begin logging
myTime = datetime.datetime.now()
LOG_H.write("%s%s%s\n" % (CODE_NAME," begin processing at ",myTime))
LOG_H.flush()

# Set environment variable to subordinate code knows to print progress messages
if VERBOSE:
    os.environ["dbPrep_VERBOSE"] = 'True'
else:
    os.environ["dbPrep_VERBOSE"] = 'False'

##### DATABASES: URLs and Files
# The following URLs etc may be modified as things change on 3rd party web sites.

# VOG Databases 
VOG_VERSION      = "99"  # Most recent version as of this writing; modify as appropriate
VOG_DOWNLOAD_URL = "http://fileshare.csb.univie.ac.at/vog/vog" + VOG_VERSION + "/"

# NCBI Databases
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

# Swissprot Database
swissprotUniprot_httpAddr = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
swissprotDBfile_gz        = "uniprot_sprot.fasta.gz"
swissprotDBfile           = "uniprot_sprot.fasta"
swissprotFastaFilename    = "swissprot.fa"

# PVOG Database
# PVOG data are available at the authors' home page and at NCBI. Choose one or the other:
#pvogHmm_httpAddr    = "http://dmk-brain.ecn.uiowa.edu/pVOGs/downloads/All/AllvogHMMprofiles.tar.gz"
#pvogAccns_httpAddr  = "http://dmk-brain.ecn.uiowa.edu/pVOGs/downloads/All/AllFalmilyProteinList.tsv"
pvogHmm_httpAddr    = "https://ftp.ncbi.nlm.nih.gov/pub/kristensen/pVOGs/downloads/All/AllvogHMMprofiles.tar.gz"
pvogAccns_httpAddr  = "https://ftp.ncbi.nlm.nih.gov/pub/kristensen/pVOGs/downloads/All/AllFamilyProteinList.tsv"

# CAZy Database http://bcb.unl.edu/dbCAN2/about.php
dbCAN2_httpAddr     = "http://bcb.unl.edu/dbCAN2/download/Databases/"
cazyProteins        = "CAZyDB.07312019.fa"
cazyFamActivities   = "CAZyDB.07312019.fam-activities.txt"
cazyFamSubfamEC     = "CAZyDB.07312019.fam.subfam.ec.txt"
cazyFamInfo         = "FamInfo.txt.04232020.xls"
cazyPrWithEC        = "CAZyDB.07312018.pr-with-ec.txt"
cazyProteins_httpAddr      = os.path.join(dbCAN2_httpAddr,cazyProteins)
cazyFamActivities_httpAddr = os.path.join(dbCAN2_httpAddr,cazyFamActivities)
cazyFamSubfamEC_httpAddr   = os.path.join(dbCAN2_httpAddr,cazyFamSubfamEC)
cazyFamInfo_httpAddr       = os.path.join(dbCAN2_httpAddr,cazyFamInfo)
cazyPrWithEC_httpAddr      = os.path.join(dbCAN2_httpAddr,cazyPrWithEC)

# DATABASES: Boolean control
BLAST               = False
NCBI_VIRUS_GENOME   = False
NCBI_VIRUS_PROTEIN  = False
REFSEQ_PROTEIN      = False
SWISSPROT           = False
NR                  = False
PHANTOME            = False 
PVOGS               = False 
PVOG_HMMS           = False
VOGS                = False
VOG_HMMS            = False
CAZY                = False

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
cazy                = ''
decision            = ''
blastPath           = '' # path to user's blast installation
cwd                 = '' # current working direcgtory
emailAddr           = '' # user's email address for ftp login

# Set up database directories
cwd                    = os.getcwd()  # current working directory
# Get directory one level up; this is where the Databases directory will go
baseDir                = Path(__file__).resolve().parents[1]  # calculate dir one dir up from cwd
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
CAZyDir                = os.path.join(dbDir,          "CAZY/")

# Additional filenames from source
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
ncbiVirusGenomes                 = "ncbiVirusGenomes.fasta"
ncbiVirusProteins                = "ncbiVirusProteins.faa"
refseqProteins                   = "refseq_protein"
swissprotProteins                = "swissprot"
nrProteins                       = "nr"
phantomeProteins                 = "Phantome_Phage_Genes.faa"
pvogProteins                     = "pVOGs.faa"
pvogProteinHmms                  = "AllvogHMMprofiles/PVOGsHmmProfilesDB.hmm"
vogGenes                         = "vog.genes.tagged.all.fa"
vogProteins                      = "vog.proteins.tagged.all.fa"
vogProteinHmms                   = "VOGsHmmProfilesDB.hmm"

# Create database directories, if they don't already exist
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
if not os.path.exists(CAZyDir):
    os.mkdir(CAZyDir)

# Calculate database path/files for user's config file
ncbiGenome_config       = os.path.join(ncbiGenomeDir,    ncbiVirusGenomes)
ncbiProtein_config      = os.path.join(ncbiProteinDir,   ncbiVirusProteins)
refseqProtein_config    = os.path.join(refseqProteinDir, refseqProteins)
swissprotProtein_config = os.path.join(swissprotDir,     swissprotProteins)
nr_config               = os.path.join(nrDir,            nrProteins)
phantome_config         = os.path.join(phantomeDir,      phantomeProteins)
pvog_config             = os.path.join(pVOGsDir,         pvogProteins)
pvogHmm_config          = os.path.join(pVOGhmmsDir,      pvogProteinHmms)
vogGenes_config         = os.path.join(VOGsDir,          vogGenes)
vogProteins_config      = os.path.join(VOGsDir,          vogProteins)
vogProteinHmms_config   = os.path.join(VOGhmmsDir,       vogProteinHmms)
cazyProteins_config     = os.path.join(CAZyDir,          cazyProteins)
cazyAnnotation_config   = os.path.join(CAZyDir,          cazyFamActivities)

# Write database path/files to file for user's benefit in constructing config file
dataFile = './' + CODE_BASE + '.lst'
DATA_H = open(dataFile,'w')
DATA_H.write("%s\n" % (ncbiGenome_config))
DATA_H.write("%s\n" % (ncbiProtein_config))
DATA_H.write("%s\n" % (refseqProtein_config))
DATA_H.write("%s\n" % (swissprotProtein_config))
DATA_H.write("%s\n" % (nr_config))
DATA_H.write("%s\n" % (phantome_config))
DATA_H.write("%s\n" % (pvog_config))
DATA_H.write("%s\n" % (pvogHmm_config))
DATA_H.write("%s\n" % (vogGenes_config))
DATA_H.write("%s\n" % (vogProteins_config))
DATA_H.write("%s\n" % (vogProteinHmms_config))
DATA_H.write("%s\n" % (cazyProteins_config))
DATA_H.write("%s\n" % (cazyAnnotation_config))
DATA_H.close()

# Pre-set download instructions; skip user input
if REMOTE:
    BLAST               = True 
    NCBI_VIRUS_GENOME   = False 
    NCBI_VIRUS_PROTEIN  = False 
    REFSEQ_PROTEIN      = False 
    SWISSPROT           = True 
    NR                  = False
    PHANTOME            = False   # Provided in distribution
    PVOGS               = False   # Provided in distribution
    PVOG_HMMS           = False 
    VOGS                = False 
    VOG_HMMS            = False 
    CAZY                = False 

# Determine download instructions via user input
elif INTERACTIVE:
    ##############################################################################
    # First, determine if user needs to download BLAST+.
    time.sleep(1)
    print ("Welcome to dbPrep_getDBs.py, a code for downloading data sets used by multiPhATE.")
    time.sleep(1)
    print ("We will be downloading data sets one at a time. You will select the data sets")
    time.sleep(1)
    print ("   that you want to download. You don't need to download all of them at once--you ")
    time.sleep(1)
    print ("   may download any or all of the them now, and you can return later to this program to ")
    time.sleep(1)
    print ("   download any data sets that you skip this time around.")
    time.sleep(1)
    print ("Let's get started.")
    time.sleep(1)
    print ("First, you need blast+ in order to install several of the databases")
    time.sleep(1)
    print ("This code does not support legacy blast. :-( ")
    time.sleep(1)
    print ("Please confirm that you have downloaded and installed blast+: type 'y' (yes) or 'n' (no)")
    blast = input()
    if re.search('Y|y|Yes|yes|YES', blast):
        BLAST = True 
        print ("That's great! If you installed blast+ globally, or if you installed it in your")
        time.sleep(1)
        print ("   conda environment, and are running this code in that environment, then we can access")
        time.sleep(1)
        print ("   it. But if not, then you will need to input the location of your blast+ program")
        time.sleep(1)
        print ("   on your computer. If so, please input the fully qualified path to blast+.")
        time.sleep(1)
        print ("   otherwise, just hit return:")
        blastPath = input()
    elif re.search('N|n|No|no|NO', blast):
        BLAST = False 
        print ("Please consult the README file for how to acquire and install BLAST+.") 
        time.sleep(1)
        print ("Note: The easiest way to install BLAST+ is within a Conda environment.")
        time.sleep(1)
        print ("Without blast+, we can still download some of the databases.")
        time.sleep(1)
        print ("Shall we continue? please respond 'y' or 'n': ")
        toContinue = input()
        if re.search('N|n|No|no|NO', toContinue):
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
    time.sleep(1)
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
    time.sleep(1)
    print ("NCBI Virus Protein database: ('y'/'n')")
    ncbi_virus_protein = input()
    if re.search('Y|y|yes|Yes|YES',ncbi_virus_protein):
        print ("Great, let's download NCBI Virus Protein")
        NCBI_VIRUS_PROTEIN = True
    elif re.search('N|n|no|No|NO',ncbi_virus_protein):
        print ("Ok, we'll skip that one")
    else:
        print ("That was not a correct response; please run this script again to download the database.")

    ##### REFSEQ PROTEIN
    time.sleep(1)
    if BLAST:
        print ("Refseq Protein database: ('y'/'n')")
        refseq_protein = input()
        if re.search('Y|y|yes|Yes|YES',refseq_protein):
            print ("Great, let's download the Refseq Protein database")
            REFSEQ_PROTEIN = True
        elif re.search('N|n|no|No|NO',refseq_protein):
            print ("Ok, we'll skip that one")
        else:
            print ("That was not a correct response; please run this script again to download the database.")

    ##### SWISSPROT
    time.sleep(1)
    if BLAST:
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
    time.sleep(1)
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
    time.sleep(1)
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
    time.sleep(1)
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
    time.sleep(1)
    if BLAST:
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

    ##### CAZY
    time.sleep(1)
    if BLAST:
        print ("CAZy database: ('y'/'n')")
        cazy = input()
        if re.search('Y|y|yes|Yes|YES',swissprot):
            print ("Great, let's download the CAZy database")
            CAZY = True
        elif re.search('N|n|no|No|NO',cazy):
            print ("Ok, we'll skip that one")
        else:
            print ("That was not a correct response; please run this script again to download the database.")

    ##### PHANTOME
    time.sleep(1)
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
    if (NCBI_VIRUS_GENOME or NCBI_VIRUS_PROTEIN or REFSEQ_PROTEIN or SWISSPROT or VOGS or VOG_HMMS or PVOG_HMMS or NR or PVOGS or CAZY or PHANTOME):

        if (NCBI_VIRUS_GENOME or NCBI_VIRUS_PROTEIN or REFSEQ_PROTEIN or SWISSPROT or VOGS or VOG_HMMS or PVOG_HMMS or NR or CAZY):
            time.sleep(1)
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
            if CAZY:
                print ("CAZy database")

        if (PVOGS or PHANTOME):
            time.sleep(1)
            print ("This is what we are going to format:")
            if PVOGS:
                print ("PVOGS")
            if PHANTOME:
                print ("PHANTOME")

        time.sleep(1)
        print ("\nType 'go' to proceed, or 'stop' to reconsider. ")
        print ("(You can always run this script again.) ")
        decision = input()
        time.sleep(1)
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
##############################################################################

# Install NCBI_VIRUS_GENOME

if NCBI_VIRUS_GENOME:
    myTime = datetime.datetime.now()
    LOG_H.write("%s%s\n" % ("NCBI VIRUS GENOME download begin processing at ",myTime))
    LOG_H.flush()
    os.chdir(ncbiGenomeDownloadDir)

    # Download directories containing virus genome fasta files
    try:
        print ("Downloading NCBI Genome fasta files.")
        print ("This may take a while...")
        command = "wget " + virGenome_httpAddr
        success = os.system(command)
    except:
        print ("WARNING: Command " + command + " failed.")

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
        print ("WARNING: Uncompression of genome fasta file was unsuccessful.")

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
        print ("WARNING: Concatenation failed.")

    # Move consolidated fasta file up one directory
    try:
        print ("Moving ",ncbiVirusGenomes," to ",ncbiGenomeDir," directory")
        command = "mv ./ncbiVirusGenomes.fasta ../."
        success = os.system(command)
    except:
        print ("WARNING: Cannot move ncbiVirusGenomes.fastaile to directory ",ncbiGenomeDir)

    # Before leaving the Download directory, remove the ncbi downloaded file
    try:
        print ("Removing the NCBI Virus Genome compressed file that we downloaded.")
        command = "rm " + virGenomeFile_tar
        success = os.system(command)
        command = "rm fileList.out"
        success = os.system(command)
    except:
        print ("WARNING: Problem removing ",ncbiVirusGenomes," or fileList.out file.")

    os.chdir(ncbiGenomeDir)  # Change to NCBI/Virus_Genome directory

    # Download the accession2taxid file
    try:
        print ("Downloading the accession2taxid file from ",accn2taxid_fileAddr)
        command = "wget " + accn2taxid_fileAddr
        success = os.system(command)
    except:
        print ("WARNING: Download of the accession2taxid file failed.")

    try:
        print ("Unpacking accession2taxid.")
        command = "gunzip " + accn2taxid_file_gz
        success = os.system(command)
        print ("accn2taxid file has been unzipped.")
    except:
        print ("WARNING: Unpacking of accession2taxid failed.")

    # Finally, format the Virus Genome database for blast
    try:
        print ("Formatting Virus Genome database for blast.")
        command = "makeblastdb -dbtype nucl -in " + ncbiVirusGenomes
        success = os.system(command)
        print ("Formatting complete.")
    except:
        print ("WARNING: Formatting of Virus Genome database for blast failed.")

    # And last but not least, remove the Download directory
    try:
        print ("Removing the Download/ directory.")
        command = "rmdir Download/"
        success = os.system(command)
    except:
        print ("WARNING: Cound not remove the Download/ directory.")

    os.chdir(cwd)  # Return to home directory
    myTime = datetime.datetime.now()
    LOG_H.write("%s%s\n" % ("NCBI Virus Genome download finish processing at ",myTime))
    LOG_H.flush()

##############################################################################
# Install NCBI VIRUS PROTEIN database

if NCBI_VIRUS_PROTEIN:

    myTime = datetime.datetime.now()
    LOG_H.write("%s%s\n" % ("NCBI Virus Protein download begin processing at ",myTime))
    LOG_H.flush()
 
    os.chdir(ncbiProteinDownloadDir)

    # Download directories containing virus peptide fasta files
    try:
        print ("Downloading NCBI Virus peptide fasta files.")
        print ("This may take a while...")
        command = "wget " + virProtein_httpAddr
        success = os.system(command)
    except:
        print ("WARNING: Command " + command + " failed.")

    # Uncompress the downloaded file
    try:
        print ("Uncompressing NCBI Virus Protein gz file.")
        command = "gunzip " + virProteinFile_gz
        success = os.system(command)
        command = "tar -xvf " + virProteinFile_tar
        success = os.system(command)
        print ("Virus protein data is now uncompressed.")
    except:
        print ("WARNING: Uncompression of file failed.")

    # Remove the DBV/ directory from download
    # This appears to hold peptides from raw Prokka gene-call predictions
    try:
        print ("Removing the DBV/ subdirectory, which we do not need.")
        command = "rm -r DBV/"
        success = os.system(command)
        print ("DBV/ has been removed.")
    except:
        print ("WARNING: DBV directory removal failed.")

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
        print ("WARNING: Listing and concatenating files failed.")

    # Move consolidated fasta file up one directory
    try:
        print ("Moving ",ncbiVirusProteins," to ",ncbiProteinDir," directory")
        command = "mv ./ncbiVirusProteins.faa ../."
        success = os.system(command)
    except:
        print ("WARNING: Cannot move ncbiVirusProteins.faa file to directory ",ncbiProteinDir)
   
    # Before leaving the Download directory, remove the ncbi downloaded file
    try:
        print ("Removing the NCBI Virus Genome compressed file that we downloaded.")
        command = "rm " + virProteinFile_tar
        success = os.system(command)
        command = "rm fileList.out"
        success = os.system(command)
        command = "rm ls.out"
        success = os.system(command)
    except:
        print ("WARNING: Problem removing ",ncbiVirusProteins," or fileList.out file.")

    os.chdir(ncbiProteinDir)  # Change to NCBI/Virus_Protein directory

    # Finally, format the Virus Protein database for blast
    try:
        print ("Formatting Virus Protein database for blast.")
        command = "makeblastdb -dbtype prot -in " + ncbiVirusProteins
        success = os.system(command)
        print ("Formatting complete.")
    except:
        print ("WARNING: Formatting of Virus Protein database for blast failed.")

    # And last but not least, remove the Download directory
    try:
        print ("Removing the Download/ directory.")
        command = "rmdir Download/"
        success = os.system(command)
    except:
        print ("WARNING: Cound not remove the Download/ directory.")

    os.chdir(cwd)  # Return to home directory
    myTime = datetime.datetime.now()
    LOG_H.write("%s%s\n" % ("NCBI Virus Protein download finish processing at ",myTime))
    LOG_H.flush()

##############################################################################
# Install REFSEQ PROTEIN database

if REFSEQ_PROTEIN:
    myTime = datetime.datetime.now()
    LOG_H.write("%s%s\n" % ("REFSEQ PROTEIN download begin processing at ",myTime))
    LOG_H.flush()
    os.chdir(refseqProteinDir)
 
    filename_root = "refseq_protein"
    try:
        print ("Downloading NCBI Refseq Protein database.")
        print ("This may take a while...")
        command = blastPath + "update_blastdb.pl" + ' ' + "refseq_protein"
        success = os.system(command)
        print ("NCBI Refseq Protein database download complete.")
    except BlastError:  
        print ("WARNING: Command " + command + " failed; please check the location of your blast executables")

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
        print ("WARNING: Error encountered in unpacking files")

    # Clean up direcgtory
    try:
        print ("Cleaning up directory")
        command = "rm ls.out"
        success = os.system(command)
        command = "rm *tar.gz*"
        success = os.system(command)
    except:
        print ("WARNING: Could not remove unneeded files.")

    os.chdir(cwd)
    myTime = datetime.datetime.now()
    LOG_H.write("%s%s\n" % ("REFSEQ PROTEIN download finish processing at ",myTime))
    LOG_H.flush()

##############################################################################
# Install SWISSPROT database

if SWISSPROT:
    myTime = datetime.datetime.now()
    LOG_H.write("%s%s\n" % ("Swissprot download begin processing at ",myTime))
    LOG_H.flush()
 
    os.chdir(swissprotDir)

    try:
        print("Downloading Swissprot fasta sequences")
        command = "wget " + swissprotUniprot_httpAddr
        success = os.system(command)
    except:
        print("WARNING: Download of Swissprot fasta sequences failed.")

    try:
        command = "gunzip " + swissprotDBfile_gz
        success = os.system(command)
    except:
        print("WARNING: Could not unzip swissprot fasta file.")

    # Rename the swissprot fasta files
    try:
        command = "mv " + swissprotDBfile + ' ' + swissprotFastaFilename
        success = os.system(command)
    except:
        print("WARNING: Cound not rename swissprot fasta file.")

    # Download the swissprot blast database (already formatted)
    try:
        print ("Downloading Swissprot blast database.")
        print ("This may take a while...")
        command = blastPath + "update_blastdb.pl" + ' ' + "swissprot"
        success = os.system(command)
        print ("Swissprot database download complete.")
    except BlastError:
        print ("WARNING: Command " + command + " failed; please check the location of your blast executables")

    print ("Unpacking files...")
    try:
        command = "ls > ls.out"
        success = os.system(command)
        ls_h = open("ls.out",'r')
        files = ls_h.read().splitlines()
        for filename in files:
            if not re.search('md5', filename) and not re.search('ls\.out', filename):
                command = "gunzip -c " + filename + " | tar xopf -"
                success = os.system(command)
        ls_h.close()
    except Exception:
        print ("WARNING: Error encountered in unpacking files")

    # Remove unneeded files
    try:
        print ("Removing unneeded files.")
        command = "rm *tar.gz*"
        success = os.system(command)
        command = "rm ls.out"
        success = os.system(command)
    except:
        print ("WARNING: Could not remove unneeded files.")

    os.chdir(cwd)
    myTime = datetime.datetime.now()
    LOG_H.write("%s%s\n" % ("Swissprot download finish processing at ",myTime))
    LOG_H.flush()

##############################################################################
# Install NR database

if NR:
    myTime = datetime.datetime.now()
    LOG_H.write("%s%s\n" % ("NR download begin processing at ",myTime))
    LOG_H.flush()
    os.chdir(nrDir)

    # Download NR database and format for blast
    try:
        print ("Downloading NR database.")
        print ("This may take a long time...")
        command = blastPath + "update_blastdb.pl --decompress nr"
        success = os.system(command)
        print ("NR database download complete.")
    except BlastError:
        print ("WARNING: Command " + command + " failed; please check the location of your blast executables")

    # Clean up
    try:
        print ("Cleaning up directory.")
        command = "rm *tar.gz*"
        success = os.system(command)
    except:
        print ("WARNING: Could not remove unneeded files.")

    myTime = datetime.datetime.now()
    LOG_H.write("%s%s\n" % ("NR download finish processing at ",myTime))
    LOG_H.flush()
 
    os.chdir(cwd)

##############################################################################
# Install VOG databases
# 
if VOGS:
    myTime = datetime.datetime.now()
    LOG_H.write("%s%s\n" % ("Vogs download begin processing at ",myTime))
    LOG_H.flush()
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
            success = os.system(command)
            command = 'gunzip ' + vogMembers_filename
            success = os.system(command)
            MEMBERS = True
        except:
            print ("WARNING: Download of vog members file unsuccessful")

        try:
            command = 'wget -O ' + vogAnnotations_filename + ' "' + VOG_DOWNLOAD_URL + vogAnnotations_filename + '"'
            success = os.system(command)
            command = 'gunzip ' + vogAnnotations_filename
            success = os.system(command)
            ANNOTATIONS = True
        except:
            print ("WARNING: Download of vog annotations file unsuccessful")

        try:
            command = 'wget -O ' + vogGenes_filename + ' "' + VOG_DOWNLOAD_URL + vogGenes_filename + '"'
            success = os.system(command)
            command = 'gunzip ' + vogGenes_filename
            success = os.system(command)
            GENES = True
        except:
            print ("WARNING: Download of vog genes fasta file unsuccessful")

        try:
            command = 'wget -O ' + vogProteins_filename + ' "' + VOG_DOWNLOAD_URL + vogProteins_filename + '"'
            success = os.system(command)
            command = 'gunzip ' + vogProteins_filename
            success = os.system(command)
            PROTEINS = True
        except:
            print ("WARNING: Download of vog proteins fasta file unsuccessful")

        try:
            command = 'wget -O ' + vogFunctionalCategories_filename + ' "' + VOG_DOWNLOAD_URL + vogFunctionalCategories_filename + '"'
            success = os.system(command)
            CATEGORIES = True
        except:
            print ("WARNING: Download of vog members file unsuccessful")

        if MEMBERS and ANNOTATIONS and GENES and PROTEINS and CATEGORIES:
            print ("VOGS database files download complete.")
            OK2FORMAT4VOG = True
        else:
            print ("WARNING: Not all data files were successfully downloaded.")
    except:
        print("WARNING: Downloading of VOG file(s) failed") 

    # Next, reformat sequence headers with VOG identifiers
    OK2FORMAT4BLAST = False
    os.chdir(cwd)
    if OK2FORMAT4VOG:
        try:
            # The dbPrep_vogTagFastas.py code will tag gene and protein sequence headers
            print ("Reformatting sequence headers with VOG identifiers.")
            print ("This may take a long time (perhap hours).")
            try:
                command = "python3 dbPrep_vogTagFastas.py " + VOGsDir 
                success = os.system(command)
                OK2FORMAT4BLAST = True
            except:
                command = "python dbPrep_vogTagFastas.py " + VOGsDir
                success = os.system(command)
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
            print ("WARNING: Command " + command + " failed; please check the location of your blast executables")

        try:
            print ("Formatting VOGs gene database for blast.")
            command = blastPath + "makeblastdb -dbtype nucl -in vog.genes.tagged.all.fa"
            success = os.system(command)
            print ("done")
        except BlastError:
            print ("WARNING: Command " + command + " failed; please check the location of your blast executables")

    os.chdir(cwd)
    myTime = datetime.datetime.now()
    LOG_H.write("%s%s\n" % ("Vogs download finish processing at ",myTime))
    LOG_H.flush()

##### Install VOG HMMS

if VOG_HMMS:
    myTime = datetime.datetime.now()
    LOG_H.write("%s%s\n" % ("Vog Hmms download begin processing at ",myTime))
    LOG_H.flush()
    os.chdir(VOGhmmsDir)

    # Download VOG Hmms, unpack, concatenate, and format for hmm search
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
            print ("WARNING: Concatenation failed.")
       
        # Format VOG HMMs
        try:
            try:
                print ("Formatting VOG HMM file.")
                command = "hmmpress VOGsHmmProfilesDB.hmm"
                success = os.system(command)
                print ("Formatting complete.")
            except:
                command = "hmmconvert VOGsHmmProfilesDB.hmm VOGsHmmProfilesDB.hmm"
                success = os.system(command)
                command = "hmmpress VOGsHmmProfilesDB.hmm"
                success = os.system(command)
                print ("Formatting complete.")
        except:
            print ("WARNING: Formatting of VOG HMM file unsuccessful.")

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
    myTime = datetime.datetime.now()
    LOG_H.write("%s%s\n" % ("Vog Hmms download finish processing at ",myTime))
    LOG_H.flush()

##############################################################################
# Prepare PVOG HMMs

if PVOG_HMMS:
    myTime = datetime.datetime.now()
    LOG_H.write("%s%s\n" % ("pVog Hmms download begin processing at ",myTime))
    LOG_H.flush()
    os.chdir(pVOGhmmsDir)

    # Download Pvog Hmms, unpack, concatenate, and format for hmm search
    HMMS = False
    try:
        print ("Downloading PVOG HMM profiles file.")
        command = 'wget "' + pvogHmm_httpAddr + '"'
        success = os.system(command)
        print ("gunzip")
        command = 'gunzip AllvogHMMprofiles.tar.gz'
        success = os.system(command)
        print ("tar")
        command = 'tar -xvf AllvogHMMprofiles.tar'
        success = os.system(command)
        HMMS = True
    except:
        print ("WARNING: Download of PVOG HMM Profiles database unsuccessful.")
  
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
            print ("WARNING: Concatenation failed.")

        # Format PVOG HMMs
        try:
            try:
                print ("Formatting PVOG HMM file.")
                command = "hmmpress PVOGsHmmProfilesDB.hmm"
                success = os.system(command)
                print ("Formatting complete.")
            except:  # Could be older format
                command = "hmmconvert PVOGsHmmProfilesDB.hmm PVOGsHmmProfilesDB.hmm"
                success = os.system(command)
                command = "hmmpress PVOGsHmmProfilesDB.hmm"
                success = os.system(command)
                print ("Formatting complete.")
        except:
            print ("WARNING: Cound not format PVOG HMM file")

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
            success = os.system(command)
            print ("Cleanup successful.")
        except:
            print ("WARNING: Could not clean up pVOGhmms dir")

    os.chdir(cwd)
    myTime = datetime.datetime.now()
    LOG_H.write("%s%s\n" % ("Pvog Hmms downloading finish processing at ",myTime))
    LOG_H.flush()

################
# Download and format CAZy sequence database and associated files

if CAZY:
    myTime = datetime.datetime.now()
    LOG_H.write("%s%s\n" % ("CAZy download begin processing at ",myTime))
    LOG_H.flush()
    os.chdir(CAZyDir)

    # Download CAZy database files
    try:
        print ("Downloading ",cazyProteins)
        command = 'wget "' + cazyProteins_httpAddr + '"'
        success = os.system(command)
    except:
        print ("WARNING: Could not download ",cazyProteins)

    try:
        print ("Downloading ",cazyFamActivities)
        command = 'wget "' + cazyFamActivities_httpAddr + '"'
        success = os.system(command)
    except:
        print ("WARNING: Could not download ",cazyFamActivities)
    """
    try:
        print ("Downloading ",cazyFamSubfamEC)
        command = 'wget "' + cazyFamSubfamEC_httpAddr + '"'
        success = os.system(command)
    except:
        print ("WARNING: Could not download ",cazyFamSubfamEC)

    try:
        print ("Downloading ",cazyFamInfo)
        command = 'wget "' + cazyFamInfo_httpAddr + '"'
        success = os.system(command)
    except:
        print ("WARNING: Could not download ",cazyFamInfo)

    try:
        print ("Downloading ",cazyPrWithEC)
        command = 'wget "' + cazyPrWithEC_httpAddr + '"'
        success = os.system(command)
    except:
        print ("WARNING: Could not download ",cazyPrWithEC)
    """
    # Format CAZy protein database for Blast
    try:
        print ("Formatting CAZy protein database for blast.")
        command = blastPath + "makeblastdb -dbtype prot -in " + cazyProteins  
        success = os.system(command)
    except:
        print ("WARNING: Could not format CAZy protein database for blast.")

    os.chdir(cwd)
    myTime = datetime.datetime.now()
    LOG_H.write("%s%s\n" % ("CAZy download finish processing at ",myTime))
    LOG_H.flush()
    os.chdir(cwd)

##############################################################################
# Run makeblastdb on the Phantome and pVOGs databases

if PHANTOME:
    myTime = datetime.datetime.now()
    LOG_H.write("%s%s\n" % ("Phantome formatting begin processing at ",myTime))
    LOG_H.flush()
 
    os.chdir(phantomeDir)
    try:
        print ("Formatting Phantome database for blast.")
        command = blastPath + "makeblastdb -dbtype prot -in Phantome_Phage_genes.faa"
        success = os.system(command)
        print ("done")
    except BlastError:
        print ("WARNING: Command " + command + " failed; please check the location of your blast executables")

    os.chdir(cwd)
    myTime = datetime.datetime.now()
    LOG_H.write("%s%s\n" % ("Phantome formatting finish processing at ",myTime))
    LOG_H.flush()

if PVOGS:
    myTime = datetime.datetime.now()
    LOG_H.write("%s%s\n" % ("pVogs formatting begin processing at ",myTime))
    LOG_H.flush()
 
    os.chdir(pVOGsDir)
    try:
        print ("Formatting pVOGs database for blast.")
        command = blastPath + "makeblastdb -dbtype prot -in pVOGs.faa"
        success = os.system(command)
        print ("done")
    except BlastError:
        print ("WARNING: Command " + command + " failed; please check the location of your blast executables")

    os.chdir(cwd)
    myTime = datetime.datetime.now()
    LOG_H.write("%s%s\n" % ("pVogs formatting finish processing at ",myTime))
    LOG_H.flush()

#############################################################################

print ("Done!")
LOG_H.write("%s%s%s\n" % (CODE_NAME," end processing at ",myTime))
LOG_H.flush()
LOG_H.close()
DATA_H.close()

##############################################################################
##############################################################################

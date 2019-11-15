##########################################################################
#
# module phate_pVOG
#
# Description:  Handles data and processing for pVOG data type
#
# Programmer:  C. E. Zhou
#
#########################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.PDF FOR DETAILS.

#DEBUG = True
DEBUG = False 
SYSTEM = 'MAC'     # Controls which method used for system call
#SYSTEM = 'LINUX'

import re
import os
import copy
import subprocess

#class pVOGrecord(object):
    
#    def __init__(self):
#        self.pVOGrecord = {
#            'summary'      : '',
#            'family'       : '',
#            'species'      : '',
#            'genomeAccn'   : '',
#            'peptideAccn'  : '',
#            'location'     : '',
#            'number'       : 0,
#            'genomeCount'  : 0,
#            'peptideCount' : 0,
#            'description'  : '',
#            }
pVOGrecord = {
    'summary'      : '',
    'family'       : '',
    'species'      : '',
    'genomeAccn'   : '',
    'peptideAccn'  : '',
    'location'     : '',
    'number'       : 0,
    'genomeCount'  : 0,
    'peptideCount' : 0,
    'description'  : '',
    }
        
class pVOG(object):
    
    def __init__(self):
        self.pVOGid         = ''  # e.g., VOG0334
        self.accessionList  = []  # list of peptide accession numbers 
        self.comment        = ''  # text from pVOG identifier line; specified peptide & genome counts 
        self.count          = 0
        self.pVOGrecordList = [] # list of pVOGrecords
        #self.pVOGrecordObj  = pVOGrecord()

class pVOGs(object):

    def __init__(self):
        self.databaseName = "pVOGs"
        self.downloadDate = "June 2017"
        self.version      = "unknown"
        self.pVOGlist     = []  # list of pVOG objects 
        self.accnList     = []  # non-redundant list of accessions; for pulling seqs from NR subset
        self.pVOGobj      = pVOG()

    def addPvogs(self,pVOGdir,logFile_h):  # Inputs the directory where pVOG library files reside
        progressCount = 0
        p_pVOG = re.compile('^#(VOG\d+[a|b]?):\shas\s(\d+)[\w\s]+(\d+)\sgenomes?')  # Note: There are VOG1984a and VOG1984b
        # pVOG accession info is tabbed, with the following columns:
        # Family Species GenomeAccn PeptideAccn Location Number Description
        p_pVOGrecord = re.compile('^[^\t]+\t')  # Checks that it's a tabbed field
        #p_pVOGrecord = re.compile('^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)') # ie, not-tabs, tab, not-tabs, tab, etc.
        #p_pVOGrecord = re.compile('^([\w\d\-\_]+)\t([\w\d\-\_\.\s]+)\t([\w\d\-\_\.]+)\t([\w\d\-\_\.]+)\t([\d\.]+)\t(\d+)\t(.*)')
        pVOG_NRlist = []  # ensures that pVOG peptide accession list is non-redundant

        # Go to pVOG directory; list files; read in pVOG lines from each file

        # First, get directory listing of files with pVOG data:  VOG number followed by number of members and members' accessions
        if SYSTEM == 'MAC':
            cwd = os.getcwd()
            os.chdir(pVOGdir)
            fileList = os.listdir(pVOGdir)
            os.chdir(cwd)
        elif SYSTEM == 'LINUX':
            command = 'ls ' + pVOGdir
            proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
            (rawlisting, err) = proc.communicate()
            listing = rawlisting.decode('utf-8')   # Python3
            fileList = listing.split('\n')
        else:
            cwd = os.getcwd()
            os.chdir(pVOGdir)
            fileList = os.listdir(pVOGdir)
            os.chdir(cwd)

        if DEBUG:
            print("fileList is", fileList)

        # From each file, extract pVOG and associated data fields
        fileList.pop(0)  # remove system message, "Terminal Saved Output".
        for fileName in fileList:
            if DEBUG:
                print("Processing fileName", fileName)
            # initialize
            pVOGid = ''; peptideCount = ''; genomeCount = ''
            family = ''; species = ''; genomeAccn = ''; peptideAccn = ''; location = ''; number = 0; description = ''
            PVOG = False

            # prepend directory
            nextFile = pVOGdir + fileName
            if DEBUG:
                print("nextFile is", nextFile)
            if os.path.isfile(nextFile):
           
                # open next pVOG library file
                file_h = open(nextFile,"r")

                # extract pVOG identifier and associated accessions
                lines = file_h.read().splitlines()
                for line in lines:
                    if DEBUG:
                        print("Processing line", line)

                    # what information is on this line?
                    match_pVOG       = re.search(p_pVOG,line)
                    match_pVOGrecord = re.search(p_pVOGrecord,line)

                    # extract the information; pVOG identifier line should be 1st
                    if match_pVOG:
                        pVOG = copy.deepcopy(self.pVOGobj)
                        pVOG.pVOGid  = match_pVOG.group(1) 
                        peptideCount = match_pVOG.group(2) 
                        genomeCount  = match_pVOG.group(3)
                        pVOG.summary = "Reported", peptideCount, "peptides from", genomeCount, "genomes"
                        PVOG = True

                    elif match_pVOGrecord:
                        fields = line.split('\t')
                        if PVOG:
                            nextPvogRecord = copy.deepcopy(pVOGrecord)
                            #nextPvogRecord['family']      = match_pVOGrecord.group(1)
                            #nextPvogRecord['species']     = match_pVOGrecord.group(2) 
                            #nextPvogRecord['genomeAccn']  = match_pVOGrecord.group(3) 
                            #nextPvogRecord['peptideAccn'] = match_pVOGrecord.group(4) 
                            #nextPvogRecord['location']    = match_pVOGrecord.group(5) 
                            #nextPvogRecord['number']      = match_pVOGrecord.group(6) 
                            #nextPvogRecord['description'] = match_pVOGrecord.group(7) 
                            if len(fields) >= 6:
                                nextPvogRecord['family']      = fields[0] 
                                nextPvogRecord['species']     = fields[1] 
                                nextPvogRecord['genomeAccn']  = fields[2] 
                                nextPvogRecord['peptideAccn'] = fields[3] 
                                nextPvogRecord['location']    = fields[4] 
                                nextPvogRecord['number']      = fields[5] 
                            if len(fields) == 7:   # Kludge:  sometimes the description field is omitted
                                nextPvogRecord['description'] = fields[6] 
                            else:
                                nextPvogRecord['description'] = "unknown"  
                            pVOG.pVOGrecordList.append(nextPvogRecord)
                            pVOG.accessionList.append(nextPvogRecord['peptideAccn'])
                            # Check for potential redundancy, an add current peptide accession to the cumulative non-redundant list accordingly
                            if nextPvogRecord['peptideAccn'] in self.accnList:
                                print("WARNING:  Redudancy encountered: peptideAccn", peptideAccn, "is already in the non-redundant list, filename", fileName)
                                logFile_h.write("%s%s%s%s\n" % ("WARNING:  Redundancy encountered: peptideAccn ",peptideAccn," is already in the non-redundant list, file ",fileName))
                            else:
                                self.accnList.append(peptideAccn)
                        else:
                            print("WARNING: Irregularity in pVOG library file:", fileName, "line", line, "pVOG", pVOG.pVOGid) 
                            logFile_h.write("%s%s%s%s%s%s\n" % ("WARNING:  Irregularity in pVOG library file ", fileName," line ",line," pVOG ",pVOG.pVOGid))

                if PVOG:
                    # Verify that the number of peptide accessions captured equals the number that were declared in file's comment 
                    if len(pVOG.pVOGrecordList) != int(peptideCount):
                        logFile_h.write("%s%s%s%s%s%s\n" % ("WARNING: Discrepancy in record list and peptide count from file ",fileName," size:",len(pVOG.pVOGrecordList)," pepCount:",peptideCount))
                    self.pVOGlist.append(pVOG)

                # clean up
                file_h.close()

                # Update progress
                progressCount += 1
                progress = progressCount % 100
                if progress == 0:
                    print("Working...", progressCount, "pVOGs have been processed")  
                    logFile_h.write("%s%s%s\n" % ("Working... ", progressCount, " pVOGs have been processed"))

    def addPvogs_old(self,pVOGlines,logFile_h):
        progressCount = 0
        p_pVOG = re.compile('^#(VOG\d+)\s+(\d+)\s+(.*)')
        pVOG_NRlist = []  # ensures that pVOG list is non-redundant

        # Create pVOG objects for each pVOG identifier; add to associated data to object
        # Then, add new pVOG object to list
        for pVOGline in pVOGlines:
            match_pVOG = re.search(p_pVOG,pVOGline) 

            # Extract associated data for this pVOG
            if match_pVOG:
                progressCount += 1
                progress = progressCount % 100
                if progress == 0:
                    print("Working...", progressCount, "pVOGs have been processed")
                    logFile_h.write("%s%s%s\n" % ("Working... ",progressCount," pVOGs have been processed."))
                pVOGid       = match_pVOG.group(1)
                count        = match_pVOG.group(2)
                memberString = match_pVOG.group(3)
                if pVOGid not in pVOG_NRlist:  # Have we seen this pVOG already? 
                    pVOG_NRlist.append(pVOGid) # proceed only if this pVOG not seen before
                    pVOG = copy.deepcopy(self.pVOGobj)  # create new pVOG object
                    pVOG.pVOGid  = pVOGid
                    pVOG.count   = count
                    memberList   = memberString.split(',')
 
                    # Extract the members (genome-accession pairs) of this pVOG
                    for member in memberList:
                        if member not in pVOG.memberList:
                            pVOG.memberList.append(member)
                            (genome,accession) = member.split('-')
                            if accession not in pVOG.accessionList:  # no duplicates allowed
                                pVOG.accessionList.append(accession) 
                            else:
                                print("WARNING: for pVOG", pVOGid, "member", member, "accession", accession, "is redundant")
                                logFile_h.write("%s%s%s%s%s%s%s\n" % ("WARNING: for pVOG ",pVOGid," member ",member," accession ",accession," is redundant "))
                        else:
                            print("WARNING: for pVOG", pVOGid, "there is a redundant member:", member)
                            logFile_h.write("%s%s%s%s\n" % ("WARNING: for pVOG ",pVOGid," there is a redundant member: ",member))
                    # Add the new pVOG object to the pVOGlist
                    self.pVOGlist.append(pVOG)

    def getPvogCount(self):
        return len(self.pVOGlist)

    def getAccessionCount(self):
        accessionCount = 0
        for pVOG in self.pVOGlist:
            accessionCount += len(pVOG.accessionList)
        return accessionCount


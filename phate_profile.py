############################################################################
# THIS MODULE IS UNDER DEVELOPMENT 
#
# Module:  phate_profile.py
# Most Recent Update: 20 December 2019
#
# This class performs hmm searches of hmm profiles or alignments against 
# various protein- or phage-related hmm profile databases.  The terminology may
# be confusing. This module uses profiles against profile databases. The "hmm"
# module run hmm codes to search sequence databases.
# 
# Programmer:  C Zhou
#
# Programmer's Notes:
# Currently working on running profile search with hmmscan. Need to install the pVOGs profile
# database on this machine; then test/fix code for running/parsing/saving profile searches.
#
# Classes and Methods:
#
############################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.PDF FOR DETAILS.

import re
import copy
import os
import subprocess
from xml.etree.ElementTree import ElementTree as ET
import phate_fastaSequence
import phate_genomeSequence
import phate_annotation
from subprocess import Popen, PIPE, STDOUT
import string

# DEBUG control
DEBUG = True
#DEBUG = False

MAX_SEQ_HITS = 100 # Upper limit

# Get environment variables (set in phate_runPipeline.py)
# Hmm searches in this module are done against fasta blast databases 

HMMER_HOME                           = os.environ["PHATE_HMMER_HOME"]
NCBI_VIRUS_GENOME_PROFILE_DB_HOME    = os.environ["PHATE_NCBI_VIRUS_GENOME_HMM_HOME"]
NCBI_VIRUS_PROTEIN_PROFILE_DB_HOME   = os.environ["PHATE_NCBI_VIRUS_PROTEIN_HMM_HOME"]
REFSEQ_PROTEIN_PROFILE_DB_HOME       = os.environ["PHATE_REFSEQ_PROTEIN_HMM_HOME"]
REFSEQ_GENE_PROFILE_DB_HOME          = os.environ["PHATE_REFSEQ_GENE_HMM_HOME"]
PVOGS_PROFILE_DB_HOME                = os.environ["PHATE_PVOGS_HMM_HOME"]
PHANTOME_PROFILE_DB_HOME             = os.environ["PHATE_PHANTOME_HMM_HOME"]
PHAGE_ENZYME_PROFILE_DB_HOME         = os.environ["PHATE_PHAGE_ENZYME_HMM_HOME"]
KEGG_VIRUS_PROFILE_DB_HOME           = os.environ["PHATE_KEGG_VIRUS_HMM_HOME"]
PFAM_PROFILE_DB_HOME                 = os.environ["PHATE_PFAM_HMM_HOME"]
SMART_PROFILE_DB_HOME                = os.environ["PHATE_SMART_HMM_HOME"]
SWISSPROT_PROFILE_DB_HOME            = os.environ["PHATE_SWISSPROT_HMM_HOME"]
UNIPROT_PROFILE_DB_HOME              = os.environ["PHATE_UNIPROT_HMM_HOME"]
NR_PROFILE_DB_HOME                   = os.environ["PHATE_NR_HMM_HOME"]
# database of pVOG headers for searching and sequences for writing fasta files (pre-alignment)
PVOG_HEADERS                         = os.environ["PHATE_PVOGS_BASE_DIR"] + "pVOGs.headers.lst"
PVOG_SEQUENCES                       = os.environ["PHATE_PVOGS_BASE_DIR"] + "pVOGs.faa"

# Verbosity and output/error capture
CLEAN_RAW_DATA                       = os.environ["PHATE_CLEAN_RAW_DATA"]
PHATE_WARNINGS                       = os.environ["PHATE_PHATE_WARNINGS"]
PHATE_MESSAGES                       = os.environ["PHATE_PHATE_MESSAGES"]
PHATE_PROGRESS                       = os.environ["PHATE_PHATE_PROGRESS"]
PHATE_ERR                            = os.environ["PHATE_PHATE_ERR"]
PHATE_OUT                            = os.environ["PHATE_PHATE_OUT"]

# Other configurables 

GENE_CALL_DIR            = ""  # set by set method, via parameter list
PROFILE_OUT_DIR          = ""  # set by set method, via parameter list
TBL  = 1
XML  = 2
LIST = 3

# templates 
annotation = phate_annotation.annotationRecord()

class multiProfile(object):

    def __init__(self):
        self.profileProgram             = 'unknown' # Select from 'hmmscan', others? 
        self.profileAnnotations         = []        # List of phate_annotation objects; blast output get temporarily stored here
        self.geneCallDir                = ""        # needs to be set
        self.profileOutDir              = ""        # needs to be set
        self.genomeProfileOutDir        = ""        # needs to be set
        self.geneProfileOutDir          = ""        # needs to be set
        self.proteinProfileOutDir       = ""        # needs to be set
        self.pVOGsOutDir                = ""        # needs to be set
        self.topHitCount                = MAX_SEQ_HITS 
        self.hmmscan                    = False     # If True, run hmmscan
        self.outputFormat               = TBL
        self.NCBI_VIRUS_PROTEIN_PROFILE = False     # Booleans control whether hmm search will be done against a given fasta blast database 
        self.REFSEQ_PROTEIN_PROFILE     = False     #
        self.REFSEQ_GENE_PROFILE        = False     #  
        self.PVOGS_PROFILE              = False     # 
        self.PHANTOME_PROFILE           = False     # 
        self.PHAGE_ENZYME_PROFILE       = False     # 
        self.KEGG_VIRUS_PROFILE         = False     # 
        self.PFAM_PROFILE               = False     # 
        self.SMART_PROFILE              = False     # 
        self.SWISSPROT_PROFILE          = False     #  
        self.UNIPROT_PROFILE            = False     # 
        self.NR_PROFILE                 = False     #  

    ##### SET AND GET PARAMETERS

    def setProfileParameters(self,paramset):
        if isinstance(paramset,dict):
            # Set by function for quality control
            if 'profileProgram' in list(paramset.keys()):
                self.profileProgram = paramset['profileProgram']
            if 'geneCallDir' in list(paramset.keys()):
                self.setGeneCallDir(paramset['geneCallDir'])
            if 'profileOutDir' in list(paramset.keys()):
                self.profileOutDir = paramset['profileOutDir']
            if 'genomeProfileOutDir' in list(paramset.keys()):
                self.genomeProfileOutDir = paramset['genomeProfileOutDir']
            if 'geneProfileOutDir' in list(paramset.keys()):
                self.geneProfileOutDir = paramset['geneProfileOutDir']
            if 'proteinProfileOutDir' in list(paramset.keys()):
                self.proteinProfileOutDir = paramset['proteinProfileOutDir']
            if 'pVOGsOutDir' in list(paramset.keys()):
                self.pVOGsOutDir = paramset['pVOGsOutDir']
            if 'hmmscan' in list(paramset.keys()):  # Only hmmscan, for now
                self.hmmscan = paramset['hmmscan'] 
                if self.hmmscan == True:
                    self.profileProgram = 'hmmscan'
            # Booleans to control execution  #*** Only pvogsHmm is currently in service
            if 'ncbiVirusGenomeHmm' in list(paramset.keys()):
                self.NCBI_VIRUS_GENOME_PROFILE = paramset["ncbiVirusGenomeHmm"]
            if 'ncbiVirusProteinHmm' in list(paramset.keys()):
                self.NCBI_VIRUS_PROTEIN_PROFILE = paramset["ncbiVirusProteinHmm"]
            if 'refseqProteinHmm' in list(paramset.keys()):
                self.REFSEQ_PROTEIN_PROFILE = paramset["refseqProteinHmm"]
            if 'refseqGeneHmm' in list(paramset.keys()):
                self.REFSEQ_GENE_PROFILE = paramset["refseqGeneHmm"]
            if 'pvogsHmm' in list(paramset.keys()):
                self.PVOGS_PROFILE = paramset["pvogsHmm"]
            if 'phantomeHmm' in list(paramset.keys()):
                self.PHANTOME_PROFILE = paramset["phantomeHmm"]
            if 'phageEnzymeHm' in list(paramset.keys()):
                self.PHAGE_ENZYME_PROFILE = paramset["phageEnzymeHmm"]
            if 'keggVirusHmm' in list(paramset.keys()):
                self.KEGG_VIRUS_PROFILE = paramset["keggVirusHmm"]
            if 'pfamHmm' in list(paramset.keys()):
                self.PFAM_PROFILE = paramset["pfamHmm"]
            if 'smartHmm' in list(paramset.keys()):
                self.SMART_PROFILE = paramset["smartHmm"]
            if 'swissprotHmm' in list(paramset.keys()):
                self.SWISSPROT_PROFILE = paramset["swissprotHmm"]
            if 'uniprotHmm' in list(paramset.keys()):
                self.UNIPROT_PROFILE = paramset["uniprotHmm"]
            if 'nrHmm' in list(paramset.keys()):
                self.NR_PROFILE = paramset["nrHmm"]

    def setTopHitCount(self,number):
        if (int(number) >= 1 and int(number) <= HIT_COUNT_MAX):
            self.topHitCount = int(number)
        else:
            if PHATE_WARNINGS == 'True':
                print("WARNING in hmm module: You may capture from 1 to", HIT_COUNT_MAX, "hits per query. If this is insufficient, you may change HIT_COUNT_MAX in phate_blast.py.")

    def setOutputFormat(self,outfmt):
        if outfmt == XML or str(outfmt) == str(XML) or outfmt.lower() == 'xml':
            self.outputFormat = XML
        elif outfmt == LIST or str(LIST) == str(LIST) or outfmt.lower() == 'list':
            self.outputFormat = LIST 
        elif outfmt == TBL or str(outfmt) == str(TBL) or outfmt.lower() == 'tbl':
            self.outputFormat = TBL 
        else:
            if PHATE_WARNINGS == 'True':
                print("WARNING in hmm module: Only acceptable output format is TBL or", TBL) 

    def setGeneCallDir(self,geneCallDir):
        self.geneCallDir = geneCallDir
        GENE_CALL_DIR = geneCallDir

    def setProfileOutDir(self,profileOutDir):
        self.profileOutDir = profileOutDir
        PROFILE_OUT_DIR = profileOutDir

    def getTopHits(self):
        return self.topHitList 

    ##### PERFORM HMM SEARCH ON PROFILE DATABASE - SINGLE SEQUENCE

    def profile1fasta(self,fasta,outfile,database,dbName): # fasta is a phate_fastaSequence.fasta object
        command = ''

        # Write fasta sequence to temporary file
        fastaFile = self.proteinProfileOutDir + "temp.fasta"
        fastaFileH = open(fastaFile,"w")
        if fasta.sequentialHeader == "unknown":  # unchanged from default
            fasta.printFasta2file(fastaFileH,"blastHeader")
        else:
            fasta.printFasta2file(fastaFileH,"sequential")  # use sequential header format to avoid special chars issue
        fastaFileH.close()

        # Construct out file names
        seqOutfile = outfile + '.hmmscan.seqout' # captures tabbed data for the sequence-level analysis
        domOutfile = outfile + '.hmmscan.domout' # captures tabbed data for the domain-level analysis
        stdOutfile = outfile + '.hmmscan.stdout' # captures data written to standard out, other than seq or dom output
 
        # Run hmm program; write output to specified file  # Note: only hmmscan is currently in service
        if self.hmmscan == True:
            command = HMMER_HOME + "hmmscan --tblout " + seqOutfile + ' --domtblout ' + domOutfile + ' ' + database + ' ' + fastaFile 

        if PHATE_OUT == 'True':
            p = Popen(command, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
            output = p.stdout.read() 
            OUT_H = open(stdOutfile, "w")
            OUT_H.write("%s\n" % (output)) 
            OUT_H.close()
        else:
            result = os.system(command)

        # For TBL parsing
        hitList = [] # contains a list of hitDataSet objects
        sequenceDataSet = {
            "hitNumber"         : 0,
            "hitSequenceName"   : "",  # ex: ??? |ref|YP_002308524.1|
            "hitDescription"    : "",  # ex: phage 21-like group II holin [Ba  # truncated
            "hitProfile"        : "",  # ex: VOG4733|ref|Yp_008126541.1|
            "hitRef"            : "",  # ex: |ref|YP_002308524.1| or |gb|ALY09357.1| 
            "hitLength"         : 0,
            "hitEvalue"         : 0.0,
            "hitScore"          : 0.0,
            "hitBias"           : 0.0,
            "domEvalue"         : 0.0, # Domain data for best 1 domain, listed on main sequence hit line
            "domScore"          : 0.0,
            "domBias"           : 0.0,
            "domExp"            : 0.0,
            "domN"              : 0,
            "hitDomainList"     : [], # list of domainDataSet objects
            }
        domainDataSet = {
            "number"            : 0,
            "score"             : 0.0,
            "c-Evalue"          : 0.0, # conditional E-value
            "i-Evalue"          : 0.0,
            "hmmFrom"           : 0,
            "hmmTo"             : 0,
            "alignFrom"         : 0,
            "alignTo"           : 0,
            "envFrom"           : 0,
            "envTo"             : 0,
            "acc"               : 0.0, 
            }

        # Capture result(s) and store as an annotation object for this fasta; Coded for output format TBL (5) 

        # Parse XML-, LIST-, or TBL-formatted hmm output #*** Is XML format available?  

        if DEBUG:
            print("self.outputFormat is", self.outputFormat)

        if self.outputFormat == XML or self.outputFormat == LIST:
            if PHATE_WARNINGS == 'True':
                print("ERROR in profile module:  only TBL format is currently being used for hmmscan output parsing")

        elif self.outputFormat == TBL:

            # Parse the sequence-level profile search output data (global hit)
            #***
            # Note:  Although data appear to be in explicitly defined fields (columns), they are not!
            #        Column locations for data values change in each file. Therefore, you cannot parse
            #        the jackhmmer output using fixed columns. Nor are the columns tab-delimited. The
            #        following parser separates data fields by splitting on >= 1 white space. However,
            #        this will fail if the subject (target) fasta header contains white space!!! It will
            #        work for the pVOGs fasta file that is prepared for use with PhATE, but will not work
            #        for most fasta input data sets. I need to explore other possible jackhmmer output
            #        formats, if there are others. XML would be great.... 
            hitCount = 0

            # Open and process the sequence hit file; note: domain hits are in a separate file
            seqOutfileH = open(seqOutfile,"r")  # Open the sequence hit file
            sLines = seqOutfileH.read().splitlines()
            for sLine in sLines:
                match_comment = re.search('^#',sLine)
                if match_comment:
                    continue
                
                else: # extract data fields and remove preceding or trailing whitespace
                    hitCount += 1
                    fields = sLine.split()  # split line on white space
                    targetName        = fields[0];   targetName = targetName.rstrip() 
                    targetAccession   = fields[1];   targetAccession = targetAccession.rstrip()
                    queryName         = fields[2];   queryName = queryName.rstrip()
                    queryAccession    = fields[3];   queryAccession = queryAccession.rstrip()
                    seqEvalue         = fields[4];   seqEvalue = seqEvalue.lstrip()
                    seqScore          = fields[5];   seqScore = seqScore.lstrip()
                    seqBias           = fields[6];   seqBias = seqBias.lstrip()
                    dom1evalue        = fields[7];   dom1evalue = dom1evalue.lstrip()
                    dom1score         = fields[8];   dom1score = dom1score.lstrip()
                    dom1bias          = fields[9];   dom1bias = dom1bias.lstrip()
                    dom1exp           = fields[10];  dom1exp = dom1exp.lstrip()
                    dom1reg           = fields[11];  dom1reg = dom1reg.lstrip()
                    dom1clu           = fields[12];  dom1clu = dom1clu.lstrip()
                    dom1ov            = fields[13];  dom1ov = dom1ov.lstrip()
                    dom1env           = fields[14];  dom1env = dom1env.lstrip()
                    dom1dom           = fields[15];  dom1dom = dom1dom.lstrip()
                    dom1rep           = fields[16];  dom1rep = dom1rep.lstrip()
                    dom1inc           = fields[17];  dom1inc = dom1inc.lstrip()
                    #targetDescription = fields[18];  targetDescription = targetDescription.rstrip()
                    targetDescription = ' '.join(fields[18:]);  targetDescription = targetDescription.rstrip()

                    if DEBUG:
                        print("DEBUG in phate_profile.py / sequence parsing")
                        print("targetAccession:", targetAccession)
                        print("hitCount:", hitCount, "targetName:", targetName, "queryName:", queryName)
                        print("queryAccession:", queryAccession, "seqEvalue:", seqEvalue, "seqScore:", seqScore, "seqBias:", seqBias)
                        print("dom1evalue:", dom1evalue, "dom1score:", dom1score, "dom1bias:", dom1bias)
                        print("dom1exp:", dom1exp, "dom1reg:", dom1reg, "dom1clu:", dom1clu, "dom1ov:", dom1ov, "dom1env:", dom1env)
                        print("dom1dom:", dom1dom, "dom1rep:", dom1rep, "dom1inc:", dom1inc)
                        print("targetDescription:", targetDescription)
                    
                    vogIDs = ''; pVOGlist = []
                    if re.search("pvog",database.lower()):
                        # Collect pVOG identifiers for this hmm/profile search hit; fasta header of target may have >= 1 pVOG identifier
                        if DEBUG:
                            print("DEBUG: targetName is", targetName)
                        pVOGlist = re.findall('VOG\d+', targetName)
                        if DEBUG:
                            print("DEBUG: VOG IDs found, pVOGlist is", pVOGlist)
                        if pVOGlist:
                            for pVOG in pVOGlist:
                                vogIDs += pVOG + ' '
                            vogIDs.rstrip()
                            if DEBUG:
                                print("DEBUG: vogIDs is", vogIDs)
                        else:
                            if DEBUG:
                                print("DEBUG: pVOGlist is empty")
                     
                    # Create new hitDataSet object and store data (note: some data may not be stored)
                    newSequenceDataSet = copy.deepcopy(sequenceDataSet)
                    newSequenceDataSet["hitNumber"]       = hitCount
                    newSequenceDataSet["hitSequenceName"] = targetName
                    newSequenceDataSet["hitDescription"]  = targetDescription
                    newSequenceDataSet["hitVOG"]          = vogIDs
                    newSequenceDataSet["hitEvalue"]       = seqEvalue
                    newSequenceDataSet["hitScore"]        = seqScore
                    newSequenceDataSet["hitBias"]         = seqBias
                    newSequenceDataSet["domEvalue"]       = dom1evalue
                    newSequenceDataSet["domScore"]        = dom1score
                    newSequenceDataSet["domBias"]         = dom1bias
                    newSequenceDataSet["domExp"]          = dom1exp
                    hitList.append(newSequenceDataSet)
                    if hitCount >= MAX_SEQ_HITS:
                        break # enough hits; we are done!

            # Collect all pVOG identifiers for this hmm search, then remove redundancy from list
            vogs = []; tempVOGlist = []; nrVOGlist = []
            for hit in hitList:
                if hit["hitVOG"]:
                    vogs = hit["hitVOG"].split(' ')
                    for vog in vogs:
                        tempVOGlist.append(vog)  # Now we have a complete (super)set of all pVOG identifiers found in hmm search 
            nrVOGlist = list(set(tempVOGlist)) # create a set, then convert back to list: presto! non-redundant list
            if DEBUG:
                print("DEBUG: nrVOGlist is", nrVOGlist)

            # Parse the domain-level hmm data
            # There can be multiple domain-level hits for a given sequence (global) hit (above)
            # Process each domain hit, then identify its corresponding sequence-level hit, and add to that object
            # Data are in explicitly defined fields (columns)
            domOutfileH = open(domOutfile,"r")
            dLines = domOutfileH.read().splitlines()
            for dLine in dLines:
                match_comment = re.search('^#',dLine)
                if match_comment:
                    continue

                else:  # extract data fields and remove preceding or trailing whitespace
                    fields = dLine.split()
                    targetName        = fields[0];  targetName = targetName.rstrip()
                    targetAccession   = fields[1];  targetAccession = targetAccession.rstrip()
                    targetLength      = fields[2];  targetLength = targetLength.lstrip()
                    queryName         = fields[3];  queryName = queryName.rstrip()
                    queryAccession    = fields[4];  queryAccession = queryAccession.rstrip()
                    queryLength       = fields[5];  queryLength = queryLength.lstrip()
                    fullSeqEvalue     = fields[6];  fullSeqEvalue = fullSeqEvalue.lstrip()
                    fullSeqScore      = fields[7];  fullSeqScore = fullSeqScore.lstrip()
                    fullSeqBias       = fields[8];  fullSeqBias = fullSeqBias.lstrip()
                    domainNumber      = fields[9];  domainNumber = domainNumber.lstrip()
                    domOf             = fields[10]; domOf = domOf.lstrip()
                    cEvalue           = fields[11]; cEvalue = cEvalue.lstrip()
                    iEvalue           = fields[12]; iEvalue = iEvalue.lstrip()
                    score             = fields[13]; score = score.lstrip()
                    bias              = fields[14]; bias = bias.lstrip()
                    hmmFrom           = fields[15]; hmmFrom = hmmFrom.lstrip()
                    hmmTo             = fields[16]; hmmTo = hmmTo.lstrip()
                    alignFrom         = fields[17]; alignFrom = alignFrom.lstrip()
                    alignTo           = fields[18]; alignTo = alignTo.lstrip()
                    envFrom           = fields[19]; envFrom = envFrom.lstrip()
                    envTo             = fields[20]; envTo = envTo.lstrip()
                    acc               = fields[21]; acc = acc.lstrip()
                    targetDescription = fields[22]; targetDescription = targetDescription.rstrip()
                    
                    if DEBUG:
                        print("DEBUG in phate_profile.py / domain parsing")
                        print("targetName:", targetName, "targetAccession:", targetAccession, "targetLength:", targetLength, "queryName:", queryName)
                        print("queryAccession:", queryAccession, "queryLength:", queryLength, "fullSeqEvalue:", fullSeqEvalue)
                        print("fullSeqScore:", fullSeqScore, "fullSeqBias:", fullSeqBias, "domainNumber:", domainNumber)
                        print("cEvalue:", cEvalue, "domOf:", domOf, "iEvalue:", iEvalue, "score:", score, "bias:", bias)
                        print("hmmFrom:", hmmFrom, "hmmTo:", hmmTo, "alignFrom:", alignFrom, "envFrom:", envFrom, "envTo:", envTo, "acc:", acc)
                        print("targetDescription:", targetDescription) 

                    # Create new domDataSet object and store data 
                    newDomainDataSet = copy.deepcopy(domainDataSet)
                    newDomainDataSet["number"]    = domainNumber
                    newDomainDataSet["score"]     = score 
                    newDomainDataSet["c-Evalue"]  = cEvalue
                    newDomainDataSet["i-Evalue"]  = iEvalue
                    newDomainDataSet["hmmFrom"]   = hmmFrom
                    newDomainDataSet["hmmTo"]     = hmmTo
                    newDomainDataSet["alignFrom"] = alignFrom
                    newDomainDataSet["alignTo"]   = alignTo
                    newDomainDataSet["envFrom"]   = envFrom
                    newDomainDataSet["envTo"]     = envTo
                    newDomainDataSet["acc"]       = acc  

                    # Figure out which seqDataSet to add this domDataSet to
                    FOUND = False
                    for hit in hitList:
                        if hit["hitSequenceName"] == targetName:
                            hit["hitDomainList"].append(newDomainDataSet)
                            FOUND = True
                    if not FOUND:
                        if PHATE_WARNINGS == 'True':
                            print("WARNING in profile module: sequence data object not found for domain object", targetName)

            # Close HMM output files
            seqOutfileH.close()
            domOutfileH.close()

            # Assuming there were hits...create an annotation object, and fill with sequence- and domain- hit data 
            if hitList:
                for hit in hitList:
                    # Extract hmm info from hitLine and stash into new annotation object
                    newAnnotation = copy.deepcopy(annotation) 
                    newAnnotation.source         = database 
                    newAnnotation.method         = self.profileProgram
                    newAnnotation.annotationType = "profile search"
                    newAnnotation.name           = hit["hitSequenceName"] # subject
                    newAnnotation.description    = hit["hitDescription"]
                    newAnnotation.start          = '1' # query start
                    newAnnotation.end            = 'N' # query end
                    newAnnotation.category       = "sequence"

                    # Extract pVOG identifier(s) from hit's header (if pVOG database)
                    pVOGlist = []; pVOGstring = ""
                    pVOGlist = hit["hitVOG"].split(' ')
                    for pVOG in pVOGlist: 
                        newAnnotation.pVOGlist.append(pVOG)
                        pVOGstring += pVOG + '_'

                    # Collapse all domain hits as into an annotation list for this sequence/global hit
                    for domain in hit["hitDomainList"]:
                        resultString  = pVOGstring                           + '|'
                        resultString += "domainNo:"    + domain["number"]    + '|'
                        resultString += "score:"       + domain["score"]     + '|'
                        resultString += "c-E:"         + domain["c-Evalue"]  + '|'
                        resultString += "i-E:"         + domain["i-Evalue"]  + '|'
                        resultString += "hmmFrom:"     + domain["hmmFrom"]   + '|'
                        resultString += "hmmTo:"       + domain["hmmTo"]     + '|'
                        resultString += "alignFrom:"   + domain["alignFrom"] + '|'
                        resultString += "alignTo:"     + domain["alignTo"]   + '|'
                        resultString += "envFrom:"     + domain["envFrom"]   + '|'
                        resultString += "envTo:"       + domain["envTo"]     + '|'
                        resultString += "acc:"         + domain["acc"]       + '|'
                        newAnnotation.annotationList.append(resultString)

                    # Add this completed annotation to growing list for this fasta
                    fasta.annotationList.append(newAnnotation)
            else:
                if PHATE_MESSAGES == 'True':
                    print("Profile module says, No Profile hit found for query", fasta.blastHeader, "against", database)    

        # Requested Hmm/Profile output format not supported
        else:
            if PHATE_WARNINGS == 'True':
                print("WARNING in profile module: Output format", self.outputFormat, "not yet supported in phate_profile.py/profile1fasta(). Use TBL format for now.")

    def runProfile(self,fastaSet,dbType="protein"): # fastaSet is a phate_fastaSequence.multiFasta object

        if DEBUG:
            print("Profile module, runProfile says: Running runProfile...")
        # Set sequence type 
        GENOME = False; GENE = False; PROTEIN = False
        if dbType.lower() == "protein" or dbType.lower == "aa" or dbType.lower == "prot" or dbType.lower == "peptide":
            PROTEIN = True
        elif dbType.lower() == "genome":
            GENOME = True
        elif dbType.lower() == "gene":
            GENE = True
        else:
            if PHATE_WARNINGS == 'True':
                print("WARNING in profile module: unrecognized database type in runProfile:", dbType)
            return
               
        # Set database variable, invoke HMM program for each fasta 
        database = ''

        if DEBUG:
            print("TESTING: GENOME, GENE, PROTEIN, are: ", GENOME, GENE, PROTEIN)
        if GENOME:
            if PHATE_WARNINGS == 'True':
                print("WARNING in profile module:  Currently genome HMM/profile search is not supported.")

        if GENE:
            if self.REFSEQ_GENE_PROFILE:
                database = REFSEQ_GENE_PROFILE_DB_HOME # same database as blast uses
                dbName   = 'refseqGene'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("Profile module says: Running Refseq gene hmm/profile search:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.geneProfileOutDir + self.profileProgram + "_refseqGene_" + str(count)
                    self.profile1fasta(fasta,outfile,database,dbName)  #*** CONTROL

        if PROTEIN:
            if self.NR_PROFILE:  
                database = NR_PROFILE_DB_HOME 
                dbName   = 'nr'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("Profile module says: Running NR hmm/profile search:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.proteinProfileOutDir + self.profileProgram + "_nrProfile_" + str(count)
                    self.profile1fasta(fasta,outfile,database,dbName)  #*** CONTROL

            if self.NCBI_VIRUS_PROTEIN_PROFILE:  
                database = NCBI_VIRUS_PROTEIN_PROFILE_DB_HOME 
                dbName   = 'ncbiVirusProteinHmm'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("Profile module says: Running NCBI_VIRUS_PROTEIN hmm/profile search:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.proteinProfileOutDir + self.profileProgram + "_ncbiVirProtProfile_" + str(count)
                    self.profile1fasta(fasta,outfile,database,dbName)  #*** CONTROL

            if self.REFSEQ_PROTEIN_PROFILE:
                database = REFSEQ_PROTEIN_PROFILE_DB_HOME 
                dbName   = 'refseqProteinHmm'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("Profile module says: Running Refseq protein hmm/profile search:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.proteinProfileOutDir + self.profileProgram + "_refseqProteinProfile_" + str(count)
                    self.profile1fasta(fasta,outfile,database,dbName)  #*** CONTROL

            if self.PHANTOME_PROFILE:
                database = PHANTOME_PROFILE_DB_HOME # same database as blast uses
                dbName   = 'phantomeProfile'
                count = 0 
                if PHATE_PROGRESS == 'True':
                    print("Profile module says: Running PHANTOME hmm/profile search:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.proteinProfileOutDir + self.profileProgram + "_phantomeProfile_" +str(count)
                    self.profile1fasta(fasta,outfile,database,dbName)

            if self.KEGG_VIRUS_PROFILE: 
                database = KEGG_VIRUS_PROFILE_DB_HOME # same database as blast uses
                dbName   = 'keggProfile'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("Profile module says: Running KEGG hmm/profile search:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.proteinProfileOutDir + self.profileProgram + "_keggProfile_" + str(count)
                    self.profile1fasta(fasta,outfile,database,dbName)

            if self.SWISSPROT_PROFILE: 
                database = SWISSPROT_PROFILE_DB_HOME # same database as blast uses
                dbName   = 'swissprotProfile'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("Profile module says: Running Swissprot hmm/profile search:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.proteinProfileOutDir + self.profileProgram + "_swissprotProfile_" + str(count)
                    self.hmm1fasta(fasta,outfile,database,dbName)   #*** CONTROL

            if self.PHAGE_ENZYME_PROFILE: 
                database = PHAGE_ENZYME_PROFILE_DB_HOME 
                dbName   = 'phageEnzProfile'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("Profile module says: Running phage enzyme hmm/profile search:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.proteinProfileOutDir + self.profileProgram + "_phageEnzProfile_" + str(count)
                    self.profile1fasta(fasta,outfile,database,dbName)   #*** CONTROL

            if self.PVOGS_PROFILE:  
                database = PVOGS_PROFILE_DB_HOME 
                pVOGseqDB = PVOG_SEQUENCES
                dbName   = 'pVOGsHmm'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("Profile module says: Running pVOGs hmm/profile search:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.proteinProfileOutDir + self.profileProgram + "_pvogHmm_" + str(count)
                    self.profile1fasta(fasta,outfile,database,dbName)

                if PHATE_PROGRESS == 'True':
                    print("Profile module says: Done!")
                    print("Profile module says: Collecting and saving pVOG sequences corresponding to hmm/profile hit(s)")

                # Next you want to create pVOG fasta group files so user can do alignments
                # You need only one "alignment" file per pVOG group that the fasta hit (under blast cutoffs)
                # Capture pVOG.faa lines
                #pVOGs_h = open(database,"r")
                pVOGs_h = open(pVOGseqDB,"r")
                pVOGlines = pVOGs_h.read().splitlines()
                pVOGs_h.close()
                if DEBUG:
                    print("There are", len(pVOGlines), "pVOG database lines to search")
                count = 0; countA = 0
                for fasta in fastaSet.fastaList:
                    pvogPrintedList = []  # keeps track of pVOGs that have already been printed for current fasta
                    count += 1 
                    countA = 0
                    for annot in fasta.annotationList:
                        if DEBUG:
                            print("DEBUG: pVOGlist is", annot.pVOGlist)
                        for pVOG in annot.pVOGlist:  # There may be multiple annotations to inspect
                            if DEBUG:
                                print("DEBUG: Processing pVOG", pVOG)
                            match_good = re.search('VOG',pVOG)
                            if match_good:
                                # Avoid redundancy in printing pVOG groups for this fasta; only once per pVOG ID that was an hmm hit
                                if pVOG not in pvogPrintedList:
                                    pvogPrintedList.append(pVOG)  # Record this pVOG identifier as "done"
                                    # create dynamic file name
                                    countA += 1 
                                    outfilePVOG = self.pVOGsOutDir + "profile_pvogGroup_" + str(count) + '_' + str(countA) + '.faa' 
                                    if DEBUG:
                                        print("Writing to outfile", outfilePVOG, "pVOG group", pVOG, "pvogPrintedList:", pvogPrintedList)
                                    # open file and write current fasta pluse each corresponding pVOG fasta
                                    outfilePVOG_h = open(outfilePVOG,'w')
                                    outfilePVOG_h.write("%c%s\n%s\n" % ('>',fasta.header,fasta.sequence)) # write the current peptide fasta,
                                    self.writePVOGsequences2file(outfilePVOG_h,pVOGlines,pVOG)                  # followed by the pVOG group
                                    outfilePVOG_h.close()
                            else:
                                if DEBUG:
                                    print("WARNING: unexpected pVOG identifier:", pVOG)        

        if CLEAN_RAW_DATA == 'True':
            if PHATE_PROGRESS == 'True':
                print("Profile module says: Removing raw data files.")
            self.cleanProfileOutDir()

    def writePVOGsequences2file(self,FILE_H,lines,pVOG):
        pVOGheader = ""; pVOGsequence = ""; GET_SEQ = False
        for i in range(0,len(lines)-1):
            nextLine = lines[i]
            match_header = re.search('>',nextLine)
            match_pVOG = re.search(pVOG,nextLine)

            if match_header:   
                if GET_SEQ:
                    if DEBUG:
                        print("Writing sequence to file for header", pVOGheader)
                    FILE_H.write("%s\n%s\n" % (pVOGheader,pVOGsequence))
                    pVOGheader = ""; pVOGsequence = ""
                    GET_SEQ = False

                if match_pVOG:
                    pVOGheader = nextLine
                    GET_SEQ = True

            elif GET_SEQ:
                pVOGsequence = pVOGsequence + nextLine

    def cleanProfileOutDir(self):  # Remove temporary files from HMM_OUT_DIR
        #command = "ls " + HMM_OUT_DIR  #*** FIX: list only files, not directories too
        #command = "ls " + self.hmmOutDir  #*** FIX: list only files, not directories too
        command = "ls -p " + self.profileOutDir + " grep -v /"  #*** should list only files, not directories too
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        (rawresult, err) = proc.communicate()
        result = rawresult.decode('utf-8')    # Python3
        if DEBUG:
            print("Result of listing profile out dir,", self.profileOutDir) 
        fileList = result.split('\n')
        for filename in fileList:
            file2delete = self.profileOutDir + filename
            #if re.search('profile',file2delete):
            match_seqout = re.search('seqout',file2delete)
            match_domout = re.search('domout',file2delete)
            match_stdout = re.search('stdout',file2delete)
            if match_seqout or match_domout or match_stdout:
                command = "rm " + file2delete
                proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
                (result, err) = proc.communicate()

    ##### PRINT METHODS

    def printParameters(self):
        print("Parameters:")
        print("   profileProgram:             ", self.profileProgram)
        print("   profileAnnotations:         ", self.profileAnnotations)
        print("   geneCallDir:                ", self.geneCallDir)
        print("   genomeProfileOutDir:        ", self.genomeProfileOutDir)
        print("   geneProfileOutDir:          ", self.geneProfileOutDir)
        print("   proteinProfileOutDir:       ", self.proteinProfileOutDir)
        print("   hmmscan:                    ", self.hmmscan)
        print("   pVOGsOutDir:                ", self.pVOGsOutDir)
        print("   NCBI_VIRUS_GENOME_PROFILE:  ", self.NCBI_VIRUS_GENOME_PROFILE)
        print("   NCBI_VIRUS_PROTEIN_PROFILE: ", self.NCBI_VIRUS_PROTEIN_PROFILE)
        print("   REFSEQ_PROTEIN_PROFILE:     ", self.REFSEQ_PROTEIN_PROFILE)
        print("   REFSEQ_GENE_PROFILE:        ", self.REFSEQ_GENE_PROFILE)
        print("   PVOGS_PROFILE:              ", self.PVOGS_PROFILE)
        print("   PHANTOME_PROFILE:           ", self.PHANTOME_PROFILE)
        print("   PHAGE_ENZYME_PROFILE:       ", self.PHAGE_ENZYME_PROFILE)
        print("   KEGG_VIRUS_PROFILE:         ", self.KEGG_VIRUS_PROFILE)
        print("   PFAM_PROFILE:               ", self.PFAM_PROFILE)
        print("   SMART_PROFILE:              ", self.SMART_PROFILE)
        print("   SWISSPROT_PROFILE:          ", self.SWISSPROT_PROFILE)
        print("   UNIPROT_PROFILE:            ", self.UNIPROT_PROFILE)
        print("   NR_PROFILE:                 ", self.NR_PROFILE)

    def printParameters2file(self,fileHandle):
        fileHandle.write("%s\n" % ("Parameters:"))
        fileHandle.write("%s\n" % ("  profile program: ",self.profileProgram))
        fileHandle.write("%s\n" % ("  profileAnnotations: ",self.profileAnnotations))
        fileHandle.write("%s\n" % ("  geneCallDir: ",self.geneCallDir))
        fileHandle.write("%s\n" % ("  profileOutDir: ",self.profileOutDir))
        fileHandle.write("%s\n" % ("  pVOGsOutDir: ",self.pVOGsOutDir))
        fileHandle.write("%s\n" % ("  NCBI_VIRUS_GENOME_PROFILE: ",self.NCBI_VIRUS_GENOME_PROFILE))
        fileHandle.write("%s\n" % ("  NCBI_VIRUS_PROTEIN_PROFILE: ",self.NCBI_VIRUS_PROTEIN_PROFILE))
        fileHandle.write("%s\n" % ("  REFSEQ_PROTEIN_PROFILE: ",self.REFSEQ_PROTEIN_PROFILE))
        fileHandle.write("%s\n" % ("  REFSEQ_GENE_PROFILE: ",self.REFSEQ_GENE_PROFILE))
        fileHandle.write("%s\n" % ("  PVOGS_PROFILE: ",self.PVOGS_PROFILE))
        fileHandle.write("%s\n" % ("  PHANTOME_PROFILE: ",self.PHANTOME_PROFILE))
        fileHandle.write("%s\n" % ("  PHAGE_ENZYME_PROFILE: ",self.PHAGE_ENZYME_PROFILE))
        fileHandle.write("%s\n" % ("  KEGG_VIRUS_PROFILE: ",self.KEGG_VIRUS_PROFILE))
        fileHandle.write("%s\n" % ("  PFAM_PROFILE: ",self.PFAM_PROFILE))
        fileHandle.write("%s\n" % ("  SMART_PROFILE: ",self.SMART_PROFILE))
        fileHandle.write("%s\n" % ("  SWISSPROT_PROFILE: ",self.SWISSPROT_PROFILE))
        fileHandle.write("%s\n" % ("  UNIPROT_PROFILE: ",self.UNIPROT_PROFILE))
        fileHandle.write("%s\n" % ("  NR_PROFILE: ",self.NR_PROFILE))

    def printAnnotations(self):  #*** Don't need this???; don't use, because code writes directoy to fasta's annotation object
        print("Profile Annotations:")
        if self.profileAnnotations:
            for annotation in self.profileAnnotations:
                print("   ", annotation.PrintAll())
        else:
            print("   There are no hmm/profile annotations")

    def printAnnotations2file(self,fileHandle):
        fileHandle.write("%s\n" % ("Profile Annotations:"))
        if self.profileAnnotations:
            for annotation in self.profileAnnotations:
                annotation.PrintAll2file(fileHandle)
        else:
            fileHandle.write("%s\n" % ("   There are no hmm/profile annotations"))

    def printAll(self):
        self.printParameters()
        self.printAnnotations()

    def printAll2file(self,fileHandle):
        self.printParameters2file(self.fileHandle)
        self.printAnnotations2file(self.fileHandle)

    ##### STATISTICS

    def calculateStatistics(self):
        pass


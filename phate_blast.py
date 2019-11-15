############################################################################
#
# Module:  phate_blast.py
#
# This class performs blast against various phage-related databases. 
#
# Programmer:  Carol Zhou
# 
# Classes and Methods:
#    multiBlast
#       setBlastParameters(dict)
#       setBlastFlavor(flavor)
#       setIdentityMin(identity)
#       setIdentitySelect(identity)
#       setEvalueMin(evalue)
#       setEvalueSelect(evalue)
#       setTopHitCount(number)
#       setTopHits
#       getTopHits
#       blast1fasta(fasta,outfile,database)
#       runBlast(fastaSet,database)
#       cleanBlastOutDir
#       printParameters
#       printParameters2file(file_h)
#       printAnnotations
#       printAll
#       printAll2file(fileH)
#       calculateStatistics  
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

# Get environment variables (set in phate_runPipeline.py)

BLAST_HOME                    = os.environ["BLAST_HOME"]
KEGG_VIRUS_BLAST_HOME         = os.environ["KEGG_VIRUS_BLAST_HOME"]
NCBI_VIRUS_BLAST_HOME         = os.environ["NCBI_VIRUS_BLAST_HOME"]
NCBI_VIRUS_PROTEIN_BLAST_HOME = os.environ["NCBI_VIRUS_PROTEIN_BLAST_HOME"]
PHANTOME_BLAST_HOME           = os.environ["PHANTOME_BLAST_HOME"]
UNIPARC_VIRUS_BLAST_HOME      = os.environ["UNIPARC_VIRUS_BLAST_HOME"]
NR_BLAST_HOME                 = os.environ["NR_BLAST_HOME"]
PVOGS_BLAST_HOME              = os.environ["PVOGS_BLAST_HOME"]
REFSEQ_PROTEIN_BLAST_HOME     = os.environ["REFSEQ_PROTEIN_BLAST_HOME"]
REFSEQ_GENE_BLAST_HOME        = os.environ["REFSEQ_GENE_BLAST_HOME"]
SWISSPROT_BLAST_HOME          = os.environ["SWISSPROT_BLAST_HOME"]

# Verbosity
CLEAN_RAW_DATA                = os.environ["CLEAN_RAW_DATA"]
PHATE_WARNINGS                = os.environ["PHATE_WARNINGS"]
PHATE_MESSAGES                = os.environ["PHATE_MESSAGES"]
PHATE_PROGRESS                = os.environ["PHATE_PROGRESS"]

# Other configurables 

GENE_CALL_DIR            = ""  # set by set method, via parameter list
BLAST_OUT_DIR            = ""  # set by set method, via parameter list
PVOGS_OUT_DIR            = ""  # set by set method, via parameter list
PVOGS_FASTA_DB_NAME      = "pVOGs.fasta"

SCORE_EDGE_MAX = 1.0
OVERHANG_MAX   = 100
HIT_COUNT_MAX  = 10  # Put a limit on the number of hits that should be processed

#DEBUG  = True 
DEBUG  = False 

# blast formats
XML = 5
LIST = 7

# templates 
annotation = phate_annotation.annotationRecord()

class multiBlast(object):

    def __init__(self):
        self.blastFlavor              = 'blastp'  # Select from 'blastp', 'blastn', 'blastx'; default = blastp
        self.identityMin              = 50        # default: Invokes blast with this value (for phage maybe set much lower)
        self.evalueMin                = 0.01      # default: Invokes blast with this value (for phage maybe set much lower)
        self.identitySelect           = 70        # default: Selects best hit if meets or exceeds this value (for phage maybe set much lower)
        self.evalueSelect             = 0.01      # default: Selects best hit if meets or exceeds this value (for phage maybe set to 10)
        self.topHitCount              = 3         # default: Number of best hits to record
        self.scoreEdge                = 0.1       # default: BLAST recommended default
        self.overhang                 = 0.1       # default: BLAST recommended default
        #self.outputFormat            = XML       # default: XML output
        self.outputFormat             = 5         # default: XML output
        self.blastAnnotations         = []        # List of phate_annotation objects; blast output get temporarily stored here
        # move to hit class:  self.topHitList     = []        # 
        self.geneCallDir              = ""        # needs to be set
        self.blastOutDir              = ""        # needs to be set
        self.pVOGsOutDir              = ""        # needs to be set
        self.NCBI_VIRUS_BLAST         = False     # assume not running this blast process, unless changed by parameter set method
        self.NCBI_VIRUS_PROTEIN_BLAST = False     # ditto 
        self.NR_BLAST                 = False     # ditto 
        self.KEGG_VIRUS_BLAST         = False     # ditto
        self.REFSEQ_PROTEIN_BLAST     = False     # ditto
        self.REFSEQ_GENE_BLAST        = False     # ditto # not yet in service
        self.PHANTOME_BLAST           = False     # ditto
        self.PVOGS_BLAST              = False     # ditto
        self.UNIPARC_BLAST            = False     # ditto # not yet in service
        self.SWISSPROT_BLAST          = False     # ditto # not yet in service

    ##### SET AND GET PARAMETERS

    def setBlastParameters(self,paramset):
        if isinstance(paramset,dict):
            if 'blastFlavor' in list(paramset.keys()):
                self.setBlastFlavor(paramset['blastFlavor'])
            if 'identityMin' in list(paramset.keys()):
                self.setIdentityMin(paramset['identityMin'])
            if 'identitySelect' in list(paramset.keys()):
                self.setIdentitySelect(paramset['identitySelect'])
            if 'evalueMin' in list(paramset.keys()):
                self.setEvalueMin(paramset['evalueMin'])
            if 'evalueSelect' in list(paramset.keys()):
                self.setEvalueSelect(paramset['evalueSelect'])
            if 'topHitCount' in list(paramset.keys()):
                self.setTopHitCount(paramset['topHitCount'])
            if 'outputFormat' in list(paramset.keys()):
                self.setOutputFormat(paramset['outputFormat'])
            if 'scoreEdge' in list(paramset.keys()):
                self.setScoreEdge(paramset['scoreEdge'])
            if 'overhang' in list(paramset.keys()):
                self.setOverhang(paramset['overhang'])
            if 'geneCallDir' in list(paramset.keys()):
                self.setGeneCallDir(paramset['geneCallDir'])
            if 'blastOutDir' in list(paramset.keys()):
                self.setBlastOutDir(paramset['blastOutDir'])
            if 'pvogsOutDir' in list(paramset.keys()):
                self.setPVOGsOutDir(paramset['pvogsOutDir'])
            if 'ncbiVirusBlast' in list(paramset.keys()):
                self.NCBI_VIRUS_BLAST = paramset["ncbiVirusBlast"]
            if 'ncbiVirusProteinBlast' in list(paramset.keys()):
                self.NCBI_VIRUS_PROTEIN_BLAST = paramset["ncbiVirusProteinBlast"]
            if 'nrBlast' in list(paramset.keys()):
                self.NR_BLAST = paramset["nrBlast"]
            if 'keggVirusBlast' in list(paramset.keys()):
                self.KEGG_VIRUS_BLAST = paramset["keggVirusBlast"]
            if 'refseqProteinBlast' in list(paramset.keys()):
                self.REFSEQ_PROTEIN_BLAST = paramset["refseqProteinBlast"]
            if 'refseqGeneBlast' in list(paramset.keys()):
                self.REFSEQ_GENE_BLAST = paramset["refseqGeneBlast"]
            if 'phantomeBlast' in list(paramset.keys()):
                self.PHANTOME_BLAST = paramset["phantomeBlast"]
            if 'pvogsBlast' in list(paramset.keys()):
                self.PVOGS_BLAST = paramset["pvogsBlast"]
            if 'uniparcBlast' in list(paramset.keys()):
                self.UNIPARC_BLAST = paramset["uniparcBlast"]
            if 'swissprotBlast' in list(paramset.keys()):
                self.SWISSPROT_BLAST = paramset["swissprotBlast"]

    def setBlastFlavor(self,flavor):
        if(flavor.lower() == 'blastp'):
            self.blastFlavor = 'blastp'
        elif(flavor.lower() == 'blastn'):
            self.blastFlavor = 'blastn'
        elif(flavor.lower() == 'blastx'):
            self.blastFlavor = 'blastx'
        elif(flavor.lower() == 'tblastx'):
            self.blastFlavor = 'tblastx'
        else:
            if PHATE_WARNINGS == 'True':
                print("WARNING in blast module: Unrecognized blast flavor:", flavor)

    def setIdentityMin(self,identity):
        if (int(identity) >= 1) and (int(identity) <= 100):
            self.identityMin = int(identity)
        else:
            if PHATE_WARNINGS == 'True':
                print("WARNING in blast module: Identity minimum should be from 1 to 100")

    def setIdentitySelect(self,identity):
        if (int(identity) >= 10) and (int(identity) <= 100):
            self.identitySelect = int(identity)
        else:
            if PHATE_WARNINGS == 'True':
                print("WARNING in blast module: Identity select should be from 10 to 100. If this is insufficient, you may change constants in phate_blast.py.")

    def setEvalueMin(self,evalue):
        if (float(evalue) >= 0.0) and (float(evalue) <= 10):
            self.evalueMin = float(evalue)
        else:
            if PHATE_WARNINGS == 'True':
                print("WARNING in blast module: Evalue minimum shold be from 0.0 to 10.0. If this is insufficient, you may change constants in phate_blast.py.")

    def setEvalueSelect(self,evalue):
        if (float(evalue) >= 0.0000001) and (float(evalue) <= 10.0):
            self.evalueSelect = float(evalue)
        else:
            if PHATE_WARNINGS == 'True':
                print("WARNING in blast module: Evalue select should be from 0.0000001 to 10.0. If this is insufficient, you may change constants in phate_blast.py.")

    def setTopHitCount(self,number):
        if (int(number) >= 1 and int(number) <= HIT_COUNT_MAX):
            self.topHitCount = int(number)
        else:
            if PHATE_WARNINGS == 'True':
                print("WARNING in blast module: You may capture from 1 to", HIT_COUNT_MAX, "hits per query. If this is insufficient, you may change HIT_COUNT_MAX in phate_blast.py.")

    def setOutputFormat(self,outfmt):
        if outfmt == XML:
            self.outputFormat = XML
        elif outfmt == LIST:
            self.outputFormat = LIST 
        else:
            if PHATE_WARNINGS == 'True':
                print("WARNING in blast module: Select from output formats", LIST, "or", XML)

    def setScoreEdge(self,scoreEdge):
        if scoreEdge > 0.0 and scoreEdge < SCORE_EDGE_MAX:
            self.scoreEdge = scoreEdge
        else:
            if PHATE_WARNINGS == 'True':
                print("WARNING in blast module: Score edge should be between 0.0 and", SCORE_EDGE_MAX, "If this is insufficient, you may change SCORE_EDGE_MAX in phate_blast.py.")

    def setOverhang(self,overhang):
        if overhang > 0 and overhang < OVERHANG_MAX:
            self.overhang = overhang
        else:
            if PHATE_WARNINGS == 'True':
                print("WARNING in blast module: Overhang should be between 0 and", OVERHANG_MAX, "If this is insufficient, you may change OVERHANG_MAX in phate_blast.py.")

    def setGeneCallDir(self,geneCallDir):
        self.geneCallDir = geneCallDir
        GENE_CALL_DIR = geneCallDir

    def setBlastOutDir(self,blastOutDir):
        self.blastOutDir = blastOutDir
        BLAST_OUT_DIR = blastOutDir

    def setPVOGsOutDir(self,pVOGsOutDir):
        self.pVOGsOutDir = pVOGsOutDir
        PVOGS_OUT_DIR = pVOGsOutDir

    def getTopHits(self):
        return self.topHitList 

    ##### PERFORM BLAST

    def blast1fasta(self,fasta,outfile,database,dbName): # fasta is a phate_fastaSequence.fasta object

        # Write fasta sequence to temporary file
        fastaFile = self.blastOutDir + "temp.fasta"
        fastaFileH = open(fastaFile,"w")
        if fasta.sequentialHeader == "unknown":  # unchanged from default
            fasta.printFasta2file(fastaFileH,"blastHeader")
        else:
            fasta.printFasta2file(fastaFileH,"sequential")  # use sequential header format to avoid special chars issue
        fastaFileH.close()

        # Run blast
        if self.blastFlavor == 'blastn':
            command = BLAST_HOME + "blastn -query " + fastaFile + " -out " + outfile + \
                " -task blastn -db " + database + " -evalue " + str(self.evalueMin) + \
                " -best_hit_score_edge " + str(self.scoreEdge) + " -best_hit_overhang " + \
                str(self.overhang) + " -outfmt " + str(self.outputFormat) + " -perc_identity " + \
                str(self.identityMin) + " -max_target_seqs " + str(self.topHitCount) 

        elif self.blastFlavor == 'blastp': # Recall: You can't specificy %identity, but can filter afterward
            command = BLAST_HOME + "blastp -query " + fastaFile + " -out " + outfile + \
                " -task blastp -db " + database + " -evalue " + str(self.evalueMin) + \
                " -best_hit_score_edge " + str(self.scoreEdge) + " -best_hit_overhang " + \
                str(self.overhang) + " -outfmt " + str(self.outputFormat) + \
                " -max_target_seqs " + str(self.topHitCount)
        else:
            if PHATE_WARNINGS == 'True':
                print("ERROR in blast module: blast flavor not currently supported: ", self.blastFlavor)
            return
        if DEBUG:
            print("command is",command)
        result = os.system(command)

        # Capture result(s) and store as an annotation object for this fasta; Coded for -outfmt 7

        # For XML parsing
        MIN_LENGTH = 100 # shortest hit length that is of sufficient coverage
        doPrint = False
        count = 0
        countGood = 0
        nrList = {}
        organismList = [] # holds the organisms identified in the deflines
        hitList = [] # contains a list of hitDataSet objects
        hitDataSet = {
            "hitNumber"    : 0,
            "hitID"        : "",
            "hitDefline"   : "",
            "hitAccession" : "",
            "gi"           : "",  # gi identifier
            "hitLength"    : 0,
            "hitHSPs"      : [], # list of hspDataSet objects
            }
        hspDataSet = {
            "hspSequence"        : "",
            "hspNumber"          : 0,
            "hspScore"           : 0,
            "hspEvalue"          : 0,
            "queryStart"         : 0,
            "queryEnd"           : 0,
            "hitStart"           : 0,
            "hitEnd"             : 0,
            "hspAlignLen"        : 0,
            "hspBitScore"        : 0.0,
            "hspIdentity"        : 0,
            "hspPositives"       : 0,
            "hspGaps"            : 0,
            "hspPercentIdentity" : 0.0,
            }
        
        p_organism = re.compile('\[[\w\d\s]+\]') # organism names are enclosed in brackets in defline
        p_gi       = re.compile('^gi\|(\d*)|')   # gi number occurs at front of hit defline; capture number only via group(0)

        # Parse from XML-formatted blast output
        if self.outputFormat == XML:

            blastDatabase = ""
            blastProgram  = ""
            blastQuery    = ""

            # Load XML tree
            tree = ET()
            if DEBUG:
                print("Attempting to parse outfile into tree", outfile)
            tree.parse(outfile)
            if DEBUG:
                print("Attempt successful?")
            root = tree.getroot()
            for child in root:
                if child.tag == 'BlastOutput_program':
                    blastProgram = child.text
                if child.tag == 'BlastOutput_db':
                    blastDatabase = child.text
                if child.tag == 'BlastOutput_query-def':  # CHECK: was this not "Blast_query-def" in other blastp outputs?
                    blastQuery = child.text

            # Find hits and extract hit data
            for hit in root.getiterator('Hit'):  
                # reset
                speciesList = []; doPrint = False
                #queryFrom = 0; queryTo = 0; span = 0; percentIdentity = 0
                querySeq = ""; subjectSeq = ""; hitDefline = ""
                # increment
                count += 1
                # capture salient data
                nextHitDataSet = copy.deepcopy(hitDataSet)

                # Identify hit
                for hitData in hit:
                    if hitData.tag == "Hit_num":
                        nextHitDataSet["hitNumber"] = hitData.text
                    if hitData.tag == "Hit_id":
                        nextHitDataSet["hitID"] = hitData.text
                    if hitData.tag == "Hit_def":
                        defline = hitData.text
                        match = re.findall(p_gi,defline)
                        gi = match[0]
                        nextHitDataSet["gi"] = gi
                        nextHitDataSet["hitDefline"] = defline
                        match = re.findall(p_organism,defline)
                        if match:
                            for m in match:
                                organismList.append(m)
                    if hitData.tag == "Hit_accession":
                        nextHitDataSet["hitAccession"] = hitData.text
                    if hitData.tag == "Hit_len":
                        nextHitDataSet["hitLength"] = hitData.text
                    if hitData.tag == "Hit_hsps":
  
                        # Capture one or more hsps
                        for hsps in hitData:
                            nextHspDataSet = copy.deepcopy(hspDataSet)
                            for hsp in hsps:
                                if hsp.tag == "Hsp_hseq":
                                     nextHspDataSet["hspSequence"]= hsp.text
                                if hsp.tag == "Hsp_num":
                                     nextHspDataSet["hspNumber"] = hsp.text
                                if hsp.tag == "Hsp_bit_score":
                                     nextHspDataSet["hspBitScore"]= hsp.text
                                if hsp.tag == "Hsp_evalue":
                                     nextHspDataSet["hspEvalue"]= hsp.text
                                if hsp.tag == "Hsp_query-from":
                                     nextHspDataSet["queryStart"]= hsp.text
                                if hsp.tag == "Hsp_query-to":
                                     nextHspDataSet["queryEnd"]= hsp.text
                                if hsp.tag == "Hsp_hit-from":
                                     nextHspDataSet["hitStart"]= hsp.text
                                if hsp.tag == "Hsp_hit-to":
                                     nextHspDataSet["hitEnd"]= hsp.text
                                if hsp.tag == "Hsp_identity":
                                     nextHspDataSet["hspIdentity"]= hsp.text
                                if hsp.tag == "Hsp_positive":
                                     nextHspDataSet["hspPositives"]= hsp.text
                                if hsp.tag == "Hsp_gaps":
                                     nextHspDataSet["hspGaps"]= hsp.text
                                if hsp.tag == "Hsp_align-len":
                                     nextHspDataSet["hspAlignLen"] = hsp.text
                                     if int(nextHspDataSet["hspAlignLen"]) > 0:
                                         nextHspDataSet["hspPercentIdentity"] = int(nextHspDataSet["hspIdentity"])*100/int(nextHspDataSet["hspAlignLen"])
                                          
                            nextHitDataSet["hitHSPs"].append(nextHspDataSet)
                hitList.append(nextHitDataSet)

                if DEBUG:
                    print("Number of hits:", len(hitList))
                    print("Writing hit data to outfile", outfile)

                # Store new blast annotation
                newAnnotation = copy.deepcopy(annotation)
                newAnnotation.source = blastDatabase
                newAnnotation.method = self.blastFlavor
                newAnnotation.annotationType = "homology"
                newAnnotation.name  = nextHitDataSet["hitDefline"]               # subject
                newAnnotation.start = nextHitDataSet["hitHSPs"][0]["queryStart"] # query start
                newAnnotation.end   = nextHitDataSet["hitHSPs"][0]["queryEnd"]   # query end
                resultString = 'identity=' + str(nextHitDataSet["hitHSPs"][0]["hspPercentIdentity"]) 
                newAnnotation.annotationList.append(resultString)
                resultString = 'alignlen=' + str(nextHitDataSet["hitHSPs"][0]["hspAlignLen"]) 
                newAnnotation.annotationList.append(resultString)
                resultString = 'evalue='   + str(nextHitDataSet["hitHSPs"][0]["hspEvalue"]) 

                # If this is a pVOGs blast result, capture the pVOG identifiers in the annotation objecta
                match_pVOG = re.search('pvog',newAnnotation.source.lower())
                if match_pVOG: 
                    pVOGidList = re.findall('VOG\d+',newAnnotation.name)
                    for pVOGid in pVOGidList:
                        if pVOGid not in newAnnotation.pVOGlist: #** cez: bug fix 05 sept 2017
                            newAnnotation.pVOGlist.append(pVOGid)
                newAnnotation.annotationList.append(resultString)
 
                # Get DBXREFs, packed into annotation object's self.description field
                if DEBUG:
                    print("TESTING: newAnnotation", newAnnotation.printAll())
                    print("TESTING: database and dbName are", database, dbName)
                newAnnotation.link2databaseIdentifiers(database,dbName) # Get DBXREFs, packed into self.description
 
                # Add this completed annotation to growing list for this fasta
                fasta.annotationList.append(newAnnotation)

            if DEBUG: 
                print("Done parsing XML tree") 

        # Parse from LIST-formatted blast output
        elif self.outputFormat == LIST:
            columns = []; hitList = []
            outfileH = open(outfile,"r")
            bLines = outfileH.read().splitlines()

            for bLine in bLines:
                match_comment = re.search('^#',bLine)
                if match_comment:
                    continue
                if fasta.sequentialHeader == "unknown":  # unchanged from default
                    match_query = re.search(fasta.blastHeader,bLine)  # Use blast header to avoid Blast's truncation issue 
                else:  # Use the code-assigned benign header to avoid special chars issue
                    match_query = re.search(fasta.sequentialHeader,bLine)  
                if match_query:
                    hitList.append(bLine) 

            if hitList:
                for hitLine in hitList:

                    # Extract blast info from hitLine and stash into new annotation object
                    newAnnotation = copy.deepcopy(annotation) 
                    newAnnotation.source = database 
                    newAnnotation.method = self.blastFlavor 
                    newAnnotation.annotationType = "homology"
                    columns = hitLine.split('\t')
                    newAnnotation.name  = columns[1] # subject
                    newAnnotation.start = columns[6] # query start
                    newAnnotation.end   = columns[7] # query end
                    resultString = 'identity=' + columns[2]
                    newAnnotation.annotationList.append(resultString)
                    resultString = 'alignlen=' + columns[3]
                    newAnnotation.annotationList.append(resultString)
                    resultString = 'evalue=' + columns[10]
                    newAnnotation.annotationList.append(resultString)

                    # Get DBXREFs, packed into annotation object's self.description field
                    newAnnotation.link2databaseIdentifiers(database,dbName) # Get DBXREFs, packed into self.description

                    # Add this completed annotation to growing list for this fasta
                    fasta.annotationList.append(newAnnotation)
            else:
                if PHATE_MESSAGES == 'True':
                    print("Blast module says: No hit found for query", fasta.blastHeader, "against", database)    
            outfileH.close 

        # Requested blast output format not supported
        else:
            if PHATE_WARNINGS == 'True':
                print("WARNING in blast module: Output format", self.outputFormat, "not yet supported in phate_blast.py/blast1fasta(). Use blast out xml or list format for now.")

    def runBlast(self,fastaSet,dbType="protein"): # fastaSet is a phate_fastaSequence.multiFasta object

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
                print("WARNING in blast module: unrecognized database type in runBlast:", dbType)
            return
               
        # Set database variable, invoke blast for each fasta 
        database = ''

        if GENOME:
            if self.NCBI_VIRUS_BLAST:
                database = NCBI_VIRUS_BLAST_HOME
                dbName   = 'ncbi'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("Blast module says: Running NCBI blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_ncbi_" + str(count)
                    self.blast1fasta(fasta,outfile,database,dbName)

        if GENE:
            if self.REFSEQ_GENE_BLAST:
                database = REFSEQ_GENE_BLAST_HOME
                dbName   = 'refseqGene'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("Blast module says: Running Refseq gene blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_refseqGene_" + str(count)
                    self.blast1fasta(fasta,outfile,database,dbName)  #*** CONTROL

        if PROTEIN:
            if self.NR_BLAST:  
                database = NR_BLAST_HOME
                dbName   = 'nr'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("Blast module says: Running NR blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_nr_" + str(count)
                    self.blast1fasta(fasta,outfile,database,dbName)  #*** CONTROL

            if self.NCBI_VIRUS_PROTEIN_BLAST:  
                database = NCBI_VIRUS_PROTEIN_BLAST_HOME
                dbName   = 'ncbiVirusProtein'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("Blast module says: Running NCBI_VIRUS_PROTEIN blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_ncbiVirProt_" + str(count)
                    self.blast1fasta(fasta,outfile,database,dbName)  #*** CONTROL

            if self.REFSEQ_PROTEIN_BLAST:
                database = REFSEQ_PROTEIN_BLAST_HOME
                dbName   = 'refseqProtein'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("Blast module says: Running Refseq protein blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_refseqProtein_" + str(count)
                    self.blast1fasta(fasta,outfile,database,dbName)  #*** CONTROL

            if self.PHANTOME_BLAST:
                database = PHANTOME_BLAST_HOME
                dbName   = 'phantome'
                count = 0 
                if PHATE_PROGRESS == 'True':
                    print("Blast module says: Running PHANTOME blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_phantome_" +str(count)
                    self.blast1fasta(fasta,outfile,database,dbName)

            if self.KEGG_VIRUS_BLAST: 
                database = KEGG_VIRUS_BLAST_HOME
                dbName   = 'kegg'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("Blast module says: Running KEGG blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_kegg_" + str(count)
                    self.blast1fasta(fasta,outfile,database,dbName)

            if self.UNIPARC_BLAST: 
                database = UNIPARC_VIRUS_BLAST_HOME
                dbName   = 'uniparc'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("Blast module says: Running UNIPARC blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_uniparc_" + str(count)
                    self.blast1fasta(fasta,outfile,database,dbName)   #*** CONTROL

            if self.SWISSPROT_BLAST: 
                database = SWISSPROT_BLAST_HOME
                dbName   = 'swissprot'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("Blast module says: Running Swissprot blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_swissprot_" + str(count)
                    self.blast1fasta(fasta,outfile,database,dbName)   #*** CONTROL

            if self.PVOGS_BLAST:  
                database = PVOGS_BLAST_HOME
                dbName   = 'pVOGs'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("Blast module says: Running pVOGs blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_pvog_" + str(count)
                    self.blast1fasta(fasta,outfile,database,dbName)

                if PHATE_PROGRESS == 'True':
                    print("Blast module says: Done!")
                    print("Blast module says: Collecting and saving pVOG sequences corresponding to blast hit(s)")

                # Next you want to create pVOG fasta group files so user can do alignments
                # You need only one "alignment" file per pVOG group that the fasta hit (under blast cutoffs)
                # Capture pVOG.faa lines
                pVOGs_h = open(database,"r")
                pVOGlines = pVOGs_h.read().splitlines()
                if DEBUG:
                    print("There are", len(pVOGlines), "pVOG database lines to search")
                count = 0; countA = 0
                for fasta in fastaSet.fastaList:
                    pvogPrintedList = []  # keeps track of pVOGs that have already been printed for current fasta
                    count += 1 
                    countA = 0
                    for annot in fasta.annotationList:
                        for pVOG in annot.pVOGlist:  # There may be multiple annotations to inspect
                            # Avoid redundancy in printing pVOG groups for this fasta; only once per pVOG ID that was a blast hit
                            if pVOG not in pvogPrintedList:
                                pvogPrintedList.append(pVOG)  # Record this pVOG identifier as "done"
                                # create dynamic file name
                                countA += 1 
                                outfilePVOG = self.pVOGsOutDir + "pvogGroup_" + str(count) + '_' + str(countA) + '.faa' 
                                if DEBUG:
                                    print("Writing to outfile", outfilePVOG, "pVOG group", pVOG, "pvogPrintedList:", pvogPrintedList)
                                # open file and write current fasta pluse each corresponding pVOG fasta
                                outfilePVOG_h = open(outfilePVOG,'w')
                                outfilePVOG_h.write("%c%s\n%s\n" % ('>',fasta.header,fasta.sequence)) # write the current peptide fasta,
                                self.writePVOGsequences2file(outfilePVOG_h,pVOGlines,pVOG)                  # followed by the pVOG group
                                outfilePVOG_h.close()

        if CLEAN_RAW_DATA == 'True':
            if PHATE_PROGRESS == 'True':
                print("Blast module says: removing raw data files.")
            self.cleanBlastOutDir()

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

    def cleanBlastOutDir(self):  # Remove temporary files from BLAST_OUT_DIR
        #command = "ls " + BLAST_OUT_DIR  #*** FIX: list only files, not directories too
        command = "ls " + self.blastOutDir  #*** FIX: list only files, not directories too
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        (rawresult, err) = proc.communicate()
        result = rawresult.decode('utf-8')
        if DEBUG:
            print("Result of listing blast out dir,", self.blastOutDir) 
        fileList = str(result).split('\n')   # Python3
        for filename in fileList:
            file2delete = self.blastOutDir + filename
            if re.search('blast',file2delete):
                command = "rm " + file2delete
                proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
                (result, err) = proc.communicate()

    ##### PRINT METHODS

    def printParameters(self):
        print("Parameters:")
        print("   blast flavor:   ", self.blastFlavor)
        print("   identity min:   ", self.identityMin)
        print("   evalue min:     ", self.evalueMin)
        print("   identity select:", self.identitySelect)
        print("   evalue select:  ", self.evalueSelect)
        print("   topHitCount:    ", self.topHitCount)
        print("   scoreEdge:      ", self.scoreEdge)
        print("   overhang:       ", self.overhang)
        print("   outputFormat:   ", self.outputFormat)

    def printParameters2file(self,fileHandle):
        fileHandle.write("%s\n" % ("Parameters:"))
        fileHandle.write("%s%s\n" % ("   blast flavor ",self.blastFlavor))
        fileHandle.write("%s%s\n" % ("   identity min ",self.identityMin))
        fileHandle.write("%s%s\n" % ("   evalue min ",self.evalueMin))
        fileHandle.write("%s%s\n" % ("   identity select: ",self.identitySelect))
        fileHandle.write("%s%s\n" % ("   evalue select: ",self.evalueSelect))
        fileHandle.write("%s%s\n" % ("   topHitCount: ",self.topHitCount))
        fileHandle.write("%s%s\n" % ("   scoreEdge: ",self.scoreEdge))
        fileHandle.write("%s%s\n" % ("   overhang ",self.overhang))
        fileHandle.write("%s%s\n" % ("   outputFormat: ",self.outputFormat))

    def printAnnotations(self):  #*** Don't need this; don't use, because code writes directoy to fasta's annotation object
        print("Annotations:")
        if self.annotations:
            for annotation in self.annotations:
                print("   ", annotation.PrintAll())
        else:
            print("   There are no annotations")

    def printAll(self):
        self.printParameters()
        #self.printAnnotations()

    def printAll2file(self,fileHandle):
        fileHandle.write("%s%s\n" % ("blast flavor:   ", self.blastFlavor))
        fileHandle.write("%s%s\n" % ("identity min:   ", self.identityMin))
        fileHandle.write("%s%s\n" % ("evalue min:     ", self.evalueMin))
        fileHandle.write("%s%s\n" % ("identity select:", self.identitySelect))
        fileHandle.write("%s%s\n" % ("evalue select:  ", self.evalueSelect))

    ##### STATISTICS

    def calculateStatistics(self):
        pass


############################################################################
#
# Module:  phate_blast.py
#
# This class performs blast against various phage-related databases. 
#
# Programmer:  Carol Zhou
#
# Last Update:  05 August 2020
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

# Constants
# Note: sometimes blast+ ignores sort based on selected output format selected
HIT_SORT_CRITERION = 3 # for sorting blast hits
    # 0 = evalue
    # 1 = bit socre
    # 2 = total score
    # 3 = percent identity
    # 4 = query coverage
HSP_SORT_CRITERION = 3 # for sorting hsp's
    # 0 = evalue
    # 1 = score
    # 2 = query start
    # 3 = percent identity
    # 4 = subject start

# Get environment variables (set in phate_runPipeline.py)

#EMBOSS_PHATE_HOME             = os.environ["PHATE_EMBOSS_PHATE_HOME"]
BLAST_HOME                    = os.environ["PHATE_BLAST_HOME"]
NCBI_VIRUS_GENOME_BLAST_HOME  = os.environ["PHATE_NCBI_VIRUS_GENOME_BLAST_HOME"]
NCBI_VIRUS_PROTEIN_BLAST_HOME = os.environ["PHATE_NCBI_VIRUS_PROTEIN_BLAST_HOME"]
REFSEQ_PROTEIN_BLAST_HOME     = os.environ["PHATE_REFSEQ_PROTEIN_BLAST_HOME"]
REFSEQ_GENE_BLAST_HOME        = os.environ["PHATE_REFSEQ_GENE_BLAST_HOME"]
PVOGS_BLAST_HOME              = os.environ["PHATE_PVOGS_BLAST_HOME"]
VOGS_BLAST_HOME               = os.environ["PHATE_VOGS_BLAST_HOME"]    #*** To be deprecated
VOG_GENE_BLAST_HOME           = os.environ["PHATE_VOG_GENE_BLAST_HOME"]
VOG_PROTEIN_BLAST_HOME        = os.environ["PHATE_VOG_PROTEIN_BLAST_HOME"]
PHANTOME_BLAST_HOME           = os.environ["PHATE_PHANTOME_BLAST_HOME"]
KEGG_VIRUS_BLAST_HOME         = os.environ["PHATE_KEGG_VIRUS_BLAST_HOME"]
SWISSPROT_BLAST_HOME          = os.environ["PHATE_SWISSPROT_BLAST_HOME"]
PHAGE_ENZYME_BLAST_HOME       = os.environ["PHATE_PHAGE_ENZYME_BLAST_HOME"]
PFAM_BLAST_HOME               = os.environ["PHATE_PFAM_BLAST_HOME"]
SMART_BLAST_HOME              = os.environ["PHATE_SWISSPROT_BLAST_HOME"]
UNIPROT_BLAST_HOME            = os.environ["PHATE_UNIPROT_BLAST_HOME"]
NR_BLAST_HOME                 = os.environ["PHATE_NR_BLAST_HOME"]
CAZY_BLAST_BASE_DIR           = os.environ["PHATE_CAZY_BASE_DIR"]
CAZY_BLAST_HOME               = os.environ["PHATE_CAZY_BLAST_HOME"]
NCBI_TAXON_DIR                = os.environ["PHATE_NCBI_TAXON_DIR"]
CUSTOM_GENOME_BLAST_HOME      = os.environ["PHATE_CUSTOM_GENOME_BLAST_HOME"]
CUSTOM_GENE_BLAST_HOME        = os.environ["PHATE_CUSTOM_GENE_BLAST_HOME"]
CUSTOM_PROTEIN_BLAST_HOME     = os.environ["PHATE_CUSTOM_PROTEIN_BLAST_HOME"]

VOGS_ANNOTATION_FILE          = os.environ["PHATE_VOGS_ANNOTATION_FILE"]

# Blast parameters
MIN_BLASTP_IDENTITY           = os.environ["PHATE_MIN_BLASTP_IDENTITY"]   # Sets a lower limit
MIN_BLASTN_IDENTITY           = os.environ["PHATE_MIN_BLASTN_IDENTITY"]   # Sets a lower limit
MAX_BLASTP_HIT_COUNT          = os.environ["PHATE_MAX_BLASTP_HIT_COUNT"]  # Sets an upper limit
MAX_BLASTN_HIT_COUNT          = os.environ["PHATE_MAX_BLASTN_HIT_COUNT"]  # Sets an upper limit
BLASTP_IDENTITY_DEFAULT       = os.environ["PHATE_BLASTP_IDENTITY_DEFAULT"]
BLASTN_IDENTITY_DEFAULT       = os.environ["PHATE_BLASTN_IDENTITY_DEFAULT"]
BLASTP_HIT_COUNT_DEFAULT      = os.environ["PHATE_BLASTP_HIT_COUNT_DEFAULT"]
BLASTN_HIT_COUNT_DEFAULT      = os.environ["PHATE_BLASTN_HIT_COUNT_DEFAULT"]
SCORE_EDGE_MAX                = float(os.environ["PHATE_SCORE_EDGE_MAX"]) 
OVERHANG_MAX                  = int(os.environ["PHATE_OVERHANG_MAX"]) 
HIT_COUNT_MAX                 = int(os.environ["PHATE_HIT_COUNT_MAX"])    # Limit set in multiPhate.py

# Verbosity
CLEAN_RAW_DATA                = os.environ["PHATE_CLEAN_RAW_DATA"]
PHATE_WARNINGS                = os.environ["PHATE_PHATE_WARNINGS"]
PHATE_MESSAGES                = os.environ["PHATE_PHATE_MESSAGES"]
PHATE_PROGRESS                = os.environ["PHATE_PHATE_PROGRESS"]

# Other configurables 

GENE_CALL_DIR            = ""  # set by set method, via parameter list
BLAST_OUT_DIR            = ""  # set by set method, via parameter list
# Pvog and Vog directories hold fasta groups; other blast processes do not generate these.
PVOGS_OUT_DIR            = ""  # set by set method, via parameter list
VOGS_OUT_DIR             = ""  # set by set method, via parameter list

DEBUG  = False 
#DEBUG  = True 

# blast formats
XML  = 5
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
        self.scoreEdge                = SCORE_EDGE_MAX # default: BLAST recommended default
        self.overhang                 = OVERHANG_MAX   # default: BLAST recommended default
        self.blastThreads             = 1         # default: run in serial
        self.outputFormat             = 5         # default: XML output
        self.blastAnnotations         = []        # List of phate_annotation objects; blast output get temporarily stored here
        # move to hit class:  self.topHitList     = []        # 
        self.geneCallDir              = ""        # needs to be set
        self.blastOutDir              = ""        # needs to be set
        self.pVOGsOutDir              = ""        # needs to be set
        self.VOGsOutDir               = ""        # needs to be set
        self.NCBI_VIRUS_GENOME_BLAST  = False     # For this and below: default is not running process; change by set method
        self.NCBI_VIRUS_PROTEIN_BLAST = False     #  
        self.NR_BLAST                 = False     #  
        self.KEGG_VIRUS_BLAST         = False     # 
        self.REFSEQ_PROTEIN_BLAST     = False     #
        self.REFSEQ_GENE_BLAST        = False     #*** To be deprecated 
        self.PHANTOME_BLAST           = False     # 
        self.PVOGS_BLAST              = False     # 
        self.VOGS_BLAST               = False     #*** To be deprecated
        self.VOG_GENE_BLAST           = False     #  
        self.VOG_PROTEIN_BLAST        = False     #  
        self.SWISSPROT_BLAST          = False     #  
        self.PHAGE_ENZYME_BLAST       = False     # not yet in service
        self.CAZY_BLAST               = False     # not yet in service
        self.CUSTOM_GENOME_BLAST      = False     #
        self.CUSTOM_GENE_BLAST        = False     #
        self.CUSTOM_PROTEIN_BLAST     = False     #
        
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
            if 'blastThreads' in list(paramset.keys()):
                self.setBlastThreads(paramset['blastThreads'])
            if 'geneCallDir' in list(paramset.keys()):
                self.setGeneCallDir(paramset['geneCallDir'])
            if 'blastOutDir' in list(paramset.keys()):
                self.setBlastOutDir(paramset['blastOutDir'])
            if 'pvogsOutDir' in list(paramset.keys()):
                self.setPVOGsOutDir(paramset['pvogsOutDir'])
            if 'vogsOutDir' in list(paramset.keys()):
                self.setVOGsOutDir(paramset['vogsOutDir'])
            if 'ncbiVirusGenomeBlast' in list(paramset.keys()):
                self.NCBI_VIRUS_GENOME_BLAST = paramset["ncbiVirusGenomeBlast"]
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
            if 'vogsBlast' in list(paramset.keys()):  #*** To be deprecated
                self.VOGS_BLAST = paramset["vogsBlast"]
            if 'vogGeneBlast' in list(paramset.keys()):
                self.VOG_GENE_BLAST = paramset["vogGeneBlast"]
            if 'vogProteinBlast' in list(paramset.keys()):
                self.VOG_PROTEIN_BLAST = paramset["vogProteinBlast"]
            if 'swissprotBlast' in list(paramset.keys()):
                self.SWISSPROT_BLAST = paramset["swissprotBlast"]
            if 'phageEnzymeBlast' in list(paramset.keys()):
                self.PHAGE_ENZYME_BLAST = paramset["phageEnzymeBlast"]
            if 'cazyBlast' in list(paramset.keys()):
                self.CAZY_BLAST = paramset["cazyBlast"]
            if 'customGenomeBlast' in list(paramset.keys()):
                self.CUSTOM_GENOME_BLAST = paramset["customGenomeBlast"]
            if 'customGeneBlast' in list(paramset.keys()):
                self.CUSTOM_GENE_BLAST = paramset["customGeneBlast"]
            if 'customProteinBlast' in list(paramset.keys()):
                self.CUSTOM_PROTEIN_BLAST = paramset["customProteinBlast"]

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
                print("phate_blast says, WARNING: Unrecognized or not-yet implemented blast flavor:", flavor)

    def setIdentityMin(self,identity):
        if (int(identity) >= 1) and (int(identity) <= 100):
            self.identityMin = int(identity)
        else:
            if PHATE_WARNINGS == 'True':
                print("phate_blast says, WARNING: Identity minimum should be from 1 to 100")

    def setIdentitySelect(self,identity):
        if (int(identity) >= 10) and (int(identity) <= 100):
            self.identitySelect = int(identity)
        else:
            if PHATE_WARNINGS == 'True':
                print("phate_blast says, WARNING: Identity select should be from 10 to 100. If this is insufficient, you may change constants in phate_blast.py.")

    def setEvalueMin(self,evalue):
        if (float(evalue) >= 0.0) and (float(evalue) <= 10):
            self.evalueMin = float(evalue)
        else:
            if PHATE_WARNINGS == 'True':
                print("phate_blast says, WARNING: Evalue minimum shold be from 0.0 to 10.0. If this is insufficient, you may change constants in phate_blast.py.")

    def setEvalueSelect(self,evalue):
        if (float(evalue) >= 0.0000001) and (float(evalue) <= 10.0):
            self.evalueSelect = float(evalue)
        else:
            if PHATE_WARNINGS == 'True':
                print("phate_blast says, WARNING: Evalue select should be from 0.0000001 to 10.0. If this is insufficient, you may change constants in phate_blast.py.")

    def setTopHitCount(self,number):
        if (int(number) >= 1 and int(number) <= HIT_COUNT_MAX):
            self.topHitCount = int(number)
        else:
            if PHATE_WARNINGS == 'True':
                print("phate_blast says, WARNING: You may capture from 1 to", HIT_COUNT_MAX, "hits per query. If this is insufficient, you may change HIT_COUNT_MAX in phate_blast.py.")

    def setOutputFormat(self,outfmt):
        if outfmt == XML:
            self.outputFormat = XML
        elif outfmt == LIST:
            self.outputFormat = LIST 
        else:
            if PHATE_WARNINGS == 'True':
                print("phate_blast says, WARNING: Select from output formats", LIST, "or", XML)

    def setScoreEdge(self,scoreEdge):
        if scoreEdge > 0.0 and scoreEdge < SCORE_EDGE_MAX:
            self.scoreEdge = scoreEdge
        else:
            if PHATE_WARNINGS == 'True':
                print("phate_blast says, WARNING: Score edge should be between 0.0 and", SCORE_EDGE_MAX, "If this is insufficient, you may change SCORE_EDGE_MAX in phate_blast.py.")

    def setOverhang(self,overhang):
        if overhang > 0 and overhang < OVERHANG_MAX:
            self.overhang = overhang
        else:
            if PHATE_WARNINGS == 'True':
                print("phate_blast says, WARNING: Overhang should be between 0 and", OVERHANG_MAX, "If this is insufficient, you may change OVERHANG_MAX in phate_blast.py.")

    def setBlastThreads(self,blastThreads):
        if int(blastThreads) >= 1:
            self.blastThreads = int(blastThreads)
        else:
            self.blastThreads = 1

    def setGeneCallDir(self,geneCallDir):
        self.geneCallDir = geneCallDir
        GENE_CALL_DIR = geneCallDir

    def setBlastOutDir(self,blastOutDir):
        self.blastOutDir = blastOutDir
        BLAST_OUT_DIR = blastOutDir

    def setPVOGsOutDir(self,pVOGsOutDir):
        self.pVOGsOutDir = pVOGsOutDir
        PVOGS_OUT_DIR = pVOGsOutDir

    def setVOGsOutDir(self,VOGsOutDir):
        self.VOGsOutDir = VOGsOutDir
        VOGS_OUT_DIR = VOGsOutDir

    def getTopHits(self):
        return self.topHitList 

    ##### PERFORM BLAST

    # Blast one fasta sequence. This method is called by runBlast 
    def blast1fasta(self,fasta,outfile,database,dbName): # fasta is a phate_fastaSequence.fasta object
        command = ""

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
                str(self.identityMin) + " -max_target_seqs " + str(self.topHitCount) + \
                " -num_threads " + str(self.blastThreads)

        elif self.blastFlavor == 'blastp': # Recall: You can't specificy %identity, but can filter afterward
            command = BLAST_HOME + "blastp -query " + fastaFile + " -out " + outfile + \
                " -task blastp -db " + database + " -evalue " + str(self.evalueMin) + \
                " -best_hit_score_edge " + str(self.scoreEdge) + " -best_hit_overhang " + \
                str(self.overhang) + " -outfmt " + str(self.outputFormat) + \
                " -max_target_seqs " + str(self.topHitCount) + \
                " -num_threads " + str(self.blastThreads)
                #" -max_target_seqs " + str(self.topHitCount) + \
                #" -sorthits " + str(HIT_SORT_CRITERION) + \
                #" -sorthsps " + str(HSP_SORT_CRITERION)
        else:
            if PHATE_WARNINGS == 'True':
                print("phate_blast says, ERROR: blast flavor not currently supported: ", self.blastFlavor)
            return
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
            tree.parse(outfile)
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

                # Store new blast annotation
                newAnnotation = copy.deepcopy(annotation)
                newAnnotation.source         = blastDatabase
                newAnnotation.method         = self.blastFlavor
                newAnnotation.annotationType = "homology"
                newAnnotation.category       = "sequence"
                newAnnotation.name           = nextHitDataSet["hitDefline"]               # subject header (possibly truncated by blast)
                newAnnotation.start          = nextHitDataSet["hitHSPs"][0]["queryStart"] # query start
                newAnnotation.end            = nextHitDataSet["hitHSPs"][0]["queryEnd"]   # query end
                resultString                 = 'identity=' + str(round(nextHitDataSet["hitHSPs"][0]["hspPercentIdentity"],2)) 
                newAnnotation.annotationList.append(resultString)
                resultString                 = 'alignlen=' + str(nextHitDataSet["hitHSPs"][0]["hspAlignLen"]) 
                newAnnotation.annotationList.append(resultString)
                resultString = 'evalue='   + str(nextHitDataSet["hitHSPs"][0]["hspEvalue"]) 
                MEETS_IDENTITY_CUTOFF = False
                if float(nextHitDataSet["hitHSPs"][0]["hspPercentIdentity"]) >= float(self.identityMin):
                    MEETS_IDENTITY_CUTOFF = True

                # CHECK THIS: code adapted assuming same structure of pVOGs vs. VOGs
                # If this is a pVOG/VOG blast result, capture the VOG identifiers in the annotation object
                match_pVOG = re.search('pvog',newAnnotation.source.lower())  # First check if it's a pVOG
                if not match_pVOG:                                           # If it's not a pVOG, then it might be a VOG
                    match_VOG  = re.search('vog', newAnnotation.source.lower())
                if match_pVOG or match_VOG: 
                    # Note that structure of pVOG annotation information differs from that of VOG.
                    # VOG hit header has VOGid(s) + proteinID; pVOG hit header has VOGid(s) + proteinID + description
                    VOGidList = re.findall('VOG\d+',newAnnotation.name)  # name holds to hit's header, which has func dscr for pVOG
                    for VOGid in VOGidList:
                        if VOGid not in newAnnotation.VOGlist:  # ensure list is non-redundant
                            newAnnotation.VOGlist.append(VOGid)
                newAnnotation.annotationList.append(resultString)
 
                # Get DBXREFs, packed into annotation object's self.description field
                newAnnotation.link2databaseIdentifiers(database,dbName) # Get DBXREFs, packed into self.description
 
                # Add this completed annotation to growing list for this fasta
                if MEETS_IDENTITY_CUTOFF:
                    fasta.annotationList.append(newAnnotation)

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
                    resultString = 'identity=' + round(columns[2],2)
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
                    print("phate_blast says, No hit found for query", fasta.blastHeader, "against", database)    
            outfileH.close 

        # Requested blast output format not supported
        else:
            if PHATE_WARNINGS == 'True':
                print("phate_blast says, WARNING: Output format", self.outputFormat, "not yet supported in phate_blast.py/blast1fasta(). Use blast out xml or list format for now.")

    # Run BLAST over a set of fasta sequences. This method calls blast1fasta for each sequence.
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
                print("phate_blast says, WARNING: unrecognized database type in runBlast:", dbType)
            return
               
        # Set database variable, invoke blast for each fasta 
        database = ''

        if GENOME:
            if self.NCBI_VIRUS_GENOME_BLAST:
                database = NCBI_VIRUS_GENOME_BLAST_HOME
                dbName = 'ncbiVirusGenome'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("phate_blast says, Running NCBI blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_ncbiVirGenome_" + str(count)
                    self.blast1fasta(fasta,outfile,database,dbName)

            elif self.CUSTOM_GENOME_BLAST:
                database = CUSTOM_GENOME_BLAST_HOME
                dbName = 'customGenome'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("phate_blast says, Running Custom Genome blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_customGenome_" + str(count)
                    self.blast1fasta(fasta,outfile,database,dbName)

        if GENE:
            if self.REFSEQ_GENE_BLAST:  #*** To be deprecated
                database = REFSEQ_GENE_BLAST_HOME
                dbName   = 'refseqGene'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("phate_blast says, Running Refseq gene blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_refseqGene_" + str(count)
                    self.blast1fasta(fasta,outfile,database,dbName)  

            elif self.VOG_GENE_BLAST:
                database = VOG_GENE_BLAST_HOME
                dbName   = 'vogGene'
                count    = 0
                if PHATE_PROGRESS == 'True':
                    print("phate_blast says, Running VOG Gene blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_vogGene_" + str(count)
                    self.blast1fasta(fasta,outfile,database,dbName) 

                if PHATE_PROGRESS == 'True':
                    print("phate_blast says, Blasting completed.")
                    print("phate_blast says, Collecting and saving VOG Gene sequences corresponding to blast hit(s)")

                # Next you want to create VOG fasta group files so user can do alignments
                # You need only one "alignment" file per VOG group that the fasta hit (under blast cutoffs)
                # Capture VOG.faa lines
                VOGs_h = open(database,"r")
                VOGlines = VOGs_h.read().splitlines()
                count = 0; countA = 0
                for fasta in fastaSet.fastaList:
                    vogPrintedList = []  # keeps track of VOGs that have already been printed for current fasta
                    count += 1 
                    countA = 0
                    for annot in fasta.annotationList:
                        for VOG in annot.VOGlist:  # There may be multiple annotations to inspect
                            # Avoid redundancy in printing VOG groups for this fasta; only once per VOG ID that was a blast hit
                            if VOG not in vogPrintedList:
                                vogPrintedList.append(VOG)  # Record this pVOG identifier as "done"
                                # create dynamic file name
                                countA += 1 
                                outfileVOG = self.VOGsOutDir + "vogGeneGroup_" + str(count) + '_' + str(countA) + '.fnt' 
                                # open file and write current fasta plus each corresponding VOG fasta
                                outfileVOG_h = open(outfileVOG,'w')
                                outfileVOG_h.write("%c%s\n%s\n" % ('>',fasta.header,fasta.sequence)) # write the current peptide fasta,
                                self.writeVOGsequences2file(outfileVOG_h,VOGlines,VOG)               # followed by the VOG group
                                outfileVOG_h.close()

            elif self.CUSTOM_GENE_BLAST:
                database = CUSTOM_GENE_BLAST_HOME
                dbName   = 'customGene'
                count    = 0
                if PHATE_PROGRESS == 'True':
                    print("phate_blast says, Running Custom Gene blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_customGene_" + str(count)
                    self.blast1fasta(fasta,outfile,database,dbName)

        if PROTEIN:
            if self.NR_BLAST:  
                database = NR_BLAST_HOME
                dbName   = 'nr'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("phate_blast says, Running NR blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_nr_" + str(count)
                    self.blast1fasta(fasta,outfile,database,dbName)  

            if self.NCBI_VIRUS_PROTEIN_BLAST:  
                database = NCBI_VIRUS_PROTEIN_BLAST_HOME
                dbName   = 'ncbiVirusProtein'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("phate_blast says, Running NCBI_VIRUS_PROTEIN blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_ncbiVirProt_" + str(count)
                    self.blast1fasta(fasta,outfile,database,dbName) 

            if self.REFSEQ_PROTEIN_BLAST:
                database = REFSEQ_PROTEIN_BLAST_HOME
                dbName   = 'refseqProtein'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("phate_blast says, Running Refseq protein blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_refseqProtein_" + str(count)
                    self.blast1fasta(fasta,outfile,database,dbName) 

            if self.PHANTOME_BLAST:
                database = PHANTOME_BLAST_HOME
                dbName   = 'phantome'
                count = 0 
                if PHATE_PROGRESS == 'True':
                    print("phate_blast says, Running PHANTOME blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_phantome_" +str(count)
                    self.blast1fasta(fasta,outfile,database,dbName)

            if self.KEGG_VIRUS_BLAST: 
                database = KEGG_VIRUS_BLAST_HOME
                dbName   = 'kegg'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("phate_blast says, Running KEGG blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_kegg_" + str(count)
                    self.blast1fasta(fasta,outfile,database,dbName)

            if self.SWISSPROT_BLAST: 
                database = SWISSPROT_BLAST_HOME
                dbName   = 'swissprot'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("phate_blast says, Running Swissprot blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_swissprot_" + str(count)
                    self.blast1fasta(fasta,outfile,database,dbName)  

            if self.PHAGE_ENZYME_BLAST: 
                database = PHAGE_ENZYME_BLAST_HOME
                dbName   = 'phageEnzyme'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("phate_blast says, Running phageEnzyme blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_phageEnz_" + str(count)
                    self.blast1fasta(fasta,outfile,database,dbName) 

            if self.PVOGS_BLAST:  
                database = PVOGS_BLAST_HOME
                dbName   = 'pVOGs'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("phate_blast says, Running pVOGs blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_pvog_" + str(count)
                    self.blast1fasta(fasta,outfile,database,dbName)

                if PHATE_PROGRESS == 'True':
                    print("phate_blast says, Blasting completed.")
                    print("phate_blast says, Collecting and saving pVOG sequences corresponding to blast hit(s)")

                # Next you want to create pVOG fasta group files so user can do alignments
                # You need only one "alignment" file per pVOG group that the fasta hit (under blast cutoffs)
                # Capture pVOG.faa lines
                pVOGs_h = open(database,"r")
                pVOGlines = pVOGs_h.read().splitlines()
                count = 0; countA = 0
                for fasta in fastaSet.fastaList:
                    pvogPrintedList = []  # keeps track of pVOGs that have already been printed for current fasta
                    count += 1 
                    countA = 0
                    for annot in fasta.annotationList:
                        for pVOG in annot.VOGlist:  # There may be multiple annotations to inspect
                            # Avoid redundancy in printing pVOG groups for this fasta; only once per pVOG ID that was a blast hit
                            if pVOG not in pvogPrintedList:
                                pvogPrintedList.append(pVOG)  # Record this pVOG identifier as "done"
                                # create dynamic file name
                                countA += 1 
                                outfilePVOG = self.pVOGsOutDir + "pvogGroup_" + str(count) + '_' + str(countA) + '.faa' 
                                # open file and write current fasta pluse each corresponding pVOG fasta
                                outfilePVOG_h = open(outfilePVOG,'w')
                                outfilePVOG_h.write("%c%s\n%s\n" % ('>',fasta.header,fasta.sequence)) # write the current peptide fasta,
                                self.writePVOGsequences2file(outfilePVOG_h,pVOGlines,pVOG)                  # followed by the pVOG group
                                outfilePVOG_h.close()

            if self.VOGS_BLAST:  #*** To be deprecated
                database = VOGS_BLAST_HOME
                dbName   = 'VOGs'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("phate_blast says, Running VOGs blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_vog_" + str(count)
                    self.blast1fasta(fasta,outfile,database,dbName)

                if PHATE_PROGRESS == 'True':
                    print("phate_blast says, Blasting completed.")
                    print("phate_blast says, Collecting and saving VOG sequences corresponding to blast hit(s)")

                # Next you want to create VOG fasta group files so user can do alignments
                # You need only one "alignment" file per VOG group that the fasta hit (under blast cutoffs)
                # Capture VOG.faa lines
                VOGs_h = open(database,"r")
                VOGlines = VOGs_h.read().splitlines()
                count = 0; countA = 0
                for fasta in fastaSet.fastaList:
                    vogPrintedList = []  # keeps track of pVOGs that have already been printed for current fasta
                    count += 1 
                    countA = 0
                    for annot in fasta.annotationList:
                        for VOG in annot.VOGlist:  # There may be multiple annotations to inspect
                            # Avoid redundancy in printing VOG groups for this fasta; only once per VOG ID that was a blast hit
                            if VOG not in vogPrintedList:
                                vogPrintedList.append(VOG)  # Record this VOG identifier as "done"
                                # create dynamic file name
                                countA += 1 
                                outfileVOG = self.VOGsOutDir + "vogGroup_" + str(count) + '_' + str(countA) + '.faa' 
                                # open file and write current fasta pluse each corresponding VOG fasta
                                outfileVOG_h = open(outfileVOG,'w')
                                outfileVOG_h.write("%c%s\n%s\n" % ('>',fasta.header,fasta.sequence)) # write the current peptide fasta,
                                self.writeVOGsequences2file(outfileVOG_h,VOGlines,VOG)               # followed by the VOG group
                                outfileVOG_h.close()

            if self.VOG_PROTEIN_BLAST:  
                database = VOG_PROTEIN_BLAST_HOME
                dbName   = 'vogProtein'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("phate_blast says, Running VOG Protein blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_vogProtein_" + str(count)
                    self.blast1fasta(fasta,outfile,database,dbName)

                if PHATE_PROGRESS == 'True':
                    print("phate_blast says, Blasting completed.")
                    print("phate_blast says, Collecting and saving VOG Protein sequences corresponding to blast hit(s)")

                # Next you want to create VOG fasta group files so user can do alignments
                # You need only one "alignment" file per VOG group that the fasta hit (under blast cutoffs)
                # Capture VOG.faa lines
                VOGs_h = open(database,"r")
                VOGlines = VOGs_h.read().splitlines()
                count = 0; countA = 0
                for fasta in fastaSet.fastaList:
                    vogPrintedList = []  # keeps track of VOGs that have already been printed for current fasta
                    count += 1 
                    countA = 0
                    for annot in fasta.annotationList:
                        for VOG in annot.VOGlist:  # There may be multiple annotations to inspect
                            # Avoid redundancy in printing VOG groups for this fasta; only once per VOG ID that was a blast hit
                            if VOG not in vogPrintedList:
                                vogPrintedList.append(VOG)  # Record this VOG identifier as "done"
                                # create dynamic file name
                                countA += 1 
                                outfileVOG = self.VOGsOutDir + "vogProtGroup_" + str(count) + '_' + str(countA) + '.faa' 
                                # open file and write current fasta pluse each corresponding VOG fasta
                                outfileVOG_h = open(outfileVOG,'w')
                                outfileVOG_h.write("%c%s\n%s\n" % ('>',fasta.header,fasta.sequence)) # write the current peptide fasta,
                                self.writeVOGsequences2file(outfileVOG_h,VOGlines,VOG)               # followed by the VOG group
                                outfileVOG_h.close()

            if self.CAZY_BLAST: 
                database = CAZY_BLAST_HOME
                dbName   = 'cazy'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("phate_blast says, Running CAZy blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_cazy_" + str(count)
                    self.blast1fasta(fasta,outfile,database,dbName) 
                if PHATE_PROGRESS == 'True':
                    print("phate_blast says, Blasting completed.")
                    print("phate_blast says, Collecting and saving CAZy EC numbers and descriptions corresponding to blast hit(s)")

            if self.CUSTOM_PROTEIN_BLAST: 
                database = CUSTOM_PROTEIN_BLAST_HOME
                dbName   = 'customProtein'
                count = 0
                if PHATE_PROGRESS == 'True':
                    print("phate_blast says, Running Custom Protein blast:", database, dbName)
                for fasta in fastaSet.fastaList:
                    count += 1
                    outfile = self.blastOutDir + self.blastFlavor + "_customProtein_" + str(count)
                    self.blast1fasta(fasta,outfile,database,dbName) 

        if CLEAN_RAW_DATA == 'True':
            self.cleanBlastOutDir()

    def writePVOGsequences2file(self,FILE_H,lines,pVOG):
        pVOGheader = ""; pVOGsequence = ""; GET_SEQ = False
        for i in range(0,len(lines)-1):
            nextLine = lines[i]
            match_header = re.search('>',nextLine)
            match_pVOG = re.search(pVOG,nextLine)

            if match_header:   
                if GET_SEQ:
                    FILE_H.write("%s\n%s\n" % (pVOGheader,pVOGsequence))
                    pVOGheader = ""; pVOGsequence = ""
                    GET_SEQ = False

                if match_pVOG:
                    pVOGheader = nextLine
                    GET_SEQ = True

            elif GET_SEQ:
                pVOGsequence = pVOGsequence + nextLine

    #*** CHECK THIS
    #*** This method is modeled after writePVOGsequences2file, pending data is structured same as pVOGs
    def writeVOGsequences2file(self,FILE_H,lines,VOG):
        VOGheader = ""; VOGsequence = ""; GET_SEQ = False
        for i in range(0,len(lines)-1):
            nextLine = lines[i]
            match_header = re.search('>',nextLine)
            match_VOG    = re.search(VOG,nextLine)

            if match_header:   
                if GET_SEQ:
                    FILE_H.write("%s\n%s\n" % (VOGheader,VOGsequence))
                    VOGheader = ""; VOGsequence = ""
                    GET_SEQ = False

                if match_VOG:
                    VOGheader = nextLine
                    GET_SEQ = True

            elif GET_SEQ:
                VOGsequence = VOGsequence + nextLine

    def cleanBlastOutDir(self):  # Remove temporary files from BLAST_OUT_DIR
        if PHATE_PROGRESS == 'True':
            print("phate_blast says, cleanBlastOutDir(): Removing raw blast output files")
        command = "ls " + self.blastOutDir
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        (rawresult, err) = proc.communicate()
        result = rawresult.decode('utf-8')
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


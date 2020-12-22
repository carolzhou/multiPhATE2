###############################################################
# Name: cgp_genomeSequence.py
#
# Programmer: Carol L. Ecale Zhou
#
# Last Update: 17 December 2020
# 
# Description:
# Module comprising data structures for organizing genome information
# Note:  EMBOSS messes with fasta headers; therefore, I am putting minimal info in the header and using only '/'
#    as a separator. EMBOSS adds "_1" to the end of a simple fasta header--can't prevent this. Using other 
#    separators wreaks havoc in terms of the header EMBOSS produces.
#
# Classes and methods:
#    genome
#        setCodeBaseDir
#        setOutputDir
#        setName
#        setGenomeType
#        setFilename
#        setSpecies
#        getCGCsubsequence
#        getSubsequence
#        getSubsequenceWithFlank
#        processGeneCalls
#        addContig
#        addContigs
#        addGene
#        addProtein
#        addAnnotation
#        countAllAnnotations
#        printGenomeData
#        printGenomeData_tab
#        printGenomeData2file_tab
#        printGenomeData2file
#        printGenomeData2file_GFF
#        printAll
#        printAll2file
#        printGenes
#        printFastas2file
#        translateGene
#        translateGenes
#        makeBlastDB
#        write2proteinSet
#        restorSlashesAfterEMBOSS
#        cleanUpAfterEmboss
#
#########################################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE GPL3 LICENSE. SEE INCLUDED FILE GPL-3.PDF FOR DETAILS.

import string
from Bio.Seq import Seq
import cgp_fastaSequence as fastaSequence
import cgp_annotation as annotation
import re, os, copy
import subprocess
import random
import time

#RANDOM_LIMIT = 10

BLAST_HOME          = os.environ["PHATE_BLAST_HOME"] 
EMBOSS_PHATE_HOME   = os.environ["PHATE_EMBOSS_PHATE_HOME"] 
CODE_BASE_DIR       = ""
OUTPUT_DIR          = ""

# Verbosity

CLEAN_RAW_DATA_STRING   = os.environ["PHATE_CLEAN_RAW_DATA"]
PHATE_PROGRESS_STRING   = os.environ["PHATE_PHATE_PROGRESS"]
PHATE_MESSAGES_STRING   = os.environ["PHATE_PHATE_MESSAGES"]
PHATE_WARNINGS_STRING   = os.environ["PHATE_PHATE_WARNINGS"]

CLEAN_RAW_DATA = False
PHATE_PROGRESS = False
PHATE_MESSAGES = False
PHATE_WARNINGS = False

if CLEAN_RAW_DATA_STRING.lower() == 'true':
    CLEAN_RAW_DATA = True
if PHATE_PROGRESS_STRING.lower() == 'true':
    PHATE_PROGRESS = True
if PHATE_MESSAGES_STRING.lower() == 'true':
    PHATE_MESSAGES = True
if PHATE_WARNINGS_STRING.lower() == 'true':
    PHATE_WARNINGS = True

#DEBUG    = True 
DEBUG    = False

# For GFF output
GFF_COMMENT = "##gff-version 3"

# Templates
annotationObj = annotation.annotationRecord()
fastaObj      = fastaSequence.fasta()

# Reverse complement
#complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB') # Python 2
#complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')     # Python 3

class genome(object):

    # Class genome uses fasta objects to represent a genome/chromosome or a plasmid and its genes and proteins.
    # The chromosome's contig set, gene set, and protein set are each represented by multi-fasta objects.

    def __init__(self):
        self.filename                = ""          # name of file that contains the genome sequence
        self.genomeType              = "unknown" # typically, "chromosome", "plasmid", "virus", "segmented", "satellite"
        self.name                    = "unknown" # common name (e.g., Drosophila)
        self.genomeName              = "unkknown" # more specific name (e.g., 'Yp phage YP123')
        self.species                 = "unknown" # Latin name (e.g., Drosophila melanogaster)
        self.annotationList          = []    # using cgp_annotation.py classes
        self.contigSet               = fastaSequence.multiFasta()
        self.geneSet                 = fastaSequence.multiFasta()
        self.proteinSet              = fastaSequence.multiFasta()
        self.contigSet.moleculeType  = 'contig'
        self.contigSet.sequenceType  = 'nt'
        self.geneSet.moleculeType    = 'gene'
        self.geneSet.sequenceType    = 'nt'
        self.proteinSet.moleculeType = 'peptide'
        self.proteinSet.sequenceType = 'aa'
        self.codeBaseDir             = ""       # needs to be set
        self.outputDir               = ""       # needs to be set

    # GET and SET
 
    def setCodeBaseDir(self,codeBaseDir):
        self.codeBaseDir = codeBaseDir
        CODE_BASE_DIR = codeBaseDir

    def setOutputDir(self,outputDir):
        self.outputDir = outputDir
        OUTPUT_DIR = outputDir

    def setName(self,name):
        self.name = name

    def setGenomeType(self,genomeType):
        self.genomeType = genomeType

    def setFilename(self,filename):
        self.filename = filename

    def setSpecies(self,species):
        self.species = species

    def getCGCsubsequence(self,start,end,strand,contig):  # Tailored for output from CGCparser.py
        subSeq = ""
        for fa in self.contigSet.fastaList:
            if fa.header == contig:
                subSeq = fa.getSubsequence(int(start)-1,int(end)) #*** ???
            else:
                if PHATE_WARNINGS == 'True':
                    print("cgp_genomeSequence says, WARNING: fa.header", fa.header, "did not match contig", contig)
        return subSeq           

    def getSubsequence(self,start,end,contig):  # Note: tailored to RAST
        subSeq = ""
        for fa in self.contigSet.fastaList:
            if contig == fa.shortHeader: # RAST truncates after 1st space
                subSeq = fa.getSubsequence(int(start)-1,int(end)) # recall: numerbing starts w/zero
        return subSeq

    def getSubsequenceWithFlank(self,start,end,contig,flank):
        subSeq = ""
        flankedStart = start - flank
        flankedEnd   = end   + flank
        if contig == "":  # take subsequence from the very first contig
            fa = self.contigSet.fastaList[0]
            subSeq = fa.getSubsequence(int(flankedStart)-1,int(flankedEnd)) # recall: numerbing starts w/zero
        else:
            for fa in self.contigSet.fastaList:
                if contig == fa.shortHeader: # RAST truncates after 1st space
                    subSeq = fa.getSubsequence(int(flankedStart)-1,int(flankedEnd)) # recall: numerbing starts w/zero
        return subSeq

    # INPUT/PROCESS SEQUENCES 

    def processGeneCalls(self,geneCallInfo,GENE_CALL_H):  # GeneCalls formatted by CGCparser.py
        # Gene calls: tabbed in 6 columns: No. \t strand(+/-) \t leftEnd \t rightEnd \t contig
        geneCaller = 'unknown'
        geneCallFile = '' 
        geneNo = ''
        strand = ''
        leftEnd = ''
        rightEnd = ''
        length = '' 
        contig = ''
        geneName = ''
        sequence = ''
        contig = ''  # header of contig on which genes occur, if not listed in gene call file

        # Gather information about these gene calls
        if isinstance(geneCallInfo,dict):
            p_digits = re.compile('\d+')
            if 'contig' in geneCallInfo:
                contig = geneCallInfo['contig']
            if 'primaryCalls' in geneCallInfo:
                gc = geneCallInfo['primaryCalls'].lower()
                match_phanotate  = re.search('phanotate',gc)
                match_PHANOTATE  = re.search('PHANOTATE',gc)
                match_rast       = re.search('rast',gc)
                match_glimmer2   = re.search('glimmer2',gc)
                match_glimmer3   = re.search('glimmer3',gc)
                match_prodigal   = re.search('prodigal',gc)
                match_geneMarkS  = re.search('genemarks',gc)
                match_custom     = re.search('custom',gc)
                match_consensus  = re.search('consensus',gc)
                match_superset   = re.search('superset',gc)
                match_commoncore = re.search('commoncore',gc)
                if match_phanotate or match_PHANOTATE:
                    geneCaller = 'phanotate'
                elif match_rast:
                    geneCaller = 'rast'
                elif match_glimmer2:
                    geneCaller = 'glimmer2' 
                elif match_glimmer3:
                    geneCaller = 'glimmer3' 
                elif match_prodigal:
                    geneCaller = 'prodigal'
                elif match_geneMarkS:
                    geneCaller = 'genemark'
                elif match_custom:
                    geneCaller = 'custom'
                elif match_consensus:
                    geneCaller = 'consensus'
                elif match_superset:
                    geneCaller = 'superset'
                elif match_commoncore:
                    geneCaller = 'commoncore'
            if 'primaryCallsPathFile' in geneCallInfo:
                geneCallFile = geneCallInfo['primaryCallsPathFile']
            else:
                if PHATE_WARNINGS == 'True':
                    print("cgp_genomeSequence says, WARNING: processGeneCalls(), no geneCall file provided")
                return (0)

        # Read gene-call lines from gene caller output file and create a new gene object
        fLines = GENE_CALL_H.read().splitlines()
        for fLine in fLines:
            match_geneCall = re.search('^\d+',fLine)
            if match_geneCall:
                # Format output by CGCparser.py is 6 columns, tabbed; final column is protein, but ignore
                columns  = fLine.split('\t')
                geneNo   = columns[0]
                strand   = columns[1]
                length   = columns[4]
                contig   = columns[5]  #*** NOTE: contig name should be in PHANOTATE gene call file, but currently not
                geneName = contig + '_' + geneCaller + '_' + geneNo + '_geneCall'

                # Remove '<' and '>' symbols in leftEnd, rightEnd, if present
                matchDigits = re.search(p_digits,columns[2])
                if matchDigits:
                    leftEnd = matchDigits.group(0)
                matchDigits = re.search(p_digits,columns[3])
                if matchDigits:
                    rightEnd = matchDigits.group(0)

                # Create new gene object and fill data
                newGene = copy.deepcopy(fastaObj)
                newGene.moleculeType = self.geneSet.moleculeType # Gene inherits from multi-fasta "parent"
                newGene.assignHeader(geneName)  # Note:  using method assignHeader establishes all header variants
                newGene.parentSequence = contig
                newGene.parentName     = contig
                newGene.number         = int(geneNo)  # as numbered in gene-call file #*** OOPS! might be a string
                newGene.order          = len(self.geneSet.fastaList) + 1 # order in which object is added
                newGene.sequenceType   = 'nt'
                newGene.geneCallFile   = geneCallFile
                newGene.geneCaller     = geneCaller
                newGene.geneCallRank   = 0  # rank matters only if comparing gene calls (not yet in service)
                # Set strand orientation (f=forward) and adjust start/end accordingly
                if strand == 'f' or strand == '+':
                    newGene.strand = '+'
                    newGene.start  = int(leftEnd)
                    newGene.end    = int(rightEnd)
                elif strand == 'r' or strand == '-':
                    newGene.strand = '-'
                    newGene.start  = int(leftEnd)
                    newGene.end    = int(rightEnd)
                else:
                    newGene.strand = 'x'
                    if PHATE_WARNINGS == 'True':
                        print("cgp_genomeSequence says, ERROR: anomalous strand setting in processGeneCalls, phate_genomeSequence module:", newGene.strand)

                # Extract gene from genome sequence
                sequence = self.getCGCsubsequence(newGene.start,newGene.end,newGene.strand,contig)
                newGene.sequence = sequence

                # Reverse complement the string if on reverse strand
                if newGene.strand == '-':
                    tempGene = Seq(newGene.sequence)
                    newGene.sequence = tempGene.reverse_complement()

                # Create new protein object and fill data
                newProtein = copy.deepcopy(fastaObj)
                newProtein.moleculeType = self.proteinSet.moleculeType   # "Inherits" moleculeType from its "parent"
                newProtein.assignHeader(geneName)   # Note: using method assignHeader establishes all header variants
                newProtein.parentSequence = newGene.parentSequence # for reporting, need to capture contig name
                newProtein.order = int(geneNo)
                newProtein.sequenceType = 'aa'
                newProtein.geneCallFile = geneCallFile
                newProtein.geneCaller = geneCaller
                newProtein.geneCallRank = 0
                newProtein.start = 1 
                newProtein.end = len(newGene.sequence)
                newProtein.parentName  = newGene.header
                newProtein.parentStart = newGene.start
                newProtein.parentEnd   = newGene.end
                newProtein.strand      = newGene.strand
                newProtein.sequence = self.translateGene(newGene)
                newProtein.annotationList = newGene.annotationList #*** Has annotation been assigned to new gene?

                # Record new gene, protein
                self.addGene(newGene)     # Add this new gene to the gene list
                self.addProtein(newProtein)

    def addContig(self,newContig):  # Adding a fasta object
        newContig.order = len(self.contigSet) 
        newContig.moleculeType = self.contigSet.moleculeType 
        self.contigSet.addFasta(newContig)

    def addContigs(self,lines):  # Reading in lines from a file
        if isinstance(lines,list):
            self.contigSet.addFastas(lines,"nt")

    def addGene(self,newGene):
        newGene.moleculeType = self.geneSet.moleculeType
        self.geneSet.addFasta(newGene)

    def addProtein(self,newProtein):
        newProtein.moleculeType = self.proteinSet.moleculeType
        self.proteinSet.addFasta(newProtein)

    # ANNOTATION

    def addAnnotation(self,newAnnot):
        self.annotationList.append(newAnnot)

    def countAllAnnotations(self):
        annotationCount = 0

        # Count total annotations at all levels
        annotationCount += len(self.annotationList)
        annotationCount += len(self.contigSet.annotationList)
        for fasta in self.contigSet.fastaList:
            annotationCount += len(fasta.annotationList)
        annotationCount += len(self.geneSet.annotationList)
        for fasta in self.geneSet.fastaList:
            annotationCount += len(fasta.annotationList)
        annotationCount += len(self.proteinSet.annotationList)
        for fasta in self.proteinSet.fastaList:
            annotationCount += len(fasta.annotationList)
        
        return annotationCount

    # PRINT METHODS

    def printGenomeData(self):
        contigCount     = len(self.contigSet.fastaList)
        geneCount       = len(self.geneSet.fastaList)
        proteinCount    = len(self.proteinSet.fastaList)
        annotationCount = len(self.annotationList)
        print(self.genomeType, "genome", self.name, self.genomeName, self.species)
        print("Number of contigs:", contigCount)
        print("Names of contigs:")
        for fa in self.contigSet.fastaList:
            print(fa.header)
        print("  gene calls:", geneCount)
        print("  proteins:", proteinCount)
        print("  annotations:", annotationCount)
        if annotationCount > 0:
            print("Annotations:")
            for annot in self.annotationList:
                print("  ", annot)

    def printGenomeData_tab(self):
        FIRST = True
        contigCount     = str(len(self.contigSet.fastaList))
        geneCount       = str(len(self.geneSet.fastaList))
        proteinCount    = str(len(self.proteinSet.fastaList))
        annotationCount = str(len(self.annotationList))

        # Print summary for this genome
        annotationCount = self.countAllAnnotations()
        #tabOut = 'Genome:' + self.filename + '\tName:' + self.name + '\tType:' + self.genomeType
        tabOut = 'Genome:' + self.genomeName + '\tName:' + self.name + '\tType:' + self.genomeType
        print(tabOut)
        tabOut = 'Contigs:' + contigCount + '\tGenes:' + geneCount + '\tProteins:' + proteinCount + '\tAnnotations:' + str(annotationCount)
        print(tabOut)

        # Print annotations for the genome sequence
        if self.annotationList:
            for annot in self.annotationList:
                if FIRST:
                    annot.printAnnotationRecord_tabHeader()
                    FIRST = False
                annot.printAnnotationRecord_tab()

        # Print annotations for each contig
        if self.contigSet:
            for annotation in self.contigSet.annotationList:
                annotation.printAnnotationRecord()
            for contig in self.contigSet.fastaList:
                contig.printAll_tab()

        # Print annotations for each gene #*** Blasting of gene sequences not yet in service
        if self.geneSet:
            for annotation in self.geneSet.annotationList:
                annotation.printAnnotationRecord()
            for gene in self.geneSet.fastaList:
                gene.printAll_tab()

        # Print annotations for the proteins in the protein set
        if self.proteinSet:
            for annotation in self.proteinSet.annotationList:
                annotation.printAnnotationRecord()
            for protein in self.proteinSet.fastaList:
                protein.printAll_tab()

    def printGenomeData2file_tab(self,FILE_HANDLE):
        FIRST = True
        contigCount     = str(len(self.contigSet.fastaList))
        geneCount       = str(len(self.geneSet.fastaList))
        proteinCount    = str(len(self.proteinSet.fastaList))
        annotationCount = str(len(self.annotationList))  #*** CHECK: why is this being calc'd 2 ways???

        # Print summary for this genome
        annotationCount = self.countAllAnnotations()
        #tabOut = 'Genome:' + self.filename + '\tName:' + self.name + '\tType:' + self.genomeType
        tabOut = 'Genome:' + self.genomeName + '\tName:' + self.name + '\tType:' + self.genomeType
        FILE_HANDLE.write("%s\n" % (tabOut))
        tabOut = 'Contigs:' + contigCount + '\tGenes:' + geneCount + '\tProteins:' + proteinCount + '\tAnnotations:' + str(annotationCount)
        FILE_HANDLE.write("%s\n" % (tabOut))

        # Print annotations for the genome sequence
        if self.annotationList:
            for annot in self.annotationList:
                if FIRST:
                    annot.printAnnotationRecord2file_tabHeader(FILE_HANDLE)
                    FIRST = False
                annot.printAnnotationRecord2file_tab(FILE_HANDLE)

        # Print annotations for each contig
        if self.contigSet:
            for annotation in self.contigSet.annotationList:
                annotation.printAnnotationRecord2file_tab(FILE_HANDLE)
            for contig in self.contigSet.fastaList:
                contig.printAll2file_tab(FILE_HANDLE)

        # Print annotations for each gene #*** Blasting of gene sequences not yet in service
        if self.geneSet:
            for annotation in self.geneSet.annotationList:
                annotation.printAnnotationRecord2file_tab(FILE_HANDLE)
            for gene in self.geneSet.fastaList:
                gene.printAll2file_tab(FILE_HANDLE)

        # Print annotations for the proteins in the protein set
        if self.proteinSet:
            for annotation in self.proteinSet.annotationList:
                annotation.printAnnotationRecord2file_tab(FILE_HANDLE)
            for protein in self.proteinSet.fastaList:
                protein.printAll2file_tab(FILE_HANDLE)

    def printGenomeData2file(self,FILE_HANDLE):  # For reporting / debugging
        contigCount = len(self.contigSet.fastaList)
        geneCount = len(self.geneSet.fastaList)
        proteinCount = len(self.proteinSet.fastaList)
        annotationCount = len(self.annotationList)
        FILE_HANDLE.write("%s%s%s%s%s" % (self.genomeType,"genome",self.name,self.species,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("Number of contigs:",contigCount,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("  gene calls:",geneCount,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("  proteins:",proteinCount,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("  annotations:",annotationCount,"\n"))
        if annotationCount > 0:
            FILE_HANDLE.write("%s" % ("Annotations:\n"))
            for annot in self.annotationList:
                FILE_HANDLE.write("%s%s%s" % ("  ",annot,"\n"))

    def printGenomeData2file_GFF(self,FILE_HANDLE):
        FILE_HANDLE.write("%s\n" % (GFF_COMMENT))
        if PHATE_PROGRESS:
            print("cgp_genomeSequence says, There are", len(self.geneSet.fastaList), "genes, and", len(self.proteinSet.fastaList), "proteins")
            print("cgp_genomeSequence says, Writing data to GFF file")
        #*** NOTE: The following code assumes that the list of proteins corresponds precisely to
        #*** the list of genes.  Oh, for the want of a pointer!!!
        for i in range(0, len(self.geneSet.fastaList)-1):
            self.geneSet.fastaList[i].printData2file_GFF(FILE_HANDLE,'gene')
            self.proteinSet.fastaList[i].printData2file_GFF(FILE_HANDLE,'CDS')

    def printAll(self):
        print("Contigs =====================================================")
        self.contigSet.printAll()
        print("Genes =======================================================")
        self.geneSet.printAll()
        print("Proteins ====================================================")
        self.proteinSet.printAll()

    def printAll2file(self,FILE_HANDLE):  # For reporting / debugging
        FILE_HANDLE.write("%s" % ("Contigs =============================================\n"))
        self.contigSet.printAll2file(FILE_HANDLE)
        FILE_HANDLE.write("%s" % ("Genes =============================================\n"))
        self.geneSet.printAll2file(FILE_HANDLE)
        FILE_HANDLE.write("%s" % ("Proteins =============================================\n"))
        self.proteinSet.printAll2file(FILE_HANDLE)

    def printGenes(self):
        print(self.filename)
        print(self.name)
        self.geneSet.printAll()

    def printFastas2file(self,kvargs): # Prints to the file holding fastas (not a report/debug file) 
        mtype = ""
        # Get arguments that were provided
        if "mtype" in list(kvargs.keys()):
            mtype = kvargs["mtype"]
        else:
            mtype = "gene"
        if "headerType" in list(kvargs.keys()):
            hdrType = kvargs["headerType"]
        else:
            hdrType = "short"
        if "filename" in list(kvargs.keys()):
            filename = kvargs["filename"]
            OPEN_FILE = open(filename,"w")
        else:
            return False
        # Get headers and sequences, depending on which set is requested
        if mtype.lower() == "gene":
            for fa in self.geneSet.fastaList:
                if hdrType.lower() == "short":
                    hdr = fa.getShortHeader()
                elif hdrType.lower() == "compound":
                    hdr = fa.getCompoundHeader()
                elif hdrType.lower() == "truncated":
                    hdr = fa.getTruncHeader()
                elif hdrType.lower() == "full":
                    hdr = fa.getHeader("full")
                else:
                    hdr = fa.getHeader("full")
                seq = fa.sequence
                if hdr and seq:
                    OPEN_FILE.write("%s\n%s\n" % (hdr,seq))
            OPEN_FILE.close()
            return True
        elif mtype.lower() == "protein":
            for fa in self.proteinSet.fastaList:
                if hdrType.lower() == "short":
                    hdr = fa.getShortHeader()
                else:
                    hdr = fa.getHeader("full")
                seq = fa.sequence
                if hdr and seq:
                    OPEN_FILE.write("%s\n%s\n" % (hdr,seq))
            OPEN_FILE.close()
            return True
        elif mtype.lower() == "contig":
            for fa in self.contigSet.fastaList:
                if hdrType.lower() == "short":
                    hdr = fa.getShortHeader()
                else:
                    hdr = fa.getHeader("full")
                seq = fa.sequence
                if hdr and seq:
                    OPEN_FILE.write("%s\n%s\n" % (hdr,seq))
            OPEN_FILE.close()
            return True
        return False   # Oops, insufficient or incorrect data provided

    # OTHER METHODS

    def translateGene(self,geneObject):  #*** This is ugly, but transeq works on files, not strings
        proteinString = ""
        geneticCode = 11
        tempGeneFile = self.outputDir + "tempGeneFile"
        tempProtFile = self.outputDir + "tempProtFile"
        GENE_H = open(tempGeneFile,"w")
        outHeader = '>' + geneObject.header
        outSequence = geneObject.sequence
        GENE_H.write("%s\n" % (outHeader)) 
        GENE_H.write("%s\n" % (outSequence)) 
        GENE_H.close()
        command = "transeq" + " -sequence " + tempGeneFile + " -outseq " + tempProtFile + " -table " + str(geneticCode) + " -frame=1"
        # OS matters; change to other system call if you get an error message on this line
        result = os.system(command)
        #result = subprocess.check_output(command,shell=True)
        PROT_H = open(tempProtFile,"r")
        fLines = PROT_H.read().splitlines()
        for i in range(1,(len(fLines))):
            proteinString += fLines[i] 
        PROT_H.close()
        return proteinString

    def translateGenes(self,kvargs):  # Clear proteinSet; translate geneSet fastas and load to proteinSet
        # Note: Recent version of EMBOSS/transeq is messing even more with headers, changing '/' to '_'
        # ...so, added code to reinstate the correct header strings (25 March 2020).
        #*** MOVE THIS CODE TO cleanupAfterEmboss
        if "geneticCode" in list(kvargs.keys()):
            geneticCode = kvargs["geneticCode"]
        else:
            geneticCode = 11 # default is bacteria
        if "geneFile" in list(kvargs.keys()):    #*** should check that files exist
            geneFile = kvargs["geneFile"]
        if "proteinFile" in list(kvargs.keys()):
            protFile = kvargs["proteinFile"]
            tempProtFile = protFile + '.temp'
        command = EMBOSS_PHATE_HOME + "transeq" + " -sequence " + geneFile + " -outseq " + tempProtFile + " -table " + str(geneticCode)
        result = os.system(command)
        temp_h = open(tempProtFile,"r")
        prot_h = open(protFile,"w")
        lines = temp_h.read().splitlines()
        for line in lines:
            match_header = re.search("^>",line)
            if match_header:
                tempHeader = re.sub('_','/',line) 
                newHeader  = re.sub('/1$','_1',tempHeader)
                prot_h.write("%s\n" % (newHeader))
            else:
                prot_h.write("%s\n" % (line))
        prot_h.close()
        temp_h.close()
        return result 
    
    # MOVING THIS METHOD TO CLASS BLAST    
    def x_makeBlastDB(self,kvargs): # Create blast DBs for contigs, genes, proteins
        if "dbType" in list(kvargs.keys()):
            databaseType = kvargs["dbType"].lower()
        else:
            databaseType = "nucl"
        if "filename" in list(kvargs.keys()):
            filename = kvargs["filename"]
        else:
            return False
        if databaseType == "nucl" or databaseType == "nt" or databaseType == "dna":
            command = "makeblastdb -in " + filename + " -dbtype nucl -logfile " + filename + ".blastdb1_nucl_log"
        elif databaseType == "prot" or databaseType == "aa" or databaseType == "protein":
            command = "makeblastdb -in " + filename + " -dbtype prot -logfile " + filename + ".blastdb1_prot_log"
        else:
            return False

        # Run makeblastdb when semiphore lock is available 
        MAKEBLASTDB_DONE = False
        count = 0
        time.sleep(random.randint(1,RANDOM_LIMIT))
        if os.environ["CGP_LOCK"] == '1':
            os.environ["CGP_LOCK"] = '0'
            myResult = os.system(command)
            os.environ["CGP_LOCK"] = '1'
            MAKEBLASTDB_DONE = True
        else:
            while os.environ["CGP_LOCK"] == '0':
                time.sleep(random.randint(1,RANDOM_LIMIT))
                count += 1
                if os.environ["CGP_LOCK"] == '1':
                    os.environ["CGP_LOCK"] = '0'
                    myResult = os.system(command)
                    os.environ["CGP_LOCK"] = '1'
                    MAKEBLASTDB_DONE = True
                if count > 100:
                    break
        if not MAKEBLASTDB_DONE:
            print("cgp_genomeSequence says, ERROR: semiphore-controlled makeblastdb was not executed!!!")
            exit(0)

        return myResult

    def write2proteinSet(self,faLines):  # input list of lines containing protein fasta sequences
        self.proteinSet.addFastas(faLines,"aa")
        if len(self.proteinSet.fastaList) > 0:
            return True
        else:
            return False 

    def restorSlashesAfterEMBOSS():
        temp_h = open(tempProtFile,"r")
        prot_h = open(protFile,"w")
        lines = temp_h.read().splitlines()
        for line in lines:
            match_header = re.search("^>",line)
            if match_header:
                tempHeader = re.sub('_','/',line) 
                newHeader  = re.sub('/1$','_1',tempHeader)
                prot_h.write("%s\n" % (newHeader))
            else:
                prot_h.write("%s\n" % (line))
        prot_h.close()

    def cleanUpAfterEMBOSS(self):  # removes the pesky characters that EMBOSS adds to fasta sequence
        self.proteinSet.removeEMBOSSpostfix()
        self.proteinSet.removeTerminalAsterisk()
        return 


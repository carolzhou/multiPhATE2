###############################################################
#
# Name: phate_genomeSequence.py
#
# Programmer: Carol L. Ecale Zhou
#
# Most recent update: 20 January 2021 
#
# Description:
# Module comprising data structures for organizing genome information
#
# Programmer's Note:  EMBOSS messes with fasta headers; therefore, I am putting minimal 
#    info in the header and using only '/'
#    as a separator. EMBOSS adds "_1" to the end of a simple fasta 
#    header--can't prevent this. Using other 
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
#        processTrnaGeneCalls
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
#        cleanUpAfterEmboss
#
#######################################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE GPL3 LICENSE. SEE INCLUDED FILE GPL-3.PDF FOR DETAILS.

import string
import phate_fastaSequence
import phate_annotation
import re, os, copy
import subprocess
 
BLAST_HOME          = os.environ["PHATE_BLAST_HOME"] 
EMBOSS_PHATE_HOME   = os.environ["PHATE_EMBOSS_PHATE_HOME"] 
CODE_BASE_DIR       = ""
OUTPUT_DIR          = ""

# Verbosity

CLEAN_RAW_DATA   = os.environ["PHATE_CLEAN_RAW_DATA"]
PHATE_WARNINGS_STRING = os.environ["PHATE_PHATE_WARNINGS"]
PHATE_MESSAGES_STRING = os.environ["PHATE_PHATE_MESSAGES"]
PHATE_PROGRESS_STRING = os.environ["PHATE_PHATE_PROGRESS"]
PHATE_WARNINGS = False
PHATE_MESSAGES = False
PHATE_PROGRESS = False
if PHATE_WARNINGS_STRING.lower() == 'true':
    PHATE_WARNINGS = True
if PHATE_MESSAGES_STRING.lower() == 'true':
    PHATE_MESSAGES = True
if PHATE_PROGRESS_STRING.lower() == 'true':
    PHATE_PROGRESS = True
# Override
#PHATE_WARNINGS = True
#PHATE_MESSAGES = True
#PHATE_PROGRESS = True

DEBUG            = False
#DEBUG           = True 

# For GFF output
GFF_COMMENT = "##gff-version 3"

# Templates
annotationObj = phate_annotation.annotationRecord()
fastaObj      = phate_fastaSequence.fasta()
trnaObj       = phate_fastaSequence.fasta()

# Reverse complement
MAKETRANS_OLD = False   # Differences in Python2 vs. Python3  (Python3.6 vs. Python 3.7???)
if MAKETRANS_OLD:
    complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')  # Python 2
else:
    complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')     # Python 3


class genome(object):

    # Class genome uses fasta objects to represent a genome/chromosome or a plasmid and its genes and proteins.
    # The chromosome's contig set, gene set, and protein set are each represented by multi-fasta objects.

    def __init__(self):
        self.filename                = ""          # name of file that contains the genome sequence
        self.genomeType              = "unknown" # typically, "chromosome", "plasmid", "virus", "segmented", "satellite"
        self.name                    = "unknown" # common name (e.g., Drosophila)
        self.genomeName              = "unkknown" # more specific name (e.g., 'Yp phage YP123')
        self.species                 = "unknown" # Latin name (e.g., Drosophila melanogaster)
        self.annotationList          = []    # using phate_annotation.py classes
        self.contigSet               = phate_fastaSequence.multiFasta()
        self.geneSet                 = phate_fastaSequence.multiFasta()
        self.proteinSet              = phate_fastaSequence.multiFasta()
        self.trnaSet                 = phate_fastaSequence.multiFasta()
        self.contigSet.moleculeType  = 'contig'
        self.contigSet.sequenceType  = 'nt'
        self.geneSet.moleculeType    = 'gene'
        self.geneSet.sequenceType    = 'nt'
        self.proteinSet.moleculeType = 'peptide'
        self.proteinSet.sequenceType = 'aa'
        self.trnaSet.moleculeType    = 'trna'
        self.trnaSet.sequenceType    = 'nt'
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

    def getCGCsubsequence(self,start,end,strand,in_contig):  # Tailored for output from CGCparser.py
        subSeq = ""
        # For safety, strip white space(s) from incoming contig
        cleanContig = in_contig.lstrip()
        contig = cleanContig.rstrip()
        for fa in self.contigSet.fastaList:
            if fa.header == contig:
                subSeq = fa.getSubsequence(int(start)-1,int(end)) #*** ???
            else:
                if PHATE_WARNINGS:
                    print("phate_genomeSequence says, WARNING: fa.header -", fa.header, "- did not match contig -", contig, "-")
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

    def processTrnaGeneCalls(self,kvargs):
        if PHATE_MESSAGES: 
            print("phate_genomeSequence says, Running processTrnaGeneCalls.")
        TRNA_GENES_FILE = ''; contig = ''; trnaNo = ''; begin = ''; end = ''; codon = ''; aCodon = ''; iBegin = ''; iEnd = ''; score = ''; annot = ''
        ANNOTATION_FOUND = False; codonAcodon = ""
        if isinstance(kvargs,dict):
            if "trnaGenesFile" in kvargs.keys():
                TRNA_GENES_FILE = kvargs["trnaGenesFile"]
        if os.path.exists(TRNA_GENES_FILE):
            TRNA_H = open(TRNA_GENES_FILE,'r')
            fLines = TRNA_H.read().splitlines()
            count = -3 
            for fLine in fLines:
                count += 1   
                ANNOTATION_FOUND = False; codonAcodon = ""
                if count > 0:  # skip 1st 3 lines of file
                    fields = fLine.split('\t')
                    if len(fields) >= 9:
                        contig = fields[0].rstrip()  # Name
                        trnaNo = fields[1]  # tRNA #
                        begin  = fields[2]  # tRNA Begin 
                        end    = fields[3]  # Bounds End 
                        codon  = fields[4]  # tRNA Type 
                        aCodon = fields[5]  # Anti Codon 
                        iBegin = fields[6]  # Intron Begin 
                        iEnd   = fields[7]  # Bounds End 
                        score  = fields[8]  # Inf Score 
                        codonAcodon = "tRNA-" + codon + '(' + aCodon + ')'
                    if len(fields) == 10:
                        annot  = fields[9]  # Note field
                        if re.search('[\w\d]',annot):  # check whether annotation field actually contains meaningful text
                            ANNOTATION_FOUND = True
                    trnaName = contig + '_trna_' + str(count)
                    newTrna = copy.deepcopy(fastaObj)
                    newTrna.moleculeType   = self.trnaSet.moleculeType
                    newTrna.assignHeader(trnaName)
                    newTrna.parentSequence = contig
                    newTrna.parentName     = contig
                    newTrna.number         = count
                    newTrna.order          = count
                    newTrna.sequenceType   = 'nt'
                    newTrna.geneCallFile   = TRNA_GENES_FILE 
                    newTrna.geneCaller     = 'trnascan-se'
                    newTrna.geneCallRank   = 0
                    newTrna.score          = score
                    if int(begin) < int(end):
                        newTrna.strand     = '+'
                        newTrna.start      = begin
                        newTrna.end        = end
                    elif int(begin) > int(end):
                        newTrna.strand     = '-'
                        newTrna.start      = end
                        newTrna.end        = begin
                    else:
                        newTrna.strand     = 'x'
                        newTrna.start      = begin
                        newTrna.end        = end
                        if PHATE_WARNINGS:
                            print("phate_genomeSequence says, WARNING/ERROR: anomalous strand setting in processGeneCalls, phate_genomeSequence module:", newTrna.strand, newTrna.contig, newTrna.number)
                    # Extract trna from genome sequence
                    sequence = self.getCGCsubsequence(newTrna.start,newTrna.end,newTrna.strand,contig)
                    newTrna.sequence = sequence

                    # Reverse complement the string if on reverse strand
                    if newTrna.strand == '-':
                        reverseComplement = newTrna.sequence.translate(complements)[::-1]
                        newTrna.sequence = reverseComplement

                    # If an annotation was detected, create an annotation object, populate it, add to newTrna.
                    newAnnot = copy.deepcopy(annotationObj) 
                    newAnnot.source           = 'phate'
                    newAnnot.method           = 'trnascan'
                    newAnnot.annotationType   = 'trna'
                    newAnnot.description      = codonAcodon 
                    newAnnot.annotationString = codonAcodon
                    if ANNOTATION_FOUND:
                        newAnnot.description      += '-' + annot
                        newAnnot.annotationString += '-' + annot
                    newTrna.annotationList.append(newAnnot)

                    # Record new trna
                    if PHATE_MESSAGES:
                        print("phate_genomeSequence says, Recording new trna: ",self.trnaSet)
                    self.addTrna(newTrna)   # Add this new trna to the trna list
            TRNA_H.close()
        return

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
            #if 'geneCaller' in geneCallInfo:
            if 'primaryCalls' in geneCallInfo:
                #gc = geneCallInfo['geneCaller'].lower()
                gc = geneCallInfo['primaryCalls'].lower()
                match_phanotate = re.search('phanotate',gc)
                match_PHANOTATE = re.search('PHANOTATE',gc)
                match_Phanotate = re.search('Phanotate',gc)
                match_rast      = re.search('rast',gc)
                match_glimmer2  = re.search('glimmer2',gc)
                match_glimmer3  = re.search('glimmer3',gc)
                match_prodigal  = re.search('prodigal',gc)
                match_geneMarkS = re.search('genemarks',gc)
                match_consensus = re.search('consensus',gc)
                match_superset  = re.search('superset',gc)
                match_custom    = re.search('custom',gc)
                if match_phanotate or match_PHANOTATE or match_Phanotate:
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
                elif match_consensus:
                    geneCaller = 'consensus'
                elif match_superset:
                    geneCaller = 'superset'
                elif match_custom:
                    geneCaller = 'custom'
            if 'primaryCallsPathFile' in geneCallInfo:
                geneCallFile = geneCallInfo['primaryCallsPathFile']
            else:
                if PHATE_WARNINGS:
                    print("phate_genomeSequence says, WARNING: processGeneCalls(), no geneCall file provided")
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
                contig   = columns[5].rstrip()  #*** NOTE: contig name should be in PHANOTATE gene call file, but currently not
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
                    if PHATE_WARNINGS:
                        print("phate_genomeSequence says, WARNING/ERROR: anomalous strand setting in processGeneCalls, phate_genomeSequence module:", newGene.strand)

                #*** BANDAID - to compensate for PHANOTATE sometimes starting gene at 0
                if newGene.start == 0:
                    newGene.start = 1

                # Extract gene from genome sequence
                sequence = self.getCGCsubsequence(newGene.start,newGene.end,newGene.strand,contig)
                newGene.sequence = sequence

                # Reverse complement the string if on reverse strand
                if newGene.strand == '-':
                    reverseComplement = newGene.sequence.translate(complements)[::-1] 
                    newGene.sequence = reverseComplement 

                # Create new protein object and fill data
                newProtein = copy.deepcopy(fastaObj)
                newProtein.moleculeType = self.proteinSet.moleculeType   # "Inherits" moleculeType from its "parent"
                newProtein.assignHeader(geneName)   # Note: using method assignHeader establishes all header variants
                newProtein.parentSequence = newGene.sequence
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

    def addTrna(self,newTrna):
        newTrna.moleculeType = self.trnaSet.moleculeType
        self.trnaSet.addFasta(newTrna)

    # ANNOTATION

    def addAnnotation(self,newAnnot):
        self.annotationList.append(newAnnot)

    def countAllAnnotations(self):
        annotationCount = 0

        # Count total annotations at all levels
        annotationCount += len(self.annotationList)            # genomic level
        annotationCount += len(self.contigSet.annotationList)  # contig set level
        for fasta in self.contigSet.fastaList:                 # contig level
            annotationCount += len(fasta.annotationList)       
        annotationCount += len(self.geneSet.annotationList)    # gene set level
        for fasta in self.geneSet.fastaList:                   # gene level
            annotationCount += len(fasta.annotationList)
        annotationCount += len(self.proteinSet.annotationList) # protein set level
        for fasta in self.proteinSet.fastaList:                # protein gene level
            annotationCount += len(fasta.annotationList)
        annotationCount += len(self.trnaSet.annotationList)    # trna set level
        for fasta in self.trnaSet.fastaList:                   # trna gene level
            annotationCount += len(fasta.annotationList)
        
        return annotationCount

    # PRINT METHODS

    def printGenomeData(self):
        contigCount     = len(self.contigSet.fastaList)
        geneCount       = len(self.geneSet.fastaList)
        proteinCount    = len(self.proteinSet.fastaList)
        annotationCount = len(self.annotationList)
        trnaCount       = len(self.trnaList)
        print(self.genomeType, "genome", self.name, self.genomeName, self.species)
        print("Number of contigs:", contigCount)
        print("Names of contigs:")
        for fa in self.contigSet.fastaList:
            print(fa.header)
        print("  gene calls: ", geneCount)
        print("  proteins:   ", proteinCount)
        print("  trnas:      ", trnaCount)
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
        trnaCount       = str(len(self.trnaSet.fastaList))
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

        # Print annotations for the trnas in the protein set
        if self.trnaSet:
            for annotation in self.trnaSet.annotationList:
                annotation.printAnnotationRecord()
            for trna in self.trnaSet.fastaList:
                trna.printAll_tab()

    def printGenomeData2file_tab(self,FILE_HANDLE):
        FIRST = True
        contigCount     = str(len(self.contigSet.fastaList))
        geneCount       = str(len(self.geneSet.fastaList))
        proteinCount    = str(len(self.proteinSet.fastaList))
        trnaCount       = str(len(self.trnaSet.fastaList))
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

        # Print annotations for the trnas in the trna set
        if self.trnaSet:
            for annotation in self.trnaSet.annotationList:
                annotation.printAnnotationRecord2file_tab(FILE_HANDLE)
            for trna in self.trnaSet.fastaList:
                trna.printAll2file_tab(FILE_HANDLE)

    def printGenomeData2file(self,FILE_HANDLE):  # For reporting / debugging
        contigCount  = len(self.contigSet.fastaList)
        geneCount    = len(self.geneSet.fastaList)
        proteinCount = len(self.proteinSet.fastaList)
        trnaCount    = len(self.trnaSet.fastaList)
        annotationCount = len(self.annotationList)
        FILE_HANDLE.write("%s%s%s%s%s" % (self.genomeType,"genome",self.name,self.species,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("Number of contigs:",contigCount,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("  gene calls:",geneCount,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("  proteins:",proteinCount,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("  trnas:",trnaCount,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("  annotations:",annotationCount,"\n"))
        if annotationCount > 0:
            FILE_HANDLE.write("%s" % ("Annotations:\n"))
            for annot in self.annotationList:
                FILE_HANDLE.write("%s%s%s" % ("  ",annot,"\n"))

    def printGenomeData2file_GFF(self,FILE_HANDLE):
        FILE_HANDLE.write("%s\n" % (GFF_COMMENT))
        if PHATE_PROGRESS:
            print("phate_genomeSequence says, There are",len(self.geneSet.fastaList),"genes,",len(self.proteinSet.fastaList),"proteins, and",len(self.trnaSet.fastaList),"trna genes")
            print("phate_genomeSequence says, Writing data to GFF file.")
        geneCount = 0; trnaCount = 0
        for contig in self.contigSet.fastaList:
            contigName = contig.header  #*** WHICH HEADER TYPE corresponds to gene and protein records?
            for i in range(0,len(self.geneSet.fastaList)):
                if contigName == self.geneSet.fastaList[i].parentName: 
                    geneCount += 1
                    self.geneSet.fastaList[i].printData2file_GFF(FILE_HANDLE,'gene',contigName)
                    # protein is in same order in respective fastaList as corresponding gene
                    self.proteinSet.fastaList[i].printData2file_GFF(FILE_HANDLE,'CDS',contigName) 
            # Once genes and proteins have been printed, print the trnas from this contig; #*** could use gene[i+1]'s left end to know when to print trna 
            for trna in self.trnaSet.fastaList:
                if contigName == trna.parentName:
                    trnaCount += 1
                    trna.printData2file_GFF(FILE_HANDLE,'trna',contigName)
        if PHATE_MESSAGES:
            print("phate_genomeSequence says, Printing of GFF file is complete:",geneCount,"genes and",trnaCount,"trnas.")
        return

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
                elif hdrType.lower() == "cgp":
                    hdr = fa.getCGPheader()
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
                elif hdrType.lower() == "cgp":
                    hdr = fa.getCGPheader()
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
        if "geneticCode" in list(kvargs.keys()):
            geneticCode = kvargs["geneticCode"]
        else:
            geneticCode = 11 # default is bacteria
        if "geneFile" in list(kvargs.keys()):    #*** should check that files exist
            geneFile = kvargs["geneFile"]
        if "proteinFile" in list(kvargs.keys()):
            protFile = kvargs["proteinFile"]
        command = EMBOSS_PHATE_HOME + "transeq" + " -sequence " + geneFile + " -outseq " + protFile + " -table " + str(geneticCode)
        result = os.system(command)
        return result 
    
    # MOVING THIS METHOD TO CLASS BLAST    
    def makeBlastDB(self,kvargs): # Create blast DBs for contigs, genes, proteins
        if "dbType" in list(kvargs.keys()):
            databaseType = kvargs["dbType"].lower()
        else:
            databaseType = "nucl"
        if "filename" in list(kvargs.keys()):
            filename = kvargs["filename"]
        else:
            return False
        if databaseType == "nucl" or databaseType == "nt" or databaseType == "dna":
            command = BLAST_HOME + "makeblastdb -in " + filename + " -dbtype nucl -logfile " + filename + ".blastdb1_nucl_log"
        elif databaseType == "prot" or databaseType == "aa" or databaseType == "protein":
            command = BLAST_HOME + "makeblastdb -in " + filename + " -dbtype prot -logfile " + filename + ".blastdb1_prot_log"
        else:
            return False
        myResult = os.system(command)
        return myResult

    def write2proteinSet(self,faLines):  # input list of lines containing protein fasta sequences
        self.proteinSet.addFastas(faLines,"aa")
        if len(self.proteinSet.fastaList) > 0:
            return True
        else:
            return False 

    def cleanUpAfterEMBOSS(self):  # removes the pesky characters that EMBOSS adds to fasta sequence
        self.proteinSet.removeEMBOSSpostfix()
        self.proteinSet.removeTerminalAsterisk()
        return 


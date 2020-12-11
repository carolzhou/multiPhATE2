###########################################################
# Name: cgp_annotation.py
#
# Programmer: Carol L. Ecale Zhou
#
# Latest update: 20 October 2020
#
# Description: Module containing classes and methods for representing annotation results from various sources 
#
# Classes and methods: 
#     annotationRecord
#         addPVOBid2list
#         getPVOGassociationList
#         enterGFFdata(gff/dict)
#         getDBXREFs
#         findInfo
#         getFigDescription
#         getPvogMembers
#         getNCBItaxonomy
#         link2databaseIdentifiers
#         returnGFFannotationRecord
#         printAnnotationRecord
#         printAnnotationRecord_tabHeader
#         printAnnotationRecord_tab
#         printAnnotationRecord2file_tabHeader
#         printAnnotationRecord2file_tab
#         printAnnotationRecord2file
#         printAll
#         printAll2file(fileH)
#         writePVOGgroups
##########################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE GPL3 LICENSE. SEE INCLUDED FILE GPL-3.PDF FOR DETAILS.

import re, os, subprocess

#DEBUG = True
DEBUG = False

p_comment = re.compile('^#')

KEGG_VIRUS_BASE_DIR   = os.environ["PHATE_KEGG_VIRUS_BASE_DIR"]
NCBI_VIRUS_BASE_DIR   = os.environ["PHATE_NCBI_VIRUS_BASE_DIR"]
PHANTOME_BASE_DIR     = os.environ["PHATE_PHANTOME_BASE_DIR"]
NCBI_TAXON_DIR        = os.environ["PHATE_NCBI_TAXON_DIR"]
PVOGS_BASE_DIR        = os.environ["PHATE_PVOGS_BASE_DIR"]

# Verbosity
CLEAN_RAW_DATA_STRING = os.environ["PHATE_CLEAN_RAW_DATA"]
PHATE_PROGRESS_STRING = os.environ["PHATE_PHATE_PROGRESS"]
PHATE_MESSAGES_STRING = os.environ["PHATE_PHATE_MESSAGES"]
PHATE_WARNINGS_STRING = os.environ["PHATE_PHATE_WARNINGS"]

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

# External links
NCBI_TAXON_LINK = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id="

# List of useful files to convert vg identifier to annotation information.
# These are located on mpath under /home/zhou4/BioRemediation/Phage/Kegg_virus/
# Use: os.environ["KEGG_VIRUS_BLAST_HOME"]
vg_enzyme            = KEGG_VIRUS_BASE_DIR + "vg_enzyme.list"
vg_ko                = KEGG_VIRUS_BASE_DIR + "vg_ko.list"
vg_ncbi_proteinid    = KEGG_VIRUS_BASE_DIR + "vg_ncbi-proteinid.list"
vg_pfam              = KEGG_VIRUS_BASE_DIR + "vg_pfam.list"
vg_tax               = KEGG_VIRUS_BASE_DIR + "vg_tax.list"
vg_uniprot           = KEGG_VIRUS_BASE_DIR + "vg_uniprot.list"
phantome_phage_genes = PHANTOME_BASE_DIR   + "Phantome_Phage_genes.faa.headers"
ncbi_taxon_lookup    = NCBI_TAXON_DIR      + "nucl_gb.accession2taxid"   
pVOGheaderFile       = PVOGS_BASE_DIR      + "pVOGs_headers.lst"

class annotationRecord(object):

    def __init__(self):
        self.source            = "unknown" # Typically RAST, LLNL, PhAnToMe, GeneMark, Glimmer, Prodigal, PHANOTATE, KEGG, NCBI, pVOGs 
        self.method            = "unknown" # Typcially RAST, PFP, PhiRAST, JGI, SDSU, BLAST, blastp, blastn, HMM, jackhmmer 
        self.annotationType    = "unknown" # gene, mRNA, polypeptide, CDS, functional, homology, hmm
        self.pVOGlist          = []        # list of pVOG identifiers (identified via blast hit)
        self.contig            = "unknown"
        self.start             = 0
        self.end               = 0
        self.strand            = 'x' 
        self.readingFrame      = 'x'
        self.identifier        = "none"
        self.locusTag          = "none"
        self.name              = "none"  # subject hit header (i.e., database identifier provider in fasta header)
        self.description       = "none"  # more information: dbxref identifiers provided via lookup-tables (see above) 
        self.annotationList    = []      # could be multiple from single source 
        self.category          = "none"  # functional categories: primary, sequence, structure, motif, etc.
        self.wraparound        = "none"  # indicates that a gene call wraps around the genome sequence as given
        self.annotationSTring = ""      # used to construct a summary of annotation(s) for GFF output

    def addPVOGid2list(self,pVOG):
        self.pVOGlist.append(pVOG)

    def getPVOGassociationList(self):
        return(self.pVOGlist)

    def enterGFFdata(self,gff):  # Input a dict object with key/values as specified
        if isinstance(gff,dict):
            self.source          = gff["source"]
            self.method          = gff["method"]
            self.annotationType  = gff["type"]
            self.contig          = gff["contig"]
            self.start           = gff["start"]
            self.end             = gff["end"]
            self.strand          = gff["strand"]
            self.readingFrame    = gff["readingFrame"]
            annotList = gff["annotation"].split(';')
            for annot in annotList:
                self.annotationList.append(annot)
            self.category        = "sequence"
            return True
        else:
            return False

    # METHODS FOR ANNOTATING FROM EXTERNAL SOURCES

    def removeRedundancy(self,inList): # Eliminate redundancy in list
        outList = []
        for i in range(len(inList)):
            item = inList.pop()
            if item not in inList:
                outList.append(item)
        outList.reverse()
        return outList

    # METHODS FOR LINKING TO FUNCTIONAL ANNOTATIONS
    # The Phantome, Kegg-virus, and Ncbi-virus database provide cryptic fasta headers
    # Need to use look-up tables to get more descriptive names for these annotations

    def getDBXREFs(self,database):
        idList = []
        dbxref = ""
        p_truncatedSearchTerm = re.compile('^([^\s]*)\s')
        match_truncate = re.search(p_truncatedSearchTerm, self.name)
        if match_truncate:
            searchTerm = match_truncate.group(1)
        else:
            searchTerm = self.name
        command = 'grep ' + searchTerm + ' ' + database
        proc = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True) 
        (rawout, err) = proc.communicate()
        out = rawout.decode('utf-8')   # Python3
        if out != "":
            lines  = out.split('\n')
            for line in lines:
                if line != "":
                    fields = line.split('\t')
                    if len(fields) > 1:
                        dbxref = fields[1]
                    else:
                        if PHATE_WARNINGS:
                            print("cgp_annotation says, WARNING: no dbxref found for", self.name, "in database", database, "given line", line)
                    idList.append(dbxref)
        return idList

    def findInfo(self,searchTerm,database):  # Searches Phantome header file for annotation information
        p_truncatedSearchTerm = re.compile('^([^\s]*)\s')
        infoLines = []
        DATABASE_H = open(database,"r")
        dLines = DATABASE_H.read().splitlines()
        match_truncate = re.search(p_truncatedSearchTerm,searchTerm)
        if match_truncate:
            truncatedSearchTerm = match_truncate.group(1)
            searchTerm = truncatedSearchTerm
        for dLine in dLines:
            match_searchTerm = re.search(searchTerm,dLine)
            if match_searchTerm:
                fields = dLine.split('\s+')
                if searchTerm == fields[0][1:]:
                    infoLines.append(dLine) 
        DATABASE_H.close()
        return infoLines

    def getFigDescription(self,database):
        dbxrefList = []; infoList = []; infoString = ""
        # Alert: Really ugly regex follows; how to make this better?
        p_phantomeInfo = re.compile('\[[\w\s\d\-\.\(\)\;\:\_\,\#]+\]')  # Gather all text enclosed within []

        # Pull figfam descriptions from phantome header file 
        match_fig = re.search('fig\|',self.name)
        if match_fig:
            searchTerm = self.name 
            figLines = self.findInfo(searchTerm,database)
            for fig in figLines:
                infoList = re.findall(p_phantomeInfo,fig)
                for info in infoList:
                    infoString += ' | ' + info
                dbxrefList.append(infoString) 
        else:
            if PHATE_WARNINGS:
                print("cgp_annotaton says, WARNING: Unexpected name encountered in phate_annotation.py/getFigDescription:", self.name) 
        return dbxrefList 

    def getPvogMembers(self,database): # database should be the pVOGs headers file, but fasta file will work (slowly)
        dbxrefList = []; infoList = []; infoString = []
        p_pVOGid = re.compile('VOG\d+?a|b')
        # self.name is the pVOGs database blast hit header
        pVOGidList = []
        pVOGidList = re.findall(p_pVOGid,self.name)  # extract pVOGid(s) from header
        if pVOGidList:
            for pVOGid in pVOGidList:
                infoLines = re.findall(pVOGid,database) 
                for line in infoLines:
                    infoString += ' | ' + line
                dbxrefList.append(infoString)
        return dbxrefList  # list(s) of pVOGs database headers with common pVOGid

    # Query a taxonomy lookup table to get taxonomy information
    def getNCBItaxonomyID(self,database):   # Database maps ncbi header to taxonomy 
        ncbiTaxonList = []; giNumber = ''; accnNumber = ''
        ncbiTaxonList = []
        p_version = re.compile('(\w+_\d+)\.\d+')
        p_taxID   = re.compile('[\d\w\_]\s+[\d\w\_\.]+\s+([\d]+)\s+[\d]+')
        if self.name != '' and self.name != 'none':
            fields = self.name.split('|')  # Hit's fasta header has several fields, delimited w/'|'
            if len(fields) > 4:
                giNumber    = fields[1]
                accnNumber  = fields[3]
                description = fields[4]
                self.description = description
                searchTermString = accnNumber
                match_version = re.search(p_version,searchTermString)
                if match_version:
                    searchTerm = match_version.group(1) 
                else:
                    searchTerm = searchTermString
                command = 'grep \"' + searchTerm + '\" ' + database
                proc = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
                out, err = proc.communicate()
                if out != '':
                    match_taxID = re.search(p_taxID,out)
                    taxonomyID = match_taxID.group(1)
                    taxonomyString = 'NCBItaxID=' + taxonomyID
                    ncbiTaxonList.append(taxonomyString)
                    ncbiTaxonLink = NCBI_TAXON_LINK + taxonomyID
                    ncbiTaxonList.append(ncbiTaxonLink)
            else:
                if PHATE_WARNINGS:
                    print("cgp_annotation says, WARNING: NCBI hit header has improper format or is missing:", self.name)
        return ncbiTaxonList

    def link2databaseIdentifiers(self,database,dbName):
        dbxrefList = []  # holds concatenation of functional annotations
        enzymeList = []; ncbiProteinList = []; taxonList = [] # hold specific annotations
        pfamList   = []; uniprotList     = []; koList    = [] # hold specific annotations
        figList    = []
        annotationString = "" # string containing all dbxref annotations found
        annotation_item = ""

        if self.name == "" or self.name == "none":
            if PHATE_WARNINGS:
                print("cgp_annotation says, WARNING: No name for identification of dbxref in phate_annotation/link2databaseIdentifiers")
            return 
        else:
            if dbName.lower() == 'kegg':
                enzymeList      = self.getDBXREFs(vg_enzyme)
                koList          = self.getDBXREFs(vg_ko)
                ncbiProteinList = self.getDBXREFs(vg_ncbi_proteinid)
                pfamList        = self.getDBXREFs(vg_pfam)
                taxList         = self.getDBXREFs(vg_tax)
                uniprotList     = self.getDBXREFs(vg_uniprot)

                # Create annotation lists that can later be parsed
                for enzyme in enzymeList:
                    dbxrefList.append(enzyme) 
                for ko in koList:
                    dbxrefList.append(ko) 
                for protein in ncbiProteinList:
                    dbxrefList.append(protein) 
                for pfam in pfamList:
                    dbxrefList.append(pfam) 
                for taxon in taxList:
                    dbxrefList.append(taxon) 
                for uniprot in uniprotList:
                    dbxrefList.append(uniprot) 
                
            elif dbName.lower() == 'phantome':
                figList = self.getFigDescription(phantome_phage_genes) 
                for fig in figList:
                    dbxrefList.append(fig)

            elif dbName.lower() == 'ncbi':
                taxonList = self.getNCBItaxonomyID(ncbi_taxon_lookup)
                for taxon in taxonList:
                    dbxrefList.append(taxon) 

            elif dbName.lower() == 'pvogs':
                pVOGlist = self.getPvogMembers(pVOGheaderFile) 

            elif dbName.lower() == 'pvogs_hmm':
                pVOGlist = self.getPvogMembers(pVOGheaderFile) 

            elif dbName.lower() == 'ncbivirusprotein':
                pass
 
            else:
                if PHATE_WARNINGS:
                    print("cgp_annotation says, WARNING: Unrecognized database:", dbName) 

        for annotation_item in dbxrefList:
            nextAnnot = ' | ' + annotation_item
            annotationString += nextAnnot
        self.description = annotationString

        return 

    # Return annotations as a semicolon-delimited string
    def returnGFFannotationRecord(self,FILE_HANDLE):
        annot = ''
        if self.annotationType == 'gene':
            annot = '(gene) ' + self.start + '/' + self.end + '/' + self.strand + ' ' + self.method 
        elif self.annotationType == 'functional':
            annot = '(func) ' + self.method + ' ' + self.description
        elif self.annotationType == 'homology':
            homologName = self.name
            newName = re.sub(';','',homologName)  # Remove embedded ';' because GFF uses this delimiter
            annot = '(homolog) ' + self.method + ' ' + newName
        elif self.annotationType == 'hmm search':
            annot = '(hmm search) ' + self.method + ' ' + self.description
        elif self.annotationType == 'profile search':
            annot = '(profile search) ' + self.method + ' ' + self.description
        elif self.annotationType == 'cds':
            annot = '(cds)' + self.method + ' ' + self.description
        elif self.annotationType == 'mrna':
            annot = '(mrna)' + self.method + ' ' + self.description
        elif self.annotationType == 'polypeptide':
            annot = '(polypeptide) ' + self.method + ' ' + self.description
        else:
            annot = '(unk type) ' + self.method + ' ' + self.description

        annotationString = ""
        for item in self.annotationList:
            annotationString += item
        annot += annotationString
        FILE_HANDLE.write("%s" % (annot))

    # PRINT METHODS

    def printAnnotationRecord(self):
        print("Annotation source:", self.source, '| Method:', self.method, '| Type:', self.annotationType)
        print("Contig:", self.contig, "| Start:", self.start, "| End:", self.end, "| Strand:", self.strand)
        print("Name:", self.name, "Description:", self.description)
        print("Annotations:", self.annotationList)

    def printAnnotationRecord_tabHeader(self):
        header = 'Source\tMethod\tType\tCategory\tStart-End/strand\tName\tDescription'
        print(header)

    def printAnnotationRecord_tab(self):
        annotationString = ""
        #print "Number of annotations:", len(self.annotationList)
        for annot in self.annotationList:
            annotationString += annot
            annotationString += ' | '
        tabLine = self.source + '\t' + self.method + '\t' + self.annotationType + '\t' + self.category + '\t'
        tabLine += str(self.start) + '-' + str(self.end) + '/' + self.strand + '\t'
        tabLine += self.name + '\t' + self.description + '\t' + annotationString
        print(tabLine)

    def printAnnotationRecord2file_tabHeader(self,FILE_HANDLE):
        header = 'Source\tMethod\tType\tCategory\tStart-End/strand\tName\tDescription'
        FILE_HANDLE.write("%s\n" % (header))

    def printAnnotationRecord2file_tab(self,FILE_HANDLE):
        annotationString = ""
        for annot in self.annotationList:
            annotationString += annot
            annotationString += ' | '
        tabLine = self.source + '\t' + self.method + '\t' + self.annotationType + '\t' + self.category + '\t'
        tabLine += str(self.start) + '-' + str(self.end) + '/' + self.strand + '\t'
        tabLine += self.name + '\t' + self.description + '\t' + annotationString
        FILE_HANDLE.write("%s\n" % (tabLine))

    def printAnnotationRecord2file(self,FILE_HANDLE):  #*** Update this
        FILE_HANDLE.write("%s%s%s" % ("Annotation source:",self.source,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("Method:",self.method,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("Contig:",self.contig,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("Annotations:",self.annotationList,"\n"))

    def printAll(self):
        print("=== Annotation record ===")
        print("Source:", self.source) 
        print("Method:", self.method) 
        print("Type:", self.annotationType)
        print("Contig:", self.contig)
        print("Start:", self.start)
        print("End:", self.end)
        print("Strand:", self.strand) 
        print("Reading Frame:", self.readingFrame)
        print("Identifier:", self.identifier)
        print("Locus Tag:", self.locusTag)
        print("Name:", self.name)
        print("Description:", self.description)
        print("Category:", self.category)
        print("Wraparound:", self.wraparound)
        print("Annotation List:")
        for annot in self.annotationList:
            print("  ", annot)
        print("Category:", self.category)
        print("========================")

    def printAll2file(self,FILE_HANDLE):
        FILE_HANDLE.write("%s" % ("Annotation record ===\n"))
        FILE_HANDLE.write("%s%s%s" % ("Source:",self.source,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("Method:",self.method,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("Type:",self.annotationType,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("Contig:",self.contig,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("Start:",self.start,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("End:",self.end,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("Strand:",self.strand,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("Reading Frame:",self.readingFrame,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("Identifier:",self.identifier,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("Locus Tag:",self.locusTag,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("Name:",self.name,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("Description:",self.description,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("Category:",self.category,"\n"))
        FILE_HANDLE.write("%s%s%s" % ("Wraparound:",self.wraparound,"\n"))
        FILE_HANDLE.write("%s" % ("Annotation List:i\n"))
        for annot in self.annotationList:
            FILE_HANDLE.write("%s%s%s" % ("  ",annot,"\n"))
        FILE_HANDLE.write("%s" % ("Paralog List:\n"))
        for paralog in self.paralogList:
            FILE_HANDLE.write("%s%s%s" % ("  ",paralog,"\n"))
        FILE_HANDLE.write("%s" % ("=======================\n"))

    def writePVOGgroups(self,FILE_HANDLE):
        for pVOG in self.pVOGlist:
            pass 

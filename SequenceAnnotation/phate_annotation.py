###########################################################
#
# Name: phate_annotation.py
#
# Programmer: Carol L. Ecale Zhou
#
# Latest update: 28 December 2020
#
# Description: Module containing classes and methods for representing annotation results from various sources 
#
# Classes and methods: 
#     annotationRecord
#         addPVOBid2list
#         getPVOGassociationList
#         enterGFFdata
#         removeRedundancy
#         getDBXREFs
#         findInfo
#         getFigDescription
#         getPvogMembers
#         findPVOGannotation
#         findVOGannotation
#         getNCBItaxonomy
#         link2databaseIdentifiers
#         getECdescription4cazy
#         printAnnotationRecord_tabHeader
#         printAnnotationRecord_tab
#         printAnnotationRecord_tagged
#         printAnnotationDescription
#         printAnnotationRecord2file_tabHeader
#         printAnnotationRecord2file_tab
#         returnGFFannotationRecord
#         printAnnotationRecord2file
#         printAll
#         printAll2file
#         writePVOGgroups
##########################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE GPL3 LICENSE. SEE INCLUDED FILE GPL-3.PDF FOR DETAILS.

import re, os, subprocess
import platform

DEBUG = False
DEBUG = True

p_comment = re.compile('^#')

#*** Replace will platform call(s)
PHATE_MAC_OSX       = os.environ["PHATE_MAC_OSX"]
if PHATE_MAC_OSX == 'True':
    MAC_OSX = True
else:
    MAC_OSX = False

KEGG_VIRUS_BASE_DIR  = os.environ["PHATE_KEGG_VIRUS_BASE_DIR"]
NCBI_VIRUS_BASE_DIR  = os.environ["PHATE_NCBI_VIRUS_BASE_DIR"]
PHANTOME_BASE_DIR    = os.environ["PHATE_PHANTOME_BASE_DIR"]
NCBI_TAXON_DIR       = os.environ["PHATE_NCBI_TAXON_DIR"]
PVOGS_BASE_DIR       = os.environ["PHATE_PVOGS_BASE_DIR"]
VOGS_BASE_DIR        = os.environ["PHATE_VOGS_BASE_DIR"]
CAZY_ANNOTATION_PATH = os.environ["PHATE_CAZY_ANNOTATION_PATH"]
VOG_ANNOTATION_FILE  = os.environ["PHATE_VOG_ANNOTATION_FILE"]
VOG_HEADERS_FILE     = os.environ["PHATE_VOG_PROTEIN_HEADER_FILE"]

# Verbosity
CLEAN_RAW_DATA = os.environ["PHATE_CLEAN_RAW_DATA"]
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
pVOGheaderFile       = PVOGS_BASE_DIR      + "pVOGs.headers.lst"
VOGheaderFile        = VOGS_BASE_DIR       + "VOGs.headers.lst"     #*** ???
VOGannotationFile    = VOGS_BASE_DIR       + "vog.annotations.tsv"  #*** ???

class annotationRecord(object):

    def __init__(self):
        self.source            = "unknown" # Typically RAST, LLNL, PhAnToMe, GeneMark, Glimmer, Prodigal, PHANOTATE, KEGG, NCBI, pVOGs 
        self.method            = "unknown" # Typcially RAST, PFP, PhiRAST, JGI, SDSU, BLAST, blastp, blastn, HMM, jackhmmer 
        self.annotationType    = "unknown" # gene, mRNA, polypeptide, CDS, functional, homology, hmm
        self.VOGlist           = []        # list of pVOG or VOG identifiers (identified via blast hit)
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
        self.annotationString = ""      # used to construct a summary of annotation(s) for GFF output

    def addVOGid2list(self,VOG):
        self.VOGlist.append(VOG)

    def getPVOGassociationList(self):
        return(self.VOGlist)

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
                            print("phate_annotation says, WARNING: no dbxref found for", self.name, "in database", database, "given line", line)
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
                print("phate_annotation says, WARNING:  Unexpected name encountered in phate_annotation.py/getFigDescription:", self.name) 
        return dbxrefList 

    # Note:  The identifiers between VOGs and pVOGs do not correspond at all!
    # pVOGs.headers.lst has the 4-digit VOGid, followed by the protein's NCBI identifier, <space>, protein description
    # These headers were generated by my code, which queries each sequence from NR and constructs the fastas,
    # because the pVOG database does not provide fasta sequences (and apparently is not updated).
    # This is a sample header from the constructed pVOG fastas:
    # >VOG0510|ref|YP_009168084.1| nuclear disruption protein [Escherichia phage AR1]
    # Furthermore, the pVOG database does not summarize annotations per VOG group (per identifier).
    # In contrast...
    # VOGs.headers.lst has the 5-digit VOGid, followed by the protein's NCBI identifier. (no protein description)
    # >VOG05026|VOG07547|100220.YP_009506734.1
    # The VOG database *does* provide the sequences (gene and protein), but my code prepends the VOG identifier(s)
    # to each fasta header. The VOG-provided fasta sequence headers contain the NCBI identifier, but no description.
    # However, the VOG identifiers are annotated with a summary functional description. The pVOG identifiers are not!
    # Thus, VOG hits should be tagged with the summary function description created by the authors of the VOG database, 
    # whereas pVOG hits must each be tagged with that specific protein's header description (taken from NR).
    # This may seem confusing: basically, pVOG and VOG database are constructed differently, and therefore need to be
    # treated as completely separate databases.
    def getVOGmembers(self,database,vogType):  # database is the pVOGs or VOGs headers file; vogType is 'pvog' or 'vog'
        dbxrefList = []; infoList = []; annotationString = '' 
        p_vog = re.compile('(VOG\d+)\|')
        vogList = re.findall(p_vog,self.name)  # extract VOG identifiers from header (blast hit header)
        for vogID in vogList:
            if vogType.lower() == 'pvog':
                vogAnnotation = self.findPVOGannotation(vogID,database)
            elif vogType.lower() == 'vogs' or vogType.lower() == 'vog':
                vogAnnotation = self.findVOGannotation(vogID,database)
            annotationString += vogAnnotation + ' | '
        dbxrefList.append(annotationString)
        return dbxrefList

    #*** This methods should not be used because: pVOG hits are single protein members of a VOG, and
    #    there is no VOG summary annotation. The pVOGs.headers.lst file contains many proteins/headers
    #    corresponding to a given VOG. Picking the 1st occurrence in the file is not a valid approach.
    #    For pVOG hits, it is necessary to least each protein hit.
    def findPVOGannotation(self,vogID,database):   # database is the pVOGs.headers.lst file
        # database lines look like this, and are tab delimited:
        # VOG00012        118     96      Xu      sp|P03795|Y28_BPT7 Protein 2.8 
        vogAnnotation = ' '
        database_h = open(database,'r')
        dLines = database_h.read().splitlines()
        database_h.close()
        for dLine in dLines:
            match_vog = re.search(vogID,dLine)
            if match_vog:
                fields = dLine.split(' ')
                words  = fields[1:]
                vogAnnotation = vogAnnotation.join(words) 
                break
        return vogAnnotation

    #*** IMPROVE CODE HERE: calling method is passing database; doesn't need to, since annotation DB is environ var now
    # This method finds the summary annotation for a given VOG.
    def findVOGannotation(self,vogID,database):   # database is the vog.annotations.tsv file
        #database_h = open(database,'r')
        database_h = open(VOG_ANNOTATION_FILE,'r')
        dLines = database_h.read().splitlines()
        database_h.close()
        for dLine in dLines:
            match_vog = re.search(vogID, dLine)
            if match_vog:
                (vogID,proteinCount,speciesCount,functionalCategory,vogAnnotation) = dLine.split('\t')
                break
        return vogAnnotation

    # Query a taxonomy lookup table to get taxonomy information
    def getNCBItaxonomyID(self,database):   # Database maps ncbi header to taxonomy 
        # First, determine what operating system this code is running on.
        command = "platform.system()"  #*** THIS DOES NOT WORK IN CODE, BUT DOES ON COMMAND LINE; WHY???
        proc = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
        result, err = proc.communicate()
        #result = os.system(command)
        resultStr = str(result,'utf-8')

        # The database is the nucl_gb.accession2taxid, and the format of each line is as follows:
        # accn<\t>accn.version<\t>taxid<\t>gi
        # taxID is in the 3rd column of the accession2taxID file, and is one or more digits
        accession = ''; accessionVersion = ''; taxID = ''; giNumber = ''
        ncbiTaxonList = []

        # The name field of the current object (self) contains the fasta header of a hit sequence
        # Current formatting is as follows:
        # >accn.version<\s>Text description of genome
        # Ex: >NC_001798.2 Human herpesvirus 2 strain HG52, complete genome 
        # self.name is the fasta header sans '>'
        # WAIT:  Now (8/8/2020) the name field looks like this:
        # gi|33438897|ref|NC_005056.1| Bacteriophage WPhi, complete genome
        # Databases keep changing !!!
        IDENTIFIERS_INDEX = 0
        ACCESSION_INDEX   = 3

        words = []; accession = ""; description = ""
        p_taxID = re.compile('^[\d\w_]+\t[\d\w_]+\t([\d]+)\t[\d]+')  # 3rd column contains taxID (see above)

        # Split self.name field to extract accn, then search accession2tax file for the taxonomy id
        if self.name != '' and self.name.lower() != 'none' and self.name.lower() != 'unknown':

            # Extract accession and sequence description from header string
            words       = self.name.split(' ')     # There is a space separating the identifiers string from the description, which may contain spaces
            description = ' '.join(words[1:])      # Join all but the first item to construct the full description string
            fullId      = words[IDENTIFIERS_INDEX] # The identifiers string is several items separated by pipe
            idWords     = fullId.split('|')        # Separate the identifiers and labels
            accession   = idWords[ACCESSION_INDEX] # The accession is the 4th item: gi|giNumber|ref|accession.version|

            # Find line in database corresponding to this accession
            try:
                #*** THIS COMMAND WORKS ON MAC OSX; I DON'T KNOW IF IT WILL WORK ON UNIX
                #*** LC_ALL=C might crash on non-OSX systems.
                command = 'LC_ALL="C" grep -a \"' + accession + '\" ' + database
                proc = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
                out, err = proc.communicate()
            except:  # Other operating systems might fail the try, drop to this block
                command = 'grep -a \"' + accession + '\" ' + database
                proc = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
                out, err = proc.communicate()

            # Apparently, grep reads as binary... so need to convert result
            newOut = out.decode('utf8')
            if newOut != '':
                (accn,accn_v,taxID,gi) = newOut.split('\t')
                taxonomyString = 'NCBItaxID=' + taxID
                ncbiTaxonList.append(taxonomyString)
                ncbiTaxonLink = NCBI_TAXON_LINK + taxID
                ncbiTaxonList.append(ncbiTaxonLink)
            else:
                if PHATE_WARNINGS:
                    print("phate_annotation says, WARNING: Unable to find database entry for accession:", accession)
        else:
            if PHATE_WARNINGS:
                print("phate_annotation says, WARNING: NCBI hit header has improper format or is missing:", self.name)
        return ncbiTaxonList

    def link2databaseIdentifiers(self,database,dbName): 
        cazyAnnotationFile = CAZY_ANNOTATION_PATH 
        vogAnnotationFile  = VOG_ANNOTATION_FILE  # for VOG hits
        vogHeadersFile     = VOG_HEADERS_FILE     # for pVOG hits
        dbxrefList = []; VOGlist = []  # holds concatenation of functional annotations
        enzymeList = []; ncbiProteinList = []; taxonList = [] # hold specific annotations
        pfamList   = []; uniprotList     = []; koList    = [] # hold specific annotations
        figList    = []
        annotationString = "" # string containing all dbxref annotations found
        annotation = ""

        if self.name == "" or self.name == "none":
            if PHATE_WARNINGS:
                print("phate_annotation says, WARNING: No name for identification of dbxref in phate_annotation/link2databaseIdentifiers")
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

            elif dbName.lower() == 'ncbivirusgenome':
                pass

            elif dbName.lower() == 'ncbivirusprotein':
                pass

            elif dbName.lower() == 'refseqvirusprotein':
                pass

            elif dbName.lower() == 'refseqprotein':
                pass

            elif dbName.lower() == 'pvogs': # self.name = hit header, which contains VOGid and NCBIid + function description
                # If >1 VOGid, no matter, as it's a single protein sequence, with membership in >1 group.
                description = ''
                if self.name:
                    words = self.name.split(' ') 
                    description = words[1:]
                VOGlist = description

            elif dbName.lower() == 'vogs':   #*** To be deprecated
                VOGlist = self.getVOGmembers(vogAnnotationFile,'vog') 

            elif dbName.lower() == 'voggene':
                VOGlist = self.getVOGmembers(vogAnnotationFile,'vog') 

            elif dbName.lower() == 'vogprotein':
                VOGlist = self.getVOGmembers(vogAnnotationFile,'vog') 

            elif dbName.lower() == 'pvogs_hmm':
                # If >1 VOGid, no matter, as it's a single protein sequence, with membership in >1 group.
                description = ''
                if self.name:
                    words = self.name.split(' ') 
                    description = words[1:]
                VOGlist = description

            elif dbName.lower() == 'vogs_hmm':
                VOGlist = self.getVOGmembers(vogAnnotationFile,'vog') 

            elif dbName.lower() == 'ncbivirusprotein':
                pass
 
            elif dbName.lower() == 'cazy':   # use dbxrefList
                dbxrefList = self.getECdescription4cazy(cazyAnnotationFile)
 
            else:
                if PHATE_WARNINGS:
                    print("phate_annotation says, WARNING: Unrecognized database:", dbName) 

        for annotation in dbxrefList:
            nextAnnot = ' | ' + annotation
            annotationString += nextAnnot

        for vog in VOGlist:
            nextAnnot = ' | ' + vog
            annotationString += nextAnnot

        self.description = annotationString

        return 

    # Parses the CAZy annotation file to capture enzyme description(s) for blast hits.
    def getECdescription4cazy(self,cazyAnnotationFile):
        # self.name is the subject hit header
        dbxrefList = []; dbxrefCode = ""  # captures codes and annotations for CAZy hits
        code = ""; codeList = []  # captures codes that are listed in the headers of CAZy hits 
        annotationCode = ""; annotationDescription = ""  # from cazy annotation file; these are captured for reporting
        # Headers look like this: ">AWI06117.1|GT2|AT46|" (with one or more codes, pipe-separated)
        # Headers can also look like this: ">AAD03276.1|GH104|4.2.2.n1|" (with an EC number following, with or w/o terminal '|') 
        # ...or even like this: ">AUH33181.1|CE14"  (without terminal pipe)
        p_header = re.compile('\|([\w\d_^\.]+)\|') # Some headers have EC numbers following CAZy code   
        CAZY_ANNOT_H = open(cazyAnnotationFile,'r')
        try:
            # Pull dbxref from header; find data line in cazyAnnotationFile; add to dbxrefList
            match_header = re.search(p_header,self.name)
            if match_header:
                dbxrefCode = match_header.group(1)
                codeList = dbxrefCode.split('|')
                ANNOT_H = open(cazyAnnotationFile,'r')
                aLines = ANNOT_H.read().splitlines()
                for code in codeList:
                    if code != "": 
                        matchCode = code + '\t'
                        for aLine in aLines:
                            match_dbxref = re.search(matchCode,aLine)
                            if match_dbxref:
                                try:
                                    (annotationCode,annotationDescription) = aLine.split('\t')
                                    if code == annotationCode:  # Need to identify the identical one (else keep looking)
                                        dbxrefList.append(annotationDescription)
                                        break 
                                    else:
                                        continue
                                except:
                                    if PHATE_WARNINGS:
                                        print("phate_annotation says, WARNING: No annotation found for dbxrefCode",dbxrefCode)
                ANNOT_H.close()
            else:
                print("phate_annotation says, WARNING: getECdescription4cazy saw an unexpected subject hit header:",self.name)
        except:
            print("phate_annotation says, WARNING: Cannot identify dbxref for subject header",self.name)
        CAZY_ANNOT_H.close()
        return dbxrefList

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

    def printAnnotationRecord_tagged(self):
        annotationString = ""
        for annot in self.annotationList:
            annotationString += annot
            annotationString += ' | '
        print("source:",self.source,"method:",self.method,"annotationType:",self.annotationType,"category:",self.category)
        print("start:",self.start,"end:",self.end,"strand:",self.strand,"name:",self.name,"description:",self.description)
        print("annotationString:",annotationString)

    def printAnnotationDescription(self):
        print("description:",self.description)

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

    # Return annotations as a semicolon-delimited string
    def returnGFFannotationRecord(self,FILE_HANDLE):
        self.annotationString = ''; annot = ''; annotationList = []
        if self.annotationType.lower() == 'gene':
            annot = '(gene) ' + self.start + '/' + self.end + '/' + self.strand + ' ' + self.method 
            annotationList.append(annot)
        elif self.annotationType.lower() == 'functional':
            annot = '(function - ' + self.method + ') ' + self.description
            annotationList.append(annot)
        elif self.annotationType.lower() == 'homology':
            homologName = self.name
            newName = re.sub(';','',homologName)  # Remove embedded ';' because GFF uses this delimiter
            description = self.description
            newDescription = re.sub(';','',description)
            annot = '(homology - ' + self.method + ') ' + newName + ' ' + newDescription
            annotationList.append(annot)
        elif self.annotationType.lower() == 'hmm search':
            annot = '(hmm search - ' + self.method + ') ' + self.name + ' ' + self.description
            annotationList.append(annot)
        elif self.annotationType.lower() == 'profile search':
            # Note: descriptions for many pVOG group members can be voluminous; best to report name (pVOGid) only
            annot = '(profile search - ' + self.method + ') ' + self.name + ' ' + self.description  
            annotationList.append(annot)
        elif self.annotationType.lower() == 'cds':
            annot = '(cds - ' + self.method + ') ' + self.description
            annotationList.append(annot)
        elif self.annotationType.lower() == 'mrna':
            annot = '(mrna - ' + self.method + ') ' + self.description
            annotationList.append(annot)
        elif self.annotationType.lower() == 'trna':
            if self.description:
                annotationList.append(self.description)
        elif self.annotationType.lower() == 'polypeptide':
            annot = '(polypeptide - ' + self.method + ') ' + self.description
            annotationList.append(annot)
        else:
            annot = '(unk type - ' + self.method + ') ' + self.description
            annotationList.append(annot)
        if len(annotationList) > 0:
            for i in range(0, len(annotationList)):
                if i > 0:
                    self.annotationString += '; ' + annotationList[i]
                else:
                    self.annotationString += annotationList[i]
        FILE_HANDLE.write("%s" % (self.annotationString))

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
        for pVOG in self.VOGlist:
            pass 

###########################################################
# Module: phate_annotation.py
# Programmer: Carol L. Ecale Zhou
#
# Description: Module containing classes and methods for representing annotation results from various sources 
#
# Classes and methods: 
#     annotationRecord
#         addPVOBid2list
#         getPVOGassociationList
#         enterGFFdata(gff/dict)
#         setPSATparameters
#         recordPSATannotations
#         updatePSATcount
#         getDBXREFs
#         findInfo
#         getFigDescription
#         getPvogMembers
#         getNCBItaxonomy
#         link2databaseIdentifiers
#         printAnnotationRecord_tabHeader
#         printAnnotationRecord_tab
#         printAnnotationRecord2file_tabHeader
#         printAnnotationRecord2file_tab
#         printAnnotationRecord2file
#         returnGFFannotationRecord
#         printAll
#         printAll2file(fileH)
#         writePVOGgroups
##########################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.PDF FOR DETAILS.

import re, os, subprocess

#DEBUG = True
DEBUG = False

p_comment = re.compile('^#')

KEGG_VIRUS_BASE_DIR = os.environ["PHATE_KEGG_VIRUS_BASE_DIR"]
NCBI_VIRUS_BASE_DIR = os.environ["PHATE_NCBI_VIRUS_BASE_DIR"]
PHANTOME_BASE_DIR   = os.environ["PHATE_PHANTOME_BASE_DIR"]
NCBI_TAXON_DIR      = os.environ["PHATE_NCBI_TAXON_DIR"]
PVOGS_BASE_DIR      = os.environ["PHATE_PVOGS_BASE_DIR"]

# Verbosity
CLEAN_RAW_DATA      = os.environ["PHATE_CLEAN_RAW_DATA"]
PHATE_WARNINGS      = os.environ["PHATE_PHATE_WARNINGS"]
PHATE_MESSAGES      = os.environ["PHATE_PHATE_MESSAGES"]
PHATE_PROGRESS      = os.environ["PHATE_PHATE_PROGRESS"]

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
        self.method            = "unknown" # Typcially RAST, PSAT, PFP, PhiRAST, JGI, SDSU, BLAST, blastp, blastn, HMM, jackhmmer 
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
        self.psat = {
            "jobID"    : "",   # PSAT job id
            "jobName"  : "",   # PSAT job name
            "fileName" : "",   # PSAT output file
            }
        self.psatOutDir = ""   # need to set
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

    # METHODS FOR ANNOTATING FROM EXTERNAL SOURCES (e.g., PSAT)

    def setPSATparameters(self,jobID,jobName,fileName,psatOutDir):
        self.psat["jobID"]    = jobID
        self.psat["jobName"]  = jobName
        self.psat["fileName"] = fileName
        self.source = "LLNL"
        self.method = "PSAT"
        self.annotationType = "functional"
        self.psatOutDir = psatOutDir
        PSAT_OUT_DIR = psatOutDir

    def removeRedundancy(self,inList): # Eliminate redundancy in list; Different PSAT annotations sources can return same annotation
        outList = []
        for i in range(len(inList)):
            item = inList.pop()
            if item not in inList:
                outList.append(item)
        outList.reverse()
        return outList

    def recordPSATannotations(self,proteinHeader,PSAT_H):  # Query PSAT annotation file for current gene's annotations 
        # Locations of annotations in PSAT annotation file
        EC_COLUMN       = 3
        EC_DESC_COLUMN  = 4
        PATHWAY_COLUMN  = 5
        INTERPRO_COLUMN = 13
        SIGNALP_COLUMN1 = 14
        SIGNALP_COLUMN2 = 15
        #TMHMM_COLUMN    = 25  # Not in service. TMHMM could be added to PSAT output in column 25

        # Patterns
        p_comment  = re.compile('^#')
        p_EC       = re.compile('\d\.\d+\.\d+\.\d+')
        p_ECdesc   = re.compile('[\w\d\.\-\_\s]+')
        p_pathway  = re.compile('ec\d{4}.*')
        p_GO       = re.compile('(MOLECULAR_FUNCTION|BIOLOGICAL_PROCESS|CELLULAR_COMPONENT),(GO:\d{7})')
        p_pfam     = re.compile('([\w\d\s\.\-\_\,]+)PFAM')
        p_smart    = re.compile('([\w\d\s\.\-\_\,]+)SMART')
        p_dbxref   = re.compile('EMPTY4NOW')
        p_signalp1 = re.compile('YES|NO')
        p_signalp2 = re.compile('\'(\d+)\-(\d+)\'')
        p_tmhmm    = re.compile('Topology=([io][io\d\-]+)\|')

        # 
        annotationString = ""
        psatAnnotationList = []
        tempList = []; columns = []

        if PHATE_MESSAGES == 'True':
            print("Annotation module says: Recording PSAT annotations.")

        ### Capture lines corresponding to the gene
        pLines = PSAT_H.read().splitlines()
        if DEBUG:
            print("There are this many pLines:", len(pLines))

        # Capture all annotation lines for this gene
        for pLine in pLines:
            matchTerm = proteinHeader + '\t'
            match_gene = re.search(matchTerm,pLine)
            match_comment = re.search(p_comment,pLine)
            if match_gene:
                tempList.append(pLine)
        if DEBUG:
            print("Protein", proteinHeader, "Temp list:", tempList)

        ### Parse each annotation line, capture annotations and format for genbank
        for line in tempList:

            if DEBUG:
                print("Processing PSAT hit line:", line)
            annotationString = ""
            EC = ""
            ecDescription = ""
            pathway = ""
            GO = ""
            pfam = ""
            smart = ""
            dbxrefID = ""
            signalp = ""; signalp1 = ""; signalp2 = ""
            tmhmm = ""
            columns = line.split('\t')

            ### Detect annotations; Note: There could be >1 in a given data line, or in a given column

            match_ec       = re.search(p_EC,       columns[EC_COLUMN])
            match_ecDesc   = re.search(p_ECdesc,   columns[EC_DESC_COLUMN])
            match_pathway  = re.search(p_pathway,  columns[PATHWAY_COLUMN])
            match_go       = re.search(p_GO,       columns[INTERPRO_COLUMN])
            match_pfam     = re.search(p_pfam,     columns[INTERPRO_COLUMN])
            match_smart    = re.search(p_smart,    columns[INTERPRO_COLUMN])
            match_dbxref   = re.search(p_dbxref,   columns[INTERPRO_COLUMN])
            match_signalp1 = re.search(p_signalp1, columns[SIGNALP_COLUMN1])
            match_signalp2 = re.search(p_signalp2, columns[SIGNALP_COLUMN2])
            #match_tmhmm    = re.search(p_tmhmm,    columns[TMHMM_COLUMN])  # Not currently in service

            if match_ec:
                EC = match_ec.group(0)
                annotationString = "EC:" + EC
                self.annotationList.append(annotationString)

            if match_ecDesc:
                ecDesc = match_ecDesc.group(0)
                annotationString = "EC description:" + ecDesc
                self.annotationList.append(annotationString)

            if match_pathway:
                pathwayString = match_pathway.group(0)
                pathways = pathwayString.split('|\s')
                for pathway in pathways:
                    annotationString = "Pathway:" + pathway + ' | '
                    self.annotationList.append(annotationString)

            if match_go:
                GO = match_go.group(0)
                category = match_go.group(1)
                goID     = match_go.group(2)
                if category == "BIOLOGICAL_PROCESS":
                    annotationString = "go_process:" + goID
                elif category == "MOLECULAR_FUNCTION":
                    annotationString = "go_function:" + goID
                elif category == "CELLULAR_COMPONENT":
                    annotationString = "go_component:" + goID
                else:
                    annotationString = "go_unknown:" + goID
                self.annotationList.append(annotationString)

            if match_pfam:
                pfam = match_pfam.group(1)
                pfam = re.sub('\,$','',pfam)
                annotationString = "pfam:" + pfam
                self.annotationList.append(annotationString)

            if match_smart:
                smart = match_smart.group(1)
                smart = re.sub('\,$','',smart)
                annotationString = "smart:" + smart
                self.annotationList.append(annotationString)

            if match_dbxref:
                dbxrefID = match_dbxref.group(0)
                for dbxref in match_dbxref:
                    annotationString = "dbxref:" + dbxrefID
                    self.annotationList.append(annotationString)

            if match_signalp1:
                signalp1 = ''; signalp2 = ''; start = 0; end = 0
                signalp1 = match_signalp1.group(0)
                if signalp1 == 'YES':
                    if match_signalp2:
                        signalp2 = match_signalp2.group(0)
                        start    = match_signalp2.group(1)
                        end      = match_signalp2.group(2)
                    else:
                        if DEBUG:
                            print("DID NOT FIND A signalp2 MATCH for YES", pLine)
                    
                annotationString = "signal_peptide:" + signalp1 + ' ' + signalp2 
                self.annotationList.append(annotationString)

            # Not in service
            #if match_tmhmm:
            #    tmhmm = match_tmhmm.group(1)
            #    print "found tmhmm:", tmhmm
            #    annotationString = str(start) + '\t' + str(end) + "\tfeature\tTMHMM\t" + tmhmm 
            #    #psatAnnotationList.append(annotationString) 

        self.annotationList = self.removeRedundancy(self.annotationList)
        self.updatePSATcount()

        return 

    def updatePSATcount(self):
        if self.annotationList == []:
            self.description = "There are no PSAT annotations"
        else:
            self.description = "Num annots:" + str(len(self.annotationList))


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
                        if PHATE_WARNINGS == 'True':
                            print("WARNING in annotation module: no dbxref found for", self.name, "in database", database, "given line", line)
                    idList.append(dbxref)
        return idList

    def findInfo(self,searchTerm,database):  # Searches Phantome header file for annotation information
        p_truncatedSearchTerm = re.compile('^([^\s]*)\s')
        infoLines = []
        DATABASE_H = open(database,"r")
        dLines = DATABASE_H.read().splitlines()
        if DEBUG:
            print("TESTING: original searchTerm is", searchTerm)
        match_truncate = re.search(p_truncatedSearchTerm,searchTerm)
        if match_truncate:
            truncatedSearchTerm = match_truncate.group(1)
            searchTerm = truncatedSearchTerm
            if DEBUG:
                print("TESTING: searchTerm was changed to", searchTerm) 
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
            if PHATE_WARNINGS == 'True':
                print("WARNING in annotation module:  Unexpected name encountered in phate_annotation.py/getFigDescription:", self.name) 
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
        if self.name == '' and self.name == 'none':
            if DEBUG:
                print("name field blank in getNCBItaxonomy") 
        else:
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
                if DEBUG:
                    print("command is", command) 
                proc = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
                out, err = proc.communicate()
                if DEBUG:
                    print("Result of grep is", out)
                if out != '':
                    match_taxID = re.search(p_taxID,out)
                    taxonomyID = match_taxID.group(1)
                    taxonomyString = 'NCBItaxID=' + taxonomyID
                    ncbiTaxonList.append(taxonomyString)
                    ncbiTaxonLink = NCBI_TAXON_LINK + taxonomyID
                    ncbiTaxonList.append(ncbiTaxonLink)
            else:
                if PHATE_WARNINGS == 'True':
                    print("WARNING: NCBI hit header has improper format or is missing:", self.name)
        if DEBUG:
            print("ncbiTaxonList is", ncbiTaxonList)
        return ncbiTaxonList

    def link2databaseIdentifiers(self,database,dbName):
        dbxrefList = []  # holds concatenation of functional annotations
        enzymeList = []; ncbiProteinList = []; taxonList = [] # hold specific annotations
        pfamList   = []; uniprotList     = []; koList    = [] # hold specific annotations
        figList    = []
        annotationString = "" # string containing all dbxref annotations found
        annotation = ""

        if self.name == "" or self.name == "none":
            if PHATE_WARNINGS == 'True':
                print("No name for identification of dbxref in phate_annotation/link2databaseIdentifiers")
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

            elif dbName.lower() == 'ncbivirusprotein':
                pass
 
            else:
                if PHATE_WARNINGS == 'True':
                    print("WARNING in annotation module: Unrecognized database:", dbName) 

        if DEBUG: 
            print("dbxrefList:", dbxrefList)

        for annotation in dbxrefList:
            nextAnnot = ' | ' + annotation
            annotationString += nextAnnot
        self.description = annotationString

        return 

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

    # Return annotations as a semicolon-delimited string
    def returnGFFannotationRecord(self,FILE_HANDLE):
        self.annotationString = ''; annot = ''; annotationList = []
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
        elif self.annotationType == 'cds':
            annot = '(cds)' + self.method + ' ' + self.description
        elif self.annotationType == 'mrna':
            annot = '(mrna)' + self.method + ' ' + self.description
        elif self.annotationType == 'polypeptide':
            annot = '(polypeptide) ' + self.method + ' ' + self.description
        else:
            annot = '(unk type) ' + self.method + ' ' + self.description
        annotationList.append(annot)

        if len(annotationList) > 0:   #*** IMPROVE THIS, there should be only one annot in list
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
        for pVOG in self.pVOGlist:
            pass 

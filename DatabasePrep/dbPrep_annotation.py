###########################################################
# Module: dbPrep_annotation.py
# Programmer: Carol L. Ecale Zhou
#
# Data of last update:  October 2016 - code being modified from CGP code base
#    03 January 2017 - modified output report in method printAnnotationRecord()
#    05 January 2017 - adding code to pull dbxrefs from local database instances
#    17 July 2018 - adapted for use in building pVOGs fasta database
#
# Module containing classes and methods for representing annotation results from various sources 
# Classes and methods: 
#     annotationRecord
#         enterGFFdata(gff/dict)
#         printAnnotationRecord
#         printAnnotationRecord2file
#         printAll
#         printAll2file(fileH)
#########################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.pdf FOR DETAILS

import re
import os
import subprocess 

CHATTY = True
DEBUG = False

p_comment = re.compile('^#')

class annotationRecord(object):

    def __init__(self):
        self.source            = "unknown" # Typically RAST, LLNL, PhAnToMe, GeneMark, Glimmer, Prodigal, PHATE, KEGG, NCBI 
        self.method            = "unknown" # Typcially RAST, PSAT, PFP, PhiRAST, JGI, SDSU, Blast, blastp, blastn 
        self.annotationType    = "unknown" # gene, mRNA, polypeptide, CDS, functional, homology
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

    def removeRedundancy(self,inList): # Eliminate redundancy in list; Different PSAT annotations sources can return same annotation
        outList = []
        for i in range(len(inList)):
            item = inList.pop()
            if item not in inList:
                outList.append(item)
        outList.reverse()
        return outList

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


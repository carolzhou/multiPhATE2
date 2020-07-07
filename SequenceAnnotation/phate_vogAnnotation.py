###########################################################################
# Module:  phate_vogAnnotation.py
# Programmer:  C. E. Zhou
#
# Latest update:  17 June 2020
#
# Classes and Methods:
# 
#
###########################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE BSD LICENSE. SEE INCLUDED FILE BSD.PDF FOR DETAILS.

import re, os, subprocess

#DEBUG = True
DEBUG = False

class vogDescription(object):

    def __init__(self):
        self.identifier             = ""
        self.description            = ""
        self.vogDescriptionFile     = ""   # path/filename to file that maps VOGid to description

    def getVOGmembers(self,database,inString,vogType):  # inString has a VOG identifier; vogType is 'pvog' or 'vog'
        # database lines look like this:
        # >VOG0510|ref|YP_009168084.1| nuclear disruption protein [Escherichia phage AR1]
        dbxrefList = []; infoList = []; annotationString = ''
        p_vog = re.compile('(VOG\d+)\|')
        vogList = re.findall(p_vog,inString)  # Extract VOG identifiers from header (blast hit header)
        for vogID in vogList:
            if vogType.lower() == 'pvog':
                vogAnnotation = self.findPVOGannotation(vogID,database)
            elif vogType.lower() == 'vog':
                vogAnnotation = self.findVOGannotation(vogID,database)
            annotationString += vogAnnotation + ' | '
        dbxrefList.append(annotationString)
        return dbxrefList

    def findPVOGannotation(self,vogID,database):  # database is the pVOGs.headers.lst file
        # database lines look like this, and are tab deliminted:
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

    def findVOGannotation(self,vogID,database):
        database_h = open(database,'r')
        dLines = database_h.read().splitlines()
        database_h.close()
        for dLine in dLines:
            match_vog = re.search(vogID, dLine)
            if match_vog:
                (vogID,proteinCount,speciesCount,functionalCategory,vogAnnotation) = dLine.split('\t')
                break
        return vogAnnotation

    def printAll(self):
        print("VOG identifier:",self.identifier)
        print("File mapping identifier to description:",self.vogDescriptionFile)
        print("VOG description:",self.description)


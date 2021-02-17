#!/usr/bin/env python3
##############################################################################################
#
# script scrutinizeCAZY.py
#
# programmer: C. Zhou
# date of last update: 12 February 2021
#
# purpose: There are inconsistencies in the CAZY database headers, and more variablility than
#   I may have anticipated. This script walks through the 1.3M headers from the CAZY sequence
#   files and determines which headers are unaccounted for based on a given regular expression.
#   The occasional errors or inconsistencies affect some entries that appear to possibly have
#   been truncated. If an entry ends without a pipe '|', it is unclear whether there has been
#   a truncation in preparation of this data source. This question should be pursued with the
#   database originators.
#
##############################################################################################

import os, re

p_header_1 = re.compile('\|([\w\d_^\.]+)\|')      # captures single ID, followed by pipe 
p_header_2 = re.compile('\|([\w\d_^\.\|]+)+\|')   # captures multiple IDs, followed by pipe 
p_header_3 = re.compile('\|([\w\d_^\.^\|]+)')     # captures single ID, not followed by pipe 

print("Scrutinizing CAZY headers")
count_yes = 0; count_no = 0; count_ids = 0; count_multiple = 0
total = 1386849
idList = []; idString = ""
HEADERS_H = open("../../Databases/CAZY/CAZyDB.07312019.headers",'r')
MISSED_H  = open("./missedHeaders.lst",'w')
hLines = HEADERS_H.read().splitlines()
for hLine in hLines:
    match_1 = re.search(p_header_1,hLine)
    match_2 = re.search(p_header_2,hLine)
    match_3 = re.search(p_header_3,hLine)
    if match_2:
        print("header:",hLine,"idString:",match_2.group(0))
        count_yes += 1
        idString = match_2.group(0) 
        idString_clean = idString.lstrip('|')
        idList   = idString_clean.split('|')
        if idList[-1] == '':
            idList.pop()
        count_IDs = len(idList)
        if count_IDs > 1:
            print("Found a set of IDs: ",idString_clean, idList)
            count_multiple += 1
    elif match_3:
        count_yes += 1
        idString = match_3.group(0)
        print("Found a single ID:",idString)
    else:
        count_no += 1
        MISSED_H.write("%s\n" % (hLine))
MISSED_H.close()
HEADERS_H.close()

percent = (count_yes / total) * 100
print("Number of headers accounted for: ",count_yes)
print("Number of missed headers: ",count_no)
print("Percent of headers accounted for: ",percent)
print("Number of headers with multiple IDs: ",count_multiple)
print("Done!")

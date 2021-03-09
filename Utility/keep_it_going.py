#!/user/bin/env python3

#######################################################################
#
# Name: keep_it_going.py
#
# Programmer:  Carol Zhou
#
# Last Update:  21 December 2020
#
# Description:  Run this script in the background to keep you ssh window
#  active, for example when downloading databases onto a remote server.
#  
#######################################################################

# This code was developed by Carol L. Ecale Zhou at Lawrence Livermore National Laboratory.
# THIS CODE IS COVERED BY THE GPL3 LICENSE. SEE INCLUDED FILE GPL-3.pdf FOR DETAILS.

import time
while True:
    time.sleep(500)
    print("Staying active...")

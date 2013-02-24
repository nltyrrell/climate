#!/usr/bin/python

import os
import glob
import sys

"""Run qsub command, and add stuff to it
"""
listtxt = sys.argv[1]

#with open(listtxt, 'r') as f:
#	dirs = f.readlines()

#for i in dirs:
outstring = 'echo '+listtxt #'qsub -v ID=' + i + ' tarballing.sh'
os.system(outstring)


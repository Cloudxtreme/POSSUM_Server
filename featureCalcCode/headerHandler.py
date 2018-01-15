#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#Copyright Chris & Young

def usage():
	print "headerHandler.py usage:"
	print "python headerHandler.py <options> <source files> "
	print "-i,--input: input a fasta format file."
	print "-o,--ouput: output a file of position-specific scoring matrix with header."
	print "-p,--prefix: input the prefix of header."
	print "-n,--number: "
	print "-h,--help: show the help information."

import re
import pandas as pd
import numpy as np
import fileinput
import sys, getopt
from os import listdir
from os.path import isfile, join
from sklearn import preprocessing

opts, args = getopt.getopt(sys.argv[1:], 'i:o:p:n:h', ['input=','output=','prefix=','number=','help'])
inputFile=""
outputFile=""
prefix=""
colNumber=""

for opt, arg in opts:
	if opt in ('-i','--input'):
		inputFile = arg
	elif opt in ('-o','--output'):
    		outputFile = arg
	elif opt in ('-p','--prefix'):
		prefix = arg
	elif opt in ('-n','--number'):
		colNumber = int(arg)
	elif opt in ('-h', '--help'):
		usage()
    		sys.exit(2)
  	else:
  		usage()
    		sys.exit(2)

df_input = pd.read_csv(inputFile, header=None)

headerRow=[]
for i in range(colNumber):
	headerRow.append(prefix+str(i))
df_input.columns=headerRow

df_input.to_csv(outputFile,index=False)

#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#Copyright Chris & Young

def usage():
	print "pssm_ac.py usage:"
	print "python pssm_ac.py <options> <source files> "
	print "-i,--input: input a fasta format file."
	print "-o,--ouput: output a file of position-specific scoring matrix with same property."
	print "-t,--type: the algorithm you select."
	print "-p,--pssmdir: input the directory of pssm files."
	print "-h,--help: show the help information."

import fileinput
import sys, getopt
from os import listdir
from os.path import isfile, join
import re
import numpy as np
from possum_ft import *

opts, args = getopt.getopt(sys.argv[1:], 'i:o:t:p:a:b:h', ['input=','output=','type=','pssmdir=','argument=','veriable=','help'])
inputFile=""
outputFile=""
algoType=""
pssmdir=""
argument=""
veriable=""

for opt, arg in opts:
	if opt in ('-i','--input'):
		inputFile = arg
    	elif opt in ('-o','--output'):
        	outputFile = arg
      	elif opt in ('-t','--type'):
        	algoType = arg
	elif opt in ('-p','--pssmdir'):
		pssmdir = arg
	elif opt in ('-a','--argument'):
		argument = int(arg)
	elif opt in ('-b','--veriable'):
		veriable = int(arg)
	elif opt in ('-h', '--help'):
		usage()
        	sys.exit(2)
      	else:
      		usage()
        	sys.exit(2)
check_head = re.compile(r'\>')  

smplist = [] #列表保存序列信息
smpcnt = 0
for line, strin in enumerate(fileinput.input(inputFile)):
	if not check_head.match(strin):
		smplist.append(strin.strip())
		smpcnt += 1
#print "smplist="
#print smplist
onlyfiles = [ f for f in listdir(pssmdir) if isfile(join(pssmdir,f)) ] #保存单独的pssm文件的名字
#print "onlyfiles="
#print onlyfiles
fastaDict = {}

for fi in onlyfiles: #循环文件名列表
	cntnt = ''
	#for line, strin in enumerate(fileinput.input(pssmdir+'/'+fi)): #循环单独的pssm文件内容
	pssmContentMatrix=readToMatrix(fileinput.input(pssmdir+'/'+fi))
	pssmContentMatrix=np.array(pssmContentMatrix)
	sequence=pssmContentMatrix[:,0]
	seqLength=len(sequence)
	for i in range(seqLength):
		cntnt+=sequence[i]
	if cntnt in fastaDict:
		print strin
		continue
	fastaDict[cntnt] = fi #字典保存{序列：文件名}

finalist = []
for smp in smplist: #循环序列列表
	#print "smp="+smp
	#print "fastaDict[smp]="+fastaDict[smp]
	finalist.append(pssmdir+'/'+fastaDict[smp])

file_out = file(outputFile,'w')

for fi in finalist:
	#print "fi=" + fi
	#pass in original matrix
	input_matrix=fileinput.input(fi)
	#output a feature vector
	#input_matrix=readToMatrix(input_matrix)
	feature_vector=calculateDescriptors(input_matrix, algoType, argument, veriable)
	np.savetxt(file_out, feature_vector, delimiter=",")



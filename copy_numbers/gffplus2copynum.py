#! /usr/bin/python

'''
Generates read-based copy number without correction for gene lengths from gff+ files
'''

import string
from sys import argv
import re

script, gffplus, avgSCOcoverage = argv

IFH = open(gffplus,'r')
lines = IFH.readlines()
avgSCOcov=float(avgSCOcoverage)

#### reading the gff+ file
geneLengths=[]
coverages=[]
geneCount=0
for line in lines:
	cols=line.split("\t")
	s=int(cols[3])
	e=int(cols[4])
	genelen=e-s
	geneLengths.append(genelen)
	geneCount+=1
	coverages.append(float(cols[11]))
	
avgGeneLen=sum(geneLengths)/len(geneLengths) 
avgCov=sum(coverages)/len(coverages)
copynum=(avgCov*geneCount)/avgSCOcov
#print avgCov, geneCount, avgSCOcov
#print "copy number:", copynum
#print "average gene length:", avgGeneLen
list2print=[str(gffplus),str(copynum), str(avgGeneLen)]
string2print="\t".join(list2print)
print string2print



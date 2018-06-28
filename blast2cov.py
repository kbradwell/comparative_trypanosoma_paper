#! /usr/bin/python

'''
Input: 
BLAST file
Output:
List of genes that are legitimately "FOUND" based on the reads and the given criteria. 
'''

import string
from sys import argv
import re
import os

script, blastOut, outGenes = argv

allBLAST = open(blastOut)
blastLines = allBLAST.readlines()

of=open(outGenes,'w') 
queryDict={}
queryLenDict={}
for line in blastLines:
	lineSt = line.strip()
	line2 = lineSt.split()
	genename = line2[0]
	readname = line2[1]
	sgene2=genename.strip()
	gsStr = line2[2].strip()
	geStr = line2[3].strip()
	gs=int(gsStr)
	ge=int(geStr)
	glen = line2[9].strip()
	if sgene2 not in queryLenDict:
		queryLenDict[sgene2]=int(glen)
	if sgene2 in queryDict:
		for pos in range(gs,ge+1):
			if pos in queryDict[sgene2]: 
				queryDict[sgene2][pos]+=1
			else:
				queryDict[sgene2][pos]=1
	else:
		queryDict[sgene2]={}
		for pos in range(gs,ge+1):
                	if pos in queryDict[sgene2]:
                        	queryDict[sgene2][pos]+=1
                	else:
                        	queryDict[sgene2][pos]=1

for g in queryLenDict:
	posCounter=0
	#print g
	for position in range(1,queryLenDict[g]):
		if position in queryDict[g]:
			if queryDict[g][position]>3:
				posCounter+=1
				#print position,"\t", queryDict[g][position]
			else:
				pass
				#print position, "\t", queryDict[g][position]
		else:
			queryDict[g][position]=0
			#print position, "\t",queryDict[g][position]
	percOver3Reads = ((float(posCounter)/float(queryLenDict[g]))*100)
	print g, "\t", percOver3Reads
	if (percOver3Reads>60):
		of.write(g)
		of.write("\n")
	
allBLAST.close()	
of.close()

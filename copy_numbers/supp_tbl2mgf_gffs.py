#! /usr/bin/python

'''
Script makes the MGF GFF files that are required for input to estimate copy num.
Input: 
<sp.>_genes_supp_tbl - a file that contains the general genes annotation table for a species
Output: 
GFF files of genes for each MGF of the selected species.
'''

import string
from sys import argv
import re

script, suppTbl = argv

supp = open(suppTbl)
lines = supp.readlines()

#### reading the supp table file

mucinList=["SMUG L", "SMUG S", "TASV", "mucin"]
mgfDict={}
for line in lines:
	mgfName="other"
	line = line.strip()
	cols = line.split("\t")
	genename = cols[0].strip()
	contig = cols[1].strip()
	start = cols[2].strip()
	end = cols[3].strip()
	strand = cols[5].strip()
	eval = cols[11].strip()
	dbAnnotation = cols[13]
	infoList=[contig,"nrBLAST","Prot",str(start),str(end),str(eval),str(strand),".",genename]
	# obtain the MGF classification
	if "sialidase" in dbAnnotation:
		# make annotation into name
		nameList=dbAnnotation.split("[")
		mgfNameJoined=" ".join(nameList[0:-1])
		mgfNameJoined=mgfNameJoined.replace('partial','')
		mgfNameJoined=mgfNameJoined.replace('putative','')
		mgfNameJoined=mgfNameJoined.replace('pseudogene','')
		mgfNameJoined=mgfNameJoined.replace(',','')
		mgfNameTemp=mgfNameJoined.strip()
		mgfName=mgfNameTemp.replace (" ", "_")
	elif ("MASP" not in dbAnnotation and any(name in dbAnnotation for name in mucinList)):
		# make annotation into name
		nameList=dbAnnotation.split("[")
		mgfNameJoined=" ".join(nameList[0:-1])
		mgfNameJoined=mgfNameJoined.replace('partial','')
		mgfNameJoined=mgfNameJoined.replace('putative','')
		mgfNameJoined=mgfNameJoined.replace('pseudogene','')
		mgfNameJoined=mgfNameJoined.replace(',','')
		mgfNameTemp=mgfNameJoined.strip()
		mgfName=mgfNameTemp.replace (" ", "_")
	elif ("MASP" in dbAnnotation or "mucin-associated surface protein" in dbAnnotation):
		mgfName="MASP"
	elif "RHS" in dbAnnotation:
		mgfName="RHS"
	elif ("GP63" in dbAnnotation or "gp63" in dbAnnotation):
		mgfName="GP63"
	elif "amastin" in dbAnnotation:
		mgfName="amastin"
	elif ("dispersed protein family protein 1" in dbAnnotation or "DGF-1" in dbAnnotation):
		mgfName="DGF-1"
	elif "beta galactofuranosyl glycosyltransferase" in dbAnnotation or "beta-galactofuranosyl transferase" in dbAnnotation:
		mgfName="beta-galac-gt"
	elif "cruzipain" in dbAnnotation:
		mgfName="cruzipain-precursor"
	elif "KMP-11" in dbAnnotation:
		mgfName="KMP-11"
	elif "syntaxin binding protein" in dbAnnotation:
		mgfName="syntaxin-BP"
	elif "target of rapamycin (TOR) kinase 1" in dbAnnotation:
		mgfName="TORkinase1"
	elif "UDP-GlcNAc" in dbAnnotation:
		mgfName="UDPgal"
	# add hit info to the mgf dict
	if mgfName in mgfDict:
		mgfDict[mgfName].append(infoList)
	else:
		mgfDict[mgfName]=[]
		mgfDict[mgfName].append(infoList)

# for each MGF, make a GFF file

for m in mgfDict:
	orgFile=str(suppTbl)
	subsetName="mgf"+"_"+m+"_"+orgFile+".gff"
	subset = open(subsetName,"w")
	for c in mgfDict[m]:
		outStr="\t".join(c)
		subset.write(outStr.strip()+"\n")
        subset.close()

print "the Subset GFF files are now created"





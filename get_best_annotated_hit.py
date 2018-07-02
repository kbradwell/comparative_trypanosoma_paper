#! /home/vnkoparde/opt/Python-2.7.2/bin/python

# Author: Katie Bradwell
# named get_best_annotated_hit_v2.py on author's server

import sys
from sys import argv
import string
import re
import os
import itertools
import subprocess

'''
script gets the top 10 gb and ref hits from nr for each called gene, performs blastdbcmd for each hit, then returns the top non hypothetical/unknown
hit for each gene along with the original BLAST result info. nrblastresults input must be sorted by genename (first column), then e-value.
'''

script, blastres, outfile = argv

top10hits={}
permittedDB=["gb","ref"]
infile = open(blastres)
lines = infile.readlines()
for line in lines:
	line = line.strip()
	cols = line.split("\t")
	genename=cols[0]
	#if "Tra_15.4_gene_00006" in genename:
	#	break
	hitID=cols[1]
	idParts=hitID.split("|")
	dbType=idParts[2]
	eval=float(cols[10])
	hitList=[hitID,line]
	if genename not in top10hits and (any(x in dbType for x in permittedDB)):
		counter=1
		top10hits[genename]=[]
		top10hits[genename].append(hitList)
	else:
		if counter<10 and (any(x in dbType for x in permittedDB)) and eval< 1e-10:
			counter+=1
			top10hits[genename].append(hitList)
	# the top 10 BLAST hits have now been obtained for each called gene

# get the descriptions for each of the top 10 genes using blastdbcmd

bestHitDesc={}
for g in top10hits:
	seen=0
	first=0
	for h in top10hits[g]:
		id=h[0]
		parts=id.split("|")
		gi=parts[1]
		#print g,gi
		dbt='prot'
		try:
			proc = subprocess.Popen(["/usr/global/blp/rmblast-1.2-ncbi-blast-2.2.23+/bin/blastdbcmd","-db","/gpfs_fs/data1/refdb/nr", "-dbtype", dbt, "-entry",gi], stdout=subprocess.PIPE)
                	#print g
                	tmp = proc.stdout.read()
                	tlList=tmp.split("\n")
                	descList=tlList[0].split(">")
                	firstDesc=descList[1]
                	#print firstDesc
                	#print "gi: ",gi
                	firstDescParts=firstDesc.split("|")
                	nonIDdesc=firstDescParts[4].strip()
                	#print nonIDdesc
                	if first==0:
                        	h.append(nonIDdesc)
                        	bestHitDesc[g]=h
                        	first=1
                	if ("hypothetical" not in nonIDdesc.lower()) and ("unknown" not in nonIDdesc.lower()) and (seen==0):
                        	if nonIDdesc not in h:
                                	h.append(nonIDdesc)
                        	bestHitDesc[g]=h
                        	seen=1
		except:
			continue


# format the output into a tab delimited file that can be used as input to create the new supp table for all genes

of=open(outfile,'w')
for g in sorted(bestHitDesc):
	hitInf=bestHitDesc[g][1].split('\t')
	percSubjLenFloat=((float(hitInf[9])-float(hitInf[8])+1)/float(hitInf[12]))*100
	percSubjLen='%.2f' % percSubjLenFloat
	outList=[hitInf[0],hitInf[1],str(hitInf[2]),str(hitInf[3]),str(hitInf[4]),str(hitInf[5]),str(hitInf[6]),str(hitInf[7]),str(hitInf[8]),str(hitInf[9]),str(hitInf[10]),str(hitInf[11]),str(hitInf[12]),str(percSubjLen),bestHitDesc[g][2]]
	outStr="\t".join(outList)
	of.write(outStr)
	of.write("\n")

			



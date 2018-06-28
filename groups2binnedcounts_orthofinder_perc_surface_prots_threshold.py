#! /usr/bin/python

'''
Input: 
OrthologousGroups.txt
Output: 
frequencies of binned number of genes in clusters e.g. for species X how many clusters and what percentage of clusters
contain 1 gene, 2 genes, 3 genes ... >10 genes, as well as percentage surface protein estimates for each cluster size bin
'''

import string
from sys import argv
import re

script, filename, lookup = argv

groupstxt = open(filename)
lines = groupstxt.readlines()
print "Your groups cluster file: %r" % filename
surface = open(lookup)
print "Your surface protein lookup file: %r" % lookup

surfaceList = [line.strip() for line in surface.readlines()]

OFH=open("surface_prots_per_cluster_threshold.txt",'w')
OFH2=open("results_with_threshold.txt",'w')

######## read OrthologousGroups.txt file

OFH.write("ORG\tNUMSURFACEPROTS\tNUMGENESINCLUSTER\tPERCSURFACEPROTS\n")
orgTotalClusters={}
speciesDict={}
for line in lines:
	line = line.strip()
	line2 = line.split(":")
	clustername = line2[0]
	genes = line2[1]
	genelist = []
	genelist = genes.strip().split(" ")
	orgList=[]
	surfaceFlag={}
	for gene in genelist:
		genenameparts = gene.split("_")
		org = genenameparts[0]
		orgList.append(org)
		if gene.strip() in surfaceList:
			if org not in surfaceFlag:
				surfaceFlag[org]=1
			else:
				surfaceFlag[org]+=1
	for o in set(orgList):
		if o not in orgTotalClusters:
			orgTotalClusters[o]=1
		else:
			orgTotalClusters[o]+=1
		ocount=orgList.count(o)
		clustercount=orgList.count(o)
		# print num of surface proteins out of all genes in each cluster
		outlist=[]
		if o in surfaceFlag:
			outlist="\t".join([str(i) for i in [o, surfaceFlag[o], ocount, "{0:.2f}".format((float(surfaceFlag[o])/float(ocount))*100)]])
			OFH.write(outlist)
			OFH.write("\n") 
		else:
			outlist="\t".join([str(i) for i in [o, "0", ocount, "0"]])
			OFH.write(outlist) 
			OFH.write("\n") 
		if ocount > 10:
			ocount=">10"
		if o not in speciesDict:
			speciesDict[o]={}
		if ocount not in speciesDict[o]:
			speciesDict[o][ocount]={}
			speciesDict[o][ocount]["numClusters"]=1
			if o in surfaceFlag:
				#print (float(surfaceFlag[o])/float(clustercount)*100)
				if (float(surfaceFlag[o])/float(clustercount)*100)>10:
					speciesDict[o][ocount]["numSurfaceProtClusters"]=1
		else:
			speciesDict[o][ocount]["numClusters"]+=1
			if o in surfaceFlag:
				if (float(surfaceFlag[o])/float(clustercount)*100)>10:
					if "numSurfaceProtClusters" not in speciesDict[o][ocount]:
						speciesDict[o][ocount]["numSurfaceProtClusters"]=1
					else:
						if (float(surfaceFlag[o])/float(clustercount)*100)>10:
							speciesDict[o][ocount]["numSurfaceProtClusters"]+=1
	

OFH2.write("BINCOUNT\tSURFACEPROTCOUNT\tTOTCOUNT\tPERCSURFACEPROT\tPERCTOTCLUSTERS\tORG\n")
for species in speciesDict:
	#print species
	for c in speciesDict[species]:
		percTot=(float(speciesDict[species][c]["numClusters"])/float(orgTotalClusters[species]))*100
		perc2print="{0:.2f}".format(percTot)
		percSurface=(float(speciesDict[species][c]["numSurfaceProtClusters"])/float(speciesDict[species][c]["numClusters"]))*100
		percSurface2print="{0:.2f}".format(percSurface)
		outres="\t".join([str(i) for i in [c, speciesDict[species][c]["numSurfaceProtClusters"],speciesDict[species][c]["numClusters"], percSurface2print, perc2print, species]])
		OFH2.write(outres)
		OFH2.write("\n") 

groupstxt.close()
surface.close()

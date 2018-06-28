import re
import sys
import os
import uuid

ssuref="/data/refdb/silva/SSURef_119_SILVA_14_07_14_opt.arb"
lsuref="/data/refdb/silva/LSURef_119_SILVA_15_07_14_opt.arb"
#

if len(sys.argv) != 7:
	print "Usage python "+sys.argv[0]+" <kingdom(bac|euk)> <out.gff> <out.fasta> <genome.fasta> <summary.txt> <sina_classification.txt> "
	sys.exit()

rnammerKingdom=sys.argv[1]
rnammerGff=sys.argv[2]
rnammerFasta=sys.argv[3]
fasta=sys.argv[4]
rnammerSummary=sys.argv[5]
sinaClassification=sys.argv[6]

rnammer_cmd="/usr/global/blp/rnammer-1.2/rnammer -S "+rnammerKingdom+" -gff "+rnammerGff+" -multi -f "+rnammerFasta+"  "+fasta+" -m tsu,ssu,lsu -T /tmp/"+str(uuid.uuid4())
print rnammer_cmd
os.system(rnammer_cmd);

gffFile=open(rnammerGff)
gfflines=gffFile.readlines()
gffFile.close()

gfflines=filter(lambda x:not x.startswith("#"),gfflines)
nrna=len(gfflines)
ncontigs=len(set(map(lambda x:x.split("\t")[0],gfflines)))
types=set(map(lambda x:x.split("\t")[0],gfflines))

tmpFasta=open(rnammerFasta)
fastaLines=tmpFasta.readlines()
tmpFasta.close()
fastaLines=filter(lambda x:x.startswith(">"),fastaLines)
name2rnatype=dict((y[0][0],y[0][1]) for y in map(lambda x:re.findall(">(\S+) /molecule=(\S+) /score=(\S+)",x),fastaLines))
name2rnascore=dict((y[0][0],y[0][2]) for y in map(lambda x:re.findall(">(\S+) /molecule=(\S+) /score=(\S+)",x),fastaLines))

outFile=open(rnammerSummary,'w')
outFile.write("rnammer output Fasta             \t\t\t:%s\n"%rnammerFasta)
outFile.write("rnammer output GFF               \t\t\t:%s\n"%rnammerGff)
outFile.write("No. of RNA features found        \t\t\t:%d\n"%nrna)
outFile.write("No. of contigs with RNA features \t\t\t:%d\n"%ncontigs)
outFile.write("Type of RNA found                \t\t\t:%s\n"%(",".join(types)))

refs=[ssuref,lsuref]
outtxts=[]
for ref in enumerate(refs):
	outarb=fasta+"_sina."+str(ref[0])+".arb"
	outtxt=fasta+"_sina."+str(ref[0])+".out"
	outlog=fasta+"_sina."+str(ref[0])+".log"
	outtxts.append(outtxt)
	sina_cmd="/usr/global/blp/bin/sina -i "+rnammerFasta+" -o "+outarb+" --ptdb "+ref[1]+" --search --search-db "+ref[1]+" --lca-fields tax_slv > "+outlog+" 2> "+outtxt
	print sina_cmd
	os.system(sina_cmd)

sinaouttxt=fasta+"_sina.out"
print "cat "+(' '.join(outtxts))+" > "+sinaouttxt
os.system("cat "+(' '.join(outtxts))+" > "+sinaouttxt)

log=open(sinaouttxt)
lines=log.readlines()
log.close()

class sina_alignment:
	def __init__(self,id):
		self.id=id
		self.rnatype=name2rnatype[id]
		self.rnascore=name2rnascore[id]
		self.cl=""
		self.len=""
	def __str__(self):
		x=self.id.split("_")
		return "%s\t%s\t%s\t%s\t%s\t%s"%(x[1],self.id,self.len,self.rnatype,self.rnascore,self.cl)

classification=dict()
for l in lines:
	if re.search('identifier:',l):
		if l.startswith('iden'):
			id=re.findall(r'^identifier:(.*)\n',l)
		else:
			id=re.findall(r'^sequence_identifier:(.*)\n',l)
		id=id[0].strip()
	cl=re.findall(r'^lca_tax_slv:(.*)\n',l)	
	if len(cl)!=0 and not cl[0].strip().startswith("Unclass"):
		classification[id]=sina_alignment(id)
		classification[id].cl=cl[0].strip()
	nuc=re.findall(r'nuc:(.*)\n',l)
	if len(nuc)!=0 and id in classification:
		classification[id].len=nuc[0].strip()

outFile.write("Number of RNA classified         \t\t\t:%d\n"%(len(classification)))
outFile.write("SINA classification output       \t\t\t:\n")
outFile.write("#SampleName\tseqID\tseqLength\trnaType\trnammerScore\tSINA\n")

sc=open(sinaClassification,'w')
for k,v in classification.iteritems():
	sc.write(str(v)+"\n")
	outFile.write(str(v)+"\n")
sc.close()
outFile.close()





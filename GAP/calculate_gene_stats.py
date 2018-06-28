#!/home/vnkoparde/opt/Python-2.7.2/bin/python
import HTSeq
import sys,os,argparse,shlex,subprocess,re

#~ Version 1.0

parser=argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description=
"""This script calculate stats for assembly and its gene calls

Author: Vishal N Koparde, Ph. D. 
Created: 120716
Modified: 120716
Version: 1.0""",
 version="1.0")
parser.add_argument('--assemblyFasta',help='input assembly fasta file',dest='afasta',required=True,metavar="<assembly.fasta>")
parser.add_argument('--genesFasta',help='input genes (nuc) fasta file',dest='gfasta',required=True,metavar="<genes.fasta>")
parser.add_argument('--genesFaa',help='input genes (aa) fasta file',dest='genesFaa',required=True,metavar="<genes.faa>")
parser.add_argument('--genesGFF',help='genes in gff format',dest='gff',required=True,metavar="<genes.gff>")
parser.add_argument('-o',help='output stats in this file',dest='outfile',required=True,metavar='<out.stats>')
args=vars(parser.parse_args())

if (not os.path.exists(args['afasta'])):
	print args['afasta']+" file not found!"
	sys.exit()

if (not os.path.exists(args['gfasta'])):
	print args['gfasta']+" file not found!"
	sys.exit()


try:
	gapdir=os.environ['GAPDIR']
except KeyError:
	gapdir="/usr/global/blp/GenomeAnnotationPipeline"

contigs=dict((s.name,s) for s in HTSeq.FastaReader(args['afasta']))
ncontigs=len(contigs)
totallen=reduce(lambda x,y:x+y,list(len(x.seq) for x in contigs.values()))
gc=sum(map(lambda x:len(re.findall('[gGcC]',x)),list(x.seq for x in contigs.values())))*100.0/(totallen-sum(map(lambda x:len(re.findall('[nN]',x)),list(x.seq for x in contigs.values()))))
halflen=0.5*totallen
sortedcontigs=sorted(contigs.items(),key=lambda x:len(x[1].seq),reverse=True)
longestname=sortedcontigs[0][0]
longestlen=len(contigs[longestname].seq)
shortestname=sortedcontigs[-1][0]
shortestlen=len(contigs[shortestname].seq)
tlen=0
for i in enumerate(sortedcontigs):
	tlen+=len(i[1][1].seq)
	if tlen > halflen:
		n50len=len(i[1][1].seq)
		n50contigs=i[0]+1
		n50avgcontiglen=tlen*1.0/n50contigs
		break


genes=dict((s.name,s) for s in HTSeq.FastaReader(args['gfasta']))
ngenes=len(genes)
totalgenelen=reduce(lambda x,y:x+y,list(len(x.seq) for x in genes.values()))
ggc=sum(map(lambda x:len(re.findall('[gGcC]',x)),list(x.seq for x in genes.values())))*100.0/(totalgenelen-sum(map(lambda x:len(re.findall('[nN]',x)),list(x.seq for x in genes.values()))))
sortedgenes=sorted(genes.items(),key=lambda x:len(x[1].seq),reverse=True)
longestgenename=sortedgenes[0][0]
longestgenelen=len(genes[longestgenename].seq)
shortestgenename=sortedgenes[-1][0]
shortestgenelen=len(genes[shortestgenename].seq)


avggenelen=totalgenelen*1.0/ngenes

gffcomp=args['gff']+".complement"

os.system("/usr/global/blp/bin/samtools faidx "+args['afasta']);
os.system("/usr/global/blp/bin/complementBed -i "+args['gff']+" -g "+args['afasta']+".fai > "+gffcomp);

gffCompFile=open(gffcomp)
gffclines=gffCompFile.readlines()
gffCompFile.close()

seqname=""
igDistances=[]
for l in gffclines:
    cols=l.split("\t")
    if cols[0]!=seqname:
        seqname=cols[0]
        continue
    else:
        igDistances.append(int(cols[2])-int(cols[1]))

avgigdist=sum(igDistances)*1.0/len(igDistances)


tmpaFasta=args['afasta'].split("/")
tmpaFasta=tmpaFasta[-1]
tmpgFasta=args['gfasta'].split("/")
tmpgFasta=tmpgFasta[-1]

outfile=open(args['outfile'],'w')

gsizemb="%10.2f"%(totallen*1.0/1000000)
gsizemb=gsizemb.lstrip()
tsizemb="%10.2f"%(totalgenelen*1.0/1000000)
tsizemb=tsizemb.lstrip()
avggenelen="%10.2f"%(avggenelen)
avggenelen=avggenelen.lstrip()
gdensity="%10.2f"%(ngenes/(totallen*1.0/1000000))
gdensity=gdensity.lstrip()
gc="%10.2f"%(gc)
gc=gc.lstrip()
ggc="%10.2f"%(ggc)
ggc=ggc.lstrip()
avgigdist="%10.2f"%(avgigdist)
avgigdist=avgigdist.lstrip()
ttog=totalgenelen*1.0/totallen
ttog="%10.5f"%(ttog)
ttog=ttog.lstrip()

outfile.write("Assembly file name         \t\t\t:%s\n"%(args['afasta']))
outfile.write("Genome Size(Kb)            \t\t\t:%d\n"%(totallen*1.0/1000))
outfile.write("Genome Size(Mb)            \t\t\t:%s\n"%(gsizemb))
outfile.write("Percent Genomic GC         \t\t\t:%s\n"%(gc))
outfile.write("Percent Transcriptomic GC  \t\t\t:%s\n"%(ggc))
outfile.write("Total number of contigs    \t\t\t:%d\n"%(ncontigs))
outfile.write("Longest contig length      \t\t\t:%d\n"%(longestlen))
outfile.write("Longest contig ID          \t\t\t:%s\n"%(longestname))
outfile.write("Shortest contig length     \t\t\t:%d\n"%(shortestlen))
outfile.write("Shortest contig ID         \t\t\t:%s\n"%(shortestname))
outfile.write("N50 length                 \t\t\t:%s\n"%(n50len))
outfile.write("N50 No. of contigs         \t\t\t:%s\n"%(n50contigs))
outfile.write("N50 avg. contig length     \t\t\t:%s\n"%(n50avgcontiglen))
outfile.write("Genes (nuc) file           \t\t\t:%s\n"%(args['gfasta']))
outfile.write("Genes (aa) file            \t\t\t:%s\n"%(args['genesFaa']))
outfile.write("Transcriptome Size(Kb)     \t\t\t:%d\n"%(totalgenelen*1.0/1000))
outfile.write("Transcriptome Size(Mb)     \t\t\t:%s\n"%(tsizemb))
outfile.write("Total number of genes      \t\t\t:%d\n"%(ngenes))
outfile.write("Average gene length        \t\t\t:%s\n"%(avggenelen))
outfile.write("Longest gene length        \t\t\t:%d\n"%(longestgenelen))
outfile.write("Longest gene ID            \t\t\t:%s\n"%(longestgenename))
outfile.write("Shortest gene length       \t\t\t:%d\n"%(shortestgenelen))
outfile.write("Shortest gene ID           \t\t\t:%s\n"%(shortestgenename))
outfile.write("Gene Density (per Mb)      \t\t\t:%s\n"%(gdensity))
outfile.write("Avg Intergenic Length      \t\t\t:%s\n"%(avgigdist))
outfile.write("Transcriptome/Genome ratio \t\t\t:%s\n"%(ttog))

outfile.close()
sys.exit()

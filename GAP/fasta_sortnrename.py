import sys,os,argparse,copy
import HTSeq


parser=argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description=
"""This script sorts sequences by length in decending order and renames them 

Author: Vishal N Koparde, Ph. D.
Created: 150202
Modified: 150202
Version: 1.0""",
 version="1.0")

parser.add_argument('-i',help='input fasta file',dest='infile',required=True,metavar="<in.fasta>")
parser.add_argument('-o',help='output fasta file',dest='outfile',required=True,metavar="<out.fasta>")
parser.add_argument('-p',help='seqid prefix',dest='p',required=True,metavar="<TCO> for TCO_00001,TCO_00002,etc.")
args=vars(parser.parse_args())


if __name__=="__main__":
	sequences = dict( (s.name, s) for s in HTSeq.FastaReader(args['infile']))
	f=open(args['outfile'],'w')
	f2=open(args['outfile']+".lookup",'w')
	ctr=0
	for id in sorted(sequences.items(), key=lambda x: len(x[1].seq),reverse=True):
		ctr+=1
		newname="%s_%05d"%(args['p'],ctr)
		f2.write("%s\t%s\t%d\n"%(id[1].name,newname,len(id[1].seq)))
		id[1].name=newname
		id[1].write_to_fasta_file(f)
	f.close()
	f2.close()

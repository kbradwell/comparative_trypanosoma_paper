#!/bin/bash
if [ $# -ne 2 ];then
echo "usage : gce_tcruzi <path to assembly fasta folder> <organism short name>"
cat << EOF
output files :
<shortName>_tblastn.blast.gz				gzipped blast output
<shortName>_tblastn.bp.gz				gzipped parsed blast output
<shortName>_tblastn_blaststats.txt			blast stats at any,25,50,75,90,99% alignment
<shortName>_tblastn.blast.gz.duplication		list level of duplications of each of the 2217 orthologs
	columns --> 	ortholog	Nhits	Max.percentAlignment	Max.Bitscore
<shortName>_tblastn.blast.gz.duplication.summary	duplication summarized
EOF
exit;
fi
cwd=`pwd`
project=$2
#project=${project}_assemblies
cd $1
mkdir gce_tcruzi
cd gce_tcruzi
for f in `ls *.gz`;do gzip -d $f;done
for f in `ls ../*.fasta`;do
if [ ! -f $f ]; then continue; fi
g=`echo $f|awk -F "/" '{print $2}'`
ln -s $f ${g%.*}_all.fasta
done
perl GenomeCompletionEvalution_without_qsub.pl $project TcruziCL_ESM_proteins.fa tblastn `pwd`
gzip -f *.bp
gzip -f *.blast
rm -f gce_tcruzi.txt
for f in `ls *.blast.gz`;do
if [ ! -f ${f}.duplication ]; then perl /home/vnkoparde/scripts/gce_duplication.pl <(zcat $f) > ${f}.duplication; fi
cat ${f}.duplication|awk '{print $2}'|sort|uniq -c|sort -k2n > ${f}.duplication.summary
done

for f in `ls ../*.fasta`;do
if [ ! -f $f ]; then continue; fi
g=`echo $f|awk -F "/" '{print $2}'`
tail -n1 ${g%.*}_${project}_tblastn_blaststats.txt >> gce_tcruzi.txt.tmp
done
cat gce_tcruzi.txt.tmp | awk '{sum=$4+$6+$8+$10+$12;printf "%s\t%d\n",$0,sum}'|sort -k14,14nr|rev|cut -f2-|rev > gce_tcruzi.txt
#sort -k12,12nr -k10,10nr -k8,8nr -k6,6nr -k4,4nr gce_tcruzi.txt.tmp > gce_tcruzi.txt
rm -f gce_tcruzi.txt.tmp
cd $cwd

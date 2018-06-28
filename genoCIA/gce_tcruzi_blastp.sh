#!/bin/bash
if [ $# -ne 2 ];then
echo "usage : gce_tcruzi_blastp <folder containing list of proteins as fasta file> <organism short name>"
exit;
fi
cwd=`pwd`
project=$2
project=${project}_proteins
cd $1
mkdir gce_tcruzi_blastp
cd gce_tcruzi_blastp
for f in `ls *.gz`;do gzip -d $f;done
for f in `ls ../*.fasta`;do
if [ ! -f $f ]; then continue; fi
g=`echo $f|awk -F "/" '{print $2}'`
ln -s $f ${g%.*}_all.fasta
done
perl GenomeCompletionEvalution_without_qsub_blastp.pl $project /gpfs_fs/data2/tol/scripts/GenomeCompletionByReferenceList/genelist_files/TcruziCL_ESM_proteins.fa blastp `pwd`
gzip *.bp
gzip *.blast
rm -f gce_tcruzi_blastp.txt
for f in `ls ../*.fasta`;do
if [ ! -f $f ]; then continue; fi
g=`echo $f|awk -F "/" '{print $2}'`
tail -n1 ${g%.*}_${project}_blastp_blaststats.txt >> gce_tcruzi_blastp.txt.tmp
done
cat gce_tcruzi_blastp.txt.tmp | awk '{sum=$4+$6+$8+$10+$12;printf "%s\t%d\n",$0,sum}'|sort -k14,14nr|rev|cut -f2-|rev > gce_tcruzi_blastp.txt
#sort -k12,12nr -k10,10nr -k8,8nr -k6,6nr -k4,4nr gce_tcruzi_blastp.txt.tmp > gce_tcruzi_blastp.txt
rm -f gce_tcruzi_blastp.txt.tmp
cd $cwd

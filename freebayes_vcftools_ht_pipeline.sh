#! /bin/bash

gffprefix=$1
assembly=$2

touch ${gffprefix}_ht_outstats

cat ../${gffprefix}.gff |cut -f 1,4,5 > ${gffprefix}.bed

samtools view -H reads2assembly.bam | sed -e 's/PG:clc_mapper/PG:clc_mapper	SM:tccl/' | samtools reheader - reads2assembly.bam > reads2assembly.reheadered.bam

samtools index reads2assembly.reheadered.bam

export PERL5LIB=/usr/global/blp/vcftools_0.1.9/perl/
echo "Total number of positions" >> ${gffprefix}_ht_outstats
cat ${gffprefix}.bed | awk '{ sum+=($3-$2)} END {print sum}' >> ${gffprefix}_ht_outstats

/usr/global/blp/freebayes/bin/freebayes --fasta-reference $assembly --targets ${gffprefix}.bed --ploidy 2 --bam reads2assembly.reheadered.bam --vcf ${gffprefix}.vcf &> freebayes.log

/usr/global/blp/bin/vcftools --vcf ${gffprefix}.vcf --remove-indels --recode --recode-INFO-all --out SNPs_only &> vcftools_snponly.log

echo "Total number of snp positions" >> ${gffprefix}_ht_outstats
grep "out of a possible" SNPs_only.log >> ${gffprefix}_ht_outstats

cat ../${gffprefix}.gff |cut -f 1,4,5,9 > ${gffprefix}.annotation

vcf-sort ${gffprefix}.annotation > ${gffprefix}.annotation.sorted

bgzip ${gffprefix}.annotation.sorted

tabix -s 1 -b 2 -e 3 ${gffprefix}.annotation.sorted.gz

cat SNPs_only.recode.vcf | vcf-annotate -a ${gffprefix}.annotation.sorted.gz -d key=INFO,ID=ANN,Number=1,Type=Integer,Description='My custom annotation' -c CHROM,FROM,TO,INFO/ANN > ${gffprefix}_annotated.vcf 
echo "Total number of HT SCO genes" >> ${gffprefix}_ht_outstats
grep "TYPE=" ${gffprefix}_annotated.vcf|cut -f 42 -d ";"| cut -f 1|sort|uniq|wc >> ${gffprefix}_ht_outstats

echo "Number of bi-tri- and tetra- alleles" >> ${gffprefix}_ht_outstats
/usr/global/blp/bin/vcftools --vcf ${gffprefix}_annotated.vcf --min-alleles 2 --max-alleles 2 > bi-allele_count
/usr/global/blp/bin/vcftools --vcf ${gffprefix}_annotated.vcf --min-alleles 3 --max-alleles 3 > tri-allele_count
/usr/global/blp/bin/vcftools --vcf ${gffprefix}_annotated.vcf --min-alleles 4 --max-alleles 4 > tetra-allele_count
for f in *allele_count;do echo $f >>  ${gffprefix}_ht_outstats && grep "out of a possible" $f >>  ${gffprefix}_ht_outstats;done

grep "TYPE=" ${gffprefix}_annotated.vcf|cut -f 42 -d ";"| cut -f 1|sort|uniq -c > ${gffprefix}_snp-counts-per-gene

grep "TYPE=" ${gffprefix}_annotated.vcf|cut -f 1 -d ";"| cut -f 2 -d "=" > ${gffprefix}_reference-allele-proportions
echo "Total number of SCO genes" >> ${gffprefix}_ht_outstats
cat ${gffprefix}.annotation|cut -f 4|sort|uniq|wc >> ${gffprefix}_ht_outstats

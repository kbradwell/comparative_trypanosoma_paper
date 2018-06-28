#!/usr/bin/perl
#Developer: Vishal Koparde, Ph.D.
#Created: 120723
#Modified:120723
#Version 1.0

use strict;
use warnings;

my ($gapdir);
if (exists $ENV{'GAPDIR'}) {
$gapdir=$ENV{'GAPDIR'};
} else {
$gapdir="/usr/global/blp/GenomeAnnotationPipeline";
}

my $gmsnFile = shift;
my $gffFile = shift;
my $shortName = shift;

my $ncontig = `grep -c FASTA $gmsnFile`;
chomp $ncontig;
$ncontig = int($ncontig);
my $ncuts = $ncontig-1;

my $cmd="/usr/bin/csplit --digits=5 --prefix=$shortName $gmsnFile /^FASTA/ {$ncuts}";
system($cmd);

my $gff = sprintf "%s%05d.gff",$shortName,1;
$cmd="cat ${shortName}00000 ${shortName}00001 > ${shortName}00000-00001";
system($cmd);
my $seqname=`grep "FASTA definition line" ${shortName}00000-00001|awk '{print \$NF}'`;
chomp $seqname;
$cmd="perl ${gapdir}/cnv_genemark2gff.pl";
$cmd.=" -i ${shortName}00000-00001";
$cmd.=" -o $gff";
$cmd.=" --seqname ${seqname}";
$cmd.=" --geneStart 1 ";
$cmd.=" --genePrefix ${shortName}_gene";
system($cmd);
$cmd="cp $gff $gffFile";
system($cmd);

my $gffEntryCount = `wc -l $gff|awk '{print \$1}'`;
chomp $gffEntryCount;
my ($currSplit,$currSplitPlus);
for (my $contig=2;$contig<=$ncontig;$contig++) {
$gff = sprintf "%s%05d.gff",$shortName,$contig;
$currSplit = sprintf "%s%05d",$shortName,$contig;
$currSplitPlus = sprintf "%s00000-%05d",$shortName,$contig;
$cmd="cat ${shortName}00000 ${currSplit} > ${currSplitPlus}";
#print $cmd."\n";
system($cmd);
$seqname=`grep "FASTA definition line" ${currSplitPlus}|awk '{print \$NF}'`;
chomp $seqname;
my $geneStart = int($gffEntryCount) + 1;
$cmd="perl ${gapdir}/cnv_genemark2gff.pl";
$cmd.=" -i ${currSplitPlus}";
$cmd.=" -o $gff";
$cmd.=" --seqname ${seqname}";
$cmd.=" --geneStart ${geneStart}";
$cmd.=" --genePrefix ${shortName}_gene";
print $cmd."\n";
system($cmd);
my $tmp = `wc -l $gff|awk '{print \$1}'`;
chomp $tmp;
$gffEntryCount += int($tmp);
$cmd="cat $gff >> $gffFile";
system($cmd);
}

unlink("${shortName}00000");
for (my $contig=1;$contig<=$ncontig;$contig++) {
$gff = sprintf "%s%05d.gff",$shortName,$contig;
$currSplit = sprintf "%s%05d",$shortName,$contig;
$currSplitPlus = sprintf "%s00000-%05d",$shortName,$contig;
unlink($gff);
unlink($currSplit);
unlink($currSplitPlus);
}

exit;

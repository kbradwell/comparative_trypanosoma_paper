#!/usr/bin/perl

#callGenesWithGlimmer
#Developer: Vishal Koparde, Ph.D.
#Created: 120710
#Modified:150203


use warnings;
use strict;
use lib qw(/usr/global/blp/perllib);
use Qsub;
use Bio::Tools::GFF;
use Bio::Tools::Glimmer;
use Getopt::Long;

Getopt::Long::Configure("bundling");

my $gapdir;
if (exists $ENV{'GAPDIR'}) {
$gapdir=$ENV{'GAPDIR'};
} else {
$gapdir="/usr/global/blp/GenomeAnnotationPipeline";
}

sub usage();

usage() if (scalar @ARGV==0);

my ($fastaSorted,$predictFile,$workingDir,$genesFasta,$genesGff,$prefix,$force,$genesFaa,$help);

GetOptions ( 'sortedFasta=s' => \$fastaSorted,
             'prefix=s' => \$prefix,
             'wd=s' => \$workingDir,
             'force' => \$force,
             'help|?' => \$help);

usage() if ($help);

die "Fasta file $fastaSorted not found!" if (! -e $fastaSorted);
die "prefix not defined!" unless (defined $prefix);
if (!defined $workingDir) {
    $workingDir=`pwd`;
    chomp $workingDir;
}
$genesFasta = $prefix."_genes.fasta";
$genesFaa = $prefix."_genes.faa";
$genesGff = $prefix."_genes.gff";
$predictFile=$prefix.".predict";

my $Glimmer_cmd1="${gapdir}/g3-from-scratch.sh $fastaSorted $prefix";
my $jname="GL1_".$prefix;
my $jnameout="GL1_".$prefix.".tmp.out";
my $job_gl1=new Qsub(name=>$jname,wd=>$workingDir,outfile=>$jnameout,cmd=>$Glimmer_cmd1);
$job_gl1->submit();

# This should generate the following files:
# $prefix.detail
# $prefix.icm
# $prefix.longorfs
# $prefix.predict
# $prefix.train

my $Glimmer_cmd2="perl ${gapdir}/predict2fastagff.pl -c $fastaSorted -l $predictFile -o $genesFasta -g $genesGff -s $prefix -p $genesFaa";
$jname="GL2_".$prefix;
$jnameout="GL2_".$prefix.".tmp.out";
my $job_gl2=new Qsub(name=>$jname,wd=>$workingDir,outfile=>$jnameout,cmd=>$Glimmer_cmd2,waitfor=>"$job_gl1->{jobid}");
$job_gl2->submit();
$job_gl2->waitForCompletion();


exit;

sub usage() 
{
print <<EOF;
Call Genes Using Glimmer

Developer : Vishal N Koparde, Ph. D.
Created   : 120710
Modified  : 150203
Version   : 1.0

options:
--sortedFasta      Fasta file sorted by size of sequences (large-to-small)
--prefix           Sequence ID prefix for the genes fasta file, generally organism short name(required)
--wd               Working Directory (Default="pwd")
--force            Overwrite all files
--help             This help
EOF
exit 1;
}

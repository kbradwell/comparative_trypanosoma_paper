#!/usr/bin/perl
#Developer: Vishal Koparde, Ph.D.
#Created: 120710
#Modified:150202
#Version 1.0 - Called by GAP

use strict;
use warnings;
use Getopt::Long;
use lib qw(/usr/global/blp/perllib);
use Qsub;
use Bio::Tools::GFF;
use Bio::Tools::Genemark;

my $gapdir;
if (exists $ENV{'GAPDIR'}) {
$gapdir=$ENV{'GAPDIR'};
} else {
$gapdir="/usr/global/blp/GenomeAnnotationPipeline";
}

sub usage();

usage() if (scalar @ARGV==0);

my ($fastaSorted,$gmsnGenesFile,$workingDir,$genesFasta,$genesGff,$prefix,$force,$genesFaa,$help);

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
$gmsnGenesFile = $prefix.".gmsn_genes"; 
$genesFasta = $prefix."_genes.fasta";
$genesFaa = $prefix."_genes.faa";
$genesGff = $prefix."_genes.gff";


die "$genesFasta exists, please use --force to overwrite!\n" if (-e $genesFasta && ! $force);

my $GeneMark_cmd1="/usr/global/blp/GeneMark/GeneMarkS/genemark_suite_linux_64/gmsuite/gmsn.pl -euk $fastaSorted -output $gmsnGenesFile";
my $jname = "GM1_".$prefix;
my $jnameout = $jname.".tmp.out";
my $job_gm1=new Qsub(name=>$jname,wd=>$workingDir,outfile=>$jnameout,cmd=>$GeneMark_cmd1);
$job_gm1->submit();

$jname = "GM2_".$prefix;
$jnameout = $jname.".tmp.out";
my $GeneMark_cmd2="perl ${gapdir}/gmsn2fastagff.pl -c $fastaSorted -l $gmsnGenesFile -o $genesFasta -g $genesGff -s $prefix -p $genesFaa";
my $job_gm2=new Qsub(name=>$jname,wd=>$workingDir,outfile=>$jnameout,cmd=>$GeneMark_cmd2,waitfor=>"$job_gm1->{jobid}");
$job_gm2->submit();

$job_gm2->waitForCompletion();

exit;

sub usage() 
{
print <<EOF;
Call Genes Using GeneMarkS

Developer : Vishal N Koparde, Ph. D.
Created   : 120710
Modified  : 150202
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

#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my $project = '';
my $gmsnGenes = '';
my $rpsCog = '';
my $rpsCogCdd = '';
my $rpsPfam = '';
my $rpsPfamCdd = '';
my $gmsnGenesGFF3 = '';
my $rpsCogGFF3 = '';
my $rpsPfamGFF3 = '';
my $bxTDDB = '';
my $bxTddbGFF3 = '';
my $asgardNr = '';
my $asgardNrGFF3 = '';
my $asgardGenesGFF3 = '';
my $asgardGenes2ContigsGFF3 = '';
my $tRna = '';
my $tRnaGFF3 = '';

if ( !GetOptions
  (
	'project=s'	=> \$project,
	'genes=s'	=> \$gmsnGenes,
	'rpscog=s'	=> \$rpsCog,
	'rpspfam=s'	=> \$rpsPfam,
	'bxttb=s'	=> \$bxTTDB,
	'nr=s'		=> \$asgardNr,
	'aggff3=s'	=> \$asgardGenesGFF3,
	'trna=s'	=> \$tRna
  )
)  { exit(1); };

if($project ne ''){
	print STDERR "No project name\n$usage\n";
}

if(-r $gmsnGenes and -r $rpsCog.cdd and -r $rpsPfam.cdd) {
	`gm2gff3.pl $project $gmsnGenes $rpsCog.cdd $rpsPfam.cdd > $projectGenes.gff3`;
}
else {
	print STDERR "No output for GeneMarkS\n";
}

if(-r $rpsCog){
	`blast2gff3.pl < $rpsCog > $rpsCog.gff3`;
}
else{
	print STDERR "No output for rpsCog\n";
}

if(-r $rpsPfam){
	`blast2gff3.pl < $rpsPfam > $rpsPfam.gff3`;
}
else{
	print STDERR "No output for rpsPfam\n";
}

if(-r $bxTTDB){
	`blast2gff3_dr.pl < $bxTTDB`;
}
else{
        print STDERR "No output for TriTrypDB\n";
}

if(-r $asgardNr){
	`asg2gff3.pl < $asgardNr > $asgardNr.gff3`;
}
else{
        print STDERR "No output for NR\n";
}

if(-r $asgardGenesGFF3 $gmsnGenes."gff3"){
	`genes2contigsGFF.pl $asgardGenesGFF3 $gmsnGenes.gff3 > $project.asgard_genes2contigs.gff3`;
}
else{
        print STDERR "No output for ASGARD annotation of genes\n";
}

if(-r $tRna){
`trnascan2gff3.pl $tRna > $tRna.gff3`;
}
else{
        print STDERR "No output for tRNAscan\n";
}


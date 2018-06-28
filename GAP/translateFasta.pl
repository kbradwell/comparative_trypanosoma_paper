#!/usr/bin/perl
#Developer: Vishal Koparde, Ph.D.
#Created: 120716
#Modified:120716
#Version 1.0

use Bio::SeqIO;
use Getopt::Long;

my ($infile,$outfile);

GetOptions(
'fna:s' => \$infile,
'faa:s' => \$outfile
);

die "Input DNA fasta filename not supplied!\n" unless (defined $infile);
die "No such file : $infile \n" unless (-e $infile);
die "Output Protein fasta filename not supplied!\n" unless (defined $outfile);

my $format = 'fasta';

my $seqin = Bio::SeqIO->new( -format => $format, -file => "$infile");
my $seqout = Bio::SeqIO->new( -format => $format, -file => ">$outfile" );

while( (my $seq = $seqin->next_seq()) ) {
my $pseq = $seq->translate();
$seqout->write_seq($pseq);
}

$seqin->close();
$seqout->close();
exit;
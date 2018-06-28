#!/usr/bin/perl -w

#  Extract genes predicted by Glimmer as .predict coordinates to .fasta sequences from contigs in .fasta format
#
#Developer: Vishal Koparde, Ph.D.
#Created: 150203
#Modified:150203
#

use strict;
use Getopt::Long;
use Bio::SeqIO;
use Bio::Seq;

my $usage = 

"
$0 -- extract genes predicted by Glimmer as .predict coordinates to .fasta and .gff from contigs in .fasta format

Usage:

$0 -c <contigs.fasta> -l <coordinates.predict> -o <genes.fasta> -g <genes.gff> -s <prefix or shortname> -p <genes.faa>

";

if ( $#ARGV == -1 ) { print $usage; exit(0); }

my $contigs = '';
my $coordinates = '';
my $genes;
my $gff;
my $prefix;
my $genesFaa;

if ( !GetOptions
	(
		'c=s' => \$contigs,
		'l=s' => \$coordinates,
		'o=s' => \$genes,
		'g=s' => \$gff,
		's=s' => \$prefix,
		'p=s' => \$genesFaa
	)
   ){
	print $usage;
	exit 1;
}
else{
	print "
INPUT
contigs:       $contigs
coordinates:   $coordinates
prefix:        $prefix
OUTPUT
genes(nuc):    $genes
genes(pro):    $genesFaa
gff:           $gff
";
}

die "File $contigs Doesn't Exist!\n" unless (-e $contigs);
die "File $coordinates Doesn't Exist!\n" unless (-e $coordinates);
die "Prefix is required!\n" unless defined $prefix;
die "Output fasta file name is required!\n" unless defined $genes;
die "Output gff file name is required!\n" unless defined $gff;
die "Output protein fasta file name is required!\n" unless defined $genesFaa;

#  LOAD ALL GENOME SEQUENCES

my $in  = Bio::SeqIO->new(-file => "$contigs", -format => 'Fasta');
my %seq;
my $sequence;
while ( $sequence = $in->next_seq() ) {
            $seq{$sequence->id}=$sequence->seq;
}
$in->close();

open LST, $coordinates;

my $out = Bio::SeqIO->new(-file => ">$genes", -format => 'Fasta');
my $out2 = Bio::SeqIO->new(-file=> ">$genesFaa", -format => 'Fasta');
open GFF, ">$gff";

my $id_current = '';
my $ctr=0;
while(my $string = <LST>){
	chomp $string;
	if($string =~ /^>(\S+)/){
	    $id_current = $1;
            next;
            #print $id_current."\n";exit;
	}
	if($string =~ /^orf/){
            $ctr++;
            my @string = split /\s+/, $string;
            my $start = $string[1];
            my $end = $string[2];
            if ($start > $end) {
                my $tmp=$start;
                $start=$end;
                $end=$tmp;
            }
	    my $offset = $start - 1;
            my $length = $end - $start + 1;
            my $strand = "+";
            $strand = "-" if $string[3] < 0;
            my $gene = substr $seq{$id_current}, $offset, $length;
	    my $geneid = sprintf "%s_gene_%05d",$prefix,$ctr;
            $sequence = Bio::Seq->new( -display_id => $geneid,-desc => "$id_current:$start-$end",-seq => $gene);
            $sequence = $sequence->revcom() if ($strand eq "-");
	    my $prot_sequence = $sequence->translate;
	    $out->write_seq($sequence);
	    $out2->write_seq($prot_sequence);
	    print GFF "$id_current\tGlimmer\tgene\t$start\t$end\t.\t$strand\t.\t$geneid\n";
	}
}

close(GFF);

close(LST);




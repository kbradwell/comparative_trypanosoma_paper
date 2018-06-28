#!/usr/bin/perl -w

#  Extract genes predicted by GeneMark as .lst coordinates to .fasta sequences from contigs in .fasta format
#
#  Sergey S. Pintus
#
#  Version 10-05-2012
#
#

use strict;
use Getopt::Long;
use Bio::SeqIO;
use Bio::Seq;

my $usage = 

"
$0 -- extract genes predicted by GeneMark as .lst coordinates to .fasta sequences from contigs in .fasta format

Usage:

$0 -c <contigs.fasta> -l <coordinates.lst> -o <genes.fasta>

";

if ( $#ARGV == -1 ) { print $usage; exit(0); }

my $contigs = '';
my $coordinates = '';
my $genes = '';

if ( !GetOptions
	(
		'c=s' => \$contigs,
		'l=s' => \$coordinates,
		'o=s' => \$genes
	)
   ){
	print $usage;
	exit 1;
}
else{
	print "
contigs:	$contigs
coordinates:	$coordinates
genes:		$genes

";
}

die "File $contigs Doesn't Exist!\n" unless (-e $contigs);
die "File $coordinates Doesn't Exist!\n" unless (-e $coordinates);

my $in  = Bio::SeqIO->new(-file => "$contigs", -format => 'Fasta');
my %seq;
my $sequence;
while ( $sequence = $in->next_seq() ) {
            $seq{$sequence->id}=$sequence->seq;
}

open LST, $coordinates;

my $out = Bio::SeqIO->new(-file => ">$genes", -format => 'Fasta');

my $id_current = '';
while(my $string = <LST>){
	chomp $string;
	if($string =~ /FASTA definition line\: ([^\s]+)/){
		$id_current = $1;
	}
	if($string =~ /\s+[\+\-]\s+/){
                $string =~ s/[\<\>]/ /g;
                $string = " ".$string;
                my @string = split /\s+/, $string;
                my $start = $string[3];
                my $end = $string[4];
				my $offset = $start - 1;
                my $length = $end - $start + 1;
                my $strand = $string[2];
        		my $gene = substr $seq{$id_current}, $offset, $length;
                $sequence = Bio::Seq->new( -display_id => "$id_current:$start-$end",-seq => $gene);
				$sequence = $sequence->revcom() if ($strand eq "-");
				$out->write_seq($sequence);
	}
}

close(LST);

exit;

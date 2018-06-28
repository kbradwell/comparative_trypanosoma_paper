
#! /usr/bin/perl -W

use lib("/home/nsheth/scripts");
use FastaReader;
use strict;


my $debug = 1;

my $fastafile=shift;

if(!defined($fastafile)) { 
    die "USAGE: perl $0 fastafile\n";
}


my ($totalseq,$totalseqlength);
my %seqHash;
my $fr = FastaReader->new();
$fr->init_file($fastafile);
my ($defn,$seq);
my ($size,$seqid);
my ($minsize,$maxsize,$avgsize,$gctotal,$gcseqcount,$gcperc);
my @contigsize;

$totalseq = $totalseqlength  = $gctotal = 0;
$maxsize = $avgsize = 0; 
$minsize = 1000000000; # Very Big Number ...
while(($defn,$seq) = $fr->raw_next()) { 
    ($seqid) = $$defn =~ />(.*)/;
	$$seq =~ s/\n//g;
#	$$seq =~ s/\r//g;
#	$$seq =~ s/\s+//g;

	$size = length($$seq);
	#print "$size\t$$defn\n";
	#if($size >= 100) { 
        if($size >= 1) { 
		$gcseqcount = calculateGCCount($$seq);
		$gctotal += $gcseqcount;
		$totalseq++;
		$totalseqlength += $size;
		if($minsize > $size) { 
		
			$minsize = $size;
			if($size == 25 ) { 
				print "$$defn\n$$seq\n";
			}
		}

		if($maxsize < $size) { 
			$maxsize = $size;
		}
		push @contigsize,$size;
	}	
}


my ($totalsizen50,$noofn50,$n50size,$midpoint,$avgn50);
$totalsizen50 = 0;
$midpoint = int($totalseqlength / 2);
my $i = $noofn50 = 0;
foreach $size (sort {$b <=> $a } @contigsize) { 
    $totalsizen50 += $size;
    if($totalsizen50 > $midpoint) { 
	$n50size = $size;
	last;
    }
    $i++;
}

$noofn50 = $i + 1;
if($noofn50 ==0) { 
    $noofn50 = 1;
}

$avgn50 = $totalsizen50 / $noofn50;


$avgsize  = int($totalseqlength / $totalseq);
$gcperc = ($gctotal /  $totalseqlength) * 100;
print "No of Sequences = $totalseq\n";
print "Total Sequence Length = $totalseqlength\n";
print "GC Total = $gctotal\n";
print "GC Perc = $gcperc %\n";
print "Avg Contig Size = $avgsize\n";
print "Min Contig Size = $minsize\n";
print "Max Contig Size = $maxsize\n";

print "Total N50 contigs = $noofn50\n";
print "N50ContigSize = $n50size\n";
print "N50Avg = $avgn50\n";

exit;


sub calculateGCCount { 
	my ($seqref) = @_;
	my $gccount = 0;
	my @gcs;
	@gcs = $seqref =~ /[GC]/ig;
	$gccount = @gcs;
	return $gccount;
}

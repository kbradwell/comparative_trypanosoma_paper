#!/usr/bin/perl -w

if ($ARGV[0] eq "") {print "\ngff2fa - extract reference sequence fragments which coordinates are indicated in gff3 file\n\nUsage $0 coordinates.gff3 reference.fasta > output.fasta\n\n";}

$gff3 = $ARGV[0];
$fasta = $ARGV[1];

open FASTA, $fasta;
while(<FASTA>){
	chomp;
	if(/^\>(\S+)/){
		$id = $1;
#		print $id."\n";
	}
	else{
		$ref{$id} .= $_;
	}
}
close FASTA;

#foreach $id (sort keys %ref){
#	print ">".$id."\n".$ref{$id}."\n";
#}

open GFF3, $gff3;
while(<GFF3>){
	chomp;
	if(/\#\#/){
		next;
	}
	@gff3 = split /\t+/;
	$id =  $gff3[0];
	$start = $gff3[3];
	$end = $gff3[4];
	$offset = $start - 1;
	$length = $end - $offset;
	$fragment = substr $ref{$id}, $offset, $length;
	$descr = $gff3[8];
	$newid = $id.':'.$start.'-'.$end.'|'.$descr;
	print '>'.$newid."\n".$fragment."\n";
}
close GFF3;


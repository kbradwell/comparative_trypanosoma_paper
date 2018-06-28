#!/usr/bin/perl -w

$cutoff = $ARGV[0];

if($ARGV[0] eq ""){
	print "\n$0: remove short sequences (eg. contigs) from fasta file\n\nUsage: $0 cutoff_length < input_fasta > output_fasta\n\nSergey S. Pintus sspintus\@vcu.edu\n\n";
	exit;
}

$write = "no";
while($string = <STDIN>){
	chomp $string;
	if($string =~ /\>(Contig|Scaffold)\d+\s+(\d+)\s+\d+/){
		$length = $2;
		if($length >= $cutoff){
			$write = "yes";
		}
		else{
			$write = "no";
		}
	}
	if($write eq "yes"){
		print $string."\n";
	}
}


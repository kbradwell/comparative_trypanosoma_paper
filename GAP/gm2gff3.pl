#!/usr/bin/perl -w

use Switch;

$project = $ARGV[0];
$gmout = $ARGV[1];

#open COG, "$project.rpsCog.cdd";
#open PFAM, "$project.rpsPfam.cdd";

open COG, $ARGV[2];
open PFAM, $ARG[3];

$i = 1;
while (defined COG and defined PFAM and $cog = <COG> and $pfam=<PFAM>){
	chomp $cog;
	switch (length $i){
		case 1	{$id = "000".$i;}
		case 2	{$id = "00".$i;}
		case 3	{$id = "0".$i;}
		else	{$id = $i;}
	}
	$id = $project.'_'.$id;

	chomp $cog;
	chomp $pfam;
	@cog = split /\s+/, $cog;
	@pfam = split /\s+/, $pfam;
	$cogID{$id} = $cog[1];
	$pfID{$id} = $pfam[1];
	$pfID{$id} =~ s/pfam/PF/;
	$i++;
}
close PFAM;
close COG;

print "##gff-version 3\n";

open GMOUT, $gmout;
while ($string = <GMOUT>){
	chomp $string;
	if ($string =~ /^FASTA definition line\:\s+([^\s]+)/){
		$seqid1 = $1;
		next;
	}
	if ($string =~ /\s[+-]\s/){
		$string = ' '.$string;
		@string = split /\s+/, $string;
#		    4        -       39580       39915          336        1  
		$id9 = $string[1];
		switch (length $id9){
			case 1	{$id9 = "000".$id9;}
			case 2	{$id9 = "00".$id9;}
			case 3	{$id9 = "0".$id9;}
		}
		if (defined $id9){
			$id9 = $project.'_'.$id9;
		}
		$id9 = "ID=".$id9.";Dbxref=COG:".$cogID{$id9}.",Pfam:".$pfID{$id9};
		$strand7 = $string[2];
		$start4 = $string[3];
		$start4 =~ s/\<//;
		$end5 = $string[4];
		$end5 =~ s/\>//;
		
		print "$seqid1\tGeneMarkS\t\.\t$start4\t$end5\t\.\t$strand7\t\.\t$id9\n";
	}
}
close GMOUT;

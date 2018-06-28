#!/usr/bin/perl -w

$threshold = 100;

$org_current = "";
while (<STDIN>){
	chomp;
	@m8 = split /\s+/;
	if (abs($m8[9] - $m8[8]) < $threshold){
		next;
	}
	if ($m8[10] > 1e-05){
		next;
	}
	$seqid1 = $m8[1];
	$source2 = $m8[0];
	$source2 =~ /^(\w+)\|([A-Za-z]+)/;
	$source2 = $2;
	$type3 = '.';
	if ($m8[8] < $m8[9]){
		$start4 = $m8[8];
		$end5 = $m8[9];
		$strand7 = '+';
	}
	else {
		$start4 = $m8[9];
		$end5 = $m8[8];
		$strand7 = '-';
	}
	
	$score6 = $m8[10];
	$phase8 = '.';
	$id9 = $m8[0];
	$id9 =~ /(\w+)\|(.+)/;
	$db = $1;
	$attrid = $2;
	$id9 = "ID=$attrid;Dbxref=$db:$attrid";
	$gff3 = $source2.".gff";
	if ($source2 ne $org_current){
		if (-e $org_current.".gff"){
			close SP;
		}
		open SP, ">>", $source2.".gff";
		print SP "##gff-version 3\n";
		$org_current = $source2;
	}
	else{
		print SP "$seqid1\t$source2\t$type3\t$start4\t$end5\t$score6\t$strand7\t$phase8\t$id9\n";
	}
}
close SP


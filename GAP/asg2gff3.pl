#!/usr/bin/perl -w


$threshold = 100;

print "##gff-version 3\n";
while (<STDIN>){
	chomp;
	@m8 = split /\s+/;

	$m8[0] =~ s/_(\d+)\-\d+$//;
	$chunkShift = $1;

	if (abs($m8[7] - $m8[6]) < $threshold){
		next;
	}
	$seqid1 = $m8[0];
	$source2 = $m8[1];
	$source2 =~ /^(\w+)\|([A-Za-z]+)/;
	$source2 = $2;
	$source2 = "ASGARD";
	$type3 = '.';
	if ($m8[6] < $m8[7]){
		$start4 = $m8[6];
		$end5 = $m8[7];
		$strand7 = '+';
	}
	else {
		$start4 = $m8[7];
		$end5 = $m8[6];
		$strand7 = '-';
	}

	$start4 += $chunkShift;
	$end5 += $chunkShift;
	
	$score6 = $m8[10];
	$phase8 = '.';
	$id9 = $m8[1];
	$id9 =~ /(\w+)\|(.+)/;
	$db = $1;
	$attrid = $2;
	$id9 = "ID=$attrid;Dbxref=$db:$attrid";
#	open SP, ">>", $source2.".gff";
	print "$seqid1\t$source2\t$type3\t$start4\t$end5\t$score6\t$strand7\t$phase8\t$id9\n";
#	close SP;
}

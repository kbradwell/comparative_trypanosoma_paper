#!/usr/bin/perl -w

while ($str = <STDIN>){
	chomp $str;
	if ($str =~ /^\#\#sequence/){
		next;
	}
	if ($str !~ /^\#\#/){
		@str = split /\s+/, $str;
		$str[0] =~ /(.+):(\d+)\-(\d+)/;
		$id = $1;
		$coordShift = ($2 - 1);
		$str[0] = $id;
		$str[3] += $coordShift; #recalculate hit's start coordinate
		$str[4] += $coordShift; #recalculate hit's end coordinate
		$str = join "\t", @str;
	}
	print $str."\n";
}
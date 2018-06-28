#!/usr/bin/env perl

#run this after ordering by column k1 then column k2, if possible (sort command)
# updated 2011-01-21, J

use strict;
use warnings;
use Getopt::Long;

my($col,$h,$inFile,$outFile);

GetOptions('c=i' => \$col, 'h' => \$h, 'i=s' => \$inFile, 'o=s' => \$outFile);

if($h) { usage(); }

my $k1 = $col ? $col - 1 : "0";
#my $k2 = 1;

my $prev1 = '';
#my $prev2 = '';

die "$inFile : File not found!\n" unless (-e $inFile);
die "Output file name required!\n" unless (defined $outFile);

open IN, "<$inFile";
open OUT, ">$outFile";

while (<IN>) {
 my @li = split(/\t/, $_);
# if($li[$k1] eq $prev1 && $li[$k2] eq $prev2) { next; }
 if($li[$k1] eq $prev1) { next; }
 print OUT $_;
 $prev1 = $li[$k1];
# $prev2 = $li[$k2];
}
close(IN);
close(OUT);

exit;


sub usage {
  print "
rem_dup
-------

Example: rem_dup -c 3 filename

Removes consecutive lines that have the same content in one column (-c).
To remove all duplicates from the whole file, sort file first by the 
column of interest (e.g. sort -k 2,2 filename) before running rem_dup.

Options:
-h\tshows this message and exits.
-c\tdefines which column to examine (default: first column).

Licensed under the GPL v. 3.
Copyright J.M.P. Alves (jmalves\@vcu.edu), 2011.
Modified by VNK 120724

";

exit;
}

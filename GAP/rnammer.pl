#!/usr/bin/perl
use strict;
use warnings;

if (scalar @ARGV != 5) {
    print "Usage: $0 <kingdom(bak|euk)> <out.gff> <out.fasta> <genome.fasta> <summary.txt>\n";
    print scalar @ARGV,"\n";
    exit;
}

my $rnammerKingdom=shift;
my $rnammerGff=shift;
my $rnammerFasta=shift;
my $fasta=shift;
my $rnammerSummary=shift;
my @chars = ("A".."Z", "a".."z");
my $randomDir;
$randomDir .= $chars[rand @chars] for 1..8;
mkdir("/tmp/$randomDir");
my $rnammer_cmd="/usr/global/blp/rnammer-1.2/rnammer -S $rnammerKingdom -gff $rnammerGff -multi -f $rnammerFasta $fasta -m tsu,ssu,lsu -T /tmp/$randomDir";
#my $rnammer_cmd="/usr/global/blp/rnammer-1.2/rnammer -S $rnammerKingdom -gff $rnammerGff -multi -f $rnammerFasta $fasta -m tsu,ssu,lsu";
print $rnammer_cmd."\n";
system($rnammer_cmd);
my $nrna=`grep -v "^#" ${rnammerGff}|wc -l`;
chomp($nrna);
my $ncontigs=`grep -v "^#" ${rnammerGff}|cut -f 1|sort|uniq|wc -l`;
chomp($ncontigs);
my @types=`grep -v "^#" ${rnammerGff}|cut -f 9|sort|uniq`;
chomp(@types);
open OUT,">$rnammerSummary";
printf OUT "No. of RNA features found        \t\t\t:%d\n",$nrna;
printf OUT "No. of contigs with RNA features \t\t\t:%d\n",$ncontigs;
printf OUT "Type of RNA found                \t\t\t:%s\n",join(",",@types);
printf OUT "rnammer Gff file                 \t\t\t:%s\n",$rnammerGff;
printf OUT "rnammer Fasta File               \t\t\t:%s\n",$rnammerFasta;
close OUT;
exit;

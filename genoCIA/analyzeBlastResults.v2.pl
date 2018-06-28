# Purpose : Analyze Blast results to count how many genes/query have atleast one hit and how many have full alignment.

use strict;
use lib("/gpfs_fs/data2/tol/scripts/GenomeCompletionByReferenceList/");
use FastaReader;

my $bpfile = shift;
my $multihit_flag = shift;
my $identity_threshold = shift;

if(! defined($bpfile)) { 
    die "USAGE: perl $0 bpfile multiple-hit flag\n";
}

if(! defined($multihit_flag)) { 
    $multihit_flag = 0;
}
my $fullhit_threshold = 90;
my $full_threshold = 99;
my $partialhit_threshold2 = 75;
my $partialhit_threshold = 50;

my $h_25hit_threshold = 25;

my $i_threshold = 25;

if(defined($identity_threshold)) { 
    $i_threshold = $identity_threshold;
}

my $DEBUG = 1; # to print individual query into  25, 50, 90, 99 % category file ...
$DEBUG=0;

my ($query,$hitStr,$totalgenes,$noofhits,$nooffullhit,$partialhit,$partialhit2,$fullhit,$identity, $h_25hit);
my @hits;
my ($query_len,$align_len,$align_perc);
my $fr = FastaReader->new();
$fr->init_file($bpfile);
$totalgenes = $noofhits = $nooffullhit = $partialhit = $partialhit2 = $fullhit = $h_25hit = 0;
while(($query,$hitStr) = $fr->raw_next()) { 
    $totalgenes++;

    if(length($$hitStr) > 100) { 
	$noofhits++;
	@hits = split(/SUB:/,$$hitStr);
	shift @hits;
	($query_len) = $$query =~ /len=(\d+)/;
	$align_len = getAlignmentLength($hits[0],$multihit_flag);
	#print "HIT - $hits[0]\n";
	$align_perc = ($align_len / $query_len) * 100;
	#print "$query_len\t$align_len\t$align_perc\n";
	if($align_perc >= $full_threshold && $identity >= $i_threshold) { 
	    $fullhit++;
	    `echo "$$query" >> query_99` if $DEBUG;
	}

	if($align_perc >= $fullhit_threshold && $identity >= $i_threshold) { 
	    $nooffullhit++;
	    `echo "$$query" >> query_90` if $DEBUG;
	}
	if($align_perc >= $partialhit_threshold2 && $identity >= $i_threshold) { 
	    $partialhit2++;
	    `echo "$$query" >> query_75` if $DEBUG;
	}
	
	if($align_perc >= $partialhit_threshold && $identity >= $i_threshold) { 
	    $partialhit++;
	    `echo "$$query" >> query_50` if $DEBUG;
	}
	if($align_perc >= $h_25hit_threshold && $identity >= $i_threshold) {                                  
            $h_25hit++;                                                                                       
	    `echo "$$query" >> query_25` if $DEBUG;
        }
    }
    
}

my @headers = ("File","Multi HSP(1=Yes,0=No)","Total Genes","Atleast Hit","Atleast Hit%",
	       "Partial Hit(>= 25%)", "Partial Hit(>= 25%) %",
	       "Partial Hit(>= 50%)", "Partial Hit(>= 50%) %",
	       "Partial Hit(>= 75%)", "Partial Hit(>= 75%) %",
	       "Partial Hit(>= 90%)", "Partial Hit(>= 90%) %",
	       "Full Hit(>= 99%)", "Full Hit(>= 99%) %"
	       );

my $DL = "\t";
my $printStr = join($DL, @headers);


print "$printStr\n";
	
		 


#print "File = $bpfile\n";
#print "Multiple Alignment Counted = $multihit_flag\n";
#print "Total Genes = $totalgenes\n";
#print "No. of Genes with atleast Hit = $noofhits";
my $noofhits_perc = sprintf("%.3f", $noofhits/$totalgenes) * 100;                                                                                               
#print " $p%\n";
#print "No. of Genes with Partial Hit (>= 25%) = $h_25hit";
my $h_25hit_perc = sprintf("%.3f", $h_25hit/$totalgenes) * 100;
#print " $p%\n";
#print "No. of Genes with Partial Hit (>= 50%) = $partialhit";
my $partialhit_perc = sprintf("%.3f", $partialhit/$totalgenes) * 100;                                                                                               
my $partialhit_perc2 = sprintf("%.3f", $partialhit2/$totalgenes) * 100;
#print " $p%\n";
#print "No. of Genes with Partial Hit(>= 90%) = $nooffullhit";
my $nooffullhit_perc = sprintf("%.3f", $nooffullhit/$totalgenes) * 100;                                                                                               
#print " $p%\n";
#print "No. of Genes with Full Hit(>= 99%) = $fullhit";
my $fullhit_perc = sprintf("%.3f", $fullhit/$totalgenes) * 100;                                                                                               
#print " $p%\n";


$printStr = join($DL, ($bpfile,$multihit_flag,$totalgenes,$noofhits,$noofhits_perc,$h_25hit,$h_25hit_perc,
		       $partialhit,$partialhit_perc,$partialhit2,$partialhit_perc2,$nooffullhit,$nooffullhit_perc,
		       $fullhit,$fullhit_perc));

print "$printStr\n";
exit;

sub getAlignmentLength { 
    my ($hitStr,$multihit_flag)  = @_;
    
    #print "HIT = $hitStr\n";
    
    my ($total_alignlen,$align_len,$i);
    $total_alignlen = 0;
    my @hits = split(/qb:/,$hitStr);
    shift @hits;
    $i = 0;
    foreach my $hit (@hits) { 
	($align_len) = $hit =~ /len:(\d+)/;
	($identity) = $hit =~ /perc:(\d+)/;
	
	if($identity >= $i_threshold) { 
	    $total_alignlen += $align_len;
	}
	if($multihit_flag && $i > 1) { 
	    last;
	}
	$i++;
    }
    
    return $total_alignlen;
}

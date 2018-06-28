#! /usr/bin/perl
# Purpose :  Calculate Genome-wide alignment between reference genome and assembled contigs. 
#            Using Nucmer.
# Nihar Sheth
# 8/16/2010

use strict;

my $projectname = shift;
my $reffile = shift;
my $queryfile = shift;

my $MIN_THRESHOLD = 50;

if(! defined($queryfile)) 
{ 
    die "USAGE: perl $0 project-name reference-fasta  query-fasta\n";
}
if (! -d "./Results")
{
	`mkdir Results`;
}
open (OFH, ">./Results/$projectname.txt") or die "Cannot open out file $projectname.txt:$!\n";
open (OFH2, ">./Results/$projectname.tabbed.txt") or die "Cannot open out file $projectname.tabbed.txt:$!\n";

my $nucmer_dir = "/home/jmalves/bin/MUMmer3.20/";
my $nucmer = "${nucmer_dir}nucmer";
my ($refsize,$querysize);

my $reflabel = "ref";
my $querylabel = "query";


my %seqHash; # keeps the record oflength of each sequence in ref and query files. 
$refsize = calculateSequenceSize($reffile,$reflabel);
$querysize = calculateSequenceSize($queryfile,$querylabel);

print OFH2 "Ref_Query\tRef File\tRef size\tQuery File\tQuerySize\t";
print OFH2 "Avg Iden\tNum of Alignments\tMax Cov\tMin Cov\tR Avg Cov\tQ Avg Cov\t";
print OFH2 "R Align Len\tQ Align Len\tR % Align Len\tQ % Align Len\t";
print OFH2 "R Iden Len\tQ Iden Len\tR % Iden Len\tQ % Iden Len\t\n";


print OFH2 "$projectname\t$reffile\t$refsize\t$queryfile\t$querysize\t";
print OFH "Ref File = $reffile\n";
print OFH "Reference size = $refsize\n";
print OFH "Query File = $queryfile\n";
print OFH "Query size = $querysize\n\n";

my $coordfile = "${projectname}.coords";
my $deltafile = "${projectname}.delta";

if(! -f $coordfile) 
{ 
    print STDERR "Doing nucmer comparison ...\n";
    `$nucmer -p $projectname $reffile $queryfile`;
    `${nucmer_dir}show-coords -r -g $deltafile > $coordfile`;
}

if (-s $coordfile)
{
	print STDERR "\nProcessing Coordinate File\n";
	processNucmerCoordinatorFile($coordfile);
	findUnAlignedRegions($coordfile,$reffile);
	print STDERR "\nGenerating graphical output\n";
	`${nucmer_dir}mummerplot -p $projectname -R $reffile -Q $queryfile --layout --color -postscript $projectname.delta`;
}
else
{
	print STDERR "No Alignments Available\n";
}

close OFH2;
close OFH;

exit;

##SUBS

sub processNucmerCoordinatorFile 
{ 
    my ($coordfile) = @_;

    open FH, "$coordfile" or die "Error in openning $coordfile \n";
    
   	<FH>;
    <FH>;
	<FH>;
    <FH>;
    <FH>;

    
    my @cols;
    my ($rs,$re,$qs,$qe,$ra,$qa,$identity,$ref,$query);

    my ($avgidentity,$noofalignments,$identitysum);
    my ($smallgap,$gapsize,$largegap,$negativegap,$previous_re);
    my ($alignment_len,$alignment_per);
    my ($q_alignment_len,$q_alignment_per);
    my ($identity_len,$identity_per);
    my ($q_identity_len,$q_identity_per);
    my ($avg_cov, $q_avg_cov);
    my ($max_cov, $min_cov);

	$max_cov = 0;
	$min_cov = 1000000000000;
    $identitysum = $noofalignments = 0;
    $alignment_len = $q_alignment_len = 0;
    $identity_len = $q_identity_len = 0;


    my ($avgidentity,$noofalignments,$identitysum);
    my ($smallgap,$gapsize,$largegap,$negativegap,$previous_re);
    my ($alignment_len,$alignment_per);
    my ($q_alignment_len,$q_alignment_per);
    my ($identity_len,$identity_per);
    my ($q_identity_len,$q_identity_per);


    $identitysum = $noofalignments = 0;
    $alignment_len = $q_alignment_len = 0;
    $identity_len = $q_identity_len = 0;


    while(<FH>) 
    { 
		$_ =~ s/\s+\|\s+/\t/g;
		$_ =~ s/^\s+//g;
		$_ =~ s/\s+/\t/g;
		($rs,$re,$qs,$qe,$ra,$qa,$identity,$ref,$query) =
			split(/\t/,$_);
		
		
		$noofalignments++;
		$identitysum += $identity;
		$alignment_len += $ra;
		$q_alignment_len += $qa;
	
			$identity_len += ($ra * ($identity / 100));
			$q_identity_len += ($qa * ($identity / 100));
			
		
			if ($max_cov < $ra)
			{
				$max_cov = $ra;
			}
		
			if ($min_cov > $ra)
			{
				$min_cov = $ra;
			}
		
    }

    $avgidentity = ($identitysum / $noofalignments);
    
    $avg_cov = ($alignment_len / $noofalignments);
    $q_avg_cov = ($q_alignment_len / $noofalignments);
    
    $alignment_per = (($alignment_len / $refsize) *100);
    $q_alignment_per = (($q_alignment_len / $querysize) *100);
    
    $identity_per = (($identity_len / $refsize) *100);
    $q_identity_per = (($q_identity_len / $querysize) *100);
    
	print OFH2 "$avgidentity\t$noofalignments\t$max_cov\t$min_cov\t$avg_cov\t$q_avg_cov\t";
	print OFH2 "$alignment_len\t$q_alignment_len\t$alignment_per\t$q_alignment_per\t";
	print OFH2 "$identity_len\t$q_identity_len\t$identity_per\t$q_identity_per\n";

    print OFH "Average Identity = $avgidentity\n\n";
    
   	print OFH "Ref Average Coverage = $avg_cov\n";
    print OFH "Query Average Coverage = $q_avg_cov\n\n";
    
    print OFH "Number of alignments = $noofalignments\n\n";
    
    print OFH "Max Coverage = $max_cov\n";
    print OFH "Min Coverage = $min_cov\n\n";
    
    print OFH "Ref Alignement Length = $alignment_len\n"; 
    print OFH "Query Alignement Length = $q_alignment_len\n"; 
    print OFH "Ref Percent Alignment Length = $alignment_per\n";
	print OFH "Query Percent Alignment Length = $q_alignment_per\n\n";

    print OFH "Ref Identity Length = $identity_len\n"; 
    print OFH "Query Identity Length = $q_identity_len\n"; 
    print OFH "Ref Percent Identity Length = $identity_per\n";
	print OFH "Query Percent Identity Length = $q_identity_per\n\n";	


    return;
}

sub findUnAlignedRegions { 
    
    my ($coordfile,$reffile) = @_;

    open FH, "$coordfile" or die "Error in openning $coordfile \n";
    
    <FH>;
    <FH>;
    <FH>;
    <FH>;
    <FH>;
    
    
    my @cols;
    my ($rs,$re,$qs,$qe,$ra,$qa,$identity,$ref,$query);
    my ($start,$end);
    my %alignHash;
    my ($seq_id,$seqlen,$label);
    my ($p_seq_id,$p_start,$p_end,$start,$end);
    $p_seq_id = "UNDEF";
    my %printSeqHash;
    while(<FH>) {
	$_ =~ s/\s+\|\s+/\t/g;
	$_ =~ s/^\s+//g;
	$_ =~ s/\s+/\t/g;
	($rs,$re,$qs,$qe,$ra,$qa,$identity,$ref,$query) =
	    split(/\t/,$_);
	$alignHash{$reflabel}{$ref} = 1;
	$alignHash{$querylabel}{$query} = 1;
	
	
	
	if($ref ne $p_seq_id) { # change of ref sequence compared to previous alignment
	    
	    if(exists $seqHash{$reflabel}{$p_seq_id}) { 
		$seqlen = $seqHash{$reflabel}{$p_seq_id}{len};
		$start = $p_end;
		$end = $seqlen - 1;

		if(($end - $start) > $MIN_THRESHOLD) { 
		    push @{${printSeqHash{$p_seq_id}}},[$start,$end];
		 #   print "$p_seq_id\t$start\t$end\n";
		}
	    }
	    $p_seq_id = $ref;
	    $p_start = $rs;
	    $p_end = $re;
	    
	    $seqlen = $seqHash{$reflabel}{$p_seq_id}{len};
	    $start = 1;
	    $end = $p_start;
	    if(($end - $start) > $MIN_THRESHOLD) { 
		push @{${printSeqHash{$p_seq_id}}},[$start,$end];
		#print "$p_seq_id\t$start\t$end\n";
	    }
	    
	} else { 

	    if($rs >= $p_end ) {
		$start = $p_end;
		$end = $rs;
		$seqlen = $seqHash{$reflabel}{$p_seq_id}{len};
		
		if(($end - $start) > $MIN_THRESHOLD) { 
		    push @{${printSeqHash{$p_seq_id}}},[$start,$end];
		    #print "$p_seq_id\t$start\t$end\n";
		} elsif(($end - $start) == 0) { 
		    print STDERR "$_\n";
		}
		$p_seq_id = $ref;
		$p_start = $rs;
		$p_end = $re;
	    } elsif($rs < $p_end && $re > $p_end ) { 
		$p_end = $re;
	    } else { 
		# alignment is totally inside previous larger alignment. Do NOTHING.
	    }
	    
	    
	}

	

    }
    
    my ($no_ref,$no_query);
    $no_ref = keys %{$alignHash{$reflabel}};
    $no_query = keys %{$alignHash{$querylabel}};
    
    print STDERR "Printing list of Ref- seqs which are unaligned \n";
    
    $label = $reflabel;

    foreach $seq_id (sort {$seqHash{$label}{$b} <=> $seqHash{$label}{$a}} 
		     keys %{$seqHash{$label}}) { 
	if(! exists $alignHash{$label}{$seq_id}) { 
	    $seqlen = $seqHash{$label}{$seq_id}{len};
	    push @{${printSeqHash{$seq_id}}},[1,$seqlen];
	    #print "$label\t$seq_id\t$seqlen\n";
	}
    }

    printUnalignedSeq($label,\%printSeqHash);

    print STDERR "Printing list of Query- seqs which are unaligned \n";
    
    $label = $querylabel;
    %printSeqHash = ();
    foreach $seq_id (sort {$seqHash{$label}{$b} <=> $seqHash{$label}{$a}} 
		     keys %{$seqHash{$label}}) { 
	if(! exists $alignHash{$label}{$seq_id}) { 
	    $seqlen = $seqHash{$label}{$seq_id}{len};
	    push @{${printSeqHash{$seq_id}}},[1,$seqlen];
	    #print "$label\t$seq_id\t$seqlen\n";
	}
    }
    printUnalignedSeq($label,\%printSeqHash);
    
    print STDERR "No. of Alignment Contigs Ref = $no_ref\n";
    print STDERR "No. of Alignment Contigs Query = $no_query\n";
    
    return;    
}

sub printUnalignedSeq { 

    my ($label,$hashRef) = @_;
    my ($seq_id,$start,$end,$seqlen,$seq,$posRef,$substr,$substr_len);
    my %printHash = %$hashRef;

    my $filename = "Unaligned_${label}_${projectname}.fa";
    open FH, ">$filename" or die "Error in openning $filename for writing\n";
    foreach $seq_id (keys %printHash) { 
	$seqlen = $seqHash{$label}{$seq_id}{len};
	$seq = $seqHash{$label}{$seq_id}{seq};
	
	foreach $posRef (@{$printHash{$seq_id}}) { 
	    ($start,$end) = @$posRef;
	    
	    if($end == $seqlen) { 
		$substr = $seq;
	    } else { 
		$substr_len = $end - $start;
		$start--;
		$substr = substr($seq,$start,$substr_len);
	    }

	    #print "FILE PRINT\t$label\t$seq_id\t$start\t$end\t$seqlen\n";

	    $substr_len=length($substr);
	    print FH ">${projectname}|${seq_id}|length=${seqlen}|start=${start}|end=${end}|len=$substr_len\n";
	    print FH "$substr\n";
	}
    }
    
    return;

}



sub calculateSequenceSize 
{ 
    my ($fastafile,$label)  = @_;

    my $size = 0;
    my ($seq,$seqlen);
    $seq = "";
    open FH, "$fastafile" or die "Error in openning $fastafile\n";
    my $seq_id = "BLANK";
    while(<FH>) 
    { 
		chomp $_;
		if($_ =~ />(.*)/) { 
		    $seq =~ s/\n//g;
		    $seq =~ s/\r//g;
		    $seqlen = length($seq);
		    $size += $seqlen;
		    #print STDERR "$seq_id\t$seqlen\n";
		    $seqHash{$label}{$seq_id}{len} = $seqlen;
		    $seqHash{$label}{$seq_id}{seq} = $seq;
		    $seq_id = $1;
		    if($seq_id =~ /(contig\d+)/) { 
			$seq_id = $1;
		    }
		    $seq = "";
		} else {
		    $seq .= $_;
		}
    }

    delete $seqHash{$label}{"BLANK"};
    #For last sequence ...
    $seq =~ s/\n//g;
    $seq =~ s/\r//g;
    $seqlen = length($seq);
    $size += $seqlen;
    $seqHash{$label}{$seq_id}{len} = $seqlen;
    $seqHash{$label}{$seq_id}{seq} = $seq;
    
    close FH;

    my $no_seq = keys %{$seqHash{$label}};
    print STDERR "$fastafile\t$no_seq\n";

    return $size;
}

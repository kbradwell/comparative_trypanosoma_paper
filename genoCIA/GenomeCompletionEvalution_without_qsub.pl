# 
use strict;

my $debug = 1;
my $FS = ",";

my $projectname = shift;
my $inputfile = shift;   # gene or protein file
my $blastprogram = shift; # tblastn 
my $input =  shift;   
#my $outdir = shift;

if(! defined($input)) { 
    die "USAGE: perl $0 projectname inputfile blasttype input file or dir\n";
}

if(-f $input) { 
    processFile($input);
} else { 
    my @files = <$input/*.fna>;
    my @files2 = <$input/*.fasta>;
    push(@files,@files2);
    
    foreach my $file (@files ) { 
	processFile($file);
    }
    
}

exit;

sub processFile { 
    my ($file) = @_;
    


    my ($DBname,$fullblastdbname,$stat_file,$blast_statfile,
	$blastout,$bpout,$command,$blasttype);
    my $blast = "/usr/global/blp/bin/blastall";
    my ($statRef,$assemblyline,$blastline,$lastele);
    my @stats;
    
    my @headers = ("Genome","No.of Contigs", "Total NT", "Total NT without gaps","GC NT", "GC %","Gap Total",
		   "GC % without Gap","Avg. Contig","Min Contig","Max Contig",
		   "Total N50 Contigs","N50ContigSize", "N50 Contig Average",
		   "Total Genes","Gens with atleast Hit", "Genes with >=50% Hit",
		   "Genes with >=90% Hit", "Genes with >=99% Hit");
    
    
    my $headerStr = join($FS,@headers);
    print "$headerStr\n";
    print "File - $file\n" if $debug;
    if($file =~ /(.*)_/ ) { 
	($DBname) = $1;
    } else { 
	($DBname) = $file =~ /(.*)\./;  
    }
	print "DB - $DBname\n" if $debug;

    $blasttype = "${projectname}_$blastprogram";
    
    $fullblastdbname = "${DBname}.nhr";
    $stat_file = "${DBname}_stats.txt";
    $blastout = "${DBname}_${blasttype}.blast";
    $bpout = "${DBname}_${blasttype}.bp";
    $blast_statfile = "${DBname}_${blasttype}_blaststats.txt";
    
    

    if(! -f $fullblastdbname) { 
	print "Formating file for Blast DB ...\n" if $debug;
	`/usr/global/blp/bin/formatdb -i $file -p F -t $DBname -n $DBname`;
    }


    if( ! -f $stat_file) { 
	print "Generating Stat file ..\n" if $debug;
	`perl calculateSequenceSize.pl $file > $stat_file`;
    }
    
    if(! -f $blastout) { 
	print "Blasting against the database ...\n" if $debug;
	#$command = "$blast   -d $DBname -i $inputfile -p $blastprogram -o $blastout -e 1e-05 -b 5 -v 5 -a 10  ";
	$command = "$blast   -d $DBname -i $inputfile -p $blastprogram -o $blastout -e 1e-05 -b 5 -v 5 -a 10 -F F ";
	print "Command - $command\n" if $debug;
	#`$command `;
	system($command);
	
	
    } 
    
    if(-f $blastout && ! -f $bpout) { 
	`perl useParse.pl $blastout `;
    }

    
    if(-f $bpout) { 
	print "Generating Blast Stat File ...\n" if $debug;
	#`perl /data2/tol/scripts/GenomeCompletionByReferenceList/analyzeBlastResults.pl $bpout > $blast_statfile `;
	`perl analyzeBlastResults.v2.pl $bpout > $blast_statfile `;
    }
    
    
    $statRef = readBlastStatFile($stat_file);
    @stats = @$statRef;
    #$lastele = pop @stats;
    #if($
    $assemblyline = join($FS,@stats);
    
    
    $statRef = readBlastStatFile($blast_statfile);
    $blastline = join($FS,@$statRef);
    
    print "${DBname}${FS}$assemblyline${FS}$blastline\n";
    
}


sub readBlastStatFile { 
    my $filename = shift;
    open FH, "$filename" or die "Error in openning file $filename \n";
    my @stats;
    my ($stat);
    while(<FH>) { 
	chomp $_;
	if($_ =~ /Gap/) { 
	    last;
	}
	$_ =~ s/\s+?\%//g;
	
	if($_ =~ /=\s+([\d\.]+)$/) { 
	    $stat = $1;
	    #print "Stat - $stat\n";
	    push @stats,$stat;
	}
    }
    
    close FH;
    return \@stats;
}

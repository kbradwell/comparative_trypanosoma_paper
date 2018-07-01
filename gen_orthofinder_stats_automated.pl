# takes groups.txt output from orthomcl and computes a summary table of numbers of groups containing 
# all the different combinations of organisms.

use warnings;
use strict;
use Data::PowerSet 'powerset';

my $inputfile = "OrthologousGroups.txt";
my $DBfile = "genes_w_db_hit";

my $outputfile = "orthomclclassif4dbcmp.txt";
my $summaryfile = "mysummary_groups.txt";

open SUMM, ">$summaryfile" or die "Error in opening $summaryfile for writing\n";

print SUMM "COMM = clusters containing specified organism and at least one other species\n";
print SUMM "MCO = multicopy paralogs, have orthologs\n";
print SUMM "MCPO = multicopy (paralogs or orthologs)\n";
print SUMM "SCO = total single copy ortholog clusters (not necessarily conserved in all species)\n";
print SUMM "SCOC = total single copy orthologs (found in every species)\n";
print SUMM "SING = singletons\n";
print SUMM "TOTF = total clusters (gene families) that the specified organism is found in\n";
print SUMM "TOTG = total number of genes\n";
print SUMM "UNIQ = specific to that species or group of spp. \n";

open IFH, "$inputfile" or die "Error in opening $inputfile\n";
open OFH, ">$outputfile" or die "Error in opening $outputfile for writing\n";

my ($line,$ss,$groupname,@cols,$j,@genenameparts,$org,$totalnumgenes,$orggene,$inputorg,%numORGgroups,%singcopyORGgrps,$key,@masterArray,@inOrgs);
my ($inputgrp,@grpmembers,$size_grpmembers,$grpmember,%seen,%seengrpmember,$foundgrpmem,@lineArray,$size_lineArray,$item);
my (%uniqgrportholog,$element,$are_equal,%commongrportholog,%onlyoutside_grps,$clustercount,%totalorggenes,%found,$onlyoutsideorg,$l);
my (%numORGgroups_genenames,$genename,$singcopyORGgrps_genename,%multicopyORGgrps,%mcoORGgrps,$count,$countgene,$geneword,$singtot);
my (%classcols,$cat,$catcount,$cat2,$catcount2,$perccat,$genes_in_cluster,%scoc);

##### get the list and combinations of organisms from OrthologousGroups.txt (singletons and orthgroups) #####

$singtot=0;

while($line = <IFH>) {
	chomp $line;
	$line =~ s/\n//g;
	$line =~ s/\r//g;
	@lineArray = ();
	if($line =~ /(.*?):\s(.*)/) {
		$groupname = $1;
		@cols = split (/\s/,$2);
		foreach $j (@cols) {
			@genenameparts = split(/\_/,$j);
			$org = $genenameparts[0];
			push @inOrgs, $org ;# creating one large array of all the org names;
		}
	}
}

close IFH;

my (@inputorgs,$total_num_orgs);
my (@grporgs,$total_num_grps);

my %orgseen = ();
@inputorgs = grep { ! $orgseen{ $_ }++ } @inOrgs;
$total_num_orgs = @inputorgs;

my $powerset = powerset( @inputorgs );
for my $p (@$powerset) {
        if (@$p){
                $ss = join("_",@$p);
                push @grporgs,$ss;
        }
}

$total_num_grps = @grporgs;

print SUMM "Number of organisms in OrthologousGroups.txt: $total_num_orgs\n";
print SUMM "Total number of groups (incl. grps containing just 1 organism) that you are comparing: $total_num_grps\n";

#### PROCESS THE FILES AND WRITE TO OUTFILES ####

open IFH, "$inputfile" or die "Error in opening $inputfile\n";

while($line = <IFH>) { 
    chomp $line;
    #$line =~ s/\n//g;
    #$line =~ s/\r//g;
    @lineArray = ();
    if($line =~ /(.*?):\s(.*)/) { 
    $groupname = $1;
    @cols = split (/\s/,$2);
    $clustercount++;
    #print @cols, "\n";
		foreach $j (@cols) {
			@genenameparts = split(/\_/,$j);
			$org = $genenameparts[0];
			$orggene = $j;
			$totalnumgenes++;
			print OFH "TOTG_",$org," ",$totalnumgenes,"\t",$orggene,"\n";
			push @masterArray, $org ;# creating one large array of all the org names;
        	push @lineArray, $org;
		} # FOREACH
	} # IF
	$genes_in_cluster = @lineArray;
	my %seen = ();
    my @unique = grep { ! $seen{ $_ }++ } @lineArray;
	$size_lineArray = @unique;
	
	####### getting the number of paralog/orthologs groups for each org - multicopy #########
	
	foreach $inputorg (@inputorgs) {
		if($line =~ /$inputorg/) {
		$numORGgroups{$inputorg}++;
			if($line =~ /(.*?):\s(.*)/) { 
			$groupname = $1;
			@cols = split (/\s/,$2);
				foreach $j (@cols) {
					@genenameparts = split(/\_/,$j);
					$org = $genenameparts[0];
					$orggene = $j;
					print OFH "TOTF_",$inputorg," ",$groupname,"\t",$orggene,"\n";
				}		
			}
		}
	}

	#       ########### singletons ###########
	
	$geneword = 'gene';
	$countgene =()= $line =~ /$geneword/g;
	if ($countgene == 1) {
	$singtot++;
	if($line =~ /(.*?):\s(.*)/) {
		$groupname = $1;
		@cols = split (/\s/,$2);
			foreach $j (@cols) {
				@genenameparts = split(/\_/,$j);
				$org = $genenameparts[0];
				$orggene = $j;
				print OFH "SING_",$org," ",$groupname,"\t",$orggene,"\n";
			}
	}
	}

	# 	########### single copy orthologs ############
	$count = 0;
	foreach $inputorg (@inputorgs) {
	$count =()= $line =~ /\s$inputorg/g;
		if ($count == 1 && $countgene > 1) {
		$singcopyORGgrps{$inputorg}++;
		if($line =~ /(.*?):\s(.*)/) { 
			$groupname = $1;
			@cols = split (/\s/,$2);
				foreach $j (@cols) {
					@genenameparts = split(/\_/,$j);
					$org = $genenameparts[0];
					$orggene = $j;
					print OFH "SCO_",$inputorg," ",$groupname,"\t",$orggene,"\n";
				}		
			}
		}
   	}
   	
   	# 	########### single copy orthologs - COMMON TO ALL ############
		foreach $inputorg (@inputorgs) {
			if ($size_lineArray == $total_num_orgs && $total_num_orgs == $genes_in_cluster) {
				$scoc{$inputorg}++;
				if($line =~ /(.*?):\s(.*)/) { 
				$groupname = $1;
				@cols = split (/\s/,$2);
					foreach $j (@cols) {
						@genenameparts = split(/\_/,$j);
						$org = $genenameparts[0];
						$orggene = $j;
						print OFH "SCOC_",$inputorg," ",$groupname,"\t",$orggene,"\n";
					}		
				}
			}
		}
   	
# 	########### multi-copy gene fams (paralogs or orthologs) ############
	$count = 0;
	foreach $inputorg (@inputorgs) {
	$count =()= $line =~ /\s$inputorg/g;
		if ($count > 1) {
		$multicopyORGgrps{$inputorg}++;
		if($line =~ /(.*?):\s(.*)/) { 
			$groupname = $1;
			@cols = split (/\s/,$2);
				foreach $j (@cols) {
					@genenameparts = split(/\_/,$j);
					$org = $genenameparts[0];
					$orggene = $j;
					print OFH "MCPO_",$inputorg," ",$groupname,"\t",$orggene,"\n";
				}		
			}
		}
   	}
   	
# 	########### multi-copy gene fams (have orthologs)  ############
	$count = 0;
	foreach $inputorg (@inputorgs) {
	$count =()= $line =~ /\s$inputorg/g;
		if ($count > 1 && $size_lineArray > 1) {
		$mcoORGgrps{$inputorg}++;
		if($line =~ /(.*?):\s(.*)/) { 
			$groupname = $1;
			@cols = split (/\s/,$2);
				foreach $j (@cols) {
					@genenameparts = split(/\_/,$j);
					$org = $genenameparts[0];
					$orggene = $j;
					print OFH "MCO_",$inputorg," ",$groupname,"\t",$orggene,"\n";
				}		
			}
		}
   	}


# ######### Orthologs within genera - unique ##########
# 	
# 	
# 	foreach $inputgrp (@grporgs) {
# 		@grpmembers = split(/_/,$inputgrp);
# 		
# 		# compare the arrays
# 		$are_equal = compare_arrays(\@unique, \@grpmembers);
# 		if ($are_equal == 1) {
# 			$uniqgrportholog{$inputgrp}++;
# 			}
# 	}
# 	
########## Orthologs within genera (use AND in reg ex) and in any other genera or org ########

	foreach $inputgrp (@grporgs) {
		@grpmembers = split(/_/,$inputgrp);
		$size_grpmembers = @grpmembers;
		my $flag = 0;
		foreach my $i (@unique) {
    		foreach my $k (@grpmembers) {
        			if ($i eq $k) {
            		$flag++;
        			}
    		}
		}
		if ($flag == $size_grpmembers && $size_lineArray > $size_grpmembers) {
		$commongrportholog{$inputgrp}++;
		if($line =~ /(.*?):\s(.*)/) { 
			$groupname = $1;
			@cols = split (/\s/,$2);
				foreach $j (@cols) {
					@genenameparts = split(/\_/,$j);
					$org = $genenameparts[0];
					$orggene = $j;
					print OFH "COMM_",$inputgrp," ",$groupname,"\t",$orggene,"\n";
				}		
			}
		}
	}
	
########## Orthologs uniq to genera/grps ########

	foreach $inputgrp (@grporgs) {
		@grpmembers = split(/_/,$inputgrp);
		$size_grpmembers = @grpmembers;
		my $flag = 0;
		foreach my $i (@unique) {
    		foreach my $k (@grpmembers) {
        			if ($i eq $k) {
            		$flag++;
        			}
    		}
		}
		if ($flag == $size_grpmembers && $size_lineArray == $size_grpmembers && $countgene > 1) {
		$uniqgrportholog{$inputgrp}++;
		if($line =~ /(.*?):\s(.*)/) { 
			$groupname = $1;
			@cols = split (/\s/,$2);
				foreach $j (@cols) {
					@genenameparts = split(/\_/,$j);
					$org = $genenameparts[0];
					$orggene = $j;
					print OFH "UNIQ_",$inputgrp," ",$groupname,"\t",$orggene,"\n";
				}		
			}
		}
	}
	

} # WHILE

    # sub compare_arrays {
#         my ($first, $second) = @_;
#         no warnings;  # silence spurious -w undef complaints
#         return 0 unless @$first == @$second;
#         for (my $i = 0; $i < @$first; $i++) {
#             return 0 if $first->[$i] ne $second->[$i];
#         }
#         return 1;
#     }

print SUMM "Total number of genes: ", $totalnumgenes, "\n";

foreach $key (sort keys %numORGgroups) {
	print SUMM "Number of groups containing ", $key,"\t",$numORGgroups{$key},"\n";
}

foreach $key (sort keys %singcopyORGgrps) {
	print SUMM "Single copy groups for ", $key,"\t",$singcopyORGgrps{$key},"\n";
}

foreach $key (sort keys %uniqgrportholog) {
	print SUMM "Groups uniquely containing ", $key,"\t",$uniqgrportholog{$key},"\n";
}

foreach $key (sort keys %commongrportholog) {
	print SUMM "Groups containing $key and other organisms ","\t",$commongrportholog{$key},"\n";
}

foreach $key (sort keys %scoc) {
	print SUMM "Single copy orthologs conserved across all species incl. $key ","\t",$scoc{$key},"\n";
}


# 	############# get total number of genes from each org - paralog or ortholog ##############


print SUMM "total number of clusters for all orgs: $clustercount\n";
print SUMM "all stats for table, Venn etc.:\n\n";

foreach $org (@masterArray) {
	foreach $inputorg (@inputorgs) {
		if($org eq $inputorg) {
		$totalorggenes{$inputorg}++;
		}
	}
}

foreach $key (sort keys %totalorggenes) {
	print SUMM "Total genes for ", $key,"\t",$totalorggenes{$key},"\n";
}

close IFH;
close OFH;

####### Running commands on command line ###########
# At this point the stats for the Venn and table have been produced but need to be fully printed and % to DB calculated
my $output;
my $command = "cat $outputfile | cut -f 1 -d '\t' |sort|uniq|cut -f 1 -d ' '|sort|uniq -c";
$output = runSystemCommand($command); # calls function that executes the command and assigns it to the variable $output.
print SUMM "$output\n";

my $classesfile = "classesfile.txt";
$command = "cat $outputfile | cut -f 1 -d '\t' |sort|uniq|cut -f 1 -d ' '|sort|uniq -c > $classesfile";
$output = runSystemCommand($command);

# if database hits file is provided

if (-e "./$DBfile") {

my $foundinDB = "foundinDB.txt";
$command = "while read a;do grep \$a $outputfile;done < $DBfile > $foundinDB";
$output = runSystemCommand($command); # calls function that executes the command and assigns it to the variable $output.

print SUMM "Stats for the number of hits to databases for each category: ";

$command = "cat $foundinDB | cut -f 1 -d '\t' |sort|uniq|cut -f 1 -d ' '|sort|uniq -c";
$output = runSystemCommand($command);
print SUMM "$output\n";

my $classesfilefound = "classesfilefound.txt";
$command = "cat $foundinDB | cut -f 1 -d '\t' |sort|uniq|cut -f 1 -d ' '|sort|uniq -c > $classesfilefound";
$output = runSystemCommand($command);

my $outfile_perc = "percentage_hits.txt";

open OFH2, ">$outfile_perc" or die "Error in opening $outfile_perc for writing\n";

open IFH4, "$classesfile" or die "Error in opening $classesfile\n";
while($line = <IFH4>) { 
    chomp $line;
    $line =~ /\s+(\d+)\s(\w+)/;
    $catcount = $1;
    $cat = $2;
    $classcols{$cat} = $catcount;
}

print SUMM "Percentages of hits to databases: \n";
open IFH5, "$classesfilefound" or die "Error in opening $classesfilefound\n";
	while($line = <IFH5>) { 
    chomp $line;
	$line =~ /\s+(\d+)\s(\w+)/;
    $catcount2 = $1;
    $cat2 = $2;
		if(exists $classcols{$cat2}) {
			$perccat = ($catcount2 / $classcols{$cat2}) * 100;
			print SUMM "$cat2 $perccat\n";
			print OFH2 "$cat2 $perccat\n";
		}
	}
	

close IFH4;
close IFH5;
close OFH2;

} # if database file is provided

close SUMM;

exit;

sub runSystemCommand { #goal of sub is to run command
    my ($cmd) = @_;
    my $output = ""; 
    print SUMM "CMD = $cmd\n";
	$output = `$cmd`;
    return $output;
}


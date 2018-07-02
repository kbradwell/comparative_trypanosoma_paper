#!/opt/bin/perl -w
#######  #######  #######  #######  
## parseBlast.pm
## Version 1.0 prefinal release
#######  #######  #######  #######  
use strict;
package Util::parseBlast;
use Symbol;
my $MTCH=0;## make zero for normal operation

##### ##### ##### ##### ##### ##### ##### 
##### new constructor.
##### ##### ##### ##### ##### ##### #####
sub new{	
 my $name = shift;
 my $class = ref($name) || $name;
 my $this = {};

 my $mtch=shift;
 if(defined($mtch)){$MTCH=$mtch;} 
 bless $this,$class;
 return $this;
}

sub read_blast{
  my $this = shift;
  my $in = shift;
  my $out = shift;

  if($MTCH==1){$this->read_blast_new($in,$out);}
  elsif($MTCH==2){$this->read_blast_old($in,$out);}
  else{$this->read_blast_old($in,$out);}
}

#########################################################################
## read_blast
#########################################################################
sub read_blast_old{
  my $this = shift;
  my $in = shift;
  my $out = shift;
  my $last_line;

  ## gets query and its length
  while(my ($cur_query,$cur_len) = $this->query_len($in,$last_line)){ 
    ## reads upto next Query=, last_line contains the next Query=
    my $report;
    $cur_query =~ s/\s*Query=\s*//;
    # print ">$cur_query len=$cur_len\n";
    $report = ">$cur_query len=$cur_len\n";
    my $P_cur_body;
    ($P_cur_body,$last_line) = $this->mtch_body($in);
    ##############################
    ## cleans up trailing garbage.
    ##############################
    #    my $rev =-1;
    my $rev =@{$P_cur_body};$rev--;
    while($P_cur_body->[$rev] !~ /Sbjct/){
      $rev--;if($rev<0){ last;}
    }
    splice(@$P_cur_body,++$rev);
    # if(@$P_cur_body ==0) { next; } 
    my $index = 0;

    ## analyse the body, get back HSP by HSP
    ###  get sub/len of HSP
    while(my ($cur_sub,$sub_len) = $this->sub_len($P_cur_body,\$index)){
      # print "$cur_sub length = $sub_len\n";
      $cur_sub =~ s/\>//g;
      $report .= "SUB:$cur_sub len=$sub_len\n";
      ## get HSP body to be studied.
      my $P_hsps = $this->hsp_body($P_cur_body,\$index);
      #      print "HSPs to be analysed\n";
      #      print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
      #      print "@$P_hsps\n";
      #      print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
      #######################################
      ## now analyse the hsps
      #######################################
      my $hsp_last_line;
      my $hindex=0;
      ## use Score to delimit the HSPs   ## get score details
      while(my $hsp_sc = $this->hsp_score($P_hsps,\$hindex)){ 
#	print "$hsp_sc\n";
# old	$hsp_sc =~ /Expect\s*=\s*([\-\w\.]+)/;my $p = $1;
##NEW  allow for Expect(2)=
	$hsp_sc =~ /Expect[^=]*\s*=\s*([\-\w\.]+)/;my $p = $1;
	my ($bts,$sc);	
	if($hsp_sc =~ /Score\s*=\s*([\-\+\w\.]+)\s*bits\s*\((\d+)\)/){
	  $bts=$1;$sc=$2;
	} elsif($hsp_sc =~ /Score\s*=\s*([\d\.]+)\s*\(([\d\.]+)\s*bits\)/){
	  ## WASHU case Score = 261 (123.1 bits)
	  $bts=$2;$sc=$1;
	} 
	my ($mtch,$allen,$perc,$pos);
	if($hsp_sc =~ /Identities\s*\=\s*(\d+)\/(\d+)\s*\(([\d\.]+)\%\)/){ 
	  $mtch=$1;$allen=$2;$perc = $3;
	}
	my ($posperc);
	if($hsp_sc =~ /Positives\s*\=\s*(\d+)\/(\d+)\s*\(([\d\.]+)\%\)/){
	    ## protein case,
	  ## WASHU Identities = 48/81 (59%), Positives = 58/81 (71%)
	  $pos=$1;$posperc = $3;
	} 

	##NEW frame is only for the tblastn case
	my ($frame);
	if($hsp_sc =~ /Frame\s*\=\s*(\S+\s*\/\s*\S+)/)
	  {$frame=$1;$frame=~s/\s+//g;}

	my ($qry,$P_qry,$sbj,$P_sbj,$mis_mat)
	  =  $this->hsp_mtch($P_hsps,\$hindex);
	my ($sb,$se) = @{$P_sbj};my ($qb,$qe) = @{$P_qry};

	
##OLD
	my $qgps = $qry =~ tr/\-/\-/;
	my $sgps = $sbj =~ tr/\-/\-/;
##NEW allow for the * in tblasn reports
	my $tqgps = $qry =~ tr/\*/\*/;
	my $tsgps = $sbj =~ tr/\*/\*/;
	$qgps += $tqgps;$sgps += $tsgps;

	$report .= "qb:$qb qe:$qe sb:$sb se:$se ";
	$report .= "sc:$sc bts:$bts perc:$perc mtch:$mtch len:$allen ";
#OLD	if(defined($pos))
	if(0){
	  if(0){ 
	   $report .= "qgps:$qgps sgps:$sgps p:$p pos:$pos posperc:$posperc\n";
	  } else { 
	    $report .= "qgps:$qgps sgps:$sgps p:$p\n";
	  } 
	} 
#NEW
	$report .= "qgps:$qgps sgps:$sgps p:$p";
	if($MTCH==2){ 
	  $report .= "\n$mis_mat";
	} 
	if(defined($pos))
	{
	  $report .= " pos:$pos posperc:$posperc";
	}  
	if(defined($frame))
	{ 
	  $report .= " frame:$frame";
	} 

	$report .="\n";

#	print "qry = $qry, @$P_qry\n";print "sbj = $sbj, @$P_sbj\n";
#	print "mis = $mis_mat\n";
#	print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
#	print $report;
#	exit;
      } 
    } 
    print $out $report;
  }
}
sub read_blast_new{
  my $this = shift;
  my $in = shift;
  my $out = shift;
  my $last_line;
  my %repHash;
  $repHash{"A"} = "a";
  $repHash{"G"} = "g";
  $repHash{"C"} = "c";
  $repHash{"T"} = "t";
  $repHash{"-"} = "-";
  $repHash{"*"} = "*";
  $repHash{"N"} = "n";
  
  ## gets query and its length
  while(my ($cur_query,$cur_len) = $this->query_len($in,$last_line)){ 
    ## reads upto next Query=, last_line contains the next Query=
    my $report;
    $cur_query =~ s/\s*Query=\s*//;
    # print ">$cur_query len=$cur_len\n";
    $report = ">$cur_query len=$cur_len\n";
    my $P_cur_body;
    ($P_cur_body,$last_line) = $this->mtch_body($in);
    ##############################
    ## cleans up trailing garbage.
    ##############################
    #    my $rev =-1;
    my $rev =@{$P_cur_body};$rev--;
    while($P_cur_body->[$rev] !~ /Sbjct/){
      $rev--;if($rev<0){ last;}
    }
    splice(@$P_cur_body,++$rev);
    # if(@$P_cur_body ==0) { next; } 
    my $index = 0;

    ## analyse the body, get back HSP by HSP
    ###  get sub/len of HSP
    while(my ($cur_sub,$sub_len) = $this->sub_len($P_cur_body,\$index)){
      # print "$cur_sub length = $sub_len\n";
      $cur_sub =~ s/\>//g;
      $report .= "SUB:$cur_sub len=$sub_len\n";
      ## get HSP body to be studied.
      my $P_hsps = $this->hsp_body($P_cur_body,\$index);
      #      print "HSPs to be analysed\n";
      #      print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
      #      print "@$P_hsps\n";
      #      print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
      #######################################
      ## now analyse the hsps
      #######################################
      my $hsp_last_line;
      my $hindex=0;
      ## use Score to delimit the HSPs   ## get score details
      while(my $hsp_sc = $this->hsp_score($P_hsps,\$hindex)){ 
#	print "$hsp_sc\n";
# old	$hsp_sc =~ /Expect\s*=\s*([\-\w\.]+)/;my $p = $1;
##NEW  allow for Expect(2)=
	$hsp_sc =~ /Expect[^=]*\s*=\s*([\-\w\.]+)/;my $p = $1;
	my ($bts,$sc);	
	if($hsp_sc =~ /Score\s*=\s*([\-\+\w\.]+)\s*bits\s*\((\d+)\)/){
	  $bts=$1;$sc=$2;
	} elsif($hsp_sc =~ /Score\s*=\s*([\d\.]+)\s*\(([\d\.]+)\s*bits\)/){
	  ## WASHU case Score = 261 (123.1 bits)
	  $bts=$2;$sc=$1;
	} 
	my ($mtch,$allen,$perc,$pos);
	if($hsp_sc =~ /Identities\s*\=\s*(\d+)\/(\d+)\s*\(([\d\.]+)\%\)/){ 
	  $mtch=$1;$allen=$2;$perc = $3;
	}
	my ($posperc);
	if($hsp_sc =~ /Positives\s*\=\s*(\d+)\/(\d+)\s*\(([\d\.]+)\%\)/){
	    ## protein case,
	  ## WASHU Identities = 48/81 (59%), Positives = 58/81 (71%)
	  $pos=$1;$posperc = $3;
	} 

	##NEW frame is only for the tblastn case
	my ($frame);
	if($hsp_sc =~ /Frame\s*\=\s*(\S+\s*\/\s*\S+)/)
	  {$frame=$1;$frame=~s/\s+//g;}

	my ($qry,$P_qry,$sbj,$P_sbj,$mis_mat)
	  =  $this->hsp_mtch($P_hsps,\$hindex);
	my ($sb,$se) = @{$P_sbj};my ($qb,$qe) = @{$P_qry};

	
##OLD
	my $qgps = $qry =~ tr/\-/\-/;
	my $sgps = $sbj =~ tr/\-/\-/;
##NEW allow for the * in tblasn reports
	my $tqgps = $qry =~ tr/\*/\*/;
	my $tsgps = $sbj =~ tr/\*/\*/;
	$qgps += $tqgps;$sgps += $tsgps;

	$report .= "qb:$qb qe:$qe sb:$sb se:$se ";
	$report .= "sc:$sc bts:$bts perc:$perc mtch:$mtch len:$allen ";
#OLD	if(defined($pos))
	if(0){
	  if(0){ 
	   $report .= "qgps:$qgps sgps:$sgps p:$p pos:$pos posperc:$posperc\n";
	  } else { 
	    $report .= "qgps:$qgps sgps:$sgps p:$p\n";
	  } 
	} 
#NEW
	$report .= "qgps:$qgps sgps:$sgps p:$p";
	$qry =~ tr/a-z/A-Z/;
	$sbj =~ tr/a-z/A-Z/;

	
       
	if($MTCH==1){ 
	    
	    my @mismtchs = ($mis_mat=~/(\s)/g);
	    my $misStr;
	    if(@mismtchs >= 1) {
		my $missize = @mismtchs;
		my $len=length $qry;
		my ($i,$count,$qc,$sc);
		$misStr="";
		for($i = 0;$i < $len;$i++) { 
		    $qc = substr($qry,$i,1);
		    $sc = substr($sbj,$i,1);
		    
		    if($qc ne $sc) {
			$qc =~ tr/A-Z/a-z/;
			$sc =~ tr/A-Z/a-z/;
			
			substr($qry,$i,1) = $qc;
			substr($sbj,$i,1) = $sc;
			$count=$i+1;
			$misStr .= "[P=$count:Q=$qc:S=$sc]";
		    }
		}
		#$report .= "\nM:$missize";
	    }
	    $report .="\nQ:$qry";
	    $report .="\nS:$sbj";
	    if(@mismtchs >= 1) {
		$report .= "\nM:$misStr";
	    }
	    
	} 
	if(defined($pos))
	{
	  $report .= " pos:$pos posperc:$posperc";
	}  
	if(defined($frame))
	{ 
	  $report .= " frame:$frame";
	} 

	$report .="\n";

#	print "qry = $qry, @$P_qry\n";print "sbj = $sbj, @$P_sbj\n";
#	print "mis = $mis_mat\n";
#	print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
#	print $report;
#	exit;
      } 
    } 
    print $out $report;
  }
} 

 
#########################################################################
#####returns hsp_score stuff
#########################################################################
sub hsp_score{
  my $this = shift;
  my $P_arr = shift;
  my $P_index = shift;
  my $line;
  ## go upto the score line
  while($$P_index<@$P_arr){
    ### for($$P_index =$$P_index;$$P_index<@$P_arr;$$P_index++){
    my $loc_line = $P_arr->[$$P_index];
    $$P_index++;
    if($loc_line =~ /Score/){$line = $loc_line;last;}
  } 
  if(!defined($line)){ return;} 
  #  $P_arr->[$$P_index];
  #  if($line !~ /Score/){die "something wrong in hsp_score";}   
  chomp($line);
  my $next_hsp_sc = $line;
#  $$P_index++;
  while($line =  $P_arr->[$$P_index]){ 
    chomp($line);
    if($line !~ /Query/){$next_hsp_sc .= " ".$line;} 
    else { last;} 
    $$P_index++;
    if($$P_index >=@$P_arr ){last;}
  }  
  return($next_hsp_sc);
}
#########################################################################
#####  called with arr,index
## returns all (HSPs,last_line);
## else returns (HSPs,undef)
#########################################################################
sub hsp_mtch{
  my $this = shift;
  my $P_arr = shift;
  my $P_index = shift;
  my ($qry,$mis_mat,$sbj,@qry,@sbj);

  ($mis_mat,$sbj,$qry) = ("","","");
  while(my $line =  $P_arr->[$$P_index]){ 
    if($line =~ /Score/){last;} 
    $$P_index++;
    if($line =~ /^\s*$/){ next;}
    if($line !~ /^(\s*Query\:\s*(\d+)\s*)([^\d][\w\-\*]+)\s*(\d+)/i){
     print "$line\n";
     die "error in Query of hsps";
    } 
    my $buf = length($1);
    if(@qry < 1){ push(@qry,$2);push(@qry,$4);} 
    else {$qry[1] = $4;}  
    $qry .= $3;

   # Query: 35  tccagcctgggggacaagagcaaaactccatctccaaaaaaaaaaagaaaagaa 88
    $line =  $P_arr->[$$P_index];  $$P_index++;
   #            |||||||||||  |||||||| ||||||| |||| ||||||||||| |||| ||
    chomp($line); #   $line =~ s/\s$//;
    if($line !~ /^\s{$buf}([\w\-\|\*\+\s]+)/i){die "error in mismat<$line>";} 
    $mis_mat .= $1; # $mis_mat =~ s/ $//;

    $line =  $P_arr->[$$P_index];  $$P_index++;
    # Sbjct: 250 tccagcctgggcaacaagagcgaaactccgtctcaaaaaaaaaaaaaaaaaaaa 303
#OLD 
## Fix the problem here and elsewhere of general mismatches getting through
## without alerting user.
#    if($line !~ /^(\s*Sbjct\:\s*(\d+)\s*)([^\d][\w\-]+)\s*(\d+)/i){
    if(0){
      my ($pre,$on,$post)=($`,$&,$');
      print "<$line>\n pre=$pre\n mtch=$on\n post=$post\n";
      print "1=$1,2=$2,3=$3,4=$4\n";
      die "bad news";
    }

    # print "sbjline=$line\n";exit;
######### NEW
    if($line !~ /^(\s*Sbjct\:\s*(\d+)\s*)([^\d][\w\-\*]*)\s*(\d+)/i)
      {die "error in Sbjct of hsps \n<$line>\n";} 
    if(@sbj < 1){ push(@sbj,$2);push(@sbj,$4);} 
    else { $sbj[1] = $4;}  
    $sbj .= $3;
    if($$P_index >= @$P_arr ){last;}
  } 
  return($qry,\@qry,$sbj,\@sbj,$mis_mat);
}
#########################################################################
##### called with current line (the guy with >blahblah
####   and index of array, returns (sub,len)
#########################################################################
sub sub_len{
  my $this = shift;
  my $P_arr = shift;
  my $P_index = shift;  
  my $line;
  ## looks for start the next HSP
  while($$P_index < @$P_arr){ 
    my $loc_line = $P_arr->[$$P_index];
    $$P_index++;
    if($loc_line =~ /^>/){$line = $loc_line;last;}
  }
  if(!defined($line)){ return ;} 
#  if($line !~ /^>/){die "something wrong in sub_len call";}   
  chomp($line);
  my $next_sub = $line;
#  $$P_index++;
  while($line =  $P_arr->[$$P_index]){ 
    chomp($line);
    if($line !~ /Length/){$next_sub .= $line;} 
    else { last;} 
    $$P_index++;
    if($$P_index > (@$P_arr-1) ){last;}
  }  
  my $next_len;
  if($line =~ /Length\s*\=\s*(\d+)/i) {$next_len = $1;} 
  else{die "something is wrong with length";}
  $$P_index++;
  return($next_sub,$next_len);
}
#########################################################################
#####  called with arr,index
## returns all (HSPs,last_line);
## else returns (HSPs,undef)
#########################################################################
sub hsp_body{
  my $this = shift;
  my $P_arr = shift;
  my $P_index = shift;
  my $line = $P_arr->[$$P_index];
  $$P_index++;
  my @hsps;
  push(@hsps,$line);
  while($line =  $P_arr->[$$P_index]){ 
    if($line =~ /^>/){last;} 
    else {push(@hsps,$line);} 
    $$P_index++;
    if($$P_index>(@$P_arr-1) ){last;}
  }  
  return(\@hsps);
}
#########################################################################
#########################################################################
##### called with current line (the guy with Query=
####   and filehandle, returns (query,len)
#########################################################################
#########################################################################
sub query_len{
  my $this = shift;
  my $fh = shift;
  my $line = shift;
  if($line !~ /^Query=/){
    while(my $line1=<$fh>){if($line1=~/^Query=/){$line=$line1;last;}}
  }   
  if(!defined($line)){return;}
  chomp($line);
  my $next_query = $line;
  while($line = <$fh>) { 
    chomp($line);
    if($line !~ /letters/){$next_query .= $line;} 
    else { last;} 
  }  
  my $next_len;
  if($line =~ /([\d\,]+)\s*letters/) {$next_len = $1;$next_len=~ s/\,//g;} 
  else{die "something is wrong with letters";}
  return($next_query,$next_len);
}
#########################################################################
#####  called with filehandle, 
## returns body,last line if it matches Query=, 
## else returns body,undef.
#########################################################################
sub mtch_body{
  my $this = shift;
  my $fh = shift;
  my @body;
  my $last_line;
  while(my $line = <$fh>) { 
    if($line =~ /^Query=/){$last_line = $line;last;} 
    else { push(@body,$line);} 
  }  
  return(\@body,$last_line);
}


1;

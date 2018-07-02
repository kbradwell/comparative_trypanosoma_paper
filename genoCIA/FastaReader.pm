#!/use/local/bin/perl -w
#######  #######  #######  #######  
## FastaReader
#######  #######  #######  #######  
use strict;
use Symbol;
use IO::File;
package FastaReader;


# my ($P_cur_seq,$P_cur_def);
# my ($cur_seq,$cur_def);
# my $cur_def;

##### ##### ##### ##### ##### ##### ##### 
##### new constructor.
##### ##### ##### ##### ##### ##### #####
sub new{	
 my $name = shift;
 my $class = ref($name) || $name;
 my $this = {};
 bless $this,$class;
 my $div =  shift; ## the demarcation used to show new records.
 if(!defined($div)){ $div=">";}
 $this->{div} = $div;
 $this->{cur_def}=undef;
 $this->{next_def} = undef;
 return $this;
}
#### ################ ############# ###############
#### mergeFiles
### assumes that all the files are in same order etc.
### inefficient use of memory, gets everything in and 
### then does the merge.
#### can be changed, later,
#### ################ ############# ###############
sub mergeFiles{
    my $this=shift;
    my $P_files=shift;
    my $outfile=shift;
    my $fnum=@{$P_files};
    my @defns;my @bodys;
    for (my $i=0;$i<$fnum;$i++) {
	$this->init_file($P_files->[$i]);
	my $cnt=0;
	while (my ($P_defn,$P_body)=$this->raw_next()) {
	    if (!defined($defns[$cnt])){
		$$P_defn=~s/len=\d+\s*$//;
		$defns[$cnt]=$$P_defn;
	    }
	    elsif($defns[$cnt] !~ /$$P_defn/)
	      {die "wrong in ".$P_files->[$i]."<$$P_defn>ne<$defns[$cnt]>";}
	    if (!defined($bodys[$cnt])){$bodys[$cnt]=$$P_body;}
	    else {$bodys[$cnt].=$$P_body;}
	    $cnt++;
	}
	$this->close_file();
    }
    my $dnum=@defns;
    open OUT,">$outfile" or die "$outfile";
    for (my $i=0;$i<$dnum;$i++) {
	print OUT "$defns[$i]\n";
	print OUT "$bodys[$i]";
    }
    close OUT;
    return 1;
}
#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### split_file
#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
sub split_file{
  my $this = shift;
  my $file=shift;
  my $unit_size=shift;

  if(!defined($file) || !defined($unit_size)){ 
    die "USage fr->split_file(filename,size,[clean =0/1],[outdir]]) where \n".
     "\tsize is number of elelments in the fragmented files created \n".
     "\tclean will result in ambiguous codes getting fixed,capitalization,\n".
     "def lines become cleaner etc.";
  }
  my $clean =shift;
  my $outdir=shift;
  if(!defined($clean)){$clean=0;}

  # if($unit_size==1){$clean=1;};

  $this->init_file($file);

  my $tmpfile = $file;
  my ($dir,$ext,$fname)=$this->getDirExt($tmpfile);

  #  my ($dir,$ext) = $tmpfile =~ /^(.+)\/[^\/]+\.([^\.]*)$/;
  #  my ($fname) = $tmpfile =~ /^.+\/([^\/]+)\.[^\.]*$/;

  ### print STDERR "fname=$fname,dir=$dir,ext=$ext\n";exit;
  ## to allow for .qual

  if(defined($outdir)){$tmpfile=~"$outdir/$fname";}
  else{$outdir=$dir;$tmpfile=~ s/\.(fa.*)$//;}

  my $filecnt=0;
  my $fcnt=0;
  $fname="${tmpfile}_$filecnt.$ext";

  ## you want to use the definition, if only one fasta per file
  if($unit_size==1){$fcnt=$unit_size+10;}
  else {open OUT,">$fname" or die ">$fname";} 

  while(my ($P_defn,$P_body)=$this->raw_next()){ 
    $fcnt++;
    if($clean==1){$this->fixDefBod($P_defn,$P_body);}
    if($fcnt>$unit_size){ 
      if($fcnt > ($unit_size+3)){ } 
      else {close OUT or die "$fname";}
      $fcnt=1;$filecnt++;
      if($unit_size==1){ 
	my ($tdefn)=$$P_defn=~/^>(\S+)/;
	$fname = "$outdir/$tdefn.$ext";
      } else{$fname="${tmpfile}_$filecnt.$ext";}
      open OUT,">$fname" or die "$fname";
    }
    print OUT "$$P_defn\n"; 
    print OUT "$$P_body"; 
  }
  close OUT or die "$fname";
  close $this->{file};
}
#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### split_file2
### does the split into a number of files, without knowing how many files
###  exist to begin with
#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
sub split_file2{
  my $this = shift;
  my $file=shift;
  my $file_num=shift;
  if(!defined($file) || !defined($file_num)){ 
    die "USage fr->split_file(filename,num_of_files,[clean =0/1]) where \n".
      "\tnum_of_files is number of fragmented files created \n".
	"\tclean will get ambiguous codes fixed, capitalization,\n".
	  "def lines become cleaner etc.";
  } 
  my $clean =shift;
  if(!defined($clean)){ $clean=0;}
  # if($unit_size==1){$clean=1;};
  $this->init_file($file);
  my $tmpfile = $file;
  #  my ($dir,$ext) = $tmpfile =~ /^(.+)\/[^\/]+\.(fa.*)$/;
  ## to allow for .qual
  #  my ($dir,$ext) = $tmpfile =~ /^(.+)\/[^\/]+\.([^\.]*)$/;

  my ($dir,$ext,$fname)=$this->getDirExt($tmpfile);

#  $tmpfile=~ s/\.(fa.*)$//;
  $tmpfile=~ s/\.$ext$//;

  my @fname;
  my @fh;
  for(my $i=0;$i<$file_num;$i++){
    $fname[$i]="${tmpfile}_$i.$ext";
    $fh[$i] = IO::File->new(">$fname[$i]") or die "$fname[$i]";
    # open $fh[$i],">$fname[$i]" or die "$fname[$i]";
  } 
  my $fcnt=0;
  while(my ($P_defn,$P_body)=$this->raw_next()){ 
    if($clean==1){$this->fixDefBod($P_defn,$P_body);}
    print {$fh[$fcnt]} "$$P_defn\n"; 
    print {$fh[$fcnt]} "$$P_body"; 
    $fcnt++;$fcnt %= $file_num;
  } 
  for(my $i=0;$i<$file_num;$i++){
    close $fh[$i] or die "$fname[$i]";
  } 
  close $this->{file};
}
##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### getDirExt
##### ##### ##### ##### ##### ##### ##### ##### ##### 
sub getDirExt{
    my $this=shift;
    my $tmpfile=shift;
    my ($ext)=$tmpfile =~ /\.([^\.]*)$/;
    my ($fname) = $tmpfile =~ /([^\/]+)\.?[^\.]*$/;
    my ($dir) = $tmpfile =~ /^(.+)\/[^\/]+\.?[^\.]*$/;
    print STDERR  "getDirExt dir=$dir, ext=$ext, fname=$fname\n";
    return ($dir,$ext,$fname);
}
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### fixDefBod
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
sub fixDefBod{
  my $this = shift;
  my $P_defn=shift;
  my $P_body=shift;
  if($$P_defn =~ /\|(\w+\.\d+)\|/){$$P_defn = ">$1";}
  elsif($$P_defn =~ /^>(\S+)/){$$P_defn = ">$1";}
  $$P_body =~ tr/a-z/A-Z/;
  $$P_body =~ tr/RYMKSWHBVD/GCCGGACGGG/; ## get rid of SNPS
}
#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### get_defns
## returns a string, 
##       if complicated defline, then >$simple_Def=>$complicatedDef in a line
##       else just >$defline in a line
#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
sub get_defns{
  my $this = shift;
  my $file=shift;
  if(!defined($file)){
    die "USage fr->get_defns(filename)--lists all defn lines returns string\n";
  } 
  $this->init_file($file);
  my $tmpfile = $file;
  my ($dir,$ext) = $tmpfile =~ /^(.+)\/[^\/]+\.(fa.*)$/;
  my $report="";
  while(my ($P_defn,$P_body)=$this->raw_next()){ 
    if($$P_defn =~ /\|(\w+\.\d+)\|/){$report .= ">$1=$$P_defn\n";}
    else{$report .="$$P_defn\n";} 
  } 
  close $this->{file};
  return $report;
} 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### close_file, initialise the file.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
sub close_file{my $this = shift; close $this->{file};} 

#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### init_file, initialise the file.
####  Give filename as stdin if you want to handle stdin.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
sub init_file{
    my $this=shift;
    my $name = shift;
    die "no input file to init_file in FastaReader" if(!defined($name));
    
    # my $file2 = main::gensym();
    my $file2;

## or else just keep the STDIN stuff, so that you can use legitimate filenames.
## useful in splitting very large files.
    if ($name=~/^stdin$/i) { ## to allow for stdin inputs, for large file handling
	$file2=\*STDIN;
    }else {
       $file2 = IO::File->new("<$name") or die "cannot open $name";
       # open $file2,"<$name" or die "cannot open $name";	
    }
    return $this->init_fileptr($file2);
}
#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### init_fileptr, initialise the file. using a filepointer instead of a name.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### 
sub init_fileptr{
  my $this=shift;
  my $file = shift;
  $this->{file} = $file;

  my $div = $this->{div};

  while($this->{next_def} = <$file>){
    if($this->{next_def} !~ /^\s*$/){last;}
    else {$this->{next_def} = undef;}
  } 
  if(!defined($this->{next_def})){return;}
  chomp($this->{next_def});
  if($this->{next_def} !~ /^$div/){
    print STDERR "error in FASTA file\n";
    close $this->{file};
    return;
  }
}
##################    ##################    ##################    
##### next -- the next definition and sequence  
##################    ##################    ##################    
sub next{
    my $this = shift;
    my ($P_cur_def,$P_cur_seq) = $this->raw_next();
    if(!defined($P_cur_seq)){return;}
    $$P_cur_seq =~ s/\s+//g; ## to take out all the spaces inside
    return($P_cur_def,$P_cur_seq);
}
##################    ##################    ##################    
##### raw_next -- the next definition and sequence(without removing any of
######             the formatting.  
##################    ##################    ##################    
sub raw_next{
  my $this = shift;
  if(!defined($this->{next_def})){return;} 
  $this->{cur_def} = $this->{next_def};
  $this->{next_def} = undef;
  my $cur_seq = "";
  my $file = $this->{file};
  my $div = $this->{div};
  while(defined(my $tmp = <$file>)){
    if($tmp=~/^$div/){
      chomp($tmp);
      $this->{next_def} = $tmp;last;            
    }
    $cur_seq .= $tmp;
  }
  if(!defined($this->{next_def})){ close $this->{file};}
  my $cur_def = $this->{cur_def};
  return(\$cur_def,\$cur_seq);
}
1;

=head1 NAME

FastaReader - reads in FASTA files.

=head1 SYNOPSIS

     my $fasta_reader = Util::FastaReader->new();

     $fasta_reader->init_file("stats.fasta");

     while(my ($defn,$P_seq) = $fasta_reader->next){  }

     next returns pointer to the sequence and the definition line 
     and returns nothing at EOF


=head1 AUTHOR

Ravi Sachidanandam, CSHL. ravi@cshl.org


=cut



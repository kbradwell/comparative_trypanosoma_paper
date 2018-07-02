
use parseBlast;
use Symbol;

my $dir = shift;
my $mtch = shift;
if(!defined($dir)){die "Usage wrapParse.pl dir_in/file_in  [1 if you want the alignment]";}

my $pb = Util::parseBlast->new($mtch);
my $in = gensym();
my $out = gensym();


my @files;
if(-d $dir){
    opendir DIR,$dir or die "cannot open dir $dir";
    @files = readdir(DIR);
    closedir DIR;	
} elsif (-f $dir){ 
  @files = ($dir); $dir="";
} else {die "Usage wrapParse.pl dir_in/file_in";}

foreach(@files){ 
  my $file_in = $_;
  print STDERR "infile = $file_in\n";
  if($file_in !~ /\..?blast.?$/ && $file_in !~ /\.out$/){next;} 
  my $file_out = $file_in; $file_out =~ s/\..?blast.?$//; $file_out .= ".bp"; 
  print STDERR "file_in = $file_in \n file_out = $file_out dir=$dir\n"; 
  
  if($dir!~/^\s*$/){ 
    $file_out=$dir."/$file_out";
    $file_in=$dir."/$file_in";
  }

  open $out,">$file_out" or die "cannot open $file_out"; 
  open $in,"$file_in" or die "cannot open $file_in";
  $pb->read_blast($in,$out);
  close $in;
  close $out;
} 
exit;

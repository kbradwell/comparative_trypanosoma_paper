use strict;
use List::Util qw/max min/;
use Bio::SearchIO; 
unless (@ARGV==1) {
	print "usage: perl $0 <blastFile>\n";
	exit;
}
my $blastFile=shift;
my $in = new Bio::SearchIO(-format => 'blast', 
                           -file   => $blastFile);
my %percentAlignment;
my %bits;
while( my $result = $in->next_result ) {
	my $qname=$result->query_name;
	my $qlen=$result->query_length;
        $percentAlignment{$qname}=();
	$bits{$qname}=();
  ## $result is a Bio::Search::Result::ResultI compliant object
my $nhits=0;
  while( my $hit = $result->next_hit ) {
	my $hname=$hit->name;
	my $alnlen=0;
	push @{$bits{$qname}},$hit->bits;
    ## $hit is a Bio::Search::Hit::HitI compliant object
    while( my $hsp = $hit->next_hsp ) {
	$alnlen+=$hsp->length('total');
        }
	my $percaln=0;
        $percaln=$alnlen*100.0/$qlen if $qlen>0;
	#printf "$qname\t$hname\t0.0\t0.0\t0.0\t%.2f\t%d\n",$percaln,$hit->length;
	#$nhits+=1 if $percaln>90;
	$nhits+=1; 
	push(@{$percentAlignment{$qname}},sprintf("%.2f",$percaln));
     }
	#printf #%d\n",$result->num_hits;
	push @{$percentAlignment{$qname}},sprintf("%.2f",0) unless defined $percentAlignment{$qname};
	push @{$bits{$qname}},sprintf("%.2f",0) unless defined $bits{$qname};
	my @pas=@{$percentAlignment{$qname}};
	my @b=@{$bits{$qname}};
	my $m1=max(@pas);
	my $m2=max(@b);
	#print "$qname\t$nhits\t@pas\t@b\n";
	print "$qname\t$nhits\t$m1\t$m2\n";
  
}
exit;

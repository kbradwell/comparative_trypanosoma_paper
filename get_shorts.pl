
use List::Util qw(max min sum); 
my $file = shift @ARGV;
open (FH, "< $file") or die "Can't open $file for read: $!";
my @a = <FH>;

sub median { $_[0]->[ @{$_[0]} / 2 ] }


my $min     = min(@a);
my $max     = max(@a);

$median=  median(\@a);

# print $median;

$align_reduction=$min/$median*100;

$file =~ s/\..*$//;
$file =~ s/aln1_clustalo_//g;

if($align_reduction<75){print $file."\n"};

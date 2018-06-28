#!/usr/bin/perl
#Developer: Vishal Koparde, Ph.D.
#Created: 120724
#Modified:120724
#Version 1.0 - Called by GAP

use strict;
use warnings;
use Getopt::Long;
use lib qw(/usr/global/blp/perllib);
use Qsub;

our ($gapdir);
if (exists $ENV{'GAPDIR'}) {
$gapdir=$ENV{'GAPDIR'};
} else {
$gapdir="/usr/global/blp/GenomeAnnotationPipeline";
}

sub usage();

usage() if (scalar @ARGV==0);

my ($fastaFaa,$db,$outFile,$workingDir,$help);

GetOptions ( 'faa=s' => \$fastaFaa,
             'db=s' => \$db,
             'out=s' => \$outFile,
             'wd=s' => \$workingDir,
             'help|?' => \$help);

usage() if ($help);

die "$fastaFaa : File not Found!\n" unless (-e $fastaFaa);
die "Output File name required!\n" unless (defined $outFile);
if (!defined $workingDir) {
    $workingDir=`pwd`;
    chomp $workingDir;
}
die "$workingDir : Folder not found!\n" unless (-d $workingDir);

my $rps_cmd="/usr/global/blp/bin/rpsblast ";
$rps_cmd.=" -i $fastaFaa";
$rps_cmd.=" -d $db";
$rps_cmd.=" -a 20";
$rps_cmd.=" -m 8";
$rps_cmd.=" -o ${outFile}.tmp";
my $jname="rps";
my $jout=$outFile.".tmp.out";
my $job_rps=new Qsub(name=>$jname,outfile=>$jout,wd=>$workingDir,cmd=>$rps_cmd,nproc=>"20");
$job_rps->submit();
$job_rps->waitForCompletion();
my $cmd="perl ${gapdir}/rem_dup.pl -i ${outFile}.tmp -o ${outFile}";
system($cmd);
$cmd="for f in `cat Cog.out|awk \'{print \$2}\'|awk -F \"|\" \'{print \$NF}\'`;do grep ^\$f /usr/global/blp/blast_db/cdd/cddid_all.tbl;done >> ${outFile}.tbl";
system($cmd);
exit;


sub usage() 
{
print <<EOF;
Run rpsblast on translated genes

Developer : Vishal N Koparde, Ph. D.
Created   : 120724
Modified  : 120724
Version   : 1.0

options:
--faa      proteins in fasta format
--db       rpsblast db
--wd       working dir (default=current directory)
--out      rpsblast output

EOF
exit 1;
}

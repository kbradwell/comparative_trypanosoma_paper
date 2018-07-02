use strict;
use warnings;
use lib qw(/usr/global/blp/perllib);
use Qsub;
use Cwd;
use Getopt::Long;
use Data::Dumper;
use Bio::DB::Sam;
use threads;
use Thread::Queue;
use Thread::Semaphore;
use threads::shared;
use Parallel::ForkManager;
use Clone 'clone';

my $BWA="/usr/global/blp/bin/bwa";
my $SAMTOOLS="/usr/global/blp/bin/samtools";
my $BEDTOOLS="/usr/global/blp/bedtools/bedtools2-master/bin/bedtools";
my $minDepth=10;
my $snpThreshold=25;

sub usage();

my ($help,$rsfile,$assembly,$gff,$ncpus,$nrbed);

usage() if (scalar @ARGV==0);

my $result = GetOptions ( "readset|rs=s" => \$rsfile,
			  "assembly|f=s" => \$assembly,
			  "gff=s" => \$gff,
			  "ncpus|p=i" => \$ncpus,
			  "nrbed=s" => \$nrbed,
			  "help|h" => \$help
			  );

usage() if (defined $help);
usage() unless defined $rsfile;
usage() unless defined $assembly;
usage() unless defined $gff;
$ncpus=10 unless defined $ncpus;

die "GFF file not found" unless -e $gff;
if (defined $nrbed) {die "No repeat bed file not found" unless -e $nrbed;}

if (-e $assembly) {
	my $found=1;
	for my $ext ('bwt','pac','ann','amb','sa') {
		$found=0 unless -e "${assembly}.${ext}";
	}
	system("$BWA index $assembly") unless $found==1;
} else {
	die "Assembly file not found!";
}

open RS, "<$rsfile";
my $cmd;
my $rsname;
my $rsfiles;
my @bams;
my @jobs;
my @files2delete;
while (my $l=<RS> ) {
	chomp $l;
	my @cols=split/\t/,$l;
	$rsname=$cols[0];
	$rsfiles=$cols[1];
	$rsfiles=~s/,/ /g;
	my $bam="${rsname}.bam";
	push @bams,$bam;
	unless ( -e $bam ) {
		my $workingDir=getcwd();
		my $outfile="${workingDir}/${rsname}.sh";
		open O,">$outfile";
		push @files2delete,$outfile;
		print O "$BWA mem -M -t $ncpus $assembly $rsfiles |$SAMTOOLS view -bS -|$SAMTOOLS sort - $rsname"; 
		close O;
		$cmd="sh $outfile";
		my $job=new Qsub(name=>$rsname,outfile=>"${rsname}.tmp.out",wd=>$workingDir,cmd=>$cmd,nproc=>$ncpus);
		$job->submit();
		push @jobs,$job;
	}
}
close RS;

for my $j (@jobs) {
	printf "Waiting for jobid %d to finish ....\n",$j->{jobid};
	$j->waitForCompletion();
}

for my $f (@files2delete) {
	unlink($f) if -e $f;
}

unless ( -e "reads2assembly.bam" ) {
$cmd="$SAMTOOLS merge reads2assembly.tmp.bam ".join(" ",@bams);
print "Merging alignments...\n";
system($cmd);
$cmd="$SAMTOOLS sort reads2assembly.tmp.bam reads2assembly && $SAMTOOLS index reads2assembly.bam && rm -f reads2assembly.tmp.bam";
system($cmd);
print "reads2assembly.bam created.\n";
}

my $outgff=$gff."+";
my $outgfftmp=$gff."+.tmp";
unless (-e $outgff) {
if (-e $outgfftmp) {unlink $outgfftmp};
open GFF,"<$gff";
my @gfflines;
while (my $l=<GFF>) {
	chomp $l;
	my %tmp = ( bam => "reads2assembly.bam", line => $l, outtmpfile => $outgfftmp);
	push @gfflines,\%tmp;
}
close GFF;


chomp(my $thread_limit = `cat /proc/cpuinfo | grep -c -P '^processor\\s+:'`);
$thread_limit-=2;
#$thread_limit=8;

 
my $pm = Parallel::ForkManager->new($thread_limit);
 
DATA_LOOP:
foreach my $tmpref (@gfflines) {
	my $pid = $pm->start and next DATA_LOOP;
	my %tmp = %$tmpref;
	my $bamfile = $tmp{bam};
	my $outfile = $tmp{outtmpfile};
	my $bam = Bio::DB::Sam->new(-bam=>$bamfile,-fasta=>$assembly);
	my $l = $tmp{line};
	my $result=checkForSNPperline($bam,$l);
	open O, ">>$outfile";
	print O "$result";
	close O;
	$bam->DESTROY;
	$pm->finish; # Terminates the child process
}
$pm->wait_all_children();

system("/bin/sort -k9 $outgfftmp -o $outgff");
generateGFFStats($outgff);
}

if (-e $nrbed) {
	open GFF, "<$outgff";
	open GFF2, ">$outgfftmp";
	my $outgff2=$gff."+.nr";
	my %gene2cov;
	my %gene2ds;
	my %gene2ns;
	while (my $l=<GFF>) {
		chomp $l;
		my @cols=split/\t/,$l;
		next unless ( scalar @cols == 12 );
		$gene2cov{$cols[8]}=pop @cols;
		$gene2ds{$cols[8]}=pop @cols;
		$gene2ns{$cols[8]}=pop @cols;
		printf GFF2 "%s\n",join("\t",@cols);
	}
	close GFF;
	close GFF2;
	system("$BEDTOOLS intersect -wa -f 1.0 -a $outgfftmp -b $nrbed > ${outgfftmp}.nr");
	open GFF3, "<${outgfftmp}.nr";
	open GFF4, ">$outgff2";
	while (my $l=<GFF3>) {
		chomp $l;
		my @cols=split/\t/,$l;
		push @cols, $gene2ns{$cols[8]};	
		push @cols, $gene2ds{$cols[8]};	
		push @cols, $gene2cov{$cols[8]};	
		printf GFF4 "%s\n",join("\t",@cols);
	}
	close GFF3;
	close GFF4;
	generateGFFStats($outgff2);
}

exit;

sub generateGFFStats {
	my ($gfffile)=@_;
	my $ngenes=0;
	my $ngenes_hetero=0;
	my $sum_ds=0.0;
	my $sum_ns=0;
	my $sum_cov=0.0;
	my $sum_cov_hetero=0.0;
	my $tlen=0;
	my $tlen_hetero=0;
	open GFF, "<$gfffile";
	while (my $l=<GFF>) {
		chomp $l;
		my @cols=split/\t/,$l;
		$ngenes+=1;
		$tlen+=$cols[4];
		$tlen-=$cols[3];
		$sum_ns+=$cols[9];
		$sum_ds+=$cols[10];
		$sum_cov+=$cols[11];
		if ($cols[9]!=0) {
			$ngenes_hetero+=1;
			$tlen_hetero+=$cols[4];
			$tlen_hetero-=$cols[3];
			$sum_cov_hetero+=$cols[11];
		}
	}
	close GFF;

	open GFFOUT, ">${gfffile}.stats";
	printf GFFOUT "Number_of_genes=%d\n",$ngenes;
	printf GFFOUT "Number_of_hetero_genes=%d\n",$ngenes_hetero;
	printf GFFOUT "Percent_hetero_genes=%.2f\n",$ngenes_hetero*100.0/$ngenes;
	printf GFFOUT "Transcriptome_length=%d\n",$tlen;
	printf GFFOUT "Transcriptome_hetero_length=%d\n",$tlen_hetero;
	printf GFFOUT "Percent_hetero_Transcriptome_length=%.2f\n",$tlen_hetero*100.0/$tlen;
	printf GFFOUT "Number_of_hetero_sites_per_gene=%.2f\n",$sum_ns/$ngenes;
	printf GFFOUT "Number_of_hetero_sites_per_hetero_gene=%.2f\n",$sum_ns/$ngenes_hetero;
	printf GFFOUT "Average_hetero_site_density_per_hetero_gene(sites_per_kb)=%.2f\n",$sum_ds/$ngenes_hetero;
	printf GFFOUT "Average_hetero_site_density_over_transcriptome(sites_per_kb)=%.2f\n",$sum_ns*1000/$tlen_hetero;
	printf GFFOUT "Average_coverage_all_genes=%.2f\n",$sum_cov/$ngenes;
	printf GFFOUT "Average_coverage_hetero_genes=%.2f\n",$sum_cov_hetero/$ngenes_hetero;
	close GFFOUT;

}

sub checkForSNPperline {
	my ($bam,$l) = (@_);
	my @cols=split/\t/,$l;
	my $contig=$cols[0];
	my $start=$cols[3];
	my $end=$cols[4];
	my ($ns,$ds,$cov)=checkForSNPs($bam,$contig,$start,$end);
	my $result=sprintf "%s\t%d\t%.2f\t%.2f\n",$l,$ns,$ds,$cov;
	return $result;
}

sub checkForSNPs {
my ($bam,$contig,$start,$end)=@_;
#my $bam = Bio::DB::Sam->new(-bam=>"reads2assembly.bam",-fasta=>$assembly);
my @SNPs;  # this will be list of SNPs
my $coverage=0;
my $snp_caller = sub {
        my ($seqid,$pos,$p) = @_;
        if ($pos>=$start && $pos<=$end) {
	my $refbase = $bam->segment($seqid,$pos,$pos)->dna;
        my ($total,$different);
	my %bases=(A=>0,C=>0,G=>0,T=>0);
	$total=0;
        for my $pileup (@$p) {
            my $b     = $pileup->alignment;
            next if $pileup->indel or $pileup->is_refskip;      # don't deal with these ;-)

            my $qbase  = substr($b->qseq,$pileup->qpos,1);
            next if $qbase =~ /[nN]/;

            my $qscore = $b->qscore->[$pileup->qpos];
            next unless $qscore > 25;

            $total++;
	    $bases{uc $qbase}++;
        }
	$coverage+=$total;
	if ($total >= $minDepth) {
		for my $b (keys %bases) {
			$bases{$b}=$bases{$b}*100.0/$total;
		}
		my @l=sort {$b <=> $a} values %bases;
		push @SNPs,"$seqid:$pos" if ( $l[1] > $snpThreshold );
	}
	}
    };
my $region="${contig}:${start}..${end}";
$bam->pileup($region,$snp_caller);
my $geneLen=($end-$start+1);
$coverage/=$geneLen;
my $nsnps=scalar @SNPs;
my $snpdensity=$nsnps*1000.0/$geneLen;
$bam->DESTROY;
undef $bam;
undef @SNPs;
return ($nsnps,$snpdensity,$coverage);
}

sub usage()
{
print <<EOF;
This script align readset to the assemblydatabase and generate statistics

Author : Vishal N Koparde, Ph. D.
Created : 150814
Modified : 150814
Version : 1.0

Required options:
--readset or --rs	readset tabular file... column1=readsetid, column2=comma-separated fastq files
--assembly or --f	assembly fasta file with uniquely identifiable seqids
--gff 			gff file with gene coordinates
Other options:
--ncpus or --p		number of processors to use for alignment, default 10
--nrbed			norepeat bed
--help			to print help
EOF
exit 1;
}

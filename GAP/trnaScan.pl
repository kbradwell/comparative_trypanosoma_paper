#!/usr/bin/perl
#Developer: Vishal Koparde, Ph.D.
#Created: 120716
#Modified:120716
#Version 1.0

use strict;
use warnings;
use Getopt::Long;
use lib qw(/usr/global/blp/perllib);
use Qsub;
use Bio::Tools::tRNAscanSE;
use Bio::Tools::GFF;
use Bio::Perl;
use Bio::SeqIO;
use Bio::Seq;


my ($project,$workingDir,$fastaSorted,$organismType,$summaryFile,$gff,$outFasta);

GetOptions(
'project:s' => \$project, # project name is shortName from GAP
'wd:s' => \$workingDir,
'fasta:s' => \$fastaSorted,
'org:i' => \$organismType,
'summary:s' => \$summaryFile,
'gff:s' => \$gff,
'outFasta:s' => \$outFasta
);

die "Project name not provided!\n" unless (defined $project);
die "Working directory not provided!\n" unless (defined $workingDir);
die "Working directory does not exist!\n" unless (-d $workingDir);
die "Fasta file name not provided!\n" unless (defined $fastaSorted);
die "Fasta file does not exist!\n" unless (-e $fastaSorted);
die "Organism type not provided!\n" unless (defined $organismType);
die "Invalid organism type, should be 1(Bacteria) or 2(Eukaryote)\n" unless ( ($organismType==1) || ($organismType==2));
die "Summary file name not provided!\n" unless (defined $summaryFile);
die "Gff3 file name not provided!\n" unless (defined $gff);
die "Output Fasta file name not provided!\n" unless (defined $outFasta);

opendir (DIR,"/usr/global/blp/tRNAscan-SE-1.23");
my @files2delete;
while (my $file=readdir(DIR)) {
	next if $file eq '.' or $file eq '..';
	next unless ( ($file =~ m/\.cm$/) or ($file =~ m/signal$/));
	my $cmd="ln -s /usr/global/blp/tRNAscan-SE-1.23/$file .";
	system($cmd);
	push @files2delete,$file;
}
closedir (DIR);

unlink("${project}_tRNA") if (-e "${project}_tRNA");
unlink("${project}_tRNA_struct") if (-e "${project}_tRNA_struct");
unlink("$summaryFile") if (-e "$summaryFile");
unlink("$gff") if (-e "$gff"); 

my $tRNAscan_cmd="/usr/global/blp/bin/tRNAscan-SE";
$tRNAscan_cmd.=" -B" if ($organismType==1);
$tRNAscan_cmd.=" -o ${project}_tRNA";
$tRNAscan_cmd.=" -f ${project}_tRNA_struct";
$tRNAscan_cmd.=" -m ${summaryFile}.details";
$tRNAscan_cmd.=" $fastaSorted";
 
my $jname="tRNA_".$project;
my $jout="tRNA_".$project.".tmp.out";
my $trna_job=new Qsub(name=>$jname,cmd=>$tRNAscan_cmd,wd=>$workingDir,outfile=>$jout);
$trna_job->submit();
$trna_job->waitForCompletion();

foreach my $file (@files2delete) {
	my $cmd="rm -f $file";
	system($cmd);
}

my $outgff = Bio::Tools::GFF->new(-file => ">$gff",-gff_version => 3);
my $trnaScanResults = Bio::Tools::tRNAscanSE->new(-file => "${project}_tRNA");

# parse the results
while( my $trna = $trnaScanResults->next_prediction() ) {
     $outgff->write_feature($trna);
     for my $trans ( $trna->get_SeqFeatures() ) {	
	 $outgff->write_feature($trans,$trans->get_SeqFeatures());
    }

}
$trnaScanResults->close();
$outgff->close();

my $start  = qr/^Summary/;
my $finish = qr/^tRNAs with introns:/;
open FILE, "<${summaryFile}.details";
open OFILE, ">$summaryFile";
while ( <FILE> ) {
    if ( /$start/ .. /$finish/ ) {
	next  if /^\s+$/;
	next  if /^Overall/;
	next  if /$finish/;
    print OFILE $_;
    }
}
print OFILE "For detailed summary see file : ${summaryFile}.details\n";
close FILE;
close OFILE;

#extract sequences to fasta file
my $out  = Bio::SeqIO->new(-file => ">$outFasta" , '-format' => 'fasta');

my @seq_object_array = read_all_sequences($fastaSorted,'fasta');
my %seq_hash;
foreach my $seq (@seq_object_array) {
    $seq_hash{$seq->id}=$seq;
}

open GFF,"<$gff";
while (my $line=<GFF>) {
	chomp $line;
	next if $line =~ /^#/;
	my @cols=split/\t/,$line;
	next unless $cols[2] eq "gene";
	my $start=$cols[3];
	my $end=$cols[4];
	my @tmp=split/;/,$cols[-1];
	my @tmp2=split/:/,$tmp[0];
	my $seqName=$tmp2[1]."_".$cols[0]."_${start}_${end}";
	my $s=$seq_hash{$cols[0]}->subseq($start,$end);
	my $seq=Bio::Seq->new(-id=>$seqName,-seq=>$s);
	$seq->display_id($seqName);
	$out->write_seq($seq);
}
close GFF;
$out->close;

exit;

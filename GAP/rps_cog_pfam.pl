#!/usr/bin/perl
#Developer: Vishal Koparde, Ph.D.
#Created: 120724
#Modified:120724
#Version 1.0 - Called by GAP
#Version 2.0 - Splits the genes into 10 smaller fasta files and qsubs them
#            - Includes Kog

use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use POSIX;
use lib qw(/usr/global/blp/perllib);
use Qsub;


sub usage();

usage() if (scalar @ARGV==0);

my ($fasta,$project,$workingDir,$summaryFile,$templateGff,$help);

GetOptions ( 'fasta=s' => \$fasta,
             'project=s' => \$project,
             'wd=s' => \$workingDir,
             'summary=s' => \$summaryFile,
	     'templateGff=s' => \$templateGff,
             'help|?' => \$help);

my $evalue=0.01;
my $praThreshold=50;
my $minAlnLen=100;

usage() if ($help);

die "$fasta : File not Found!\n" unless (-e $fasta);
if (!defined $workingDir) {
    $workingDir=`pwd`;
    chomp $workingDir;
}
die "$workingDir : Folder not found!\n" unless (-d $workingDir);

my $infile=Bio::SeqIO->new(-file=>"<$fasta",-format=>'fasta');
my $nseq=`grep -c ">" $fasta`;
chomp($nseq);
my $nseqPerFile=POSIX::ceil($nseq/10);
my ($i,$j,$newfasta);
for ($i=0;$i<10;$i++) {
    $j=$i+1;
    $newfasta=$fasta."_".$j;
    my $outfile=Bio::SeqIO->new(-file=>">$newfasta",-format=>'fasta');
    my $seqCtr=0;
    while(my $seq=$infile->next_seq()) {
	$outfile->write_seq($seq);
	$seqCtr++;
	last if $seqCtr==$nseqPerFile;	
    }
    $outfile->close();
}

my @jobs;

for (my $i=0;$i<10;$i++) {
    my $j=$i+1;
    my $newfasta=$fasta."_".$j;
    
    my $outFile1=$project.".COG.rps_".$j;
    my $rps_cmd1="/usr/global/blp/bin/rpsblast ";
    $rps_cmd1.=" -i $newfasta";
    $rps_cmd1.=" -d /usr/global/blp/blast_db/cdd/Cog";
    $rps_cmd1.=" -p F";
    $rps_cmd1.=" -e $evalue";
    $rps_cmd1.=" -m 8";
    $rps_cmd1.=" -o ${outFile1}";
    my $jname="rps1_".$i."_".$project;
    my $jout=$jname.".tmp.out";
    my $job_rps1=new Qsub(name=>$jname,outfile=>$jout,wd=>$workingDir,cmd=>$rps_cmd1,memory=>"10G");
    $job_rps1->submit();
    push(@jobs,$job_rps1);
    
    my $outFile2=$project.".PFAM.rps_".$j;
    my $rps_cmd2="/usr/global/blp/bin/rpsblast ";
    $rps_cmd2.=" -i $newfasta";
    $rps_cmd2.=" -d /usr/global/blp/blast_db/cdd/Pfam";
    $rps_cmd2.=" -p F";
    $rps_cmd2.=" -e $evalue";
    $rps_cmd2.=" -m 8";
    $rps_cmd2.=" -o ${outFile2}";
    $jname="rps2_".$i."_".$project;
    $jout=$jname.".tmp.out";
    my $job_rps2=new Qsub(name=>$jname,outfile=>$jout,wd=>$workingDir,cmd=>$rps_cmd2,memory=>"10G");
    $job_rps2->submit();
    push(@jobs,$job_rps2);
    
    my $outFile3=$project.".KOG.rps_".$j;
    my $rps_cmd3="/usr/global/blp/bin/rpsblast ";
    $rps_cmd3.=" -i $newfasta";
    $rps_cmd3.=" -d /usr/global/blp/blast_db/cdd/Kog";
    $rps_cmd3.=" -p F";
    $rps_cmd3.=" -e $evalue";
    $rps_cmd3.=" -m 8";
    $rps_cmd3.=" -o ${outFile3}";
    $jname="rps3_".$i."_".$project;
    $jout=$jname.".tmp.out";
    my $job_rps3=new Qsub(name=>$jname,outfile=>$jout,wd=>$workingDir,cmd=>$rps_cmd3,memory=>"10G");
    $job_rps3->submit();
    push(@jobs,$job_rps3);
    
}

foreach my $job (@jobs) {
    $job->waitForCompletion();
}

my $o1=$project.".COG.rps";
my $o2=$project.".PFAM.rps";
my $o3=$project.".KOG.rps";
unlink($o1) if (-e $o1);
unlink($o2) if (-e $o2);
unlink($o3) if (-e $o3);


for ($i=0;$i<10;$i++) {
    $j=$i+1;
    my $outFile1=$o1."_".$j;
    my $outFile2=$o2."_".$j;
    my $outFile3=$o3."_".$j;
    system("cat ${outFile1} >> ${o1}");
    system("cat ${outFile2} >> ${o2}");
    system("cat ${outFile3} >> ${o3}");
    $newfasta=$fasta."_".$j;
    unlink($newfasta);
    unlink($outFile1);
    unlink($outFile2);
    unlink($outFile3);
}

my $cddAllFile="/usr/global/blp/blast_db/cdd/cddid_all.tbl";
open CDD,"<$cddAllFile";
my %cdHash;
my $pssmID;
while (my $cddline=<CDD>) {
    chomp($cddline);
    my @columns=split /\t/,$cddline;
    next if ($columns[0] =~ /^0/ );
    $pssmID=shift @columns;
    $cdHash{$pssmID}=join "|",@columns;
}
close(CDD);

my @rpsOutFiles=($o1,$o2,$o3);
foreach my $rpsOutFile (@rpsOutFiles) {
    my $rpsOutFileSorted=$rpsOutFile.".sorted";
    my $rpsOutFileWithDesc=$rpsOutFile.".desc";
    my $evalueMin=1e-4;
    system("sort -k1,1 -k11,11g ${rpsOutFile} > ${rpsOutFileSorted}");
    rename $rpsOutFileSorted,$rpsOutFile;
    open RPSOUT,"<$rpsOutFile";
    open RPSDESC,">$rpsOutFileWithDesc";

    while (my $line=<RPSOUT>) {
        chomp($line);
	my @cols=split /\t/,$line;
	my $alnLen=$cols[3];
	my @pssm_cols=split /\|/,$cols[1];
	$pssmID=$pssm_cols[-1];
	my @cddcols=split /\|/,$cdHash{$pssmID};
	my $cdLen=$cddcols[-1];
	my $perc=100.0*$alnLen/$cdLen;
	if ($alnLen > $minAlnLen or $perc > $praThreshold) {
	    print RPSDESC "$line\t$cdHash{$pssmID}\n";        
	}
    }

    close(RPSOUT);
    close(RPSDESC);
    system("/usr/bin/pigz -p4 $rpsOutFile");``
}


my %rpsOutFilesHash;
$rpsOutFilesHash{"COG"}=${o1}.".desc";
$rpsOutFilesHash{"PFAM"}=${o2}.".desc";
$rpsOutFilesHash{"KOG"}=${o3}.".desc";

if (defined $templateGff) {

    open F1, "<$templateGff";
    my %templateGffHash;
    while (<F1>) {
        chomp;
        my @tmpLine = split /\t/;
        $templateGffHash{$tmpLine[8]}=$_;
    }
    close(F1);

    foreach my $cd (keys %rpsOutFilesHash) {
	my $cdgff=$project.".".$cd.".gff";
	open F2,"<$rpsOutFilesHash{$cd}";
	open O2, ">$cdgff";
	while (<F2>) {
	    chomp;
	    my @cdLine = split /\t/;
	    die "Incorrect template gff file or gene nomenclature not consistent!\n" unless (exists $templateGffHash{$cdLine[0]});
	    my @gffLine = split /\t/,$templateGffHash{$cdLine[0]};
	    $gffLine[1] = "RPSBLAST";
	    $gffLine[2] = $cd;
	    my $start = $cdLine[6];
	    my $end = $cdLine[7];
	    if ($start > $end) {
	    	$gffLine[6]="-";
		$start = $cdLine[7];
		$end = $cdLine[6];
	    } else {
		$gffLine[6]="+";
	    }
	    $gffLine[3] = $start;
	    $gffLine[4] = $end;
	    $gffLine[5] = $cdLine[2];
	    $gffLine[8] = $cdLine[12];
	    my $newGffLine = join ("\t",@gffLine);
	    printf O2 "%s\n",$newGffLine;
	}
	close(F2);
        close(O2);
    }
}

if (defined $summaryFile) {
my $ncog=`wc -l ${o1}.desc|awk '{print \$1}'`;
my $ncoguniq=`cut -f 1 ${o1}.desc|uniq|wc -l|awk '{print \$1}'`;
chomp $ncog;
my $npfam=`wc -l ${o2}.desc|awk '{print \$1}'`;
my $npfamuniq=`cut -f 1 ${o2}.desc|uniq|wc -l|awk '{print \$1}'`;
chomp $npfam;
my $nkog=`wc -l ${o3}.desc|awk '{print \$1}'`;
my $nkoguniq=`cut -f 1 ${o3}.desc|uniq|wc -l|awk '{print \$1}'`;
chomp $nkog;
open SF,">$summaryFile";
print SF "Number of COGs hit              \t\t\t:$ncog\n";
print SF "No. of genes with 1 or more COG \t\t\t:$ncoguniq";
print SF "GffFile                         \t\t\t:${project}.COG.gff\n" if (defined $templateGff);
print SF "For details see files           \t\t\t:${o1}.desc\n\n";
print SF "Number of PFAMs hit             \t\t\t:$npfam\n";
print SF "No. of genes with 1 or more PFAM\t\t\t:$npfamuniq";
print SF "GffFile                         \t\t\t:${project}.PFAM.gff\n" if (defined $templateGff);
print SF "For details see files           \t\t\t:${o2}.desc\n\n";
print SF "Number of KOGs hit              \t\t\t:$nkog\n";
print SF "No. of genes with 1 or more KOG \t\t\t:$nkoguniq";
print SF "GffFile                         \t\t\t:${project}.KOG.gff\n" if (defined $templateGff);
print SF "For details see files           \t\t\t:${o3}.desc\n\n";
close(SF);
}
exit;


sub usage() 
{
print <<EOF;
Run rpsblast on translated genes against COG and PFAM

Developer : Vishal N Koparde, Ph. D.
Created   : 120724
Modified  : 120809
Version   : 2.0

options:
--fasta    genes in fasta format
--project  project name
--wd       working dir (default=current directory)
--summary  summary file
--templateGff   gff file to use as template (optional)

EOF
exit 1;
}
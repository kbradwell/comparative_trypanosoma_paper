#!/usr/bin/perl
#Developer: Vishal Koparde, Ph.D.
#Created: 120709
#Modified:120709
#Version 1.0

use strict;
use warnings;
use Getopt::Long;
use lib qw(/usr/global/blp/perllib);
use Qsub;
use File::Basename;
use File::Copy;
use Cwd;

Getopt::Long::Configure ("bundling");

our ($gapdir,$revision);
if (exists $ENV{'GAPDIR'}) {
$gapdir=$ENV{'GAPDIR'};
} else {
$gapdir="/usr/global/blp/GenomeAnnotationPipeline";
}

our $python27="/home/vnkoparde/opt/Python-2.7.2/bin/python";
our $asgardDB="/gpfs_fs/data1/refdb/asgardDB";

our $cwd;
$cwd=`pwd`;
chomp($cwd);
chdir($gapdir);
#$revision=`/usr/bin/svn info|grep Revision|awk '{print \$2}'`;
#$revision=`/usr/bin/git log -1 --format="%cD"`; # this only gives last commit date ... not really revision like svn
if (-r "${gapdir}/version") {
	$revision=`cat ${gapdir}/version`;chomp $revision;
} else {
	$revision="unknown";
}
chdir($cwd);
chomp $revision;

our $commandLineArguments=join(" ",@ARGV);

our ($date_started,$date_ended);
$date_started=`date`;
chomp($date_started);

sub usage();

usage() if (scalar @ARGV==0);

our ($all, $callGenes, $geneCaller, $organismType, $fastaFile, $alienHunter, $outdir, $shortName, $force);
our ($trnaScan, $rpsblast, $blastx, $blastxdb, $rnammer, $doAsgard);
our $asgard = "/usr/global/blp/bin/asgard";

GetOptions ( 'all' => \$all,
             'call_genes' => \$callGenes,
             'gene_caller=i' => \$geneCaller,
             'rpsblast'=> \$rpsblast,
             'trna_scan' => \$trnaScan,
             'org=i' => \$organismType,
             'fasta=s' => \$fastaFile,
             'outdir=s' => \$outdir,
             'short_name=s' => \$shortName,
             'blastx' => \$blastx,
             'blastxdb=s' => \$blastxdb,
	     'rnammer' => \$rnammer,
	     'asgard' => \$doAsgard,
             'ah' => \$alienHunter,
             'force' => \$force);

die "Assembly fasta file required!\n" unless(defined $fastaFile);
die "File  $fastaFile Doesn't Exist!\n" unless (-e $fastaFile);

$geneCaller=2 unless (defined $geneCaller); #default geneCaller is GeneMark
die "Invalid GeneCaller!\n" if ($geneCaller!=1 && $geneCaller!=2); 

$organismType=2 unless (defined $organismType);
die "Invalid Organism Type!\n" if ($organismType!=1 && $organismType!=2);

die "Short organism name is required!\n" unless(defined $shortName);

if (defined $all) {
    print "will do all\n";
    $callGenes=1;
    $trnaScan=1;
    $rpsblast=1;
    $blastx=1;
    $rnammer=1;
    $doAsgard=1;
}

if (defined $blastx) {
    $blastxdb="/gpfs_fs/data1/refdb/nr" unless (defined $blastxdb);
}

our ($fastaSorted, $workingDir);
die "Output folder is required!\n" unless defined $outdir;
chomp $outdir;

$workingDir=$outdir;
die "$workingDir is not an absolute path! Complete path to output folder is required!\n" unless ( $workingDir =~ /^\// );
chop($workingDir) if ( $workingDir =~ /\/$/ );
if (! -d "${workingDir}") {
	mkdir("${workingDir}") || die "Unable to create folder ${workingDir}\n";
}

my ($fastaFileName,$fastaFilePath,$fastaFileExt);
($fastaFileName,$fastaFilePath,$fastaFileExt) = fileparse($fastaFile,qr"\..[^.]*$");
my $fastatmp = "${workingDir}/${fastaFileName}.fasta.tmp";
$fastaSorted = "${workingDir}/${fastaFileName}.sorted.fasta";
open CMD, ">${workingDir}/gap.commands";

# SORT AND CHANGE SEQIDS IN FASTA FILE

if ((defined $force) or (! -e $fastaSorted)) {
    my $cmd="$python27 ${gapdir}/fasta_sortnrename.py -i $fastaFile -o $fastaSorted -p \"${shortName}_contig\"";
    print CMD "$cmd\n";
    system($cmd);
}    
our ($genesFasta, $genesFaa, $genesGff, %jobs, %summaryFiles, $geneCallLastJob, $statusFile);
$statusFile = "${workingDir}/${shortName}.STATUS";
$genesFasta = "${workingDir}/${shortName}_genes.fasta";
$genesFaa = "${workingDir}/${shortName}_genes.faa";
$genesGff = "${workingDir}/${shortName}_genes.gff";

our ($job_gl, $job_gm, $job_pgc, $job_gcs);
our ($job_trna);
our ($job_rps);
our ($job_ah);


$summaryFiles{"GeneCalling"}="${workingDir}/${shortName}_GC.summary";
$summaryFiles{"tRNAScan"}="${workingDir}/${shortName}_tRNA.summary";
$summaryFiles{"rpsBlast"}="${workingDir}/${shortName}_rpsblast.summary";
$summaryFiles{"Blastx_contigs"}="${workingDir}/blastxC/${shortName}_blastx_contigs.summary";
$summaryFiles{"Blastx_genes"}="${workingDir}/blastxG/${shortName}_blastx_genes.summary";
$summaryFiles{"rnammer"}="${workingDir}/rnammer/${shortName}.RNAMMER.summary";




#           _ _  ____                      
#  ___ __ _| | |/ ___| ___ _ __   ___  ___ 
# / __/ _` | | | |  _ / _ \ '_ \ / _ \/ __|
#| (_| (_| | | | |_| |  __/ | | |  __/\__ \
# \___\__,_|_|_|\____|\___|_| |_|\___||___/
# 
if (defined $callGenes)  {
    print "will call genes\n";
	my $gcSummaryFile=$summaryFiles{"GeneCalling"};
	my $waitForJobID;
    if ((! -e "${genesFasta}") or (! -e "${genesGff}") or (! -e "${genesFaa}")or (defined $force)) {
	unlink ("${gcSummaryFile}") if (-e "${gcSummaryFile}");
	my $genesFastaTmp = $genesFasta.".tmp.fasta";
        if ($geneCaller == 2) {
            #my $gmsn_genes_file="${workingDir}/${shortName}.gmsn_genes";
            my $GeneMark_cmd="${gapdir}/callGenesWithGeneMark.pl";
            $GeneMark_cmd.=" --sortedFasta $fastaSorted";
            $GeneMark_cmd.=" --prefix $shortName";
            $GeneMark_cmd.=" --wd $workingDir";
	    $GeneMark_cmd.=" --force" if (defined $force);
            my $jname = "GM_".$shortName;
	    print CMD "$GeneMark_cmd\n";
            $job_gm=new Qsub(name=>$jname,wd=>$workingDir,outfile=>"GeneMark.tmp.out",cmd=>$GeneMark_cmd);
            $job_gm->submit();
	    $waitForJobID=$job_gm->{jobid};
#            push @jobs,$job_gm;
        }
	if ($geneCaller == 1) {
	    #my $coordsFile="${workingDir}/${shortName}.glimmercoords";
	    my $Glimmer_cmd="/usr/bin/perl ${gapdir}/callGenesWithGlimmer.pl";
            $Glimmer_cmd.=" --sortedFasta $fastaSorted";
            $Glimmer_cmd.=" --prefix $shortName";
            $Glimmer_cmd.=" --wd $workingDir";
	    $Glimmer_cmd.=" --force" if (defined $force);
	    my $jname="GL_".$shortName;
	    print CMD "$Glimmer_cmd\n";
	    $job_gl=new Qsub(name=>$jname,wd=>$workingDir,outfile=>"Glimmer.tmp.out",cmd=>$Glimmer_cmd);
	    $job_gl->submit();
	    $waitForJobID=$job_gl->{jobid};
#	    push @jobs,$job_gl;
	}
	#if ((! -e "${genesFasta}") or (! -e "${genesGff}") or (! -e "${genesFaa}")) {
	#    print STDERR "Gene calling failed ! $genesFasta or $genesGff or $genesFaa is (are) missing.\n";
	#    exit;
	#}
    }
    
    if ((! -e "${gcSummaryFile}") or (defined $force)) {
	my $jname="GCS_".$shortName;
	my $geneStats_cmd="$python27 ${gapdir}/calculate_gene_stats.py --assemblyFasta $fastaSorted --genesFasta $genesFasta --genesFaa $genesFaa --genesGFF $genesGff -o $gcSummaryFile";
	print CMD "$geneStats_cmd\n";
	if (defined $waitForJobID) {
		$job_gcs=new Qsub(name=>$jname,wd=>$workingDir,outfile=>"GeneStats.tmp.out",cmd=>$geneStats_cmd,waitfor=>$waitForJobID);
	} else {
		$job_gcs=new Qsub(name=>$jname,wd=>$workingDir,outfile=>"GeneStats.tmp.out",cmd=>$geneStats_cmd);
	}
	$job_gcs->submit();
	$jobs{"GeneCalling"}=$job_gcs;
    }
}
#$geneCallLastJob->waitForCompletion();


# _   ____  _   _    _    ____                  
#| |_|  _ \| \ | |  / \  / ___|  ___ __ _ _ __  
#| __| |_) |  \| | / _ \ \___ \ / __/ _` | '_ \ 
#| |_|  _ <| |\  |/ ___ \ ___) | (_| (_| | | | |
# \__|_| \_\_| \_/_/   \_\____/ \___\__,_|_| |_|
# 
if (defined $trnaScan) {
    print "will run trnaScan-SE\n";
    my $trnaSummaryFile=$summaryFiles{"tRNAScan"};
    my $trnaGffFile="${workingDir}/${shortName}_tRNA.gff3";
    my $trnaFastaFile="${workingDir}/${shortName}_tRNA.fasta";
    if (((! -e "$trnaSummaryFile") or (! -e "$trnaGffFile")) or (! -e "$trnaFastaFile") or (defined $force)) {
	my $trnaScan_cmd="/usr/bin/perl ${gapdir}/trnaScan.pl";
    	$trnaScan_cmd.=" --project $shortName";
    	$trnaScan_cmd.=" --wd $workingDir";
    	$trnaScan_cmd.=" --fasta $fastaSorted";
    	$trnaScan_cmd.=" --org $organismType";
    	$trnaScan_cmd.=" --summary $trnaSummaryFile";
    	$trnaScan_cmd.=" --gff $trnaGffFile";
	$trnaScan_cmd.=" --outFasta $trnaFastaFile";
    	print CMD "$trnaScan_cmd\n";
    	$job_trna=new Qsub(name=>"tRNA_".$shortName,outfile=>"tRNA_".$shortName.".tmp.out",wd=>$workingDir,cmd=>$trnaScan_cmd);
    	$job_trna->submit();
    	$jobs{"tRNAScan"}=$job_trna;
    }
}
#                ____  _           _   
# _ __ _ __  ___| __ )| | __ _ ___| |_ 
#| '__| '_ \/ __|  _ \| |/ _` / __| __|
#| |  | |_) \__ \ |_) | | (_| \__ \ |_ 
#|_|  | .__/|___/____/|_|\__,_|___/\__|
# 
if (defined $rpsblast) {
    print "will run rbsblast\n";
    die "Have you called Genes!\n Could not find ${genesFasta} or ${genesGff}\n" if ((! $callGenes) and ((! -e "${genesFasta}") or (! -e "${genesGff}")));
    if ( (!-e "${workingDir}/${shortName}.COG.rps.gz") or
     (!-e "${workingDir}/${shortName}.PFAM.rps.gz") or
     (!-e "${workingDir}/${shortName}.KOG.rps.gz") or
     (!-e "${workingDir}/${shortName}.COG.gff") or
     (!-e "${workingDir}/${shortName}.PFAM.rps") or
     (!-e "${workingDir}/${shortName}.KOG.rps") or
     (!-e "${workingDir}/${shortName}_rpsblast.summary") or (defined $force)) {
	my $rps_cmd="/usr/bin/perl ${gapdir}/rps_cog_pfam.pl ";
        $rps_cmd.=" --project ${shortName}";
	$rps_cmd.=" --wd ${workingDir}";
	$rps_cmd.=" --fasta ${genesFasta}";
	$rps_cmd.=" --templateGff ${genesGff}";
	$rps_cmd.=" --summary ${workingDir}/${shortName}_rpsblast.summary";
	print CMD "$rps_cmd\n";
	$job_rps=new Qsub(name=>"RPS_".$shortName,outfile=>"RPS_".$shortName."tmp.out",wd=>$workingDir,cmd=>$rps_cmd);
	$job_rps->{waitfor}=($jobs{"GeneCalling"})->{jobid} if ($callGenes);
	$job_rps->submit();
	$jobs{"rpsBlast"}=$job_rps;
    }
}
# _     _           _        
#| |__ | | __ _ ___| |___  __
#| '_ \| |/ _` / __| __\ \/ /
#| |_) | | (_| \__ \ |_ >  < 
#|_.__/|_|\__,_|___/\__/_/\_\
# 
if (defined $blastx) {
    print "will run blastx\n";

    ($fastaFileName,$fastaFilePath,$fastaFileExt) = fileparse($fastaSorted,qr"\..[^.]*$");
    my $blastxContigsOutFile="${workingDir}/blastxC/blastres_${fastaFileName}.fasta";
    my $blastxContigsInFasta="${workingDir}/blastxC/${fastaFileName}.fasta";
    ($fastaFileName,$fastaFilePath,$fastaFileExt) = fileparse($genesFasta,qr"\..[^.]*$");
    my $blastxGenesOutFile="${workingDir}/blastxG/blastres_${fastaFileName}.fasta";
    my $blastxGenesInFasta="${workingDir}/blastxG/${fastaFileName}.fasta";

    my $blastxContigsGff;
    my $blastxGenesGff;
    if ($blastxdb eq "/data/refdb/nr") {
        $blastxContigsGff="${workingDir}/blastxC/blastx_contigs_nr.gff";
        $blastxGenesGff="${workingDir}/blastxG/blastx_genes_nr.gff";
    } else {
	$blastxContigsGff="${workingDir}/blastxC/blastx_contigs.gff";
	$blastxGenesGff="${workingDir}/blastxG/blastx_genes.gff";
    }

    mkdir("${workingDir}/blastxC") if (! -d "${workingDir}/blastxC");
    mkdir("${workingDir}/blastxG") if (! -d "${workingDir}/blastxG");

    if ((defined $force) or ((! -e $blastxContigsOutFile) or (! -e $blastxContigsGff))) {
	unlink("$blastxContigsOutFile") if (! -e "$blastxContigsOutFile");
	unlink("$blastxContigsGff") if (! -e "$blastxContigsGff");
	symlink("$fastaSorted","$blastxContigsInFasta");
	my $tmpSummaryBlastxContigFile=$summaryFiles{"Blastx_contigs"};
	my $blastxC_cmd="/usr/bin/perl ${gapdir}/blastx_contigs_asgard.pl ";
	$blastxC_cmd.=" --fasta $blastxContigsInFasta";
	$blastxC_cmd.=" --project ${shortName}";
	$blastxC_cmd.=" --db ${blastxdb}";
	$blastxC_cmd.=" --wd ${workingDir}/blastxC";
	$blastxC_cmd.=" --summary ${tmpSummaryBlastxContigFile}";
    	$blastxC_cmd.=" -gff ${blastxContigsGff}";
	print CMD "${blastxC_cmd}\n";
	my $job_blastx=new Qsub(name=>"BxC_".$shortName,outFile=>"BxC.tmp.out",wd=>"${workingDir}/blastxC",cmd=>$blastxC_cmd);
	$job_blastx->submit();
	$jobs{"BlastxContigs"}=$job_blastx;
    }
    if ((defined $force) or ((! -e $blastxGenesOutFile) or (! -e $blastxGenesGff))) {
	
	my $tmpSummaryBlastxGenesFile=$summaryFiles{"Blastx_genes"};
	if ( (! -e $genesFasta) and (! exists $jobs{"GeneCalling"}) ) {
		print "Genes file does not exist! Call Genes first!\n";
	} 
	symlink("$genesFasta","$blastxGenesInFasta");
	my $blastxG_cmd="/usr/bin/perl ${gapdir}/blastx_genes_asgard.pl ";
	$blastxG_cmd.=" --fasta $blastxGenesInFasta";
	$blastxG_cmd.=" --project ${shortName}";
	$blastxG_cmd.=" --db ${blastxdb}";
	$blastxG_cmd.=" --wd ${workingDir}/blastxG";
	$blastxG_cmd.=" --summary ${tmpSummaryBlastxGenesFile}";
    	$blastxG_cmd.=" -gff ${blastxGenesGff}";
	my $job_blastx;
        if (exists $jobs{"GeneCalling"}) {
		print "\n\n$blastxG_cmd\n\n";
		print "\n\nwaiting!\n\n";
        	my $waitforid=$jobs{"GeneCalling"}->{jobid};
		print CMD "$blastxG_cmd\n";
		$job_blastx=new Qsub(name=>"BxG_".$shortName,outFile=>"BxG.tmp.out",wd=>"${workingDir}/blastxG",cmd=>$blastxG_cmd,waitfor=>$waitforid);
		$job_blastx->submit();
		$jobs{"BlastxGenes"}=$job_blastx;
	} else {
		if ( -e $genesFasta ) {
			print "\n\n$blastxG_cmd\n\n";
			print "\n\nnot waiting and found $genesFasta\n\n";
			print CMD "$blastxG_cmd\n";
			$job_blastx=new Qsub(name=>"BxG_".$shortName,outFile=>"BxG.tmp.out",wd=>"${workingDir}/blastxG",cmd=>$blastxG_cmd);
			$job_blastx->submit();
			$jobs{"BlastxGenes"}=$job_blastx;
		} else {
			print "Blastx will not be performed on genes as genes fasta does not exist and is not being created!\n";
		}
	}
    }
}

#                               _ 
#  __ _ ___  __ _  __ _ _ __ __| |
# / _` / __|/ _` |/ _` | '__/ _` |
#| (_| \__ \ (_| | (_| | | | (_| |
# \__,_|___/\__, |\__,_|_|  \__,_|
#           |___/

if (defined $doAsgard) {
    print "will run asgard\n";
    mkdir("${workingDir}/asgard_contigs") if (! -d "${workingDir}/asgard_contigs");
    ($fastaFileName,$fastaFilePath,$fastaFileExt) = fileparse($fastaSorted,qr"\..[^.]*$");
    my $asgardContigsInputFasta="${workingDir}/asgard_contigs/${fastaFileName}.fasta";
    my $asgardContigsInputFastaNoPath="${fastaFileName}.fasta";
    system("rm -f $asgardContigsInputFasta") if -e $asgardContigsInputFasta;
    symlink("$fastaSorted","$asgardContigsInputFasta") or die "SymLink of $fastaSorted into ${workingDir}/asgard_contigs failed: $!";
    #copy("${fastaSorted}","${workingDir}/asgard_contigs/${shortName}_contigs.fasta") or die "Copy failed: $!";
    my $asgard_cmd1="$asgard -i $asgardContigsInputFastaNoPath -p blastx -d ${asgardDB}/UniRef100 -d ${asgardDB}/KEGG -f ${asgardDB}/uniref100.fasta.gz -f ${asgardDB}/genes.pep.gz -l ${asgardDB}/";
    print CMD "$asgard_cmd1\n";
    my $job_asgard1=new Qsub(name=>"ASG_C".$shortName,outFile=>"ASG_C.tmp.out",wd=>"${workingDir}/asgard_contigs",cmd=>$asgard_cmd1);
    $job_asgard1->submit();
    $jobs{"AsgardContigs"}=$job_asgard1;
    mkdir("${workingDir}/asgard_genes") if (! -d "${workingDir}/asgard_genes");
    ($fastaFileName,$fastaFilePath,$fastaFileExt) = fileparse($genesFasta,qr"\..[^.]*$");
    my $asgardGenesInputFasta="${workingDir}/asgard_genes/${fastaFileName}.fasta";
    my $asgardGenesInputFastaNoPath="${fastaFileName}.fasta";
    system("rm -f $asgardGenesInputFasta") if -e $asgardGenesInputFasta;
    symlink("$genesFasta","$asgardGenesInputFasta") or die "SymLink of $genesFasta into ${workingDir}/asgard_genes failed: $!";;
    #copy("${genesFasta}","${workingDir}/asgard_contigs/${shortName}_genes.fasta") or die "Copy failed: $!";
    my $asgard_cmd2="$asgard -i $asgardGenesInputFastaNoPath -p blastx -d ${asgardDB}/UniRef100 -d ${asgardDB}/KEGG -f ${asgardDB}/uniref100.fasta.gz -f ${asgardDB}/genes.pep.gz -l ${asgardDB}/";
    print CMD "$asgard_cmd2\n";
    my $job_asgard2;
	if (exists $jobs{"GeneCalling"}) {
    my $waitforid=$jobs{"GeneCalling"}->{jobid};
    $job_asgard2=new Qsub(name=>"ASG_G".$shortName,outFile=>"ASG_G.tmp.out",wd=>"${workingDir}/asgard_genes",cmd=>$asgard_cmd2,waitfor=>$waitforid);
    $job_asgard2->submit();
    $jobs{"AsgardGenes"}=$job_asgard2;
	} else {
		if ( -e $genesFasta ) {
    $job_asgard2=new Qsub(name=>"ASG_G".$shortName,outFile=>"ASG_G.tmp.out",wd=>"${workingDir}/asgard_genes",cmd=>$asgard_cmd2);
    $job_asgard2->submit();
    $jobs{"AsgardGenes"}=$job_asgard2;
		} else {
    print "ASGARD will not be performed on genes as genes fasta does not exist and is not being created!\n";
		}
	}
}

#                                                
# _ __ _ __   __ _ _ __ ___  _ __ ___   ___ _ __ 
#| '__| '_ \ / _` | '_ ` _ \| '_ ` _ \ / _ \ '__|
#| |  | | | | (_| | | | | | | | | | | |  __/ |   
#|_|  |_| |_|\__,_|_| |_| |_|_| |_| |_|\___|_|   
# 

if (defined $rnammer) {
    print "will run rnammer\n";
    my $rnammerGff="${workingDir}/rnammer/${shortName}_rna.gff";
    my $rnammerFasta="${workingDir}/rnammer/${shortName}_rna.fasta";
    my $rnammerSINA="${workingDir}/rnammer/${shortName}_rna_sina_classification.txt";
    my $rnammerSummary=$summaryFiles{"rnammer"};
    my $rnammerKingdom;
    
    $rnammerKingdom="bac" if ($organismType==1);
    $rnammerKingdom="euk" if ($organismType==2);
    if ((! -e $rnammerGff) or (! -e $rnammerFasta) or (! -e $rnammerSummary) or (defined $force)) {
	system("rm -rf ${workingDir}/rnammer") if (-d "${workingDir}/rnammer");
	mkdir("${workingDir}/rnammer");
	($fastaFileName,$fastaFilePath,$fastaFileExt) = fileparse($fastaSorted,qr"\..[^.]*$");
	my $rnammerInputFasta="${workingDir}/rnammer/${fastaFileName}.fasta";
	symlink("$fastaSorted","$rnammerInputFasta") or die "SymLink of $fastaSorted into ${workingDir}/rnammer failed: $!";
	my $rnammer_cmd="$python27 ${gapdir}/rnammer_w_sina_classification.py $rnammerKingdom $rnammerGff $rnammerFasta $rnammerInputFasta $rnammerSummary $rnammerSINA";
	print CMD "$rnammer_cmd\n";
	my $job_rnammer=new Qsub(name=>"rnammer_".$shortName,outFile=>"rnammer.tmp.out",wd=>"${workingDir}/rnammer",cmd=>$rnammer_cmd,memory=>"10G",nproc=>"4");
	$job_rnammer->submit();
	$jobs{"rnammer"}=$job_rnammer;
    }
}

#       _ _            _   _             _            
#  __ _| (_) ___ _ __ | | | |_   _ _ __ | |_ ___ _ __ 
# / _` | | |/ _ \ '_ \| |_| | | | | '_ \| __/ _ \ '__|
#| (_| | | |  __/ | | |  _  | |_| | | | | ||  __/ |   
# \__,_|_|_|\___|_| |_|_| |_|\__,_|_| |_|\__\___|_|   
# 
if (defined $alienHunter) {
    print "will do alien hunter\n";
    my $alienFile = "${shortName}.alien";
    if ((! -e "${workingDir}/alien_hunter/${alienFile}") or (defined $force)) {
	system("rm -rf ${workingDir}/alien_hunter") if (-d "${workingDir}/alien_hunter");
	mkdir("${workingDir}/alien_hunter");
    	#system("cp -r /usr/global/blp/alien_hunter ${workingDir}/alien_hunter");
	($fastaFileName,$fastaFilePath,$fastaFileExt) = fileparse($fastaSorted,qr"\..[^.]*$");
	my $alienHunterInputFasta="${workingDir}/alien_hunter/${fastaFileName}.fasta";
    	copy("$fastaSorted","$alienHunterInputFasta");
    	my $AlienHunter_cmd="/usr/global/blp/bin/alien_hunter $alienHunterInputFasta $alienFile";
    	my $jname = "AH_".$shortName;
   	print CMD "AlienHunter_cmd\n";
    	$job_ah=new Qsub(name=>$jname,wd=>"${workingDir}/alien_hunter",outfile=>"AlienHunter.tmp.out",cmd=>$AlienHunter_cmd);
    	$job_ah->submit();
	$jobs{"AlienHunter"}=$job_ah; 
    }   
}

while (checkAllJobStatus(\%jobs) ne "Completed") {
    printJobStatus(\%jobs,$statusFile);
    sleep(10);
}
printJobStatus(\%jobs,$statusFile);

$date_ended=`date`;
chomp($date_ended);
#open SF,">",$summaryFile;

my $summaryFile="${workingDir}/GAP.summary";
system("cat ${gapdir}/figlets/GenomeAnnotationPipeline.fig > ${summaryFile}");
open SF, ">>${summaryFile}";
print SF "Version:\t${revision}\n";
print SF "Start:\t${date_started}\n";
print SF "End:\t${date_ended}\n";
print SF "Command Line Arguments:\t${commandLineArguments}\n";
close SF;
my $tmpSummary;
if (exists $summaryFiles{"GeneCalling"}) {
    $tmpSummary=$summaryFiles{"GeneCalling"};
    if (-e $tmpSummary) {
	system ("cat ${gapdir}/figlets/GeneCalling.fig >> ${summaryFile}");
	system ("cat ${tmpSummary} >> ${summaryFile}")
    };
}
if (exists $summaryFiles{"tRNAScan"}) {
    $tmpSummary=$summaryFiles{"tRNAScan"};
    if (-e $tmpSummary){
        system ("cat ${gapdir}/figlets/tRNAScan.fig >> ${summaryFile}");
        system ("cat ${tmpSummary} >> ${summaryFile}");
    }
}
if (exists $summaryFiles{"rpsBlast"}) {
    $tmpSummary=$summaryFiles{"rpsBlast"};
    if (-e $tmpSummary) {
        system ("cat ${gapdir}/figlets/rpsBlast.fig >> ${summaryFile}");
        system ("cat ${tmpSummary} >> ${summaryFile}") ;
    }
}
if (exists $summaryFiles{"BlastxContigs"}) {
    $tmpSummary=$summaryFiles{"BlastxContigs"};
    my $tmpSummary2=$summaryFiles{"BlastxGenes"};
    if (-e $tmpSummary) {
        system ("cat ${gapdir}/figlets/Blastx.fig >> ${summaryFile}");
        system ("cat ${tmpSummary} >> ${summaryFile}") ;
    }
    if (-e $tmpSummary2) {
	system ("cat ${tmpSummary2} >> ${summaryFile}");
    }
}
if (exists $summaryFiles{"rnammer"}) {
    $tmpSummary=$summaryFiles{"rnammer"};
    if (-e $tmpSummary) {
        system ("cat ${gapdir}/figlets/rnammer.fig >> ${summaryFile}");
        system ("cat ${tmpSummary} >> ${summaryFile}") ;
    }
}




exit;











sub printJobStatus {
my ($jobs_hash,$statusFile) = @_;
my %jobs=%$jobs_hash;
open SFILE,">", $statusFile;
my $date=`date`;
print SFILE "\n$date";
print "\n$date";
#foreach my $jobname (sort { ($jobs{$a}->{jobid}) <=> ($jobs{$b}->{jobid}) } keys %jobs) {
foreach my $jobname (keys %jobs) {
	my $job=$jobs{$jobname};
	my $jobid=$job->{jobid};
	my $jobStatus=$job->status();
	print SFILE "$jobname\t$jobid\t$jobStatus\n";
	print "$jobname\t$jobid\t$jobStatus\n";
}
close(SFILE);
return "Completed";
}

sub checkAllJobStatus {
my ($jobs_hash) = @_;
my %jobs = %$jobs_hash;
foreach my $jobname (keys %jobs) {
	my $job=$jobs{$jobname};
	return "Running" if ($job->status() ne "Completed");
}
return "Completed";
}

sub usage() 
{
print <<EOF;
Genome Annotation Pipeline

Developer : Vishal N Koparde, Ph. D.
Created   : 120709
Modified  : 150403
EOF
print "Version   : ",$revision."\n";
print <<EOF;

options:
--outdir         Output directory path (REQUIRED)
--fasta          Assembly in fasta format (REQUIRED)
--short_name     Short Name for the organism (REQUIRED)
--org            Type of organism, i.e., 1=Bacteria 2=Eukaryote(default)
--call_genes     Call genes using Glimmer or GeneMark
--gene_caller    Specify gene caller, i.e., 1=Glimmer 2=GeneMark(default)
--trna_scan      Run trnaScan-SE
--rpsblast       rpsblast against COG/PFAM (requires called genes)
--blastx         blast assembly against protein database (NR by default)
--blastxdb       full path of protein database (NR by default)
--rnammer        find RNAs using rnammer-1.2
--asgard         perform metabolic reconstruction using ASGARD
--ah             Run Alien Hunter (Warning:Will take long time)
--all            Run everything, except alien hunter
--force          Overwrite files
EOF
exit 1;
}

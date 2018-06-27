# comparative_trypanosoma_paper
Custom scripts for recreating some of the analyses from the trypanosoma genomic comparisons work.

[1. Genome Completion Pipeline](#1-genome-completion-pipeline)  
[2. Genome Annotation Pipeline](#2-genome-annotation-pipeline)  
[3. Edit Alignments for Phylogeny](#3-edit-alignments-for-phylogeny)  
[4. Percent Identity Calculation](#4-percent-identity-calculation)  

## 1. Genome Completion Pipeline
## 2. Genome Annotation Pipeline

Script name: gap.pl

Description: Pipeline for genome annotation.

Usage: `perl gap.pl`

Input  
--outdir         Output directory path (REQUIRED)
--fasta          Assembly in fasta format (REQUIRED)
--short_name     Short Name for the organism (REQUIRED)
--org            Type of organism, i.e., 1=Bacteria 2=Eukaryote(default)

Output  
Results for the following options are available:
--call_genes     Call genes using Glimmer or GeneMark
--trna_scan      Run trnaScan-SE
--rpsblast       rpsblast against COG/PFAM (requires called genes)
--blastx         blast assembly against protein database (NR by default)
--rnammer        find RNAs using rnammer-1.2
--asgard         perform metabolic reconstruction using ASGARD

Script name: get_best_annotated_hit.py (separate but complementary script to the Genome Annotation Pipeline)

Description: Searches a database of the best annotation available.

Usage: `python get_best_annotated_hit.py <nr_blastp_results> <outfile> > <logfile>`

Input  
- nr BLASTp results

Output  
- Best hits for each gene

## 3. Edit Alignments for Phylogeny
## 4. Percent Identity Calculation

Script name: get_perc_id_general.py

Description: Script for calculating percent identity from aligned nucleotide or amino acid sequences.

Usage: `python get_perc_id_general.py <fastafile>`

Input  
- fastafile: alignment in fasta format

Output  
- percent identity



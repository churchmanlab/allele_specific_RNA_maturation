#!/bin/bash

# This script will call a Python script to reformat the VCF file for use with the DRIMSeq tuQTL tool, for genes in which allelic splicing order could be analyzed

module load gcc/6.2.0
module load python/3.7.4
source ~/myenv_python3/bin/activate

baseDir="/path/to/projects/LCLs/directRNAseq"

vcf1="$baseDir/analysis_new/allSamples/LCLs_studied_cell_lines_all_SNPs.genes_splicing_order.extended.vcf" # all SNPs
vcf2="$baseDir/analysis_new/allSamples/LCLs_studied_cell_lines_het_SNPs_phased.genes_splicing_order.extended.vcf" # het SNPs phased by LORALS
out_vcf="$baseDir/analysis_new/allSamples/LCLs_studied_cell_lines.genes_splicing_order.extended.formatted_for_drimseq.txt"

my_script="$baseDir/scripts/23_process_vcf_for_drimseq.py"


# Command to run
python $my_script $vcf2 $vcf1 $out_vcf

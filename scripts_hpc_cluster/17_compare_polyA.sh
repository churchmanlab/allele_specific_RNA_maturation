#!/bin/bash

# This script will call a Python script that will compare poly(A) tail length and read end position between alleles

sample="GM18861_chr_rep1"


###################

# CONFIG
baseDir="/path/to/projects/LCLs/directRNAseq"
outDir="$baseDir/dorado/$sample"

# Input files
nano_df="$baseDir/dorado/$sample/$sample.polyA_dorado.txt"
readmap="$baseDir/analysis_new/$sample/${sample}_read_mapping_per_gene.RefSeq.txt"
allelemap="$baseDir/hapcut2/post_dorado_merged_vcf/$sample/${sample}_read_to_allele_mapping.txt"
read_ends="$baseDir/analysis_new/$sample/${sample}_read_ends.both_alleles.txt"

# Annotation files
hg38_intron_df="/path/to/reference_libraries/annotations/hg38/hg38_all_intron_features.txt"

# Output files
out_pA="$outDir/${sample}_polyA_comparison_between_alleles.txt"
out_summary="$outDir/${sample}_polyA_summary_measures.txt"
sig_read_pA="$outDir/${sample}_polyA_comparison_between_alleles.individual_reads.txt"
out_read_ends="$outDir/${sample}_read_ends_polyA_comparison.txt"
sig_read_ends="$outDir/${sample}_read_ends_polyA_comparison.individual_reads.txt"

# Python script
compare_py="$baseDir/scripts/17a_compare_polyA_tails.py"
summary_py="$baseDir/scripts/17b_summary_measures_polyA_tails.py"

######################

module load gcc/6.2.0
module load samtools/1.9
module load bedtools/2.27.1
module load python/3.7.4
source ~/myenv_python3/bin/activate

# Create output directory if non existent
if [ ! -d $outDir ]; then
    mkdir -p $outDir
fi

# Command
# Compare poly(A) tail lengths and reads end positions between alleles
python $compare_py $nano_df $readmap $allelemap $hg38_intron_df $out_pA $sig_read_pA $read_ends $out_read_ends $sig_read_ends
# Get summary measures for poly(A) tail length (e.g. average, CV) independent of allele
python $summary_py $nano_df $readmap $hg38_intron_df $out_summary



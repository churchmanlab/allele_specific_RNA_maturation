#!/bin/bash

# This script will call a Python script to compare poly(A) tail length between alleles as a function of the read splicing status

sample="GM18861_chr_rep1"

###################

# CONFIG
baseDir="/path/to/projects/LCLs/directRNAseq"
outDir="$baseDir/dorado/$sample"

# Input files
nano_df="$baseDir/dorado/$sample/$sample.polyA_dorado.txt"
readmap="$baseDir/analysis_new/$sample/${sample}_read_mapping_per_gene.RefSeq.txt"
allelemap="$baseDir/hapcut2/post_dorado_merged_vcf/$sample/${sample}_read_to_allele_mapping.txt"
multi_introns_df="$baseDir/analysis_new/$sample/${sample}_hg38_multi_introns_isoforms_df.withSKP.RefSeq_all.txt"

# Annotation files
hg38_intron_df="/path/to/reference_libraries/annotations/hg38/hg38_all_intron_features.txt"

# Output files
out_splice_pA="$outDir/${sample}_splicing_polyA_comparison.txt"
sig_read_splice_pA="$outDir/${sample}_splicing_polyA_comparison.individual_reads.txt"

# Python script
compare_py="$baseDir/scripts/26_compare_polyA_tails_splicing.py"

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
python $compare_py $nano_df $readmap $allelemap $multi_introns_df $hg38_intron_df $out_splice_pA $sig_read_splice_pA

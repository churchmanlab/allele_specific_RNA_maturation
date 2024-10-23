#!/bin/bash

# This script will merge multiple BAM files from the same individual


# CONFIG

sample1="GM19152_chr_rep1"
sample2="GM19152_chr_rep2"
merged_sample="GM19152_chr_merged"

# File locations
baseDir="/path/to/projects/LCLs/directRNAseq/alignment/dorado"
bam1="$baseDir/$sample1/minimap2/${sample1}_minimap2_uniq_sort.bam"
bam2="$baseDir/$sample2/minimap2/${sample2}_minimap2_uniq_sort.bam"


#Output
outDir="$baseDir/$merged_sample/minimap2"
out="$outDir/${merged_sample}_minimap2_uniq_sort"

#####################################################################

module load bedtools/2.26.0
module load samtools/1.3.1

#######################################
# Create output directory if non existent
if [ ! -d $outDir ]; then
    mkdir -p $outDir
fi


# Command to run
samtools merge $out.bam $bam1 $bam2
samtools index ${out}.bam


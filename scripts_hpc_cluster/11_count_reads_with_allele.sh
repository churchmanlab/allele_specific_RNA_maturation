#!/bin/bash

# This script will count the total number of reads per sample assigned to each parental allele

sample="GM18861_chr_rep1"

###################

# CONFIG
baseDir="/path/to/projects/LCLs/directRNAseq"
outDir="$baseDir/hapcut2/post_dorado_merged_vcf/$sample"

# Input files
bed="$baseDir/lorals/alignment_post_dorado_merged_vcf/$sample/${sample}.dorado_reads_aln_sorted.merged.bam"
allelemap="$baseDir/hapcut2/post_dorado_merged_vcf/$sample/${sample}_read_to_allele_mapping.txt"

# Output files
out_stats="$outDir/${sample}_read_to_allele_stats.txt"


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

# Command to run
alleles="M P"
for i in $alleles; do
  echo $i `grep $i $allelemap | wc -l | awk '{print $1}'` >> $out_stats
done

#!/bin/bash

# This script will split a BAM file into two BAM files that each contain only the reads mapping to each parental allele


# CONFIG

sample="GM18861_chr_rep1"

# File locations
baseDir="/path/to/projects/LCLs/directRNAseq/lorals/alignment_post_dorado_merged_vcf"
bam="$baseDir/$sample/${sample}.dorado_reads_aln_sorted.merged.bam"

#Output
outDir="/path/to/projects/LCLs/directRNAseq/hapcut2/post_dorado_merged_vcf/${sample}"
out_mat="$baseDir/$sample/${sample}_RNAtoDNA_reads_aln_sorted.maternal.bam"
out_pat="$baseDir/$sample/${sample}_RNAtoDNA_reads_aln_sorted.paternal.bam"


# END OF CONFIG
#####################################################################

module load bedtools/2.26.0
module load samtools/1.13

#######################################
# Create output directory if non existent
if [ ! -d $outDir ]; then
    mkdir -p $outDir
fi


# Command to run
# Extract reads mapping to each parental allele
awk '$2 == "M" {print $1}' $outDir/${sample}_read_to_allele_mapping.txt > $outDir/${sample}_read_to_allele_mapping.maternal_only.txt
awk '$2 == "P" {print $1}' $outDir/${sample}_read_to_allele_mapping.txt > $outDir/${sample}_read_to_allele_mapping.paternal_only.txt

# Extract corresponding reads from BAM file and index
samtools view -b -N $outDir/${sample}_read_to_allele_mapping.maternal_only.txt $bam > $out_mat
samtools index $out_mat
samtools view -b -N $outDir/${sample}_read_to_allele_mapping.paternal_only.txt $bam > $out_pat 
samtools index $out_pat

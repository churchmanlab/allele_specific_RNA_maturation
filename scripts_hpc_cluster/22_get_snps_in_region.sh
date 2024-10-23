#!/bin/bash

# This script will extract all SNPs located in genes for which we were able to analyze allelic splicing order, for downstream analyses

module load gcc/6.2.0
module load python/3.7.4
module load samtools/1.3.1
module load bcftools/1.13
module load bedtools/2.26.0

# CONFIG
baseDir="/path/to/projects/LCLs/directRNAseq"
vcf1="$baseDir/lorals/vcf/LCLs_studied_cell_lines_all_SNPs.vcf.gz" # all SNPs
vcf2="$baseDir/lorals/new_vcf_dorado/LCLs_studied_cell_lines_het_SNPs_phased.vcf.gz" # het SNPs from phased LORALS VCF

# Annotation files
genes_of_interest="$baseDir/from_local/LCL_genes_splicing_order_analysis_pairs_and_triplets_2024-08-12.txt"
genes_bed="/path/to/reference_libraries/annotations/hg38/NCBI_RefSeq_hg38_genes_parsed.bed"

# Output files
outDir="$baseDir/analysis_new/allSamples"
mkdir -p $outDir
genes_bed_sub="$outDir/LCL_genes_splicing_order_analysis_pairs_and_triplets_2024-08-12"
out1="$outDir/LCLs_studied_cell_lines_all_SNPs.genes_splicing_order"
out2="$outDir/LCLs_studied_cell_lines_het_SNPs_phased.genes_splicing_order"


# Make BED file with the genes
grep -f $genes_of_interest $genes_bed > $genes_bed_sub.bed

# Extend regions by 1000 bp on each side
awk 'BEGIN{OFS="\t"}{print $1,$2-1000,$3+1000,$4,$5,$6}' $genes_bed_sub.bed > $genes_bed_sub.extended.bed

# Subset VCF files
bcftools view -H -O v -o $out1.extended.vcf.gz -R $genes_bed_sub.extended.bed $vcf1
bcftools view -H -O v -o $out2.extended.vcf.gz -R $genes_bed_sub.extended.bed $vcf2

# Get the gene to SNP correspondence
zcat $out1.extended.vcf.gz | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+1,$3}' | bedtools intersect -wo -a stdin -b $genes_bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$8}' > $out1.bed
zcat $out2.extended.vcf.gz | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+1,$3}' | bedtools intersect -wo -a stdin -b $genes_bed | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$8}' > $out2.bed

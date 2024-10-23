#!/bin/bash

#This script uses the LORALS script hap_aligner.sh to perform haplotype-aware alignment of long-read RNA-seq data

sample="GM18861_chr_rep1" # prefix of sequencing files for this sample
id="GM18861_chr_merged" # prefix of VCF

# Load necessary packages and virtual environment
module load gcc/6.2.0
module load samtools/1.15.1
module load bedtools/2.27.1
module load python/3.7.4
module load bcftools/1.13
module load gatk/4.1.9.0
source ~/myenv_python3/bin/activate

# CONFIG
baseDir="/path/to/projects/LCLs/directRNAseq"
fastq_dir="$baseDir/dorado/$sample"
fastq="$fastq_dir/${sample}.dorado.fastq"
fasta_path="$baseDir/lorals/new_vcf_dorado/$id/${id}_minimap2_uniq_sort"
outDir="/path/to/projects/LCLs/directRNAseq/lorals/alignment_post_dorado_merged_vcf/$sample"

# Software
hap_aligner="/path/to/software/lorals/scripts/hap_aligner.sh"

# Create output directory if non existent
if [ ! -d $outDir ]; then
    mkdir -p $outDir
fi


# Make sym links to find the correct fasta files
ln -s ${fasta_path}.hap1.fa ${fasta_path}_hap1.fa
ln -s ${fasta_path}.hap2.fa ${fasta_path}_hap2.fa

# Haplotype-aware alignment
$hap_aligner -f $fastq --reference $fasta_path -o $outDir

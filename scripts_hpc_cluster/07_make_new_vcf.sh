#!/bin/bash

# This script takes advantage of the LORALS make_new_vcf_sh script to correct phased genotypes based on aligned long-read RNA-seq data

id="19152" # suffix that goes after "NA" or "GM" in LCL name

# Load necessary packages and virtual environment
module load gcc/6.2.0
module load samtools/1.15.1
module load bedtools/2.27.1
module load python/3.7.4
module load bcftools/1.13
module load gatk/4.1.9.0
module load htslib/1.9.0
source ~/myenv_python3/bin/activate

# CONFIG
ref_genome="/path/to/genomes/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
baseDir="/path/to/projects/LCLs/directRNAseq"
bam="$baseDir/alignment/dorado/GM${id}_chr_merged/minimap2/GM${id}_chr_merged_minimap2_uniq_sort.bam" # BAM file from regular alignment (not haplotype-aware) to reference genome with minimap2
input_vcf="$baseDir/lorals/vcf/NA${id}_het.vcf" # phased VCF for cell line of interest
outDir="$baseDir/lorals/new_vcf_dorado/GM${id}_chr_merged"

# Software
make_vcf="/path/to/software/lorals/scripts/make_new_vcf.sh"

# Create output directory if non existent
if [ ! -d $outDir ]; then
    mkdir -p $outDir
fi

# Command to run
$make_vcf -b $bam -G $ref_genome -V $input_vcf -o $outDir


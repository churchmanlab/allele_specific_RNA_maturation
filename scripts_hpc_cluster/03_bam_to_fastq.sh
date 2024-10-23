#!/bin/bash


# This script will convert the BAM file obtained from Dorado to a fastq file

sample="GM18861_chr_rep1"

#####################################################################
# CONFIG

# File locations
baseDir="/path/to/projects/LCLs/directRNAseq/dorado"

# Input
bam="$baseDir/$sample/$sample.bam"

#Output
out="$baseDir/$sample/$sample.dorado.fastq"

# END OF CONFIG
#####################################################################

module load bedtools/2.26.0
module load samtools/1.3.1

#######################################
# Create output directory if non existent
if [ ! -d $outDir ]; then
    mkdir -p $outDir
fi


# Command to run
samtools fastq $bam > $out


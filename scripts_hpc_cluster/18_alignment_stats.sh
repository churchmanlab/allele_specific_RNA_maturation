#!/bin/bash

# This script will output alignment statistics (e.g. read length, percent match)

sample="GM18861_chr_rep1"

###################

# CONFIG
baseDir="/path/to/projects/LCLs/directRNAseq"
outDir="$baseDir/analysis_new/$sample"

# Input files
fastq="$baseDir/dorado/$sample/$sample.dorado.fastq"
bam="$baseDir/lorals/alignment_post_dorado_merged_vcf/$sample/${sample}.dorado_reads_aln_sorted.merged.bam"

# Output files
stats="$baseDir/analysis_new/$sample/${sample}_global_seq_stats.txt" # overall measures
fastq_stats="$baseDir/analysis_new/$sample/${sample}_read_length_raw.txt" # per-read measures from fastq
bam_stats="$baseDir/analysis_new/$sample/${sample}_stats_aligned.txt" # per-read measures from BAM


# Python scripts
script_stats="$baseDir/scripts/18_get_stats.py"

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
python $script_stats $fastq $bam $stats $fastq_stats $bam_stats

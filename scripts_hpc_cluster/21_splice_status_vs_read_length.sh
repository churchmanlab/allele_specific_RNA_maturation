#!/bin/bash

# This script will call a Python script to calculate the median read length as a function of the read splicing status

sample="GM18861_chr_rep1"

###################

# CONFIG
baseDir="/path/to/projects/LCLs/directRNAseq"
outDir="$baseDir/analysis_new/$sample"

# Input files
bam_stats="$baseDir/analysis_new/$sample/${sample}_stats_aligned.txt"
splice_df="$baseDir/analysis_new/$sample/${sample}_hg38_read_splicing_status.RefSeq_all.txt"

# Output files
out_df="$baseDir/analysis_new/$sample/${sample}_median_length_per_splice_status.txt"

# Python scripts
script_stats="$baseDir/scripts/21_read_splicing_status_vs_length.py"

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
python $script_stats $splice_df $bam_stats $out_df

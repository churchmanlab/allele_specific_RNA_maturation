#!/bin/bash

# This script will call a Python script to calculate the percent match per gene from the LORALS alignment and the nascent transcriptome alignment

sample="GM18861_chr_rep1"

###################

# CONFIG
baseDir="/path/to/projects/LCLs/directRNAseq"
outDir="$baseDir/analysis_new/$sample"

# Input files
stats_df1="$baseDir/analysis_new/$sample/${sample}_stats_aligned.HLA_nascent_transcriptome.txt"
stats_df2="$baseDir/analysis_new/$sample/${sample}_stats_aligned.txt"
readmap="$baseDir/analysis_new/$sample/${sample}_read_mapping_per_gene.RefSeq.txt"
allelemap="$baseDir/hapcut2/post_dorado_merged_vcf/$sample/${sample}_read_to_allele_mapping.txt"

# Output files
match_df1="$baseDir/analysis_new/$sample/${sample}_stats_aligned.HLA_nascent_transcriptome.match_per_gene.txt"
match_df2="$baseDir/analysis_new/$sample/${sample}_stats_aligned.lorals_alignment.match_per_gene.txt"

# Python scripts
script_stats="$baseDir/scripts/20_get_percent_match_per_gene.py"

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
python $script_stats $stats_df1 $readmap $allelemap $match_df1
python $script_stats $stats_df2 $readmap $allelemap $match_df2

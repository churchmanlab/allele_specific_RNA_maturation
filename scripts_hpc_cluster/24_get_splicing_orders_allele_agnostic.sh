#!/bin/bash


# This script will call a Python script that will compute allele-agnostic splicing orders

sample="GM18861_chr_rep1"

###################

# CONFIG
baseDir="/path/to/projects/LCLs/directRNAseq"
outDir="$baseDir/analysis_new/$sample"

# Input files
chr_splicing_info="$baseDir/analysis_new/$sample/${sample}_hg38_splicing_info_per_intron.withSKP.RefSeq.txt"
chr_multi_intron_df="$baseDir/analysis_new/$sample/${sample}_hg38_multi_introns_isoforms_counts.withSKP.RefSeq_all.txt"

# Output files
splicing_paths="$outDir/${sample}_splicing_order_paths_alleles_combined.all_introns.min10reads.txt"

# Annotation files
annDir="/path/to/reference_libraries/annotations/hg38"
gene_names_df="$annDir/hg38_UCSC_refGene_names.txt"
intron_df="$annDir/NCBI_RefSeq_hg38_introns_parsed.bed"

# Python scripts
script_splice_paths="$baseDir/scripts/24_get_splicing_orders_allele_agnostic.py"

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

python $script_splice_paths $gene_names_df $intron_df $chr_multi_intron_df $chr_splicing_info $splicing_paths


#!/bin/bash

# This script will call a Python script that will analyze allele-specific AS in GTEx sQTL genes

sample1="GM18861_chr_rep1"
sample2="GM18861_chr_rep2"
id="18861"

###################

# CONFIG
baseDir="/path/to/projects/LCLs/directRNAseq"
outDir="$baseDir/analysis_new/GM${id}_chr_merged"

# Input files
mat_bam1="$baseDir/lorals/alignment_post_dorado_merged_vcf/$sample1/${sample1}_RNAtoDNA_reads_aln_sorted.maternal.bam"
pat_bam1="$baseDir/lorals/alignment_post_dorado_merged_vcf/$sample1/${sample1}_RNAtoDNA_reads_aln_sorted.paternal.bam"
mat_bam2="$baseDir/lorals/alignment_post_dorado_merged_vcf/$sample2/${sample2}_RNAtoDNA_reads_aln_sorted.maternal.bam"
pat_bam2="$baseDir/lorals/alignment_post_dorado_merged_vcf/$sample2/${sample2}_RNAtoDNA_reads_aln_sorted.paternal.bam"
snp_df="$baseDir/hapcut2/post_dorado_merged_vcf/$sample1/NA${id}_het.ann.canon.snpsift.vcf"

# Output files
splicing_info="$outDir/GM${id}_chr_merged_splicing_info_per_intron.withSKP.sGenes_Gtex.v8.txt"
splicing_test="$outDir/GM${id}_chr_merged_AS_test.sGenes_Gtex.v8.txt"

# Annotation files
intron_bed_file="/path/to/projects/LCLs/directRNAseq/from_local/Cells_EBV-transformed_lymphocytes.v8.sgenes.GTEX.txt"

# Python scripts
script_splice="$baseDir/scripts/25_AS_for_sQTL_genes.py"

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
python $script_splice $mat_bam1 $mat_bam2 $pat_bam1 $pat_bam2 $snp_df $intron_bed_file $splicing_info $splicing_test


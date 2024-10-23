#!/bin/bash

# This script will assign reads to genes and then call a Python script to count the number of reads per gene and per allele

sample="GM18861_chr_rep1"

###################

# CONFIG
baseDir="/path/to/projects/LCLs/directRNAseq"
outDir="$baseDir/analysis_new/$sample"

# Input files
bed="$baseDir/lorals/alignment_post_dorado_merged_vcf/$sample/${sample}.dorado_reads_aln_sorted.merged.bed"
allelemap="$baseDir/hapcut2/post_dorado_merged_vcf/$sample/${sample}_read_to_allele_mapping.txt"

# Output files
out_cov="$outDir/${sample}_read_mapping_per_gene"
out_counts="$outDir/${sample}_read_counts_per_gene_and_allele"
out_counts2="$outDir/${sample}_read_counts_per_gene"

# Annotation files
ncbi_bed_file="/path/to/reference_libraries/annotations/hg38/NCBI_RefSeq_hg38_genes_parsed.bed"

# Counting script
count_py="$baseDir/scripts/14_gene_counts.py"

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


# Get read to gene correspondance with bedtools intersect, keeping only the read name, gene and overlap columns
bedtools intersect -s -F 0.5 -wo -a $ncbi_bed_file -b $bed | awk 'BEGIN{OFS="\t"} {print $4,$16,$19}' > $out_cov.RefSeq.txt

# Counts reads per gene and allele, or per gene independent of allele
python $count_py $out_cov.RefSeq.txt $allelemap $out_counts.RefSeq.txt $out_counts2.RefSeq.txt


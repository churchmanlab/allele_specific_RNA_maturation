#!/bin/bash


#This script will call a python script that will compute splicing order for each allele for groups of 3 introns

sample="GM18861_chr_rep1"

###################

# CONFIG
baseDir="/path/to/projects/LCLs/directRNAseq"
outDir="$baseDir/analysis_new/$sample"

# Input files
chr_multi_intron_df="$baseDir/analysis_new/$sample/${sample}_hg38_multi_introns_isoforms_df.withSKP.RefSeq_all.txt"
allele_map="$baseDir/hapcut2/post_dorado_merged_vcf/$sample/${sample}_read_to_allele_mapping.txt"

# Output files
splicing_paths="$outDir/${sample}_splicing_order_paths_per_allele.hac.all_introns.min10reads.filterND.stringent.txt"
interm_counts="$outDir/${sample}_interm_counts_per_allele.hac.all_introns.min10reads.filterND.stringent.txt"

# Annotation files
annDir="/path/to/reference_libraries/annotations/hg38"
gene_names_df="$annDir/hg38_UCSC_refGene_names.txt"
intron_df="$annDir/NCBI_RefSeq_hg38_introns_parsed.bed"

# Python scripts
script_splice_paths="$baseDir/scripts/12_get_splicing_orders_per_allele.py"

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

# Compute splicing orders
python $script_splice_paths $gene_names_df $intron_df $chr_multi_intron_df $allele_map $splicing_paths

# Extract intermediate counts only and remove duplicates
awk 'BEGIN{OFS="\t"} {print $18,$19,$16,$17,$7,$3,$5}' $splicing_paths | grep -v "NO_NO_NO" | grep -v "YES_YES_YES" | awk '!seen[$0]++' > $interm_counts


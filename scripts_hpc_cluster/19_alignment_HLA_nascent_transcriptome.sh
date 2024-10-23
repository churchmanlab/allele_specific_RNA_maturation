#!/bin/bash

# This script will create a personalized nascent transcriptome for HLA genes for each cell line, align reads to it and compute splicing orders from this alignment

sample="GM18861_chr_rep1"
id=$(echo "$sample" | grep -oP 'GM\K.{5}') #suffix after GM


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
fastq="$baseDir/dorado/$sample/$sample.dorado.fastq"
fasta_path="$baseDir/lorals/new_vcf_dorado/GM${id}_chr_merged/GM${id}_chr_merged_minimap2_uniq_sort"
bam="$baseDir/lorals/alignment_post_dorado_merged_vcf/$sample/${sample}.dorado_reads_aln_sorted.merged.bam"
allele_map="$baseDir/hapcut2/post_dorado_merged_vcf/$sample/${sample}_read_to_allele_mapping.txt"
outDir="/path/to/projects/LCLs/directRNAseq/HLA_nascent_transcriptome/post_dorado_merged_vcf/$sample"

# Annotation files
HLA_bed="$baseDir/scripts/HLA_genes.bed"
annDir="/path/to/reference_libraries/annotations/hg38"
gene_names_df="$annDir/hg38_UCSC_refGene_names.txt"
hg38_gene_df="$annDir/from_Heather/NCBI_RefSeq_hg38_merge_parsed_sortByNameCoord.bed"

# Output files
sub_bam="$outDir/${sample}_pass_RNAtoDNA_reads_aln_sorted.merged.HLA_only.bam"
sub_fq="$outDir/${sample}_pass_RNAtoDNA_reads_aln_sorted.merged.HLA_only.fastq"
splicing_paths="$outDir/${sample}_splicing_order_paths_per_allele.hac.HLA_transcriptome.min10reads.filterND.txt"

# Scripts
python_script="$baseDir/scripts/19a_make_HLA_nascent_transcriptome_fasta.py"
hap_aligner="/path/to/software/lorals/scripts/hap_aligner.transcriptome.sh"
order_script="$baseDir/scripts/19b_get_splicing_paths_HLA_nascent_transcriptome.py"

# Create output directory if non existent
if [ ! -d $outDir ]; then
    mkdir -p $outDir
fi


# 1) Extract reads mapping to HLA genes and convert to fastq
samtools view -hb $bam --region-file $HLA_bed > $sub_bam
bedtools bamtofastq -i $sub_bam -fq $sub_fq

# 2) Make allele-specific transcriptome reference for each sample
python $python_script $gene_names_df $hg38_gene_df ${fasta_path}_hap1.fa $outDir/NA${id}_HLA_genes_hap1.fa
python $python_script $gene_names_df $hg38_gene_df ${fasta_path}_hap2.fa $outDir/NA${id}_HLA_genes_hap2.fa

# 3) Align reads from 1) to reference from 2) without the splice option
$hap_aligner -f $sub_fq --reference $outDir/NA${id}_HLA_genes -o $outDir

# 4) Convert BAM to BED
HLA_bam="$outDir/${sample}_pass_RNAtoDNA_reads_aln_sorted.merged.HLA_only_reads_aln_sorted.merged.bam"
HLA_bed="$outDir/${sample}_pass_RNAtoDNA_reads_aln_sorted.merged.HLA_only_reads_aln_sorted.merged.bed"
bedtools bamtobed -cigar -tag NM -i $HLA_bam > $HLA_bed

# 5) Compute splicing orders
python $order_script $HLA_bed $allele_map $hg38_gene_df $gene_names_df $splicing_paths

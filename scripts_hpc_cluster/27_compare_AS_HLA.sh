#!/bin/bash


# This script will call a Python script that will compare AS in HLA genes between alleles


sample1="GM19144_chr_rep1"
sample2="GM19144_chr_rep2"
merged_sample="GM19144_chr_merged"


###################

# CONFIG
baseDir="/path/to/projects/LCLs/directRNAseq"
outDir="$baseDir/analysis_new/${merged_sample}"

# Input files
allelemap1="$baseDir/hapcut2/post_dorado_merged_vcf/$sample1/${sample1}_read_to_allele_mapping.txt"
allelemap2="$baseDir/hapcut2/post_dorado_merged_vcf/$sample2/${sample2}_read_to_allele_mapping.txt"

splice_df1="$baseDir/analysis_new/$sample1/${sample1}_hg38_splicing_info_per_intron.withSKP.RefSeq.txt"
splice_df2="$baseDir/analysis_new/$sample2/${sample2}_hg38_splicing_info_per_intron.withSKP.RefSeq.txt"
path_df1="$baseDir/analysis_new/$sample1/${sample1}_splicing_order_paths_per_allele.hac.all_introns.min10reads.filterND.stringent.txt"
path_df2="$baseDir/analysis_new/$sample2/${sample2}_splicing_order_paths_per_allele.hac.all_introns.min10reads.filterND.stringent.txt"

# Annotation files
hg38_intron_df="/path/to/reference_libraries/annotations/hg38/hg38_all_intron_features.txt"

# Output files
out_AS="$outDir/${merged_sample}_AS_events.HLA-C.txt"

# Python script
AS_py="$baseDir/scripts/27_identify_AS_introns.py"

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
python $AS_py $splice_df1 $splice_df2 $path_df1 $path_df2 $allelemap1 $allelemap2 $hg38_intron_df $out_AS


#!/bin/bash


# This script will call a Python script to compute splicing order for pairs of consecutive introns

sample="GM18861_chr_rep1"

###################

# CONFIG
baseDir="/path/to/projects/LCLs/directRNAseq"
outDir="$baseDir/analysis_new/$sample"

# Input and output files
splicing_dict="$outDir/${sample}_hg38_splicing_dictionary.withSKP.RefSeq"
intron_pairs_df="$outDir/${sample}_hg38_intron_pairs_df"
allele_df="$baseDir/hapcut2/post_dorado_merged_vcf/$sample/${sample}_read_to_allele_mapping.txt"
intron_pairs_counts="$outDir/${sample}_hg38_intron_pairs_counts.txt"

# Annotation files
intron_bed_file="/path/to/reference_libraries/annotations/hg38/NCBI_RefSeq_hg38_introns_parsed_bedtool_for_intron_pairs.bed"
hg38_introns="/path/to/reference_libraries/annotations/hg38/hg38_all_intron_features.txt"
gene_names_df="/path/to/reference_libraries/annotations/hg38/hg38_UCSC_refGene_names.txt"

# Python scripts
script_py="$baseDir/scripts/13a_get_intron_pairs_df.split_by_chrom.py"
script_py2="$baseDir/scripts/13b_process_intron_pairs_df.py"

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


# Compute the splicing status of intron pairs
# Get splicing statuses per chromosome, then concatenate
cd $outDir
echo -e "read\tgene\tint1_count\tint2_count\tstrand\tint1_splice\tint2_splice" > ${intron_pairs_df}_all.txt
for i in *RefSeq_chr*.npy ; do
  chrom=$(echo $i | cut -f9- -d_ | cut -f1 -d.) 
  echo $chrom
  python $script_py $i ${intron_pairs_df}_${chrom}.txt
  cat ${intron_pairs_df}_${chrom}.txt >> ${intron_pairs_df}_all.txt
  rm ${intron_pairs_df}_${chrom}.txt
done

# Match to allele and get counts
python $script_py2 ${intron_pairs_df}_all.txt $allele_df $hg38_introns $gene_names_df $intron_pairs_counts

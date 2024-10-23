#!/bin/bash


#This script will call a python script that determines the splicing status of individual introns and reads from nanopore sequencing data

sample="GM18861_chr_rep1"

###################

# CONFIG
baseDir="/path/to/projects/LCLs/directRNAseq"
outDir="$baseDir/analysis_new/$sample"

# Input files
bam="$baseDir/lorals/alignment_post_dorado_merged_vcf/$sample/${sample}.dorado_reads_aln_sorted.merged.bam"
bed="$baseDir/lorals/alignment_post_dorado_merged_vcf/$sample/${sample}.dorado_reads_aln_sorted.merged.bed"

# Output files
splicing_info="$outDir/${sample}_hg38_splicing_info_per_intron.withSKP.RefSeq.txt"
splicing_dict="$outDir/${sample}_hg38_splicing_dictionary.withSKP.RefSeq"
multi_introns_df="$outDir/${sample}_hg38_multi_introns_isoforms_df.withSKP.RefSeq"
multi_introns_counts="$outDir/${sample}_hg38_multi_introns_isoforms_counts.withSKP.RefSeq"
splicing_status="$outDir/${sample}_hg38_read_splicing_status.RefSeq"
splicing_status_counts="$outDir/${sample}_hg38_read_splicing_status.counts.RefSeq"

# Annotation files
intron_bed_file="/path/to/reference_libraries/annotations/hg38/NCBI_RefSeq_hg38_introns_parsed_bedtool_for_intron_pairs.bed"


# Python scripts
script_splice_dict="$baseDir/scripts/09a_splice_dict_dataset_per_chrom.py"
script_multi_introns="$baseDir/scripts/09b_get_multi_introns_isoforms_df.per_chrom.py"
script_splice_stats="$baseDir/scripts/09c_read_splicing_stats.per_chrom.py"

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



# Get bed file from BAM file and split into multiple lines while both mates are still there (after sorting it by read name)
echo "Converting BAM file to BED..."
bedtools bamtobed -cigar -tag NM -i $bam > $bed

# Intersect with introns
echo "Intersecting with annotated introns..."
bedtools intersect -s -a $bed -b $intron_bed_file -wo > $outDir/$sample.intersect_introns.bed

# Make splicing info and splicing dictionary datasets
python $script_splice_dict $outDir/$sample.intersect_introns.bed $splicing_info $splicing_dict

# Get splicing statuses per chromosome, then concatenate
cd $outDir
echo -e "read\tgene\tstrand\tintron_numbers\tsplice_status" > ${multi_introns_df}_all.txt
echo -e "gene\tstrand\tintron_numbers\tsplice_status\tcount" > ${multi_introns_counts}_all.txt
echo -e "read\tgene\tsplicing_status" > ${splicing_status}_all.txt
echo -e "all_spliced\tall_unspliced\tintermediate" > ${splicing_status_counts}_all.txt
for i in *RefSeq_chr*.npy ; do 
  chrom=$(echo $i | cut -f9- -d_ | cut -f1 -d.)
  echo $chrom
  python $script_multi_introns $i ${multi_introns_df}_${chrom}.txt ${multi_introns_counts}_${chrom}.txt
  cat ${multi_introns_df}_${chrom}.txt >> ${multi_introns_df}_all.txt
  cat ${multi_introns_counts}_${chrom}.txt >> ${multi_introns_counts}_all.txt
  rm ${multi_introns_df}_${chrom}.txt
  rm ${multi_introns_counts}_${chrom}.txt
  python $script_splice_stats $i ${splicing_status}_${chrom}.txt ${splicing_status_counts}_${chrom}.txt
  cat ${splicing_status}_${chrom}.txt >> ${splicing_status}_all.txt
  cat ${splicing_status_counts}_${chrom}.txt >> ${splicing_status_counts}_all.txt
  rm ${splicing_status}_${chrom}.txt
  rm ${splicing_status_counts}_${chrom}.txt
done

rm $outDir/$sample.intersect_introns.bed


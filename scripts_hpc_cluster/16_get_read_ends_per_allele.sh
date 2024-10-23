#!/bin/bash


# This script will call a Python script that extracts the read end position of each read

sample="GM18861_chr_rep1"

###################

# File locations
baseDir="/path/to/projects/LCLs/directRNAseq"

#Input: the two parental BAM files and the read to gene mapping
bam_mat="$baseDir/lorals/alignment_post_dorado_merged_vcf/$sample/${sample}_RNAtoDNA_reads_aln_sorted.maternal.bam"
bam_pat="$baseDir/lorals/alignment_post_dorado_merged_vcf/$sample/${sample}_RNAtoDNA_reads_aln_sorted.paternal.bam"
readmap="$baseDir/analysis_new/$sample/${sample}_read_mapping_per_gene.RefSeq.txt"

#Output
out_mat="$baseDir/analysis_new/$sample/${sample}_read_ends.maternal.txt"
out_pat="$baseDir/analysis_new/$sample/${sample}_read_ends.paternal.txt"
out_both="$baseDir/analysis_new/$sample/${sample}_read_ends.both_alleles.txt"

# Python scripts
python_script="$baseDir/scripts/16_get_read_ends.py"

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
echo -e "chrom\tstart\tend\treadname\tscore\tstrand\tgene\toverlap\tallele" > $out_both
python $python_script $bam_mat $out_mat $readmap 'M'
python $python_script $bam_pat $out_pat $readmap 'P'
cat $out_mat $out_pat >> $out_both

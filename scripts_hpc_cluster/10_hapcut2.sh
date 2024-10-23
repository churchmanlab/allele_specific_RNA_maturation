#!/bin/bash

#This script will assign reads to alleles using HAPCUT2 and then call a python script adapted from https://github.com/nanopore-wgs-consortium/NA12878/blob/master/nanopore-human-transcriptome/scripts/ase.py to parse the output


##########################################################################################

id="18861" # suffix that goes after "NA" or "GM" in LCL name
sample_path="GM18861_chr_merged/GM18861_chr_merged" # prefix of VCF created with the merged samples from the same cell line
sample="GM18861_chr_rep1" # replicate being analyzed


# Load necessary packages and virtual environment
module load gcc/6.2.0
module load samtools/1.15.1
module load bedtools/2.27.1
module load python/3.7.4
module load bcftools/1.13
module load gatk/4.1.9.0
module load htslib/1.9.0
source ~/myenv_python3/bin/activate


# CONFIG
ref_genome="/path/to/reference_libraries/genomes/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
baseDir="/path/to/projects/LCLs/directRNAseq"
bam="$baseDir/lorals/alignment_post_dorado_merged_vcf/$sample/${sample}.dorado_reads_aln_sorted.merged.bam"
input_vcf="$baseDir/lorals/new_vcf_dorado/${sample_path}_minimap2_uniq_sort.vcf"
outDir="/path/to/projects/LCLs/directRNAseq/hapcut2/post_dorado_merged_vcf/$sample"

# Software
extractHAIRS="/home/kc248/ld.test/HapCUT2-v.1.3.3/build/extractHAIRS"
HapCUT22VCF_python="/path/to/software/HAPCUT2VCF_script.py"
HAPCUT2="/home/kc248/ld.test/HapCUT2-v.1.3.3/build/HAPCUT2"
snpEff="/path/to/software/snpEff"
python_script="$baseDir/scripts/read_to_allele_mapping.py"

# Create output directory if non existent
if [ ! -d $outDir ]; then
    mkdir -p $outDir
fi


# Command to run
$extractHAIRS --nf 1 --ont 1 --bam $bam --VCF $input_vcf --out $outDir/${sample}_het.fullgenome.extractHairs --ref $ref_genome
java -Xmx12g -jar $snpEff/snpEff.jar -c $snpEff/snpEff.config -d -v -canon -no-downstream -no-intergenic -no-upstream GRCh38.86 $input_vcf > $outDir/NA${id}_het.ann.canon.vcf
java -jar $snpEff/SnpSift.jar extractFields -s "," -e "." $outDir/NA${id}_het.ann.canon.vcf CHROM POS REF ALT "ANN[*].GENE" "ANN[*].EFFECT" "GEN[*].GT" > $outDir/NA${id}_het.ann.canon.snpsift.vcf
python $python_script $outDir/NA${id}_het.ann.canon.snpsift.vcf $outDir/${sample}_het.fullgenome.extractHairs $outDir/${sample}_ase.txt $outDir/${sample}_read_to_allele_mapping.txt


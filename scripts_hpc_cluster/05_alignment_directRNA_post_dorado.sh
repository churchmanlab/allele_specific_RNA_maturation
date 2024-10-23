#!/bin/bash

# This script will align reads to the human genome with minimap2 (first alignment, not haplotype-aware)

module load gcc/6.2.0
module load samtools/1.3.1

sample="GM19152_chr_rep1"

#############################################

# Run minimap2 alignment (with ONT default parameters) for all samples

baseDir="/path/to/projects/LCLs/directRNAseq"
programDir="/path/to/software/minimap2"
reference="/path/to/reference_libraries/genomes/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.and.sacCer3_pool1.fa"
readDir="$baseDir/dorado/$sample"
reads="$readDir/$sample.dorado.fastq"
outDir="$baseDir/alignment/dorado/$sample/minimap2"


# Create output directory if non existent
if [ ! -d $outDir ]; then
    mkdir -p $outDir
fi

# Align with minimap2
${programDir}/minimap2 -ax splice -uf -k14 ${reference} ${reads} > ${outDir}/${sample}_minimap2.sam
samtools sort ${outDir}/${sample}_minimap2.sam -o ${outDir}/${sample}_minimap2_sort.bam
samtools index ${outDir}/${sample}_minimap2_sort.bam


###############################################3


# Extract unique reads

for f in $sample
do
    Dir=$baseDir/alignment/dorado/${f}/minimap2

    echo "Doing file "$f
    cd ${Dir}

    grep ^@ ${f}_minimap2.sam > headers.sam
    grep -v ^@  ${f}_minimap2.sam | sort > alignment.sam
    awk '$3!="*" {print}' alignment.sam | cut -f 1 | sort | uniq -c | awk '$1==1 {print $2}' > uniq_names.txt
    grep -F -f uniq_names.txt alignment.sam > uniq_alignment.sam
    cat headers.sam uniq_alignment.sam > headers_uniq_temp.sam
    samtools view -bT ${reference} headers_uniq_temp.sam > headers_uniq_temp.bam
    samtools sort headers_uniq_temp.bam -o ${f}_minimap2_uniq_sort.bam
    samtools index ${f}_minimap2_uniq_sort.bam
    rm headers.sam
    rm alignment.sam
    rm uniq_names.txt
    rm uniq_alignment.sam
    rm headers_uniq_temp.sam
    rm headers_uniq_temp.bam
done


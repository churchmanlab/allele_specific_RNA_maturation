#!/bin/bash


# This script will extract the poly(A) tail length for each read from the BAM file output from Dorado

#####################################################################

sample="GM18861_chr_rep1"

#####################################################################
# CONFIG

# File locations
baseDir="/path/to/projects/LCLs/directRNAseq/dorado"

# Input
bam="$baseDir/$sample/$sample.bam"

#Output
out="$baseDir/$sample/$sample.polyA_dorado.txt"

# END OF CONFIG
#####################################################################

module load bedtools/2.26.0
module load samtools/1.3.1

#######################################
# Create output directory if non existent
if [ ! -d $outDir ]; then
    mkdir -p $outDir
fi


# Command to run
echo -e "readname\tpolya_length" > $out
samtools view $bam | awk 'BEGIN{OFS="\t"} {if ($NF ~ /pt:i/) {print $1, $NF}}' | sed 's/pt:i://g' >> $out


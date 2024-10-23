#!/bin/bash

# This script will perform high accuracy basecalling with Dorado


sample="GM19223_cyto_rep1"

module load gcc/9.2.0 
module load cuda/11.7


# CONFIG
baseDir="/path/to/projects/LCLs/directRNAseq"
inDir="$baseDir/startingFiles/$sample/pod5"
outDir="$baseDir/dorado/$sample"

#Software
dorado="/path/to/software/dorado-0.7.0-linux-x64/bin/dorado"
rna002_model="/path/to/software/dorado-0.7.0-linux-x64/models/rna002_70bps_hac@v3"

# Create output directory if non existent
if [ ! -d $outDir ]; then
    mkdir -p $outDir
fi


# Basecall with dorado
# For SQK-RNA004:
#$dorado basecaller hac $inDir/ --estimate-poly-a > $outDir/$sample.bam

# For SQK-RNA002:
$dorado basecaller $rna002_model $inDir/ --estimate-poly-a > $outDir/$sample.bam



#!/bin/bash


# This script will convert fast5 files to pod5 format for basecalling with Dorado


sample="GM19223_cyto_rep1"

module load python/3.10.11
source ~/py310_pod5_env/bin/activate

# CONFIG
baseDir="/path/to/projects/LCLs/directRNAseq"
inDir="$baseDir/startingFiles/$sample/fast5_pass"
outDir="$baseDir/startingFiles/$sample/pod5"

mkdir -p $outDir

# Convert fast5 to pod5
pod5 convert fast5 $inDir/*.fast5 --output $outDir/ --one-to-one $inDir/

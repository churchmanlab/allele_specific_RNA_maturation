"""

Author: Karine Choquet

Date: January 24, 2022

This script will analyze basic statistics from direct RNA-seq datasets

"""



import sys
import numpy as np
import pandas as pd
import pysam
import pybedtools

# CONFIG
fastq = open(sys.argv[1], 'r')
my_iBAM = pysam.Samfile(sys.argv[2], 'rb')


##########

# Functions

# Plot read length
# get read length from fastq file
def get_read_length_df(fastq):
    count=0
    readLength_list = []

    for line in fastq:

        count+=1

        if (count%4==1):
            name = line.split('@')[1].split(' ')[0]

        if (count%4==2):
            seq = line
            read_length = len(line)

        if (count%4==0):
            readLength_list.append([name, read_length])

    readLength_df = pd.DataFrame(readLength_list)
    readLength_df.columns = ['name','length']

    return readLength_df



# get alignment stats information from bamFile
def get_alignment_stats(iBAM):

    alignment_stats = []

    # read in Bam file line by line
    for read in iBAM:

        # set variables to cigar string before summing
        M = 0    # number of alignment matches
        I = 0    # number of alignment insertions
        D = 0    # number of alignment deletions
        N = 0    # number of skipped nucleotdies from splicing
        S = 0    # number of softclipped nucleotides
        H = 0    # number of hardclipped nucleotides
        P = 0
        E = 0
        X = 0
        B = 0

        # only record data if read is aligned 
        if(read.rname!=-1):

            name = read.qname     # get read name  
            cigar = read.cigar    # get cigar string for read 

            # get NM which details edit distance to reference genome
            for i in range(0,len(read.tags)):
                if(read.tags[i][0]=='NM'):
                    NM=read.tags[i][1]

            # sum up the number of matches for each representative 
            # letter in the cigar string
            for i in range(0,len(cigar)):
                if(cigar[i][0]==0):
                    M = M + cigar[i][1]
                if(cigar[i][0]==1):
                    I = I + cigar[i][1]
                if(cigar[i][0]==2):
                    D = D + cigar[i][1]
                if(cigar[i][0]==3):
                    N = N + cigar[i][1]
                if(cigar[i][0]==4):
                    S = S + cigar[i][1]
                if(cigar[i][0]==5):
                    H = H + cigar[i][1]
                if(cigar[i][0]==6):
                    P = P + cigar[i][1]
                if(cigar[i][0]==7):
                    E = E + cigar[i][1]
                if(cigar[i][0]==8):
                    X = X + cigar[i][1]
                if(cigar[i][0]==9):
                    B = B + cigar[i][1]

            readLength = M + I + S
            alignLength = M + D
            align_percent = float(M) / float(readLength) * 100.0
            match_count = M + I - NM       
            match_percent = float(match_count) / float(M + I) * 100.0
            alignment_stats.append([read.qname,readLength,alignLength,align_percent,match_count,match_percent])
    alignment_stats_df = pd.DataFrame(alignment_stats)
    alignment_stats_df.columns = ['read','read_length','aligned_length','align_percent','match_count','match_percent']

    iBAM.close()
    
    return alignment_stats_df



############

# Get statistics for raw reads
fastq_stats = get_read_length_df(fastq)
total_reads = len(fastq_stats)

# Get statistics for aligned reads
aln_stats_df  = get_alignment_stats(my_iBAM)
uniq_reads = len(aln_stats_df)
fraction_uniq_align = uniq_reads/total_reads
median_read_length = aln_stats_df['read_length'].median()

# Output global stats to file
stats_file = open(sys.argv[-3], 'w')
stats_file.write("total reads: " + str(total_reads) + '\n' + "uniquely mapped reads: " + str(uniq_reads) + '\n' + "fraction uniquely mapped: " + str(fraction_uniq_align) + '\n' + "median_read_length: " + str(median_read_length) + '\n')
stats_file.close()

# Output per-read statistics
fastq_stats.to_csv(sys.argv[-2], sep="\t", header=True, index=False)
aln_stats_df.to_csv(sys.argv[-1], sep="\t", header=True, index=False)

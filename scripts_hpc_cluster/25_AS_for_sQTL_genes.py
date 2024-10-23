"""
This script will analyze allele-specific AS in GTEx sQTL genes

"""

import sys
import numpy as np
import pandas as pd
import pysam
from collections import Counter

import scipy
from scipy import stats

from statsmodels.stats.multitest import multipletests

import re
import math

import pybedtools
from pybedtools import BedTool

#####################################

# CONFIG

# Allele-specific BAM files
mat_bam_1 = pybedtools.BedTool(sys.argv[1])
mat_bam_2 = pybedtools.BedTool(sys.argv[2])
pat_bam_1 = pybedtools.BedTool(sys.argv[3])
pat_bam_2 = pybedtools.BedTool(sys.argv[4])

# Per-sample  heterozygous SNPs
het_df = pd.read_table(sys.argv[5])

# Annotation file with sGenes
gtex_df = pd.read_table(sys.argv[6])

# Get relevant columns from gtex file
gtex_df_sub = gtex_df.copy()[['phenotype_id','gene_name','chr','variant_pos','rs_id_dbSNP151_GRCh38p7','ref','alt','strand']]
gtex_df_sub['chr'] = gtex_df_sub['chr'].str[3:].astype(int)

# Merge with het SNPs
merge_df = het_df.merge(gtex_df_sub, left_on=['CHROM','POS','REF','ALT'], right_on=['chr','variant_pos','ref','alt'])

# Make BED file with corresponding introns
merge_df['chrom'] = merge_df['phenotype_id'].str.split(":").str[0].str[3:]
merge_df['intron_start'] = merge_df['phenotype_id'].str.split(":").str[1]
merge_df['intron_end'] = merge_df['phenotype_id'].str.split(":").str[2]
merge_df['cluster'] = merge_df['phenotype_id'].str.split(":").str[3]
merge_df['gene_id'] = merge_df['phenotype_id'].str.split(":").str[4]
gtex_bed_df = merge_df.copy()[['chrom','intron_start','intron_end','gene_id','cluster','strand']]

# Convert to bed
gtex_bed = pybedtools.BedTool.from_dataframe(gtex_bed_df)


######################################


# function to create a dataframe with reads that span 3'SS positions
def get_intron_intersect(introns_df, bam_file):
    # get reads that span 3' splice sites and convert to a dataframe
    bedFile = bam_file.bam_to_bed(cigar=True, tag='NM') # convert bam file to bed file, keep cigar string and NM (edit distance) tag
    intersect = bedFile.intersect(introns_df, s=True, wo=True) # intersect reads from bam file with 3' splice site coordinates
    df = intersect.to_dataframe(names=['chr_aln', 'start_aln', 'end_aln', 'name_aln', 'qual_aln', \
                                           'strand_aln', 'cigar_aln', 'chr_intron', 'start_intron', \
                                           'end_intron', 'name_gene', 'cluster', 'strand_gene', 'count'], \
                               dtype={"chr_aln": str, "start_aln": int, "end_aln": int, \
                                     "name_aln": str, "qual_aln": int, "strand_aln": str, \
                                     "cigar_aln": str, "chr_intron": str, "start_intron": int, \
                                     "end_intron": int, "name_gene": str, \
                                     "cluster": str,"strand_gene": str, "count": int}) # convert to a dataframe
    return df


# function to create a dataframe with splicing information for
# every read that spans an intron in the dataset
def get_splicing_info(intersect_df, min_overlap):
   
    df = intersect_df 

    # prepare a list for splice calls
    spliceCalls = []

    # set variables for parsing the cigar string
    pattern = re.compile('([MIDNSHPX=])')
    Consumes_Query = ["M", "I", "S", "=", "X"]
    Consumes_Reference = ["M", "D", "N", "=", "X"]    

    # loop through all read-intron intersects
    for i in range(0,df.shape[0]):


        # ignore reads that do not overlap intron by minimum threshold
        if (df['count'].iloc[i] < min_overlap):
            continue

        # record the start and ends of reads
        # record the start and ends of reads
        aln_start = df['start_aln'].iloc[i] # record the start of the read
        aln_end = df['end_aln'].iloc[i] # record the end of the read
        intron_start = df['start_intron'].iloc[i] # record the end of the intron
        intron_end = df['end_intron'].iloc[i] # record the end of the intron
        intron_chr = df['chr_intron'].iloc[i]
        aln_strand = df['strand_aln'].iloc[i]
        aln_name = df['name_aln'].iloc[i]
        gene_strand = df['strand_gene'].iloc[i]
        cigar_aln = df['cigar_aln'].iloc[i]
        name_gene = df['name_gene'].iloc[i]
        cluster = df['cluster'].iloc[i]

        # parse cigar string into a list of tuples for easy parsing
        Sep_Values = Sep_Values = pattern.split(cigar_aln)[:-1]
        CigarPairs = list((Sep_Values[n:n+2] for n in range(0, len(Sep_Values), 2)))  

        # set up variables for measuring the length of cigar string operators
        CigarOp_counts = {'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0, 'X': 0}
        start_counts = {'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0, 'X': 0}
        end_counts = {'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0, 'X': 0}
        intron_counts = {'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0, 'X': 0}
        currentloc = aln_start

        # go through list of cigar strings and grab splicing information
        for cigar_Entry in CigarPairs:

            op_Length = int(cigar_Entry[0]) # get length of cigar operator
            cigarOp = cigar_Entry[1] # get type of cigar operator  
            CigarOp_counts[cigarOp] += op_Length # add the cigar operator length to the counts dictionary
            cigarOp_start=currentloc # get the starting coordinate of the cigar operator

            if (cigarOp in Consumes_Reference):
                currentloc=currentloc+op_Length # add the cigar operator length to the current location coordinate 

            cigarOp_end=currentloc # get the ending coordinate of the cigar operator

            # gather information if the portion of the cigar string spans the designated intron start
            if (cigarOp_start<intron_start-min_overlap and cigarOp_end>=intron_start-min_overlap):
                if (cigarOp_end>=intron_start+min_overlap):
                    count=min_overlap*2
                else:
                    count=cigarOp_end-(intron_start-min_overlap)+1
                start_counts[cigarOp] += count # add the cigar operator length to the counts dictionary

            elif (cigarOp_start>=intron_start-min_overlap and cigarOp_end<intron_start+min_overlap):
                count=op_Length
                start_counts[cigarOp] += count # add the cigar operator length to the counts dictionary       

            elif (cigarOp_start<intron_start+min_overlap and cigarOp_end>=intron_start+min_overlap):
                if (cigarOp_start<=intron_start-min_overlap):
                    count=min_overlap*2
                else:
                    count=(intron_start+min_overlap)-cigarOp_start-1
                start_counts[cigarOp] += count # add the cigar operator length to the counts dictionary 

            # gather information if the portion of the cigar string is within the intron
            if (cigarOp_start<intron_start and cigarOp_end>=intron_start):
                if (cigarOp_end>=intron_end):
                    count=intron_end-intron_start
                else:
                    count=cigarOp_end-intron_start
                intron_counts[cigarOp] += count # add the cigar operator length to the counts dictionary

            elif (cigarOp_start>=intron_start and cigarOp_end<intron_end):
                count=op_Length
                intron_counts[cigarOp] += count # add the cigar operator length to the counts dictionary

            elif (cigarOp_start<intron_end and cigarOp_end>=intron_end):
                if (cigarOp_start<=intron_start):
                    count=intron_end-intron_start
                else:
                    count=intron_end-cigarOp_start
                intron_counts[cigarOp] += count # add the cigar operator length to the counts dictionary 

            # gather information if the portion of the cigar string spans the designated intron end
            if (cigarOp_start<intron_end-min_overlap and cigarOp_end>=intron_end-min_overlap):
                if (cigarOp_end>=intron_end+min_overlap):
                    count=min_overlap*2
                else:
                    count=cigarOp_end-(intron_end-min_overlap)
                end_counts[cigarOp] += count # add the cigar operator length to the counts dictionary

            elif (cigarOp_start>=intron_end-min_overlap and cigarOp_end<intron_end+min_overlap):
                count=op_Length
                end_counts[cigarOp] += count # add the cigar operator length to the counts dictionary

            elif (cigarOp_start<intron_end+min_overlap and cigarOp_end>=intron_end+min_overlap):
                if (cigarOp_start<=intron_end-min_overlap):
                    count=min_overlap*2
                else:
                    count=(intron_end+min_overlap)-cigarOp_start
                end_counts[cigarOp] += count # add the cigar operator length to the counts dictionary 

        # get length of the aligned portion of this read from cigar string
        aligned_read_length = CigarOp_counts['M']+CigarOp_counts['D']

       # get 5'SS and 3'SS counts as determined by gene strand
        strand = gene_strand
        if (strand == '+'):
            aln_start = aln_start  # record the start of the read
            aln_end = aln_end # record the end of the read
            intron_5SS_counts = start_counts # record the cigar string counts over the 5'SS
            intron_3SS_counts = end_counts # record the cigar string counts over the 3'SS
            read_overlap = intron_end - (aln_end - aligned_read_length + min_overlap)
            
        if (strand == '-'):
            aln_start = aln_end # record the start of the read
            aln_end = aln_start  # record the end of the read
            intron_5SS_counts = end_counts # record the cigar string counts over the 5'SS
            intron_3SS_counts = start_counts # record the cigar string counts over the 3'SS  
            read_overlap = (aln_end + aligned_read_length - min_overlap) - intron_start
            
        # annotate splicing status based on CIGAR string information around splice sites
        splice='UND'

        if (intron_5SS_counts['N']==0 and intron_3SS_counts['N']==0):
            if (intron_3SS_counts['M']+intron_3SS_counts['D']==min_overlap*2):
                if (intron_3SS_counts['M']>min_overlap):
                    splice = 'NO'

        if (intron_5SS_counts['N']>0 and intron_5SS_counts['N']<min_overlap*2):
            if (intron_3SS_counts['N']>0 and intron_3SS_counts['N']<min_overlap*2):
                splice = 'YES'
            if (intron_3SS_counts['N']==min_overlap*2):
                splice = 'SKP3SS'

        if (intron_5SS_counts['M']==0 and intron_3SS_counts['M']==0):
            if (intron_3SS_counts['N']==min_overlap*2):
                    splice = 'SKP'

        if (intron_5SS_counts['N']==min_overlap*2):
            if (intron_3SS_counts['N']>0 and intron_3SS_counts['N']<min_overlap*2):
                splice = 'SKP5SS'


        # annotate splicing status based on CIGAR string information within the intron 
        if (splice == 'YES'):
            if (float(intron_end-intron_start) > 0.0):
                ratio = float(intron_counts['N'])/float(intron_end-intron_start)
                difference = abs(intron_counts['N']-(intron_end-intron_start))

                # if read is spliced, between 90-100% of the intron has to be spliced 
                # and no more than 100 nucleotides within the intron can be matching the intron sequence
                if( ratio < 0.8 or ratio > 1.1 or difference > 50):
                    splice='UND'
            if (float(intron_end-intron_start) == 0.0):
                splice='UND'

        if (splice == 'NO'):
            if (float(intron_end-intron_start) > 0.0):
                ratio = float(intron_counts['M'])/(float(intron_counts['M'])+float(intron_counts['N'])+float(intron_counts['D'])+1)

                # if read is unspliced, at least 75% of the read has to match (CIGAR=M) the intron sequence
                if(intron_counts['M'] < min_overlap/2 or ratio < 0.75):
                    splice='UND'
            
            if (float(intron_end-intron_start) == 0.0):
                splice='UND'

        if (splice == 'SKP_3SS') or (splice == 'SKP_5SS'):
            if (float(intron_end-intron_start) > 0.0):
                ratio = float(intron_counts['N'])/float(intron_end-intron_start)
                difference = abs(intron_counts['N']-(intron_end-intron_start))

                # if read is spliced, between 90-110% of the intron has to be spliced
                # and no more than 100 nucleotides within the intron can be matching the intron sequence
                if( ratio < 0.8 or ratio > 1.1 or difference > 50): # changed this for DDX39A
                    splice='UND'
            if (float(intron_end-intron_start) == 0.0):
                splice='UND'


        # save read, intron, and splicing information
        spliceCalls.append([aln_name,intron_chr,intron_start,intron_end,gene_strand,name_gene,cluster,read_overlap,splice])

    spliceCalls_df = pd.DataFrame(spliceCalls)
    spliceCalls_df.columns = ["read_name","chrom","intron_start","intron_end","strand","gene_name","cluster","read_overlap","splice_status"]

    return spliceCalls_df



##############################################

# set all variables for analysis
min_overlap = 25

# intersect reads and sIntrons
intersect_mat_1 = get_intron_intersect(gtex_bed, mat_bam_1)
intersect_mat_2 = get_intron_intersect(gtex_bed, mat_bam_2)
intersect_pat_1 = get_intron_intersect(gtex_bed, pat_bam_1)
intersect_pat_2 = get_intron_intersect(gtex_bed, pat_bam_2)

# combine replicates
intersect_mat = pd.concat([intersect_mat_1,intersect_mat_2]).reset_index(drop=True)
intersect_pat = pd.concat([intersect_pat_1,intersect_pat_2]).reset_index(drop=True)

# get splicing information for every read that spans an intron
print("getting splicing info...")
splice_info_mat = get_splicing_info(intersect_mat,min_overlap)
splice_info_pat = get_splicing_info(intersect_pat,min_overlap)

# apply some filters and sort the dataframe
splice_info_mat['allele'] = 'M'
splice_info_pat['allele'] = 'P'
splice_info = pd.concat([splice_info_mat,splice_info_pat]).reset_index(drop=True)
splice_info = splice_info.sort_values(by=['chrom','intron_start','intron_end']).reset_index(drop=True)
splice_info.to_csv(sys.argv[-2], sep='\t', index=False, header=True)


# count the total number of reads per cluster
t = 20
splice_counts = pd.DataFrame(splice_info.groupby(['gene_name','cluster','allele','splice_status'])['read_name'].count()).reset_index().rename(columns={'gene_name':'gene'})

# Filter introns with a sufficient number of YES reads on each allele
splice_counts_M = splice_counts[splice_counts['allele']=='M'].reset_index(drop=True)
splice_counts_M_piv = splice_counts_M.pivot_table(index=['gene','cluster'], columns='splice_status', values='read_name').fillna(0).reset_index()
splice_counts_M_piv['ALT'] = splice_counts_M_piv['SKP'] + splice_counts_M_piv['SKP3SS'] + splice_counts_M_piv['SKP5SS'] + splice_counts_M_piv['UND']
splice_counts_M_piv_sub = splice_counts_M_piv[splice_counts_M_piv['YES']+splice_counts_M_piv['ALT']>=t].reset_index(drop=True)

splice_counts_P = splice_counts[splice_counts['allele']=='P'].reset_index(drop=True)
splice_counts_P_piv = splice_counts_P.pivot_table(index=['gene','cluster'], columns='splice_status', values='read_name').fillna(0).reset_index()
splice_counts_P_piv['ALT'] = splice_counts_P_piv['SKP'] + splice_counts_P_piv['SKP3SS'] + splice_counts_P_piv['SKP5SS'] + splice_counts_P_piv['UND']
splice_counts_P_piv_sub = splice_counts_P_piv[splice_counts_P_piv['YES']+splice_counts_P_piv['ALT']>=t].reset_index(drop=True)

# Merge the two alleles and filter for a sufficient number of ALT reads on at least one allele
fields = ['gene','cluster','YES','ALT']
splice_counts_merged = pd.merge(splice_counts_M_piv_sub[fields],splice_counts_P_piv_sub[fields], on=['gene','cluster'])

# Test for differential transcript usage between alleles
def fisher_exact(YES_M, ALT_M, YES_P, ALT_P):
    pvalue = scipy.stats.fisher_exact(np.array([[YES_M,ALT_M],[YES_P,ALT_P]]))[1]
    return(pvalue)

splice_counts_merged['pvalue'] = splice_counts_merged.apply(lambda row: fisher_exact(row.YES_x, row.ALT_x, row.YES_y, row.ALT_y),axis=1)
splice_counts_merged['FDR'] = multipletests(splice_counts_merged['pvalue'], alpha=0.05, method='fdr_bh')[1]
splice_counts_merged['freq_ALT_M'] = splice_counts_merged['ALT_x'] / (splice_counts_merged['YES_x']+splice_counts_merged['ALT_x'])
splice_counts_merged['freq_ALT_P'] = splice_counts_merged['ALT_y'] / (splice_counts_merged['YES_y']+splice_counts_merged['ALT_y'])
splice_counts_merged['delta'] = splice_counts_merged['freq_ALT_M'] - splice_counts_merged['freq_ALT_P']

# Merge results with annotations
splice_results = splice_counts_merged.merge(merge_df, on='cluster', how='left')

# Write results to file
splice_results.to_csv(sys.argv[-1], sep='\t', index=False, header=True)




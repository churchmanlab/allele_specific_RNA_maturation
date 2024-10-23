"""

Author: Karine Choquet

Date: May 9, 2024

This script will count the number of reads for each splicing order and allele for intron pairs


"""

import sys
import numpy as np
import pandas as pd
import pysam
from collections import Counter

import re
import math

import pybedtools
from pybedtools import BedTool

# CONFIG

pairs_df = pd.read_table(sys.argv[1])
allele_df = pd.read_table(sys.argv[2])

# Merge intron pairs with allele
pairs_allele_df = pd.merge(pairs_df,allele_df, left_on='read', right_on='readname')

# Get intron info
hg38_intron_df = pd.read_table(sys.argv[3])
gene_names_df = pd.read_table(sys.argv[4],header=None)
gene_names_df.columns = ['gene_name','gene_id']
hg38_intron_df['gene_id'] = hg38_intron_df['gene'].str.split("\\.").str[0]
hg38_intron_df = hg38_intron_df.merge(gene_names_df, on='gene_id')
hg38_intron_coord = hg38_intron_df.copy()[['chrom','start','end','gene','gene_name','count']]

################

def remove_duplicate_introns(pairs_df, hg38_intron_coord):
        
    # Remove duplicates intron pairs because they belong to different transcripts
    intron_groups = pairs_df[['gene','int1_count','int2_count']].drop_duplicates().reset_index(drop=True)

    intron_groups = intron_groups.merge(hg38_intron_coord, left_on=['gene','int1_count'], right_on=['gene','count']).rename(columns={'chrom':'chrom_1', 'start':'start_1', 'end':'end_1'})
    intron_groups = intron_groups.merge(hg38_intron_coord, left_on=['gene_name','gene','int2_count'], right_on=['gene_name','gene','count']).rename(columns={'chrom':'chrom_2', 'start':'start_2', 'end':'end_2'})

    intron_groups_nodup = intron_groups.sort_values(by=['gene','int1_count','int2_count']).drop_duplicates(subset=['chrom_1','start_1','end_1','chrom_2','start_2','end_2']).reset_index(drop=True)
    
    # Merge back with splicing orders
    pairs_df_nodup = pairs_df.merge(intron_groups_nodup, on=['gene','int1_count','int2_count'])
    
    return(pairs_df_nodup)



def get_splicing_order(int1_splice, int2_splice, strand):
    
    if strand == "+":
        if int1_splice == "YES" and int2_splice == "NO":
            splice_order = "in_order"
        elif int1_splice == "NO" and int2_splice == "YES":
            splice_order = "out_order"
    
    elif strand == "-":
        if int1_splice == "YES" and int2_splice == "NO":
            splice_order = "out_order"
        elif int1_splice == "NO" and int2_splice == "YES":
            splice_order = "in_order"
            
    return(splice_order)

def reset_intron_number(intron_number, strand, total_introns):
    
    if strand == "+":
        new_number = intron_number + 1
        
    elif strand == "-":
        new_number = total_introns - intron_number
        
    return(new_number)

def filter_intron_pairs(pairs_allele_df):
    # Select intermediate pairs
    interm_pairs = pairs_allele_df[(pairs_allele_df['int1_splice']!=pairs_allele_df['int2_splice']) & (pairs_allele_df['int1_splice'].isin(['YES','NO'])) & (pairs_allele_df['int2_splice'].isin(['YES','NO']))].reset_index(drop=True)

    # Get total number of introns per gene
    intron_numbers = hg38_intron_df[['gene','intron_total']].drop_duplicates().reset_index()

    # Add correct intron numbers
    interm_pairs = interm_pairs.merge(intron_numbers, on='gene')
    interm_pairs['int1_count_new'] = interm_pairs.apply(lambda row: reset_intron_number(row.int1_count, row.strand, row.intron_total), axis=1)
    interm_pairs['int2_count_new'] = interm_pairs.apply(lambda row: reset_intron_number(row.int2_count, row.strand, row.intron_total), axis=1)
    interm_pairs = interm_pairs.drop(columns=['int1_count','int2_count']).rename(columns={'int1_count_new':'int1_count','int2_count_new':'int2_count'})

    # Compute splicing order for each intron pair/read combination
    interm_pairs['splice_order'] = interm_pairs.apply(lambda row: get_splicing_order(row.int1_splice, row.int2_splice, row.strand),axis=1)

    # Count the number of reads with each order, per intron pair and allele
    interm_counts = pd.DataFrame(interm_pairs.groupby(['gene','int1_count','int2_count','splice_order','allele'])['read'].count()).reset_index()

    # Count the total number of reads per allele to filter for intron pairs with sufficient coverage
    interm_counts_tot = pd.DataFrame(interm_pairs.groupby(['gene','int1_count','int2_count','allele'])['read'].count()).reset_index()
    interm_counts_tot_piv = interm_counts_tot.pivot_table(index=['gene','int1_count','int2_count'],columns='allele', values='read').reset_index().fillna(0)
    interm_counts_tot_piv.columns = ['gene','int1_count','int2_count','ND','M','P']
    interm_counts_tot_sub = interm_counts_tot_piv[(interm_counts_tot_piv['M']>=10) & (interm_counts_tot_piv['P']>=10) & (interm_counts_tot_piv['M']+interm_counts_tot_piv['P']>=2*interm_counts_tot_piv['ND'])].reset_index(drop=True)

    # Add gene name to facilitate spotchecking
    interm_counts['gene_id'] = interm_counts['gene'].str.split("\\.").str[0]
    interm_counts = interm_counts.merge(gene_names_df, on='gene_id')

    # Merge with intron pairs that have sufficient coverage
    interm_counts = interm_counts.merge(interm_counts_tot_sub[['gene','int1_count','int2_count']], on=['gene','int1_count','int2_count'])
    interm_counts = interm_counts[interm_counts['allele'].isin(['M','P'])].reset_index(drop=True)

    return(interm_counts)



##########################

# Remove duplicate intron pairs within the same read
pairs_nodup = remove_duplicate_introns(pairs_allele_df, hg38_intron_coord)


# Compute splicing order for pairs that pass certain filters (see function)
pairs_counts = filter_intron_pairs(pairs_nodup)

pairs_counts.to_csv(sys.argv[-1], header=True, index=False, sep="\t")

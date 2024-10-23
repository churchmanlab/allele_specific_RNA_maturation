"""

This script will compare poly(A) tail length distributions and read end distributions between alleles

"""


import sys
import pandas as pd
from statsmodels.stats.multitest import multipletests
import scipy
from scipy import stats
import numpy as np
import math
from collections import Counter


# Load files
nano_df = pd.read_table(sys.argv[1])
read_map = pd.read_table(sys.argv[2], header=None)
read_map.columns = ['gene','readname','overlap']
allele_map = pd.read_table(sys.argv[3])

# Annotation files
hg38_intron_df = pd.read_table(sys.argv[4])
hg38_intron_df = hg38_intron_df[['gene','intron_total','strand']].drop_duplicates().reset_index(drop=True)

# Filter and merge files
df = nano_df.merge(read_map, on='readname').merge(allele_map, on='readname')

# Get genes with a minimum reads
t=20
cov_df = pd.DataFrame(df.groupby(['gene','allele'])['readname'].count()).reset_index().pivot_table(index='gene', columns='allele', values='readname').reset_index().fillna(0)
cov_df = cov_df[(cov_df['M']>t) & (cov_df['P']>t) & (cov_df['M']>2*cov_df['.']) & (cov_df['P']>2*cov_df['.'])].reset_index(drop=True)
genes_cov = cov_df['gene'].drop_duplicates().tolist()


# Test statistical difference between poly(A) tail length 
results_list = []
    
for gene in genes_cov:
        
    sub_df = df[df['gene']==gene].reset_index(drop=True)
    pA_M = sub_df[sub_df['allele']=='M']['polya_length']
    pA_P = sub_df[sub_df['allele']=='P']['polya_length']
            
    if len(pA_M) < t or len(pA_P) < t:
        continue
            
    M_mean = pA_M.mean()
    P_mean = pA_P.mean()
    M_CV = pA_M.std() / M_mean
    P_CV = pA_P.std() / P_mean
        
    pvalue_wilcox = scipy.stats.ranksums(pA_M,pA_P)[1]
    results_list.append([gene, len(pA_M), len(pA_P), M_mean, P_mean, M_CV, P_CV, pvalue_wilcox])
            
results_df = pd.DataFrame(results_list)
results_df.columns = ['gene','M_count','P_count','M_mean','P_mean','M_CV','P_CV','pvalue']
results_df['FDR'] = multipletests(results_df['pvalue'], alpha=0.05, method='fdr_bh')[1]

sig_pA_genes = results_df[results_df['FDR']<0.05]['gene'].drop_duplicates().tolist()
read_sig_df = df[df['gene'].isin(sig_pA_genes)].reset_index(drop=True)

    
# Write results to file
results_df.to_csv(sys.argv[5], sep="\t", header=True, index=False)
read_sig_df.to_csv(sys.argv[6], sep="\t", header=True, index=False)

# Determine whether there are differences in read end positions
read_ends = pd.read_table(sys.argv[7])
read_ends_pA = df.merge(read_ends, on=['readname','gene','allele'])

# Test statistical difference between read end distribution
results_list_read_ends = []

for gene in genes_cov:

    sub_df = read_ends_pA[read_ends_pA['gene']==gene].reset_index(drop=True)
    read_ends_M = sub_df[sub_df['allele']=='M']['start']
    read_ends_P = sub_df[sub_df['allele']=='P']['start']

    if len(read_ends_M) < t or len(read_ends_P) < t:
        continue

    M_mean = read_ends_M.mean()
    P_mean = read_ends_P.mean()

    pvalue_wilcox = scipy.stats.ranksums(read_ends_M,read_ends_P)[1]
    results_list_read_ends.append([gene, len(read_ends_M), len(read_ends_P), M_mean, P_mean, pvalue_wilcox])

results_df_read_ends = pd.DataFrame(results_list_read_ends)
results_df_read_ends.columns = ['gene','M_count','P_count','M_mean','P_mean','pvalue']
results_df_read_ends['FDR'] = multipletests(results_df_read_ends['pvalue'], alpha=0.05, method='fdr_bh')[1]

results_df_read_ends.to_csv(sys.argv[8], sep="\t", header=True, index=False)
# Get read ends from genes that are significant
read_ends_sig_genes = results_df_read_ends[results_df_read_ends['FDR']<0.05]['gene'].drop_duplicates().tolist()
read_ends_sig = read_ends_pA[read_ends_pA['gene'].isin(read_ends_sig_genes)].reset_index(drop=True)
read_ends_sig.to_csv(sys.argv[9], sep="\t", header=True, index=False)

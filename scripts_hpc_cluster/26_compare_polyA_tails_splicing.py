'''

This script will compare poly(A) tail length between alleles as a function of the read splicing status

'''


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
multi_introns_df = pd.read_table(sys.argv[4])

# Annotation files
hg38_intron_df = pd.read_table(sys.argv[5])
hg38_intron_df = hg38_intron_df[['gene','intron_total','strand']].drop_duplicates().reset_index(drop=True)

# Filter and merge files
df = nano_df.merge(read_map, on='readname').merge(allele_map, on='readname').merge(multi_introns_df, left_on=['readname','gene'], right_on=['read','gene'])

# Add an overall splicing status
df['splice_status_all'] = 'other'
df.loc[(df['splice_status'].str.contains("NO")) & (df['splice_status'].str.contains("YES")), 'splice_status_all'] = 'partially_spliced'
df.loc[(~df['splice_status'].str.contains("NO")) & (df['splice_status'].str.contains("YES")), 'splice_status_all'] = 'all_spliced'
df = df[df['splice_status_all'].isin(['partially_spliced','all_spliced'])].reset_index(drop=True)

# Get genes with a minimum reads
t=20
cov_df = pd.DataFrame(df.groupby(['gene','allele','splice_status_all'])['readname'].count()).reset_index().pivot_table(index=['gene','splice_status_all'], columns='allele', values='readname').reset_index().fillna(0)
cov_df = cov_df[(cov_df['M']>t) & (cov_df['P']>t) & (cov_df['M']>2*cov_df['.']) & (cov_df['P']>2*cov_df['.'])].reset_index(drop=True)
genes_cov = cov_df[['gene','splice_status_all']].drop_duplicates().reset_index(drop=True)


# Test statistical difference between poly(A) tail length 
results_list = []
    
for line in range(len(genes_cov)):

    gene = genes_cov.loc[line]['gene']
    splice_status = genes_cov.loc[line]['splice_status_all']
        
    sub_df = df[df['gene']==gene].reset_index(drop=True)
    pA_M = sub_df[sub_df['allele']=='M']['polya_length']
    pA_P = sub_df[sub_df['allele']=='P']['polya_length']
            
    if len(pA_M) < t or len(pA_P) < t:
        continue
            
    M_mean = pA_M.mean()
    P_mean = pA_P.mean()
        
    pvalue_wilcox = scipy.stats.ranksums(pA_M,pA_P)[1]
    results_list.append([gene, splice_status, len(pA_M), len(pA_P), M_mean, P_mean, pvalue_wilcox])
            
results_df = pd.DataFrame(results_list)
results_df.columns = ['gene','splice_status','M_count','P_count','M_mean','P_mean','pvalue']
results_df['FDR'] = multipletests(results_df['pvalue'], alpha=0.05, method='fdr_bh')[1]

sig_pA_genes = results_df[results_df['FDR']<0.05]['gene'].drop_duplicates().tolist()
read_sig_df = df[df['gene'].isin(sig_pA_genes)].reset_index(drop=True)

    
# Write results to file
results_df.to_csv(sys.argv[-2], sep="\t", header=True, index=False)
read_sig_df.to_csv(sys.argv[-1], sep="\t", header=True, index=False)

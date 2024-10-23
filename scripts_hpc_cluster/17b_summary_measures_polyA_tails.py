'''

This script will output summary measures for poly(A) tail length per gene, without allele information

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

# Annotation files
hg38_intron_df = pd.read_table(sys.argv[3])
hg38_intron_df = hg38_intron_df[['gene','intron_total','strand']].drop_duplicates().reset_index(drop=True)

# Filter and merge files
df = nano_df.merge(read_map, on='readname')

# Get genes with a minimum reads
t=20
cov_df = pd.DataFrame(df.groupby(['gene'])['readname'].count()).reset_index()
cov_df = cov_df[(cov_df['readname']>t)].reset_index(drop=True)
genes_cov = cov_df['gene'].drop_duplicates().tolist()


# Get summary measures for polyA tail length per gene
results_list = []
    
for gene in genes_cov:
        
    sub_df = df[df['gene']==gene].reset_index(drop=True)
    pA_count = len(sub_df['polya_length'])
    pA_avg = sub_df['polya_length'].mean()
    pA_std = sub_df['polya_length'].std()
    pA_CV = pA_std / pA_avg
            
    results_list.append([gene, pA_count, pA_avg, pA_std, pA_CV])
            
results_df = pd.DataFrame(results_list)
results_df.columns = ['gene','count','pA_avg','pA_std','pA_CV']
    
# Write results to file
results_df.to_csv(sys.argv[-1], sep="\t", header=True, index=False)

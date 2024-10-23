'''

This script will compare AS between alleles for HLA genes, using intron-centric splicing statuses

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
splice_info1 = pd.read_table(sys.argv[1])
splice_info2 = pd.read_table(sys.argv[2])
splice_info = pd.concat([splice_info1,splice_info2]).reset_index(drop=True)
path_df1 = pd.read_table(sys.argv[3])
path_df2 = pd.read_table(sys.argv[4])
path_df = pd.concat([path_df1,path_df2]).reset_index(drop=True)
allele_map1 = pd.read_table(sys.argv[5])
allele_map2 = pd.read_table(sys.argv[6])
allele_map = pd.concat([allele_map1,allele_map2]).reset_index(drop=True)


# Annotation files
hg38_intron_df = pd.read_table(sys.argv[7])
hg38_intron_df = hg38_intron_df[['gene','intron_total','strand']].drop_duplicates().reset_index(drop=True)

# Merge splice_info with annotation file and with allele
splice_info = splice_info.merge(hg38_intron_df, left_on=['strand','gene_name'], right_on=['strand','gene']).merge(allele_map, left_on='read_name', right_on='readname')


# Fix intron number
def fix_intron_number(intron_count, strand, intron_total):
    if strand == "+":
        intron_number = intron_count + 1
    elif strand == "-":
        intron_number = intron_total - intron_count
    return(intron_number)

splice_info['intron_number'] = splice_info.apply(lambda row: fix_intron_number(row.intron_count, row.strand, row.intron_total),axis=1)

# Get regions that were analyzed for splicing order
regions_df = path_df[['gene','gene_name','analyzed_introns']].drop_duplicates().reset_index(drop=True)
regions_df = regions_df.merge(hg38_intron_df, on='gene')
regions_df = regions_df[regions_df['gene_name']=='HLA-C'].reset_index(drop=True)


results_list = []
for region in range(len(regions_df)):
    gene = regions_df.loc[region]['gene']
    gene_name = regions_df.loc[region]['gene_name']
    analyzed_introns = regions_df.loc[region]['analyzed_introns']
 
    # Separate the 3 analyzed introns:
    intron_list = [int(a) for a in analyzed_introns.split("_")]

    # For each intron, compare SKP and YES
    for i in intron_list:
        sub_intron_df = splice_info[(splice_info['gene_name'] == gene) & (splice_info['intron_number']==i) & (splice_info['splice_status']!='NO')].reset_index(drop=True)
        sub_intron_df_grp = pd.DataFrame(sub_intron_df.groupby(['splice_status','allele','gene','intron_number'])['read_name'].count()).reset_index()
        yes_df = sub_intron_df_grp[sub_intron_df_grp['splice_status'] == 'YES'].reset_index(drop=True).pivot_table(index=['gene','intron_number'], columns='allele', values='read_name').reset_index().fillna(0)
        if len(yes_df) > 0 and 'M' in yes_df.columns and 'P' in yes_df.columns:
            M_yes = int(yes_df.loc[0]['M'])
            P_yes = int(yes_df.loc[0]['P'])
            other_df = sub_intron_df_grp[sub_intron_df_grp['splice_status'] != 'YES'].reset_index(drop=True).pivot_table(index=['gene','intron_number','splice_status'], columns='allele', values='read_name').reset_index().fillna(0)
            if len(other_df) > 0 and 'M' in other_df.columns and 'P' in other_df.columns:
                for j in range(len(other_df)):
                    splice_status = other_df.loc[j]['splice_status']
                    M_other = int(other_df.loc[j]['M'])
                    P_other = int(other_df.loc[j]['P'])
                    if (M_yes >= 10 and P_other >=10) or (P_yes >= 10 and M_other >= 10):
                        pvalue = scipy.stats.fisher_exact(np.array([[M_yes,M_other],[P_yes,P_other]]))[1]
                        results_list.append([gene_name, gene, i, splice_status, M_yes, M_other, P_yes, P_other, pvalue])

results_df = pd.DataFrame(results_list)
results_df.columns = ['gene_name', 'gene','intron_number','splice_status','count_M_YES', 'count_M_other', 'count_P_YES', 'count_P_other', 'pvalue']
results_df_nodup = results_df.drop_duplicates().reset_index(drop=True)
results_df_nodup['FDR'] = multipletests(results_df_nodup['pvalue'], alpha=0.05, method='fdr_bh')[1]

results_df_nodup.to_csv(sys.argv[-1], header=True, index=False, sep="\t")

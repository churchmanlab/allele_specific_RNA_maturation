'''

# This script will calculate the percent match from the alignment per gene

'''



import pandas as pd
import numpy as np
import sys

# Import stats file
stats_df = pd.read_table(sys.argv[1])

# Import read map
readmap = pd.read_table(sys.argv[2], header=None)
readmap.columns = ['gene','read','overlap']

# Import allele map
allele_df = pd.read_table(sys.argv[3])

# Merge the three together
both_df = stats_df.merge(readmap, on='read').merge(allele_df, left_on='read', right_on='readname')

# Add a category for assigned and unassigned alleles
both_df['category'] = 'unassigned'
both_df.loc[(both_df['allele'] == 'P') |(both_df['allele'] == 'M'),'category'] = 'assigned' 

# Get the median % match per gene
med_df = pd.DataFrame(both_df.groupby(['gene','category'])['match_percent'].median()).reset_index()

# Write to file
med_df.to_csv(sys.argv[4], index=False, header=True, sep="\t")

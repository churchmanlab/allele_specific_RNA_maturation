'''

This script will calculate the median read length per read splicing status

'''


import pandas as pd
import numpy as np
import sys


# Load splicing dictionary
splicing_status = pd.read_table(sys.argv[1])
stats_df = pd.read_table(sys.argv[2])

# Remove duplicates
splicing_status_nodup = splicing_status.drop_duplicates(subset=['read']).reset_index(drop=True)

# Merge
both_df = splicing_status_nodup.merge(stats_df, on="read")

# Calculate median length by read splicing status
med_df = pd.DataFrame(both_df.groupby('splicing_status')['read_length'].median()).reset_index()

# Write to file
med_df.to_csv(sys.argv[-1], sep="\t", header=True, index=False)

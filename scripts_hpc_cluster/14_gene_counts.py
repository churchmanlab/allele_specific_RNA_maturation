'''

This script will merge read to gene assignments with read to allele assignments, and count the number of reads per gene and allele, or per gene only

'''


import sys
import pandas as pd

# Read to gene assignments
readmap = pd.read_table(sys.argv[1], header=None)
readmap.columns = ['gene','readname','overlap']

# Read to allele assignments
allele_map = pd.read_table(sys.argv[2])

# Merge the two
both_map = pd.merge(readmap, allele_map, on="readname") 

# Count the number of reads per gene and allele
counts_df = pd.DataFrame(both_map.groupby(['gene','allele'])['readname'].count()).reset_index()


# Count all reads independently of allele
counts_df_bis = pd.DataFrame(readmap.groupby(['gene'])['readname'].count()).reset_index()

counts_df.to_csv(sys.argv[3], index=False, header=True, sep="\t")
counts_df_bis.to_csv(sys.argv[4], index=False, header=True, sep="\t")

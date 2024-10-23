'''

This script will extract the position of the 3'-end of each read

'''



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

bamFile = pybedtools.BedTool(sys.argv[1])
readmap= pd.read_table(sys.argv[3], header=None)
readmap.columns = ['gene','readname','overlap']
allele = sys.argv[4]

def get_read_end_bedtool(bamFile):

    bedFile = bamFile.bam_to_bed()
    bedFile_df = bedFile.to_dataframe(low_memory=False)
        
    read_end = []
        
    for i in range(0,len(bedFile_df)):

        chrom = str(bedFile_df['chrom'].iloc[i])
        start = bedFile_df['start'].iloc[i]
        end = bedFile_df['end'].iloc[i]
        read = bedFile_df['name'].iloc[i]
        score = bedFile_df['score'].iloc[i]
        strand = bedFile_df['strand'].iloc[i]

        if (strand == "-"):
            pos_1 = start
            pos_2 = start + 1

        if (strand == "+"):
            pos_1 = end - 1
            pos_2 = end

        read_end.append([chrom,pos_1,pos_2,read,score,strand])

    read_end_df = pd.DataFrame(read_end)
    read_end_df.columns = ['chrom','start','end','readname','score','strand']
    return read_end_df


# Apply function
read_ends_df = get_read_end_bedtool(bamFile)
read_ends_df = read_ends_df.merge(readmap, on='readname')
read_ends_df['allele'] = allele

read_ends_df.to_csv(sys.argv[2], sep="\t", index=False, header=False)


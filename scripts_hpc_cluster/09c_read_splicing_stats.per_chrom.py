"""

This script will assign an overall splicing status to each read: all_spliced, partially_spliced, all_unspliced

"""


import pandas as pd
import numpy as np
import sys
from collections import Counter


# Load splicing dictionary
splice_dict = np.load(sys.argv[1], encoding='latin1', allow_pickle=True).item()

def get_read_splicing_counts(read_junctions):
    
    read_property = []

    # loop through all reads in the dictionary
    for read in read_junctions.keys():

        # make a set for all intron pairs within a read
        # this will avoid duplicate pairs being called due to alternative splicing
        uniq_splice_pattern = set()

        # loop through all genes that has introns that a read maps to
        for gene in read_junctions[read].keys():

            # only go through genes that have 2 or more introns
            if (len(read_junctions[read][gene]) > 1 ):

                # characterize the number of spliced and unspliced introns in the read
                splice_status = [row[6] for row in read_junctions[read][gene]]
                splice_status_join = '_'.join(splice_status)
                status_count = Counter(splice_status)

                # only process the file if intron pattern hasn't been seen previously
                # for a gene that this read aligns to
                if (splice_status_join not in uniq_splice_pattern):
                    uniq_splice_pattern.add(splice_status_join)

                    spliced_count = status_count['YES']
                    unspliced_count = status_count['NO']

                    if (unspliced_count>1) and (spliced_count==0):
                        read_property.append([read, gene, 'all_unspliced'])

                    if (unspliced_count==0) and (spliced_count>1):
                        read_property.append([read, gene, 'all_spliced'])

                    if (unspliced_count>=1) and (spliced_count>=1):
                        read_property.append([read, gene, 'intermediate'])

    read_splicing_df = pd.DataFrame(read_property)
    read_splicing_df.columns = ['read','gene','splicing_status']        

    count = []
    spliced_count = len(read_splicing_df[read_splicing_df['splicing_status']=='all_spliced'])
    unspliced_count = len(read_splicing_df[read_splicing_df['splicing_status']=='all_unspliced'])
    intermediate_count = len(read_splicing_df[read_splicing_df['splicing_status']=='intermediate'])
    count.append([spliced_count, unspliced_count, intermediate_count])

    count_df = pd.DataFrame(count)
    count_df.columns = ['all_spliced','all_unspliced','intermediate']

    return read_splicing_df, count_df


splicing_status, splicing_counts = get_read_splicing_counts(splice_dict)

splicing_status.to_csv(sys.argv[2], sep="\t", header=False, index=False)
splicing_counts.to_csv(sys.argv[3], sep="\t", header=False, index=False)

"""

Author: Karine Choquet

Date: Feb 16, 2021

This script will determine the splicing order for each pair of consecutive intron spanned by a read

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

#####################################

# CONFIG


# Splicing dictionary for that chromosome
splicing_dict_file = sys.argv[1]

# Output files
out_df = sys.argv[2]

######################################

def get_intron_pairs_df(read_junctions):
    intron_pairs = []

    # loop through all reads in the dictionary
    for read in read_junctions.keys():

        # make a set for all intron pairs within a read
        # this will avoid duplicate pairs being called due to alternative splicing
        uniq_pairs = set()
        uniq_splice_pattern = set()
        
        # loop through all genes that has introns that a read maps to
        for gene in read_junctions[read].keys():

            # only go through genes that have 2 or more introns
            if (len(read_junctions[read][gene]) > 1 ):

                # characterize the number of spliced and unspliced introns in the read
                splice_status = [row[6] for row in read_junctions[read][gene]]
                splice_status_join = '_'.join(splice_status)
                status_count = Counter(splice_status)
                 
                spliced_count = status_count['YES']
                unspliced_count = status_count['NO']

                # build a dataframe of introns in the gene that map to this read
                # and are capable of being sequenced if the read has no splicing
                read_introns_df = pd.DataFrame(read_junctions[read][gene])
                read_introns_df.columns = ['chrom','start','end','intron_count','strand','read_overlap','splice_status']
                #read_introns_df = read_introns_df[read_introns_df['read_overlap'] > 0].sort_values('intron_count').reset_index(drop=True) # not applicable to polyA+ RNA
                read_introns_df = read_introns_df.sort_values('intron_count').reset_index(drop=True)

                # loop through introns that read maps to and find pairs
                prev_intron_count = -2    # counter for the start becuase no intron should have a negative count

                for i in range(len(read_introns_df)):
                    intron_count = read_introns_df.iloc[i]['intron_count']
                    intron_chrom = str(read_introns_df.iloc[i]['chrom'])
                    intron_start = str(read_introns_df.iloc[i]['start'])
                    intron_end = str(read_introns_df.iloc[i]['end'])
                    intron_strand = read_introns_df.iloc[i]['strand']
                    intron_splice = read_introns_df.iloc[i]['splice_status']
                    intron_coord = intron_chrom+'_'+intron_start+'_'+intron_end
                    intron_count = read_introns_df.iloc[i]['intron_count']

                    # if intron counts are sequential (one follows the next)
                    # it is a true intron pair (i.e. neighboring introns)
                    if (intron_count - prev_intron_count == 1):
                        intron_pair_coord = prev_intron_coord+'_'+intron_coord

                        # record information about the read pair only if the coordinates of this
                        # intron pair have not yet been seen
                        uniq_pairs.add(intron_pair_coord)
                        prev_intron_start = prev_intron_coord.split('_')[1]
                        prev_intron_end = prev_intron_coord.split('_')[2]

                        # append intron pair coordinate and splicing information to a list
                        intron_pairs.append([read,gene,prev_intron_count,intron_count,intron_strand,
                                                    prev_intron_splice, intron_splice])

                    # save information about this intron for the next pair
                    prev_intron_count = intron_count
                    prev_intron_coord = intron_coord
                    prev_intron_splice = intron_splice

    intron_pairs_df = pd.DataFrame(intron_pairs)
    intron_pairs_df.columns = ['read','gene','int1_count','int2_count','strand','int1_splice','int2_splice']        

    return intron_pairs_df


##############################################

# set all variables for analysis
min_overlap = 25

splice_dictionary = np.load(splicing_dict_file, encoding='latin1', allow_pickle=True).item()
# Get multi-introns dataset
intron_pairs_df = get_intron_pairs_df(splice_dictionary)
intron_pairs_df.to_csv(out_df, sep='\t', index=False, header=False)

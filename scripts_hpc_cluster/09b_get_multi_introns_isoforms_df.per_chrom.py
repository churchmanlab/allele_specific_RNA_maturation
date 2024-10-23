"""

Author: Karine Choquet

Date: Feb 16, 2021

This script will concatenate the splicing status of all introns within a transcript for each read

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
out_multi_introns_df = sys.argv[2]
out_multi_introns_counts = sys.argv[3]

######################################


def get_all_multi_introns(read_junctions, intron_min):

    multi_introns_list = []

    for read in read_junctions.keys():

        # make a set for all intron pairs within a read
        # this will avoid duplicate pairs being called due to alternative splicing
        uniq_splice_pattern = set()

        # loop through all genes that has introns that a read maps to
        for gene in read_junctions[read].keys():

            # only go through genes that have 2 or more introns
            if len(read_junctions[read][gene]) >= intron_min:

                # characterize the number of spliced and unspliced introns in the read
                strand = [row[4] for row in read_junctions[read][gene]][0]
                splice_status = [row[6] for row in read_junctions[read][gene]]
                intron_numbers = [row[3] for row in read_junctions[read][gene]]
                splice_status = [row[6] for row in read_junctions[read][gene]]
                intron_numbers = [row[3] for row in read_junctions[read][gene]]
                #if strand == "+":    
                #    splice_status = [row[6] for row in read_junctions[read][gene]]
                #    intron_numbers = [row[3] for row in read_junctions[read][gene]]
                #elif strand == "-":
                #    splice_status = reversed([row[6] for row in read_junctions[read][gene]])
                #    intron_numbers = reversed([row[3] for row in read_junctions[read][gene]])
                splice_status_join = '_'.join(splice_status)
                intron_numbers_join = '_'.join([str(i) for i in intron_numbers])
                status_count = Counter(splice_status)
                multi_introns_list.append([read,gene,strand,intron_numbers_join,splice_status_join])
                
    multi_introns_df = pd.DataFrame(multi_introns_list)
    
    if len(multi_introns_df)>0:
        multi_introns_df.columns = ['read','gene','strand','intron_numbers','splice_status']
        #multi_introns_counts = pd.DataFrame(multi_introns_df.groupby(['gene','strand','intron_numbers','splice_status'])['read'].count()).reset_index()
       	#multi_introns_counts.columns = ['gene','intron_numbers','splice_status','count']        
    else:
        multi_introns_df.columns = pd.DataFrame(['read','gene','strand','intron_numbers','splice_status'])

    return multi_introns_df


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
                 
                # only process the file if intron pattern hasn't been seen previously
                # for a gene that this read aligns to
                if (splice_status_join not in uniq_splice_pattern):
                    uniq_splice_pattern.add(splice_status_join)

                    spliced_count = status_count['YES']
                    unspliced_count = status_count['NO']

                    # build a dataframe of introns in the gene that map to this read
                    # and are capable of being sequenced if the read has no splicing
                    read_introns_df = pd.DataFrame(read_junctions[read][gene])
                    read_introns_df.columns = ['chrom','start','end','intron_count','strand','read_overlap','splice_status']
                    read_introns_df = read_introns_df[read_introns_df['read_overlap'] > 0].sort_values('intron_count').reset_index(drop=True)

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

                        # if intron counts are sequential (one follows the next)
                        # it is a true intron pair (i.e. neighboring introns)
                        if (intron_count - prev_intron_count == 1):
                            intron_pair_coord = prev_intron_coord+'_'+intron_coord

                            # record information about the read pair only if the coordinates of this
                            # intron pair have not yet been seen
                            if (intron_pair_coord not in uniq_pairs): 
                                uniq_pairs.add(intron_pair_coord)
                                prev_intron_start = prev_intron_coord.split('_')[1]
                                prev_intron_end = prev_intron_coord.split('_')[2]

                                # append intron pair coordinate and splicing information to a list
                                intron_pairs.append([read,intron_chrom,prev_intron_start,prev_intron_end,
                                                     int(intron_start),int(intron_end),intron_strand,
                                                    prev_intron_splice, intron_splice])

                        # save information about this intron for the next pair
                        prev_intron_count = intron_count
                        prev_intron_coord = intron_coord
                        prev_intron_splice = intron_splice

    intron_pairs_df = pd.DataFrame(intron_pairs)
    intron_pairs_df.columns = ['read','chrom','int1_start','int1_end','int2_start','int2_end','strand','int1_splice','int2_splice']        

    return intron_pairs_df


##############################################

# set all variables for analysis
min_overlap = 25

splice_dictionary = np.load(splicing_dict_file, encoding='latin1', allow_pickle=True).item()
# Get multi-introns dataset
multi_introns_df = get_all_multi_introns(splice_dictionary, 2)
multi_introns_df.to_csv(out_multi_introns_df, sep='\t', index=False, header=False)


# Summarize and count
if len(multi_introns_df)>0:
    multi_introns_counts = pd.DataFrame(multi_introns_df.groupby(['gene','strand','intron_numbers','splice_status'])['read'].count()).reset_index()
    multi_introns_counts.columns = ['gene','strand','intron_numbers','splice_status','count']
    multi_introns_counts.to_csv(out_multi_introns_counts, sep='\t', index=False, header=False)


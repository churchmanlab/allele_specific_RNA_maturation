"""

Author: Karine Choquet

Date: January 29, 2022

This script will identify possible splicing orders on each allele and assign them a score based on the isoform counts from direct chromatin RNA sequencing

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

import itertools
from more_itertools import consecutive_groups

from interruptingcow import timeout

#############

min_reads = 10

gene_names_df = pd.read_table(sys.argv[1])
gene_names_df.columns = ['gene_name','gene_id']
intron_df = pd.read_table(sys.argv[2], header=None, names=['chrom','start','end','intron_name','score','strand'], dtype={'chrom':str, 'start':int, 'end':int, 'intron_name':str, 'score':str, 'strand':str})
intron_df['gene_id'] = intron_df['intron_name'].str.split("\\.").str[0]
intron_total_df = pd.DataFrame(intron_df.groupby('gene_id')['intron_name'].count()).rename(columns={'intron_name':'intron_total'})
gene_names_df = gene_names_df.merge(intron_total_df, on='gene_id')


# Load multi_intron_df
multi_intron_df = pd.read_table(sys.argv[3])
multi_intron_df['gene_id'] = multi_intron_df['gene'].str.split("\\.").str[0]
multi_intron_df = multi_intron_df.merge(gene_names_df, on=['gene_id'])

# Get allele mapping
allele_map_df = pd.read_table(sys.argv[4])

# Merge with multi_intron_df
multi_intron_df = multi_intron_df.merge(allele_map_df, left_on='read', right_on='readname').replace('.','ND')

# Count number of reads per splicing status per allele
multi_intron_counts = pd.DataFrame(multi_intron_df.groupby(['gene','gene_name','strand','intron_total','intron_numbers','splice_status','allele'])['readname'].count()).reset_index().rename(columns={'readname':'count'})
allele_counts = pd.DataFrame(multi_intron_counts.groupby(['gene','allele'])['count'].sum()).pivot_table(index='gene', columns='allele', values='count').reset_index()
allele_counts_sub = allele_counts[((allele_counts['M'] >= min_reads) & (allele_counts['P'] >= min_reads))].reset_index(drop=True)
transcript_list = allele_counts_sub['gene'].drop_duplicates().tolist()

####################

def get_paths(pattern_dict, spliced_counts_dict, introns_of_interest_list, introns_of_interest_fix):
    
    # Initiate path_dict
    path_dict = {}
    for n in range(len(introns_of_interest_list)):
        path_dict[n] = {}

    # For each combination of isoforms between levels (e.g. number of splicing events),
    # retrieve the frequencies and calculate scores
    # Below, the following scores are calculated:
    # freq: number of reads for that isoform / total number of reads for that level
    # sub_path_score : freq for that isoform * freq for the isoform from which it is derived
    # path_score: freq for that isoform * path_score from the previous isoform (so multiplicative frequencies)
    # full_path_score: the path_score for the last level of that path
    for levels in itertools.product(pattern_dict.keys(), pattern_dict.keys()):
        level1 = levels[0]
        level2 = levels[1]
        if level1 != level2: # we don't want to compare isoforms from the same level

            # Iterate through each pair of isoforms that are from different levels
            for pair in itertools.product(pattern_dict[level1].keys(), pattern_dict[level2].keys()):
                pattern1 = pair[0].split("_")
                pattern2 = pair[1].split("_")

                # Get the splicing level of the isoform
                unspliced_introns_pattern1 = len([i for i, x in enumerate(pattern1) if x == "NO"])
                unspliced_introns_pattern2 = len([i for i, x in enumerate(pattern2) if x == "NO"])
                spliced_introns_pattern2 = Counter(pattern2)['YES'] + Counter(pattern2)['SKP']

                # Define the level that will be used below
                level = level1

                # retrieve the positions of the difference(s) between the two patterns
                diff_index_list = [i for i, x in enumerate(pattern2) if pattern1[i]!=x]

                # If pattern2 has one more spliced intron than pattern1:
                if len(diff_index_list) == 1:
                    diff_index = diff_index_list[0]
                    if pattern1[diff_index] == "NO" and pattern2[diff_index] == "YES":
                        count_pattern1 = pattern_dict[level1][pair[0]]
                        count_pattern2 = pattern_dict[level2][pair[1]]

                        if level == 0: # this means we are starting from the unspliced isoform
                            # define a new splicing path
                            new_intron_spliced = str(introns_of_interest_list[diff_index])
                            path_name = new_intron_spliced + "->"
                            freq = count_pattern2 / spliced_counts_dict[level+1]
                            path_score = freq
                            sub_path_score = freq
                            path_dict[level][path_name] = [pattern1, pattern2, pair[0], pair[1], count_pattern1, count_pattern2, "spliced", level, new_intron_spliced, spliced_counts_dict[level+1], freq, sub_path_score, path_score]

                        elif level > 0 and level < len(introns_of_interest_list)-1: # this means we are at an intermediate isoform
                            for k in path_dict[level-1].keys():
                                if path_dict[level-1][k][1] == pattern1:
                                    new_intron_spliced = str(introns_of_interest_list[diff_index])
                                    if unspliced_introns_pattern2 == 0:
                                        path_name = str(k) + new_intron_spliced
                                    elif unspliced_introns_pattern2 > 0:
                                        path_name = str(k) + new_intron_spliced + "->"
                                    freq = count_pattern2 / spliced_counts_dict[level+1]
                                    path_score = path_dict[level-1][k][-1] * freq
                                    sub_path_score = path_dict[level-1][k][-3] * freq # only the frequency from the previous level and this level
                                    path_dict[level][path_name] = [pattern1, pattern2, pair[0], pair[1], count_pattern1, count_pattern2, "spliced", level, new_intron_spliced, spliced_counts_dict[level+1], freq, sub_path_score, path_score]

                        elif level == len(introns_of_interest_list)-1: # this means we are at the fully spliced isoform
                            for k in path_dict[level-1].keys():
                                if path_dict[level-1][k][1] == pattern1:
                                    new_intron_spliced = str(introns_of_interest_list[diff_index])
                                    path_name = str(k) + new_intron_spliced
                                    freq = 1
                                    path_score = path_dict[level-1][k][-1] * freq
                                    sub_path_score = path_dict[level-1][k][-3] * freq
                                    path_dict[level][path_name] = [pattern1, pattern2, pair[0], pair[1], count_pattern1, count_pattern2, "spliced", level, new_intron_spliced, spliced_counts_dict[level+1], freq, sub_path_score, path_score]


    # Now loop through the dictionary to match each pair with the possible final paths
    # The final level contains only final paths, so first retrieve those
    try:
        final_level = list(path_dict.keys())[-1]
        final_level_df = pd.DataFrame.from_dict(path_dict[final_level], orient='index').reset_index()
        final_level_df.columns = ['path_name','pattern1_list','pattern2_list','pattern1','pattern2','count_pattern1','count_pattern2','event_type','level','new_intron_spliced','total_counts_level','freq','sub_path_score','path_score']
        final_level_df = final_level_df.drop(['pattern1_list','pattern2_list'],axis=1)
        final_level_df['full_path'] = final_level_df['path_name']
        final_df = final_level_df.copy()

        # Iterate through each of the levels to match the partial paths with all the possible final paths
        for lev in list(reversed(list(path_dict.keys())))[:-1]:
            # For the two final levels, merge the second last with the last to retrieve the final path and add the final path score
            if lev == final_level:
                df1 = pd.DataFrame.from_dict(path_dict[lev], orient='index').reset_index()
                df2 = pd.DataFrame.from_dict(path_dict[lev-1], orient='index').reset_index()
                df1.columns = ['path_name','pattern1_list','pattern2_list','pattern1','pattern2','count_pattern1','count_pattern2','event_type','level','new_intron_spliced','total_counts_level','freq','sub_path_score','path_score']
                df2.columns = ['path_name','pattern1_list','pattern2_list','pattern1','pattern2','count_pattern1','count_pattern2','event_type','level','new_intron_spliced','total_counts_level','freq','sub_path_score','path_score']

                fields = ['pattern1','path_name']

                new_df = df2.merge(df1[fields], left_on='pattern2', right_on='pattern1', how='left')
                new_df = new_df.rename(columns={'path_name_x':'path_name','path_name_y':'full_path','pattern1_x':'pattern1'}).drop(['pattern1_y','pattern1_list','pattern2_list'], axis=1)

                new_df = new_df.fillna(0)

                # If full_path is null, that means that it wasn't present in the level above and therefore the full path
                # is likely to be the current path, so replace it that way
                new_df_sub1 = new_df[new_df['full_path']==0].reset_index(drop=True)
                new_df_sub2 = new_df[new_df['full_path']!=0].reset_index(drop=True)

                new_df_sub1['full_path'] = new_df_sub1['path_name']

                new_df = pd.concat([new_df_sub1, new_df_sub2]).reset_index(drop=True)

                final_df = pd.concat([final_df, new_df]).reset_index(drop=True)

             # For any previous levels, repeat the previous steps
            elif (lev - 1) >= 0:
                df1 = new_df.copy()
                df2 = pd.DataFrame.from_dict(path_dict[lev-1], orient='index').reset_index()

                df2.columns = ['path_name','pattern1_list','pattern2_list','pattern1','pattern2','count_pattern1','count_pattern2','event_type','level','new_intron_spliced','total_counts_level','freq','sub_path_score','path_score']

                fields = ['pattern1','full_path']

                new_df = df2.merge(df1[fields], left_on='pattern2', right_on='pattern1', how='left')
                new_df = new_df.rename(columns={'pattern1_x':'pattern1'}).drop(['pattern1_y','pattern1_list','pattern2_list'], axis=1)

                new_df = new_df.fillna(0)

                new_df_sub1 = new_df[new_df['full_path']==0].reset_index(drop=True)
                new_df_sub2 = new_df[new_df['full_path']!=0].reset_index(drop=True)

                new_df_sub1['full_path'] = new_df_sub1['path_name']

                new_df = pd.concat([new_df_sub1, new_df_sub2]).reset_index(drop=True)

                final_df = pd.concat([final_df, new_df]).reset_index(drop=True)
        # Now make sure that the start of the full path is the same as the path name, since merging on patterns as above
        # will give rows where that is not the case
        final_final_df = final_df[final_df.apply(lambda row: row.full_path.startswith(row.path_name), axis=1)].drop_duplicates().sort_values(by='level').reset_index(drop=True)


        # Get the final score for the path and express it so that the total score of all final isoforms is 1
        last_isos = final_final_df[~final_final_df['path_name'].str.endswith("->")][['full_path','path_score']].drop_duplicates().reset_index(drop=True)
        #last_isos2 = final_final_df[(final_final_df['path_name']==final_final_df['full_path']) & (final_final_df['full_path'].str.endswith("->"))][['full_path','path_score']].drop_duplicates().reset_index(drop=True)
        #last_isos = pd.concat([last_isos1, last_isos2]).reset_index(drop=True)
        last_isos['full_path_score'] = last_isos['path_score'] / last_isos['path_score'].sum()
        last_isos = last_isos.drop('path_score',axis=1).sort_values(by='full_path_score', ascending=False).reset_index(drop=True)
        last_isos['rank'] = last_isos.index + 1

        # Modify the levels and rows where the unspliced isoform is pattern2
        final_final_df['level'] = final_final_df['level'] + 1

        NO_df = final_final_df[final_final_df['level']==1].reset_index(drop=True)
        NO_df['pattern2'] = NO_df['pattern1']
        NO_df['count_pattern2'] = NO_df['count_pattern1']
        NO_df['level'] = 0
        NO_df['new_intron_spliced'] = 0

        final_final_df = pd.concat([final_final_df,NO_df]).sort_values(by='level').reset_index(drop=True)

        final_final_df = final_final_df.merge(last_isos, on='full_path')

        final_final_df['analyzed_introns'] = introns_of_interest_fix

        return(final_final_df)
            
    except ValueError:
        pass
    

def get_score_per_path_non_consec(multi_introns_df, total_introns_of_interest, introns_of_interest_list_tmp, n_introns_gene, strand):

    # Make a dictionary with the patterns and the counts and another with the number of introns spliced and the counts    
    pattern_dict = {}
    results_list = []
    spliced_counts_dict = {}

    if strand == "-":
        introns_of_interest_list = [str(n_introns_gene - int(i)) for i in introns_of_interest_list_tmp]
        introns_of_interest_fix = "_".join(introns_of_interest_list)
    elif strand == "+":
        introns_of_interest_list = [str(int(i)+1) for i in introns_of_interest_list_tmp]
        introns_of_interest_fix = "_".join(introns_of_interest_list)
        
    # Initiate pattern_dict and spliceed_counts_dict
    alleles = ['M','P','ND']
    for allele in alleles:
        pattern_dict[allele] = {}
        for n in range(len(introns_of_interest_list)+1):
            pattern_dict[allele][n] = {}
        spliced_counts_dict[allele] = {}

    # Iterate to get counts for isoforms
    for row in range(len(multi_introns_df)):
        gene_name = multi_introns_df.loc[row]['gene_name']
        splice_status_temp = multi_introns_df.loc[row]['splice_status']
        allele = multi_introns_df.loc[row]['allele']
        intron_numbers_list_tmp = multi_introns_df.loc[row]['intron_numbers'].split("_")
        if strand == "-":
            intron_numbers_list = [str(n_introns_gene - int(i)) for i in intron_numbers_list_tmp]
            intron_numbers = "_".join(intron_numbers_list)
        elif strand == "+":
            intron_numbers_list = [str(int(i)+1) for i in intron_numbers_list_tmp]
            intron_numbers = "_".join(intron_numbers_list)

        # Determine if all introns of interest are present in the splice status
        common_introns = [a for a in introns_of_interest_list if a in intron_numbers_list]

        if len(common_introns) == total_introns_of_interest:
            introns_of_interest_pos = [i for i, x in enumerate(intron_numbers_list) if x in introns_of_interest_list]
            splice_status_list_temp1 = splice_status_temp.split("_")
            splice_status_list_temp = [splice_status_list_temp1[a] for a in introns_of_interest_pos]
            splice_status_list = ["SKP" if "SKP" in a else a for a in splice_status_list_temp]
            splice_status = "_".join(splice_status_list)
            pattern_count = multi_introns_df.loc[row]['count']
            skipped_count = Counter(splice_status_list)['SKP']
            spliced_count = Counter(splice_status_list)['YES']
            unspliced_count = Counter(splice_status_list)['NO']
            undetermined_count = Counter(splice_status_list)['UND']
            
            if skipped_count == 0 and undetermined_count == 0: # no skipped or undetermined introns
                level = spliced_count
                if skipped_count < total_introns_of_interest:
                    if splice_status not in pattern_dict[allele][level].keys():
                        pattern_dict[allele][level][splice_status] = pattern_count
                    elif splice_status in pattern_dict[allele][level].keys():
                        pattern_dict[allele][level][splice_status] += pattern_count
                    if level not in spliced_counts_dict[allele].keys():
                        spliced_counts_dict[allele][level] = pattern_count
                    else:
                        spliced_counts_dict[allele][level] += pattern_count


    for allele in alleles:
        if len(pattern_dict[allele][0]) == 0:
            unspliced_iso = "_".join(["NO" for i in range(len(introns_of_interest_list))])
            pattern_dict[allele][0][unspliced_iso] = 0
            
    # Filter for a certain number of reads at each intermediate level
    level_threshold = 10
    good_levels = {}
    for allele in ['M','P']:
        if allele == "M":
            other_allele = "P"
        elif allele == "P":
            other_allele = "M"
        good_levels[allele] = []
        for level in sorted(list(spliced_counts_dict[allele].keys()))[1:-1]: # exclude the first and the final levels (all unspliced or all spliced)
            try:
                if spliced_counts_dict[allele][level] >= level_threshold:
                    try:
                        if spliced_counts_dict[allele][level] > 2*spliced_counts_dict['ND'][level]:
                            good_levels[allele].append(level)
                    except KeyError:
                        good_levels[allele].append(level)
            except KeyError:
                pass

    print(good_levels)
    mega_path_list = []
    if len(good_levels['M']) == total_introns_of_interest - 1 and len(good_levels['P']) == total_introns_of_interest - 1: # there is at least one read at each intermediate level
        for allele in ['M','P']:
            path_df = get_paths(pattern_dict[allele], spliced_counts_dict[allele], introns_of_interest_list, introns_of_interest_fix)
            if len(path_df) > 0:
                path_df['allele'] = allele
                path_df['gene_name'] = gene_name
                mega_path_list.append(path_df)
    
    if len(mega_path_list) == 2:
        mega_path_df = pd.concat(mega_path_list)
        
        return(mega_path_df)
    
    else:
        pass


################################

results_df_list = []


for transcript in transcript_list:
    
    gene_results_list = []
    print(transcript)
   
    # Retrieve necessary information for that transcript
    multi_introns_df = multi_intron_counts[multi_intron_counts['gene']==transcript].reset_index(drop=True)
    n_introns_gene = multi_intron_counts[multi_intron_counts['gene']==transcript].reset_index(drop=True)['intron_total'].drop_duplicates().tolist()[0]
    strand = multi_intron_counts[multi_intron_counts['gene']==transcript].reset_index(drop=True)['strand'].drop_duplicates().tolist()[0]
    introns_of_interest = range(0,n_introns_gene+1)
    total_introns_of_interest = n_introns_gene
    
    # Iterate through groups of 3-6 introns in both replicates to retrieve the longest possible paths
    for i in range(3,4):
        for x in range(total_introns_of_interest-i+1):
            introns_of_interest_sub = [introns_of_interest[a] for a in range(x,x+i)]
            # Compute splicing order paths for the introns of interest
            results_df = get_score_per_path_non_consec(multi_introns_df, i, introns_of_interest_sub, n_introns_gene, strand)

    
            if results_df is not None:
                results_df['gene'] = transcript
                results_df['n_analyzed_introns'] = i
                gene_results_list.append(results_df)
                
    if len(gene_results_list)>0:
        gene_results_df = pd.concat(gene_results_list).reset_index(drop=True)
        results_df_list.append(gene_results_df)
    
   
final_results_df = pd.concat(results_df_list).reset_index(drop=True)

final_results_df.to_csv(sys.argv[-1], sep="\t", header=True, index=False)


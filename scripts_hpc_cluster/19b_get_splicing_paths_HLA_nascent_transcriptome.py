'''

This script will compute splicing order from HLA nascent transcriptome alignments

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

import itertools
from more_itertools import consecutive_groups

# Load BED file from transcriptome alignment and merge with allele map
HLA_reads_bed = pd.read_table(sys.argv[1], header=None)
HLA_reads_bed.columns = ['transcript','start','end','readname','score','strand','cigar']
allele_map = pd.read_table(sys.argv[2])
HLA_reads_bed = HLA_reads_bed.merge(allele_map, on='readname')

# Load exon and intron annotations
hg38_gene_df = pd.read_table(sys.argv[3], header=None)
gene_names_df = pd.read_table(sys.argv[4], header=None)
hg38_gene_df.columns = ['chrom','start','end','gene','feature','strand']
gene_names_df.columns = ['gene_name', 'gene_id']
hg38_gene_df['gene_id'] = hg38_gene_df['gene'].str.split("\\.").str[0]
hg38_gene_df = hg38_gene_df.merge(gene_names_df, on='gene_id')

# Get coordinates for introns and exons within nascent isoforms
def coord_intermediate_isos(gene_name):
    
    gene_df = hg38_gene_df[hg38_gene_df['gene_name']==gene_name].reset_index(drop=True)
    
    # Make a dictionary with all the exon and intron coordinates
    exon_dict = {}
    intron_dict = {}
    
    # Make a dictionary that will contain the intermediate isoforms and their sequences
    iso_dict = {}
    
    for i in range(len(gene_df)):
        feature_tmp = gene_df.loc[i]['feature']
        chrom = gene_df.loc[i]['chrom']
        start = gene_df.loc[i]['start']
        end = gene_df.loc[i]['end']
        strand = gene_df.loc[i]['strand']
        
        if feature_tmp != "gene":
            feature = feature_tmp.split("_")[0]
            feature_number = int(feature_tmp.split("_")[1])
            
            if feature == "exon":
                exon_dict[feature_number] = [start+1,end]
                
            elif feature == "intron":
                intron_dict[feature_number] = [start+1,end]
                
    # Get all possible combinations of retained introns
    n_introns = len(intron_dict.keys())
    
    coord_list = []
    
    # Add fully spliced isoform
    new_end = 0
    for pos in range(1,n_introns+2):
        l_feature = exon_dict[pos][1] - exon_dict[pos][0]
        new_start = new_end + 1
        new_end = new_start + l_feature
        exon_name = "exon" + '_' + str(pos)
        iso_name = gene_name + ":" + "all_spliced"
        coord_list.append([iso_name,exon_name,new_start,new_end])
    
    for level in range(1, n_introns+1):
        for isoform in itertools.combinations(range(1,n_introns+1), level):
            # give isoform a name - for easier integration into downstream analyses, use transcript intron numbers
            # not genomic intron numbers
            if strand == "+":
                iso_join = ",".join([str(a) for a in isoform])
                iso_name = gene_name + ":" + iso_join
            elif strand == "-":
                iso_join = ",".join([str(n_introns-a+1) for a in isoform])
                iso_name = gene_name + ":" + iso_join
            # initiate variables
            new_end = 0
            # Retrieve the sequence for the exons and introns (since this is by genomic position,
            # it will not be stranded and will correspond to the plus strand)
            for pos in range(1,n_introns+1):
                # add every exon in order
                l_feature = exon_dict[pos][1] - exon_dict[pos][0]
                new_start = new_end + 1
                new_end = new_start + l_feature
                exon_name = "exon" + '_' + str(pos)
                coord_list.append([iso_name,exon_name,new_start,new_end])
                # add intron only if it's still present in the isoform
                if pos in isoform:
                    l_feature = intron_dict[pos][1] - intron_dict[pos][0]
                    new_start = new_end + 1
                    new_end = new_start + l_feature
                    intron_name = "intron" + '_' + str(pos)
                    coord_list.append([iso_name,intron_name,new_start,new_end])
            # add last exon
            l_feature = exon_dict[n_introns+1][1] - exon_dict[n_introns+1][0]
            new_start = new_end + 1
            new_end = new_start + l_feature
            exon_name = "exon" + '_' + str(n_introns+1)
            coord_list.append([iso_name,exon_name,new_start,new_end])
    
    coord_df = pd.DataFrame(coord_list)
    coord_df.columns = ['isoform','feature','start','end']
                
    return(coord_df)



def get_read_start_overlap(transcript, aln_start, coord_df):
    
    iso_coord = coord_df[coord_df['isoform']==transcript].reset_index(drop=True)
    try:
        start_feature = iso_coord[(iso_coord['start']<=aln_start) & (iso_coord['end']>=aln_start)]['feature'].tolist()[0]
        return(start_feature)
    except IndexError:
        if aln_start == 0:
            return("exon_1")
        else:
            return("error")
        
def get_read_end_overlap(transcript, aln_end, coord_df):
    
    iso_coord = coord_df[coord_df['isoform']==transcript].reset_index(drop=True)
    try:
        end_feature = iso_coord[(iso_coord['start']<=aln_end) & (iso_coord['end']>=aln_end)]['feature'].tolist()[0]
        return(end_feature)
    except IndexError:
        return("error")

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
        introns_of_interest_list = [n_introns_gene - int(i) for i in introns_of_interest_list_tmp]
        introns_of_interest_fix = "_".join([str(a) for a in introns_of_interest_list])
    elif strand == "+":
        introns_of_interest_list = [int(i)+1 for i in introns_of_interest_list_tmp]
        introns_of_interest_fix = "_".join([str(a) for a in introns_of_interest_list])

    # Initiate pattern_dict and spliceed_counts_dict
    alleles = ['M','P','ND']
    for allele in alleles:
        pattern_dict[allele] = {}
        for n in range(len(introns_of_interest_list)+1):
            pattern_dict[allele][n] = {}
        spliced_counts_dict[allele] = {}

    # Iterate to get counts for isoforms
    for row in range(len(multi_introns_df)):
        transcript = multi_introns_df.loc[row]['transcript']
        gene_name = transcript.split(":")[0]
        allele = multi_introns_df.loc[row]['allele']
        retained_introns = transcript.split(":")[1].split(",")
        start_feature = multi_introns_df.loc[row]['start_feature'].split("_")[0]
        start_pos = int(multi_introns_df.loc[row]['start_feature'].split("_")[1])
        end_feature = multi_introns_df.loc[row]['end_feature'].split("_")[0]
        end_pos = int(multi_introns_df.loc[row]['end_feature'].split("_")[1])
        
        # Get the introns covered by the read
        if start_feature == "exon" and end_feature == "exon":
            intron_numbers_list_tmp = [x for x in range(start_pos,end_pos)]
        elif start_feature == "intron" and end_feature == "exon":
            intron_numbers_list_tmp = [x for x in range(start_pos+1,end_pos)]
        elif start_feature == "exon" and end_feature == "intron":
            intron_numbers_list_tmp = [x for x in range(start_pos,end_pos)]
        elif start_feature == "intron" and end_feature == "intron":
            intron_numbers_list_tmp = [x for x in range(start_pos+1,end_pos)]
        
        # Keep only reads covering more than 2 introns
        if len(intron_numbers_list_tmp) > 2:
            if strand == "-":
                intron_numbers_list = [n_introns_gene - int(i)+1 for i in intron_numbers_list_tmp]
                intron_numbers = "_".join([str(a) for a in intron_numbers_list])
            elif strand == "+":
                intron_numbers_list = [int(i)+1 for i in intron_numbers_list_tmp]
                intron_numbers = "_".join([str(a) for a in intron_numbers_list])
        
            # Retrieve splicing status of the read
            if retained_introns[0] != "all_spliced":
                retained_introns_bis = [int(a) for a in retained_introns]
                splice_status_all = []
                for i in intron_numbers_list:
                    if i in retained_introns_bis:
                        splice_status_all.append("NO")
                    else:
                        splice_status_all.append("YES")
            

                # Determine if all introns of interest are present in the splice status
                common_introns = [a for a in introns_of_interest_list if a in intron_numbers_list]
                if len(common_introns) == total_introns_of_interest:
                    introns_of_interest_pos = [i for i, x in enumerate(intron_numbers_list) if x in introns_of_interest_list]
                    splice_status_list = [splice_status_all[a] for a in introns_of_interest_pos]
                    splice_status = "_".join(splice_status_list)
                    pattern_count = multi_introns_df.loc[row]['count']
                    spliced_count = Counter(splice_status_list)['YES']
                    unspliced_count = Counter(splice_status_list)['NO']
                    level = spliced_count
                    if splice_status not in pattern_dict[allele][level].keys():
                        pattern_dict[allele][level][splice_status] = pattern_count
                    elif splice_status in pattern_dict[allele][level].keys():
                        pattern_dict[allele][level][splice_status] += pattern_count
                    if level not in spliced_counts_dict[allele].keys():
                        spliced_counts_dict[allele][level] = pattern_count
                    else:
                        spliced_counts_dict[allele][level] += pattern_count
            elif retained_introns[0] == "all_spliced":
                pattern_count = multi_introns_df.loc[row]['count']
                for allele in alleles:
                    pattern_dict[allele][3]["YES_YES_YES"] = pattern_count
                    spliced_counts_dict[allele][3] = pattern_count

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

    mega_path_list = []
    if len(good_levels['M']) == total_introns_of_interest - 1 and len(good_levels['P']) == total_introns_of_interest - 1: # there is at least one read at each intermediate level
        for allele in ['M','P']:
            path_df = get_paths(pattern_dict[allele], spliced_counts_dict[allele], introns_of_interest_list, introns_of_interest_fix)
            if len(path_df) > 0:
                path_df['allele'] = allele
                #path_df['gene_name'] = gene_name
                mega_path_list.append(path_df)

    if len(mega_path_list) == 2:
        mega_path_df = pd.concat(mega_path_list)

        return(mega_path_df) 

    else:
        pass


#######################


# Combine the three HLA genes of interest
coord_HLAA = coord_intermediate_isos("HLA-A")
coord_HLAB = coord_intermediate_isos("HLA-B")
coord_HLAC = coord_intermediate_isos("HLA-C")

coord_HLA_df = pd.concat([coord_HLAA,coord_HLAB,coord_HLAC]).reset_index(drop=True)


# Count reads
HLA_reads_bed['start_feature'] = HLA_reads_bed.apply(lambda row: get_read_start_overlap(row.transcript, row.start, coord_HLA_df),axis=1)
HLA_reads_bed['end_feature'] = HLA_reads_bed.apply(lambda row: get_read_end_overlap(row.transcript, row.end, coord_HLA_df),axis=1)
HLA_reads_counts = pd.DataFrame(HLA_reads_bed.groupby(['transcript','allele','start_feature','end_feature'])['readname'].count()).reset_index().rename(columns={'readname':'count'}).replace('.','ND')



# Compute paths 
strand_dict = {'HLA-A':'+','HLA-B':'-','HLA-C':'-'}
transcript_list = ['HLA-A','HLA-B','HLA-C']

total_introns_of_interest = 7
introns_of_interest = range(0,total_introns_of_interest+1)

results_df_list = []

for transcript in transcript_list:

    gene_results_list = []
    print(transcript)

    # Retrieve necessary information for that transcript
    multi_introns_df = HLA_reads_counts[HLA_reads_counts['transcript'].str.startswith(transcript)].reset_index(drop=True)
    n_introns_gene = total_introns_of_interest
    strand = strand_dict[transcript]

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


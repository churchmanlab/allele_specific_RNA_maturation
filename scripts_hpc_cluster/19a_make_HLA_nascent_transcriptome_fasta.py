"""

Date: November 15, 2022

Author: Karine Choquet

This script will generate an allele-specific fasta reference transcriptome for all possible intermediate isoforms in HLA genes


"""

import pandas as pd
import sys
import pybedtools
from pybedtools import BedTool
import itertools

# Input files
gene_names_df = pd.read_table(sys.argv[1], header=None)
gene_names_df.columns = ['gene_name','gene_id']

hg38_gene_df = pd.read_table(sys.argv[2], header=None)
hg38_gene_df.columns = ['chrom','start','end','gene','feature','strand']
hg38_gene_df['gene_id'] = hg38_gene_df['gene'].str.split("\\.").str[0]
hg38_gene_df = hg38_gene_df.merge(gene_names_df, on='gene_id')

hg38_fasta = pybedtools.example_filename(sys.argv[3])

def fasta_intermediate_isos(gene_name):
    
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
                exon_dict[feature_number] = chrom + ":" + str(start+1) + "-" + str(end)
                
            elif feature == "intron":
                intron_dict[feature_number] = chrom + ":" + str(start+1) + "-" + str(end)
                
    # Get all possible combinations of retained introns
    n_introns = len(intron_dict.keys())
    
    for level in range(1, n_introns+1):
        for isoform in itertools.combinations(range(1,n_introns+1), level):
            seq_list = []
            # Retrieve the sequence for the exons and introns (since this is by genomic position,
            # it will not be stranded and will correspond to the plus strand)
            for pos in range(1,n_introns+1):
                # add every exon in order
                seq_list.append(pybedtools.bedtool.BedTool.seq(exon_dict[pos], hg38_fasta))
                # add intron only if it's still present in the isoform
                if pos in isoform:
                    seq_list.append(pybedtools.bedtool.BedTool.seq(intron_dict[pos], hg38_fasta))
            # add last exon
            seq_list.append(pybedtools.bedtool.BedTool.seq(exon_dict[n_introns+1], hg38_fasta))
            # give isoform a name - for easier integration into downstream analyses, use transcript intron numbers
            # not genomic intron numbers
            if strand == "+":
                iso_join = ",".join([str(a) for a in isoform])
            elif strand == "-":
                iso_join = ",".join([str(n_introns-a+1) for a in isoform])
            iso_name = gene_name + ":" + iso_join
            # join all sequences from one isoform
            iso_dict[iso_name] = "".join(seq_list)
    # Add the fully spliced isoform
    seq_list = []
    for pos in range(1,n_introns+2):
        # add every exon in order
        seq_list.append(pybedtools.bedtool.BedTool.seq(exon_dict[pos], hg38_fasta))
    if strand == "+":
        iso_join = ",".join([str(a) for a in isoform])
    elif strand == "-":
        iso_join = ",".join([str(n_introns-a+1) for a in isoform])
    iso_name = gene_name + ":" + "all_spliced"
    iso_dict[iso_name] = "".join(seq_list)
                
    return(iso_dict)



# Apply function
HLA_A_dict = fasta_intermediate_isos("HLA-A")
HLA_B_dict = fasta_intermediate_isos("HLA-B")
HLA_C_dict = fasta_intermediate_isos("HLA-C")

# Write to file
out_fasta = open(sys.argv[4], 'w')

for k in HLA_A_dict.keys():
    out_fasta.write(">" + k + '\n' + HLA_A_dict[k] + '\n')

for k in HLA_B_dict.keys():
    out_fasta.write(">" + k + '\n' + HLA_B_dict[k] + '\n')

for k in HLA_C_dict.keys():
    out_fasta.write(">" + k + '\n' + HLA_C_dict[k] + '\n')

out_fasta.close()

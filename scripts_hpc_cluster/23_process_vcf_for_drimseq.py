'''

This script will reformat VCF files to have one file containing genotypes from all SNPs in genes for which splicing order could be computed.
Genotypes for heterozygous SNPs in the initial VCF file are replaced by those that were phased by LORALS

'''


import pandas as pd
import sys

# Load filtered vcf
het_vcf = pd.read_table(sys.argv[1])
hom_vcf = pd.read_table(sys.argv[2])

# Note that GM19102 and GM18870 columns are reversed in the two VCF files
hom_vcf.columns = ["CHROM","START","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","GM19102","GM19099","GM18870","GM19223","GM18501","GM18510","GM18853","GM19144","GM19152","GM18522","GM18861","GM19209"]
het_vcf.columns = ["CHROM","START","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","GM18870","GM19099","GM19102","GM19223","GM18501","GM18510","GM18853","GM19144","GM19152","GM18522","GM18861","GM19209"]

# Reformat VCF files
hom_vcf_bis = hom_vcf.melt(id_vars=['CHROM','START','ID','REF','ALT'], value_vars=["GM19102","GM19099","GM18870","GM19223","GM18501","GM18510","GM18853","GM19144","GM19152","GM18522","GM18861","GM19209"], var_name='cell_line', value_name='genotype')
hom_vcf_bis = hom_vcf_bis[hom_vcf_bis['genotype'].isin(["0|0","1|1"])].reset_index(drop=True)

het_vcf_bis = het_vcf.melt(id_vars=['CHROM','START','ID','REF','ALT'], value_vars=["GM19102","GM19099","GM18870","GM19223","GM18501","GM18510","GM18853","GM19144","GM19152","GM18522","GM18861","GM19209"], var_name='cell_line', value_name='genotype')
het_vcf_bis = het_vcf_bis[het_vcf_bis['genotype'].isin(["0|1","1|0"])].reset_index(drop=True)

# Combine the two VCF files
vcf = pd.concat([hom_vcf_bis,het_vcf_bis]).reset_index(drop=True)

# Split the genotype into two alleles
vcf['P'] = vcf['genotype'].str.split("|").str[0]
vcf['M'] = vcf['genotype'].str.split("|").str[1]

# Reformat, then make one column per allele and convert genotype 1 to 2
vcf_m = vcf.melt(id_vars=['CHROM','START','ID','REF','ALT','cell_line'], value_vars=['P','M'], var_name='allele', value_name='genotype2')
vcf_m['genotype2'] = vcf_m['genotype2'].str.replace("1","2")
vcf_m['allele_name'] = vcf_m['cell_line'] + '_' + vcf_m['allele']
vcf_m['genotype2'] = vcf_m['genotype2'].astype(int)

# Pivot table, replace NA by -1
vcf_piv = vcf_m.pivot_table(index=['CHROM','START','ID','REF','ALT'], columns='allele_name', values='genotype2', fill_value=-1).reset_index()

vcf_piv.to_csv(sys.argv[-1], index=False, sep="\t", header=True)



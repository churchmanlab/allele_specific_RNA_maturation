library(DRIMSeq)
library(GenomicRanges)
library(BiocParallel)
library(ggplot2)

set.seed(123)

outdir <- "results/drimseq_tuQTL/LCL_merged_interm_counts_per_allele_stringent_version_window500/"
n_thread <- 5
counts_file <- "../data/drimseq_tuqtl_input/LCLs_merged_interm_counts_per_allele_filtered_min_2_sig_cell_lines_3_isoforms.tsv"
genes_file <- "../data/drimseq_tuqtl_input/NCBI_RefSeq_hg38_genes_parsed_uniq_id.bed"
genotypes_file <- "../data/drimseq_tuqtl_input/LCL_genotypes.tsv"
window <- 500
min_gene_expr <- 10
min_samps_gene_expr <- 10
min_feature_expr <- 0
min_samps_feature_expr <- 0
minor_allele_freq <- 2

bpparam <- BiocParallel::MulticoreParam(workers = n_thread)


if (!base::file.exists(outdir)) {
    dir.create(outdir)
}

# Counts
counts <- read.delim(counts_file)

# Samples
samples <- data.frame(sample_id = colnames(counts)[-c(1, 2)])

# Gene GRanges. It must contain gene names when calling names()
genes_bed <- read.delim(file = genes_file)
colnames(genes_bed) <- c("chr", "start", "end", "gene_id", "score", "strand")
gene_ranges <- makeGRangesFromDataFrame(genes_bed)
mcols(gene_ranges)$name <- genes_bed$gene_id
names(gene_ranges) <- genes_bed$gene_id

# Genotypes
genotypes <- read.delim(genotypes_file)

# SNP GRanges. It must contain SNP names when calling names()
snp_ranges <- makeGRangesFromDataFrame(genotypes)
names(snp_ranges) <- genotypes$snp_id

# dmSQTLdata
d <- dmSQTLdata(
    counts = counts, gene_ranges = gene_ranges,
    genotypes = genotypes, snp_ranges = snp_ranges,
    samples = samples, window = window
)

d <- dmFilter(d,
    min_samps_gene_expr = 10, min_samps_feature_expr = 0,
    minor_allele_freq = 2, min_gene_expr = 10, min_feature_expr = 0
)


ggsave(
    paste0(outdir, "/transcripts_per_gene.png"),
    plotData(d, plot_type = "features")
)

ggsave(
    paste0(outdir, "/snps_per_gene.png"),
    plotData(d, plot_type = "snps")
)

ggsave(
    paste0(outdir, "/blocks_per_gene.png"),
    plotData(d, plot_type = "blocks")
)

d <- dmPrecision(d, BPPARAM = bpparam)
ggsave(
    paste0(outdir, "/precision_by_gene_expression.png"),
    plotPrecision(d)
)

d <- dmFit(d, BPPARAM = bpparam)

d <- dmTest(d, BPPARAM = bpparam)
ggsave(
    paste0(outdir, "/pvalues.png"),
    plotPValues(d)
)

res <- results(d)

write.table(res, paste0(outdir, "/drimseq_tuQTL.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)

filtered_res <- na.omit(res)
filtered_res <- filtered_res[filtered_res$adj_pvalue < 0.05, ]

write.table(filtered_res, paste0(outdir, "/drimseq_tuQTL_filtered.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
)

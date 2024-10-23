library(biomaRt)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
attributes <- c("ensembl_gene_id", "ensembl_transcript_id", "refseq_mrna", "refseq_ncrna", "external_gene_name", "gene_biotype")
ensembl_refseq <- getBM(attributes = attributes, mart = ensembl)

write.table(ensembl_refseq, "data/ensembl2refseq.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE)


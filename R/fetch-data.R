load("data/mdd.RNAseq.rda")
object.size(mdd.RNAseq)
object.size(mdd.RNAseq$rnaSeq)

# randomly select 500 genes
set.seed(1)
mdd.RNAseq.small <- mdd.RNAseq
n_genes <- ncol(mdd.RNAseq$rnaSeq)
mdd.RNAseq.small$rnaSeq <- mdd.RNAseq.small$rnaSeq[, sample.int(n_genes, 500)]
usethis::use_data(mdd.RNAseq.small, overwrite = TRUE)

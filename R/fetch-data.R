# load("data/mdd.RNAseq.rda")
# 
# # randomly select 500 genes
# set.seed(1)
# mdd.RNAseq.small <- mdd.RNAseq
# n_genes <- ncol(mdd.RNAseq$rnaSeq)
# mdd.RNAseq.small$rnaSeq <- mdd.RNAseq.small$rnaSeq[, sample.int(n_genes, 500)]
# 
# library(dplyr)
# set.seed(1618)
# ##### simulate case-control interaction effect data
# case.control.3sets <- privateEC::createSimulation(
#   num.samples = 300, # 100 samples in train/holdout/test
#   num.variables = 100, # 100 features
#   pct.signals = 0.1, # pct functional features
#   label = "class", # tells simulator to do case/control and adds this colname
#   bias = 0.4, # moderate effect size
#   pct.train = 1 / 3,
#   pct.holdout = 1 / 3,
#   pct.validation = 1 / 3,
#   sim.type = "interactionErdos", # or mainEffect
#   save.file = NULL,
#   verbose = FALSE
# )
# 
# qtrait.3sets <- privateEC::createSimulation(
#   num.samples = 300,
#   num.variables = 100,
#   pct.signals = 0.1, # pct functional features
#   label = "qtrait", # quantitative trait, adds this colname
#   bias = 0.6, # moderate effect size
#   pct.train = 1 / 3,
#   pct.holdout = 1 / 3,
#   pct.validation = 1 / 3,
#   sim.type = "mainEffect",
#   save.file = NULL,
#   verbose = FALSE
# )
# 
# usethis::use_data(mdd.RNAseq.small, case.control.3sets, qtrait.3sets, overwrite = TRUE)

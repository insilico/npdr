load("data-raw/mdd.RNAseq.rda")

# randomly select 500 genes
set.seed(1)
mdd.RNAseq.small <- mdd.RNAseq
n_genes <- ncol(mdd.RNAseq$rnaSeq)
mdd.RNAseq.small$rnaSeq <- mdd.RNAseq.small$rnaSeq[, sample.int(n_genes, 500)]

library(dplyr)
set.seed(1618)
##### simulate case-control interaction effect data
case.control.3sets <- privateEC::createSimulation(
  num.samples = 300, # 100 samples in train/holdout/test
  num.variables = 100, # 100 features
  pct.signals = 0.1, # pct functional features
  label = "class", # tells simulator to do case/control and adds this colname
  bias = 0.4, # moderate effect size
  pct.train = 1 / 3,
  pct.holdout = 1 / 3,
  pct.validation = 1 / 3,
  sim.type = "interactionErdos", # or mainEffect
  save.file = NULL,
  verbose = FALSE
)

qtrait.3sets <- privateEC::createSimulation(
  num.samples = 300,
  num.variables = 100,
  pct.signals = 0.1, # pct functional features
  label = "qtrait", # quantitative trait, adds this colname
  bias = 0.6, # moderate effect size
  pct.train = 1 / 3,
  pct.holdout = 1 / 3,
  pct.validation = 1 / 3,
  sim.type = "mainEffect",
  save.file = NULL,
  verbose = FALSE
)

predictors.mat <- case.control.3sets$train
predictors.mat <- predictors.mat[, names(predictors.mat) != "class"]

usethis::use_data(mdd.RNAseq.small, case.control.3sets, qtrait.3sets, predictors.mat, overwrite = TRUE)

### correlation data
# p <- 100
# attr.idx.list <- list()
# for (i in 1:p) {
#   lo.idx <- (i - 1) * (p - 1) + 1
#   hi.idx <- i * (p - 1)
#   attr.idx.list[[i]] <- c(lo.idx:hi.idx)
# }
# 
# m <- 100
# p <- 100
# 
# myfisher <- function(rho) {
#   z <- 0.5 * log((1 + rho) / (1 - rho))
#   return(z)
# }
# 
# stretch_mat <- function(M) {
#   mat <- numeric()
#   for (k in 1:nrow(M)) {
#     mat <- c(mat, M[k, -k])
#   }
#   return(mat)
# }
# 
# data.mat <- matrix(rep(0, length = m * (p * (p - 1))), nrow = (p * (p - 1)), ncol = m)
# for (k in 1:m) {
#   # print(k)
#   X <- matrix(runif(p^2, min = -1, max = 1), ncol = p)
#   cov <- X %*% t(X)
#   mycorr <- cov2cor(cov)
#   zcorr <- apply(matrix(stretch_mat(mycorr), ncol = 1), 1, myfisher)
#   data.mat[, k] <- matrix(zcorr, nrow = (p * (p - 1)), ncol = 1)
# }
# 
# var.names <- c(
#   paste("main", 1:10, sep = ""),
#   paste("var", 1:90, sep = "")
# )
# 
# class.vec <- c(rep(1, length = round(m / 2)), rep(0, length = round(m / 2)))
# case.control.data <- t(data.mat)
# # case.control.data <- cbind(case.control.data, class.vec)
# colnames(case.control.data) <- c((1:(ncol(case.control.data))))
# row.names(case.control.data) <- c(
#   paste("case", 1:floor(m / 2), sep = ""),
#   paste("ctrl", 1:floor(m / 2), sep = "")
# )
# 
# # m=100 subjects
# # p=100 attributes, p*(p-1)=9,900 correlation predictors
# npdr.cc.results <- npdr(class.vec, case.control.data,
#                         regression.type = "binomial",
#                         attr.diff.type = "correlation-data",
#                         nbd.method = "relieff", nbd.metric = "manhattan", msurf.sd.frac = .5, k = 0,
#                         neighbor.sampling = "none",
#                         dopar.nn = FALSE, padj.method = "bonferroni", verbose = TRUE,
#                         corr.attr.names = var.names
# )

library(privateEC)
library(broom)
library(dplyr)

## npdr install
# library(devtools)
# install_github("insilico/npdr")
library(npdr)

#==================== Simulate Data with pEC ===================================#
##### simulate case-control interaction effect data 
n.samples <- 300     # 100 samples in train/holdout/test
n.variables <- 100   # 100 features
label <- "class" # tells simulator to do case/control and adds this colname
#type <- "interactionErdos" # or mainEffect
type <-"mainEffect"
bias <- 0.6          # moderate effect size
pct.signals <- 0.1   # pct functional features
verbose <- FALSE
case.control.3sets <- createSimulation(num.samples = n.samples,
                                       num.variables = n.variables,
                                       pct.signals = pct.signals,
                                       label = label,
                                       bias = bias,
                                       pct.train = 1/3,
                                       pct.holdout = 1/3,
                                       pct.validation = 1/3,
                                       sim.type = type,
                                       save.file = NULL,
                                       verbose = verbose)
# combine train and holdout into 200 samples x 100 attributes
# ignore validation set
case.control.data <- rbind(case.control.3sets$train,case.control.3sets$holdout)
#validation.data <- data.sets$validation
n.samples.case.control <- dim(case.control.data)[1]
pheno.case.control <- as.factor(case.control.data[,"class"])
functional.case.control <- case.control.3sets$signal.names # functional attributes

##### Univariate logistic regression
# linear regression on all predictors, fdr adjust, check functional hits
# standardized beta and p-value
# npdr utulity function
univariate.cc.results <- uniReg(outcome="class", dataset=case.control.data, regression.type="binomial")
univariate.cc.results[1:10,]
# don't expect any less than .05 for interaction simulations
univariate.cc.results[univariate.cc.results[,"p.adj"]<.05,]

##### Run npdr
npdr.cc.results <- npdr("class", case.control.data, regression.type="binomial", attr.diff.type="numeric-abs",
                        nbd.method="multisurf", nbd.metric = "manhattan", msurf.sd.frac=.5, k=0,
                        neighbor.sampling="none", dopar.nn = F,
                        padj.method="bonferroni", verbose=T)
npdr.cc.results[npdr.cc.results$pval.adj<.05,] # pval.adj, first column

npdr.cc.results <- npdr("class", case.control.data, regression.type="binomial", attr.diff.type="numeric-abs",
                        nbd.method="relieff", nbd.metric = "manhattan", msurf.sd.frac=.5, k=0,
                        neighbor.sampling="none", separate.hitmiss.nbds=T,
                        dopar.nn = F, padj.method="bonferroni", verbose=T)
# attributes with npdr adjusted p-value less than .05 
npdr.cc.results[npdr.cc.results$pval.adj<.05,] # pval.adj, first column
selected.features <- npdr.cc.results %>% filter(pval.adj<.05) %>% pull(att)

## Neighbors
cc.attrs <- case.control.data[,colnames(case.control.data)!="class"]
cc.pheno <- as.numeric(as.character(case.control.data[,colnames(case.control.data)=="class"]))
cc.nbpairs.idx <- nearestNeighbors(cc.attrs, 
                                          nb.method="multisurf", nb.metric="manhattan", 
                                          sd.frac = .5, k=0)
cc.nbpairs.idx <- nearestNeighborsSeparateHitMiss(cc.attrs, cc.pheno, 
                                   nb.method="multisurf", nb.metric="manhattan", 
                                   sd.frac = .5, k=0)
Ri.phenos <- cc.pheno[cc.nbpairs.idx$Ri_idx]
NN.phenos <- cc.pheno[cc.nbpairs.idx$NN_idx]
hitmiss.vec <- npdrDiff(Ri.phenos, NN.phenos, diff.type="match-mismatch")
pct.misses <- sum(hitmiss.vec)/length(hitmiss.vec)
pct.hits <- 1-sum(hitmiss.vec)/length(hitmiss.vec)
pct.misses
pct.hits

#================= Test Nearest Neighbor Classification ================================#

# 100 train, 100 test
test_results <- npdrLearner(train.outcome="class", train.data=case.control.3sets$train, 
                            test.outcome="class", test.data=case.control.3sets$holdout,
                            nbd.method = "relieff",  # relieff, surf, or multisurf
                            nbd.metric = "manhattan", 
                            msurf.sd.frac = 0.5, dopar.nn=F,
                            knn=0) # k=0 uses multisurf k estimate
test_acc <- test_results$accuracy
cat("\n Train Accuracy [",test_acc,"]\n")

# 200 train, 100 test
test_results <- npdrLearner(train.outcome="class", train.data=case.control.data, 
                            test.outcome="class", test.data=case.control.3sets$validation,
                            nbd.method = "relieff", 
                            nbd.metric = "manhattan", 
                            msurf.sd.frac = 0.5, 
                            knn=0) # k=0 uses multisurf k estimate
test_acc <- test_results$accuracy
cat("\n Train Accuracy [",test_acc,"]\n")

#=======================consenses nested CV from pEC library ===========================#
# 100 train, 100 test plus a validation set
cncv.case.control <- consensus_nestedCV(train.ds = case.control.data, 
                                        validation.ds =  case.control.3sets$validation, 
                                        label = "class",
                                        method.model = "classification",
                                        is.simulated = TRUE,
                                        ncv_folds = c(10, 10),
                                        param.tune = FALSE,
                                        learning_method = "rf", 
                                        importance.algorithm = "ReliefFequalK",
                                        relief.k.method = "k_half_sigma",     # surf k
                                        num_tree = 500,
                                        verbose = F)

cat("\n Train Accuracy [",cncv.case.control$cv.acc,"]\n")
cat("\n Validation Accuracy [",cncv.case.control$Validation,"]\n")
cat("\n Selected Features \n [",cncv.case.control$Features,"]\n")
#cat("\n Elapsed Time [",cncv.case.control$Elapsed,"]\n")
#cat(detectionStats(functional.case.control, cncv.case.control$Features)$report)

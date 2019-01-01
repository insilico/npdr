library(privateEC)
library(broom)

## glmSTIR install
library(devtools)
install_github("insilico/glmSTIR")
library(glmSTIR)

##### simulate case-control interaction effect data 
n.samples <- 300     # 100 samples in train/holdout/test
n.variables <- 100   # 100 features
label <- "class" # tells simulator to do case/control and adds this colname
type <- "interactionErdos" # or mainEffect
bias <- 0.4          # moderate effect size
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
# glmSTIR utulity function
univariate.results <- univariateRegression(outcome="class", dataset=case.control.data, regression.type="glm")
univariate.results[1:10,]
# none will be less than .05 for interaction simulations
#univariate.results[univariate.results[,"p.adj"]<.05,]

##### Run glmSTIR
glm.stir.cc.results <- glmSTIR("class", case.control.data, regression.type="glm", attr.diff.type="numeric-abs",
                            nbd.method="multisurf", nbd.metric = "manhattan", sd.frac=.5, 
                            fdr.method="bonferroni", verbose=T)
# attributes with glmSTIR adjusted p-value less than .05 
glm.stir.cc.results[glm.stir.cc.results$pval.adj<.05,] # pval.adj, first column
# attributes with glmSTIR raw/nominal p-value less than .05
#rownames(glm.stir.cc.results)[glm.stir.cc.results$pval.attr<.05] # pval.attr, second column

# functional attribute detection stats
glm.stir.cc.positives <- row.names(glm.stir.cc.results[glm.stir.cc.results[,1]<.05,]) # p.adj<.05
glm.stir.cc.detect.stats <- detectionStats(functional.case.control, glm.stir.cc.positives)
cat(glm.stir.cc.detect.stats$report)

##### original (pseudo t-test) STIR
# impression is that glmSTIR gives same resutls as original t-STIR
#install_github("insilico/stir")
library(stir)
# stir interface requires splitting phenotype and predictor matrix, 
# and requires finding the neighborhood separately.
predictors.cc.mat <- case.control.data[, - which(colnames(case.control.data) == "class")]
case.control.data[, "class"] <- as.factor(case.control.data[, "class"]) 
pheno.case.control <- case.control.data[, "class"]
neighbor.idx.observed <- find.neighbors(predictors.cc.mat, pheno.case.control, k = 0, method = "multisurf")
results.list <- stir(predictors.cc.mat, neighbor.idx.observed, k = 0, metric = "manhattan", method = "multisurf")
t_sorted_multisurf <- results.list$STIR_T[, -3]  # remove cohen-d
colnames(t_sorted_multisurf) <- paste(c("t.stat", "t.pval", "t.pval.adj"), "stir", sep=".")
(t_sorted_multisurf[t_sorted_multisurf[,3]<.05,])

# functional attribute detection stats
tstat_stir.detect.stats <- detectionStats(functional.case.control, 
                                          row.names(t_sorted_multisurf[t_sorted_multisurf[,3]<.05,]))
cat(tstat_stir.detect.stats$report)

##### CORElearn ReliefF with surf fixed k
# fixed k with theoretical surf value
library(CORElearn)
core.learn.case.control <- CORElearn::attrEval("class", data = case.control.data,
                                      estimator = "ReliefFequalK",
                                      costMatrix = NULL,
                                      outputNumericSplits=FALSE,
                                      kNearestEqual = knnSURF(n.samples.case.control,.5))
core.learn.case.control.order <- order(core.learn.case.control, decreasing = T)
t(t(core.learn.case.control[core.learn.case.control.order[1:20]]))

# functional attribute detection stats
core.learn.detect.stats <- detectionStats(functional.case.control, 
                                          names(core.learn.case.control)[core.learn.case.control.order[1:20]])
cat(core.learn.detect.stats$report)


##### Consensus Nested Cross Validation with ReliefF with surf fixed k
# selects features and learns classification model.

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
cat("\n Elapsed Time [",cncv.case.control$Elapsed,"]\n")
cat(detectionStats(functional.case.control, cncv.case.control$Features)$report)

##### Regular Nested Cross Validation with ReliefF with surf fixed k
# selects features and learns classification model.

rncv.case.control <- regular_nestedCV(train.ds = case.control.data, 
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

cat("\n Train Accuracy [",rncv.case.control$cv.acc,"]\n")
cat("\n Validation Accuracy [",rncv.case.control$Validation,"]\n")
cat("\n Selected Features \n [",rncv.case.control$Features,"]\n")
cat("\n Elapsed Time [",rncv.case.control$Elapsed,"]\n")
cat(detectionStats(functional.case.control, rncv.case.control$Features)$report)

##### GLMnet (penalized regression) comparison. Don't expect good performance for interaction models. 

library(glmnet)
# cc short for case-control
predictors.cc.mat <- case.control.data[, - which(colnames(case.control.data) == "class")]
pheno.case.control <- case.control.data[, "class"]

glmnet.cc.model<-cv.glmnet(as.matrix(predictors.cc.mat), pheno.case.control, alpha=.1, family="binomial", type.measure="class")
glmnet.cc.coeffs<-predict(glmnet.cc.model,type="coefficients")
#glmnet.cc.coeffs  # maybe 3 is most important, Excess kurtosis
model.cc.terms <- colnames(predictors.cc.mat)  # glmnet includes an intercept but we are going to ignore
nonzero.glmnet.cc.coeffs <- model.terms[glmnet.cc.coeffs@i[which(glmnet.cc.coeffs@i!=0)]] # skip intercept if there, 0-based counting
nonzero.glmnet.cc.coeffs

library(privateEC)
library(broom)

## glmSTIR install
library(devtools)
install_github("insilico/glmSTIR")
library(glmSTIR)

##### simulate case-control interaction effect data 
n.samples <- 300     # 100 samples in train/holdout/test
n.variables <- 100   # 100 features
label <- "class"
#method.model <- "regression"
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
lr.results <- univariateLogisticRegression(outcome="class", dataset=case.control.data)
lr.results[1:10,]
# none will be less than .05 for interaction simulations
#lr.results[lr.results[,"p.adj"]<.05,]

##### Run glmSTIR
glm.stir.results.df <- glmSTIR("class", case.control.data, regression.type="glm", 
                            nbd.method="multisurf", nbd.metric = "manhattan", 
                            attr.diff.type="numeric-abs", covar.diff.type="numeric-abs", 
                            sd.frac=.5, fdr.method="bonferroni")
glm.stir.results.df[glm.stir.results.df[,1]<.05,]

# functional attribute detection stats
glm.stir.positives <- row.names(glm.stir.results.df[glm.stir.results.df[,1]<.05,]) # p.adj<.05
glm.stir.detect.stats <- detectionStats(functional.case.control, glm.stir.positives)
cat(glm.stir.detect.stats$report)

##### original STIR (pseudo t-test)
# t-test STIR
#install_github("insilico/stir")
library(stir)
# stir interface requires splitting phenotype and predictor matrix, 
# and requires finding the neighborhood separately.
predictors.mat <- case.control.data[, - which(colnames(case.control.data) == "class")]
case.control.data[, "class"] <- as.factor(case.control.data[, "class"]) 
pheno.case.control <- case.control.data[, "class"]
neighbor.idx.observed <- find.neighbors(predictors.mat, pheno.case.control, k = 0, method = "multisurf")
results.list <- stir(predictors.mat, neighbor.idx.observed, k = k, metric = "manhattan", method = "multisurf")
t_sorted_multisurf <- results.list$STIR_T[, -3]  # remove cohen-d
colnames(t_sorted_multisurf) <- paste(c("t.stat", "t.pval", "t.pval.adj"), "stir", sep=".")
(t_sorted_multisurf[t_sorted_multisurf[,3]<.05,])

# functional attribute detection stats
tstat_stir.detect.stats <- detectionStats(functional.case.control, 
                                          row.names(t_sorted_multisurf[t_sorted_multisurf[,3]<.05,]))
cat(tstat_stir.detect.stats$report)

##### CORElearn ReliefF with surf fixed k

# fixed k with theoretical surf value
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
k.surf.fn <- function(m,f) {floor((m-1)*(1-erf(f/sqrt(2)))/2)}
library(CORElearn)
core.learn.case.control <- CORElearn::attrEval("class", data = case.control.data,
                                      estimator = "ReliefFequalK",
                                      costMatrix = NULL,
                                      outputNumericSplits=FALSE,
                                      kNearestEqual = k.surf.fn(n.samples.case.control,.5))
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

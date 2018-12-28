library(privateEC)
library(broom)

## glmSTIR install
library(devtools)
install_github("insilico/glmSTIR")
library(glmSTIR)

##### simulate case-control interaction effect data 
n.samples <- 300     # 100 samples in train/holdout/test
n.variables <- 100   # 100 features
label <- "qtrait"   # tells simulator to do quantitative trait and adds this colname
type <- "mainEffect"
bias <- 0.6          # moderate effect size
pct.signals <- 0.1   # pct functional features
verbose <- FALSE
qtrait.3sets <- createSimulation(num.samples = n.samples,
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
qtrait.data <- rbind(qtrait.3sets$train,qtrait.3sets$holdout)
#validation.data <- data.sets$validation
n.samples.qtrait <- dim(qtrait.data)[1]
pheno.qtrait <- qtrait.data[,"qtrait"]
functional.qtrait <- qtrait.3sets$signal.names # functional attributes

##### Univariate logistic regression
# linear regression on all predictors, fdr adjust, check functional hits
# standardized beta and p-value
# glmSTIR utulity function
lr.results <- univariateLogisticRegression(outcome="qtrait", dataset=qtrait.data)
lr.results[1:10,]
# none will be less than .05 for interaction simulations
#lr.results[lr.results[,"p.adj"]<.05,]

##### Run glmSTIR
glm.stir.results.df <- glmSTIR("class", qtrait.data, regression.type="glm", 
                            nbd.method="multisurf", nbd.metric = "manhattan", 
                            attr.diff.type="numeric-abs", covar.diff.type="numeric-abs", 
                            sd.frac=.5, fdr.method="bonferroni")
glm.stir.results.df[glm.stir.results.df[,1]<.05,]

# functional attribute detection stats
glm.stir.positives <- row.names(glm.stir.results.df[glm.stir.results.df[,1]<.05,]) # p.adj<.05
glm.stir.detect.stats <- detectionStats(functional.qtrait, glm.stir.positives)
cat(glm.stir.detect.stats$report)

##### original (pseudo t-test) STIR
#install_github("insilico/stir")
library(stir)
# stir interface requires splitting phenotype and predictor matrix, 
# and requires finding the neighborhood separately.
predictors.mat <- qtrait.data[, - which(colnames(qtrait.data) == "class")]
qtrait.data[, "class"] <- as.factor(qtrait.data[, "class"]) 
pheno.qtrait <- qtrait.data[, "class"]
neighbor.idx.observed <- find.neighbors(predictors.mat, pheno.qtrait, k = 0, method = "multisurf")
results.list <- stir(predictors.mat, neighbor.idx.observed, k = k, metric = "manhattan", method = "multisurf")
t_sorted_multisurf <- results.list$STIR_T[, -3]  # remove cohen-d
colnames(t_sorted_multisurf) <- paste(c("t.stat", "t.pval", "t.pval.adj"), "stir", sep=".")
(t_sorted_multisurf[t_sorted_multisurf[,3]<.05,])

# functional attribute detection stats
tstat_stir.detect.stats <- detectionStats(functional.qtrait, 
                                          row.names(t_sorted_multisurf[t_sorted_multisurf[,3]<.05,]))
cat(tstat_stir.detect.stats$report)

##### CORElearn ReliefF with surf fixed k

# fixed k with theoretical surf value
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
k.surf.fn <- function(m,f) {floor((m-1)*(1-erf(f/sqrt(2)))/2)}
library(CORElearn)
core.learn.qtrait <- CORElearn::attrEval("class", data = qtrait.data,
                                      estimator = "ReliefFequalK",
                                      costMatrix = NULL,
                                      outputNumericSplits=FALSE,
                                      kNearestEqual = k.surf.fn(n.samples.qtrait,.5))
core.learn.qtrait.order <- order(core.learn.qtrait, decreasing = T)
t(t(core.learn.qtrait[core.learn.qtrait.order[1:20]]))

# functional attribute detection stats
core.learn.detect.stats <- detectionStats(functional.qtrait, 
                                          names(core.learn.qtrait)[core.learn.qtrait.order[1:20]])
cat(core.learn.detect.stats$report)


##### Consensus Nested Cross Validation with ReliefF with surf fixed k
# selects features and learns classification model.

cncv.qtrait <- consensus_nestedCV(train.ds = qtrait.data, 
                                  validation.ds =  qtrait.3sets$validation, 
                                  label = "class",
                                  method.model = "classification",
                                  is.simulated = TRUE,
                                  ncv_folds = c(10, 10),
                                  param.tune = FALSE,
                                  learning_method = "rf", 
                                  importance.algorithm = "RReliefFequalK",
                                  relief.k.method = "k_half_sigma",     # surf k
                                  num_tree = 500,
                                  verbose = F)

cat("\n Train Accuracy [",cncv.qtrait$cv.acc,"]\n")
cat("\n Validation Accuracy [",cncv.qtrait$Validation,"]\n")
cat("\n Selected Features \n [",cncv.qtrait$Features,"]\n")
cat("\n Elapsed Time [",cncv.qtrait$Elapsed,"]\n")
cat(detectionStats(functional.qtrait, cncv.qtrait$Features)$report)

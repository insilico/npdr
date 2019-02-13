library(privateEC)
library(broom)

## glmSTIR install
library(devtools)
install_github("insilico/npdr")
library(npdr)

##### simulate case-control interaction effect data 
n.samples <- 300     # 100 samples in train/holdout/test
n.variables <- 100   # 100 features
label <- "class" # tells simulator to do case/control and adds this colname
type <- "interactionErdos" # or mainEffect
#type <-"mainEffect"
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
univariate.cc.results <- univariateRegression(outcome="class", dataset=case.control.data, regression.type="glm")
#univariate.cc.results[1:10,]
# don't expect any less than .05 for interaction simulations
univariate.cc.results[univariate.cc.results[,"p.adj"]<.05,]

##### Run glmSTIR
glm.stir.cc.results <- glmSTIR("class", case.control.data, regression.type="glm", attr.diff.type="numeric-abs",
                            nbd.method="multisurf", nbd.metric = "manhattan", msurf.sd.frac=.5, 
                            fdr.method="bonferroni", verbose=T)
# attributes with glmSTIR adjusted p-value less than .05 
glm.stir.cc.results[glm.stir.cc.results$pval.adj<.05,] # pval.adj, first column
# attributes with glmSTIR raw/nominal p-value less than .05
#rownames(glm.stir.cc.results)[glm.stir.cc.results$pval.attr<.05] # pval.attr, second column

# functional attribute detection stats
glm.stir.cc.positives <- row.names(glm.stir.cc.results[glm.stir.cc.results[,1]<.05,]) # p.adj<.05
glm.stir.cc.detect.stats <- detectionStats(functional.case.control, glm.stir.cc.positives)
cat(glm.stir.cc.detect.stats$report)

### Compare univariate and glmSTIR
univ.log10.df <- data.frame(vars=rownames(univariate.cc.results),univ.log10=-log10(univariate.cc.results[,"pval"]))
glmstir.cc.log10.df <- data.frame(vars=rownames(glm.stir.cc.results),glmstir.cc.log10=-log10(glm.stir.cc.results$pval.attr))

univ.cc.pcutoff <- max(-log10(univariate.cc.results[,"pval"]))
#univ.pcutoff <- -log10(t_sorted_multisurf$t.pval.stir[which(t_sorted_multisurf$t.pval.adj.stir>.05)[1]-1])
glmstir.cc.pcutoff <- -log10(glm.stir.cc.results$pval.attr[which(glm.stir.cc.results$pval.adj>.05)[1]-1])

library(ggplot2)
test.cc.df <- merge(univ.log10.df,glmstir.cc.log10.df)
functional <- factor(c(rep("Func",length(functional.case.control)),rep("Non-Func",n.variables-length(functional.case.control))))
ggplot(test.cc.df, aes(x=univ.log10,y=glmstir.cc.log10)) + geom_point(aes(colour = functional), size=4) +
  theme(text = element_text(size = 20)) +
  #geom_vline(xintercept=univ.cc.pcutoff, linetype="dashed") +
  geom_hline(yintercept=glmstir.cc.pcutoff, linetype="dashed") +
  xlab("Logistic Regression -log10(P)") + ylab("NPDR -log10(P)") 

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

### Compare STIR and glmSTIR
stir.log10.df <- data.frame(vars=rownames(t_sorted_multisurf),stir.log10=-log10(t_sorted_multisurf$t.pval.stir))
glmstir.log10.df <- data.frame(vars=rownames(glm.stir.cc.results),glmstir.log10=-log10(glm.stir.cc.results$pval.attr))

stir.pcutoff <- -log10(t_sorted_multisurf$t.pval.stir[which(t_sorted_multisurf$t.pval.adj.stir>.05)[1]-1])
glmstir.pcutoff <- -log10(glm.stir.cc.results$pval.attr[which(glm.stir.cc.results$pval.adj>.05)[1]-1])

library(ggplot2)
test.df <- merge(stir.log10.df,glmstir.log10.df)
functional <- factor(c(rep("Func",length(functional.case.control)),rep("Non-Func",n.variables-length(functional.case.control))))
ggplot(test.df, aes(x=stir.log10,y=glmstir.log10)) + geom_point(aes(colour = functional), size=4) +
  theme(text = element_text(size = 20)) +
  geom_vline(xintercept=stir.pcutoff, linetype="dashed") +
  geom_hline(yintercept=glmstir.pcutoff, linetype="dashed") +
  xlab("STIR -log10(P)") + ylab("NPDR -log10(P)") 

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
#core.learn.detect.stats <- detectionStats(functional.case.control, 
#                                          names(core.learn.case.control)[core.learn.case.control.order[1:20]])

arbitrary.cc.threshold <- .0072
core.learn.cc.detect <- detectionStats(functional.case.control, 
                                           names(core.learn.case.control)[core.learn.case.control>arbitrary.cc.threshold])
cat(core.learn.cc.detect$report)

### Compare corelearn and glmSTIR
corelearn.cc.df <- data.frame(vars=names(core.learn.case.control),rrelief=core.learn.case.control)
glmstir.cc.beta.df <- data.frame(vars=rownames(glm.stir.cc.results),glmstir.beta=(glm.stir.cc.results$beta.attr))

corelearn.cc.cutoff <- arbitrary.cc.threshold
glmstir.cc.pcutoff <- (glm.stir.cc.results$beta.attr[which(glm.stir.cc.results$pval.adj>.05)[1]-1])


library(ggplot2)
test.cc.df <- merge(corelearn.cc.df,glmstir.cc.beta.df)
functional <- factor(c(rep("Func",length(functional.case.control)),rep("Non-Func",n.variables-length(functional.case.control))))
ggplot(test.cc.df, aes(x=rrelief,y=glmstir.beta)) + geom_point(aes(colour = functional), size=4) +
  theme(text = element_text(size = 20)) +
  geom_vline(xintercept=corelearn.cc.cutoff, linetype="dashed") +
  geom_hline(yintercept=glmstir.cc.pcutoff, linetype="dashed") +
  xlab("Relief Score") + ylab("NPDR Coefficient") 

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

glmnet.cc.model<-cv.glmnet(as.matrix(predictors.cc.mat), pheno.case.control, alpha=1, family="binomial", type.measure="class")
#glmnet.cc.coeffs<-as.matrix(predict(glmnet.cc.model,type="coefficients"))
#glmnet.cc.coeffs  # maybe 3 is most important, Excess kurtosis
#model.cc.terms <- colnames(predictors.cc.mat)  # glmnet includes an intercept but we are going to ignore
# skip intercept if there, 0-based counting
#nonzero.glmnet.cc.coeffs <- model.cc.terms[glmnet.cc.coeffs@i[which(glmnet.cc.coeffs@i!=0)]] 
# finds some interactions maybe because of the co-expression in the simulation
#nonzero.glmnet.cc.coeffs
#nonzero.glmnet.cc.mask <- abs(glmnet.cc.coeffs[,1])>.05  
#as.matrix(glmnet.cc.coeffs[nonzero.glmnet.cc.mask],ncol=1)

# do not expect correct associations for interaction simulations. 
glmnet.cc.coeffs<-as.matrix(predict(glmnet.cc.model,type="coefficients"))
row.names(glmnet.cc.coeffs) <- c("intercept", colnames(predictors.cc.mat))  # add variable names to results
glmnet.cc.sorted<-as.matrix(glmnet.cc.coeffs[order(abs(glmnet.cc.coeffs),decreasing = T),],ncol=1) # sort
glmnet.cc.sorted[abs(glmnet.cc.sorted)>0,]

#####
## EXPERIMENTAL Needs penalty for ordinal regression (not part of glmnet)
##### Run glmnetSTIR, penalized glmSTIR
glmnetSTIR.cc.results <- glmSTIR("class", case.control.data, regression.type="glmnet", attr.diff.type="numeric-abs",
                               nbd.method="multisurf", nbd.metric = "manhattan", msurf.sd.frac=.5,
                               glmnet.alpha=1, glmnet.family="binomial",
                               fdr.method="bonferroni", verbose=T)
# attributes with glmSTIR adjusted p-value less than .05 
glmnetSTIR.cc.results.mat <- as.matrix(glmnetSTIR.cc.results)
# .05 regression coefficient threshold is arbitrary
# not sure why glment did not force zeros
# Negative coefficients mean irrelevant attributes for Relief scores.
# However, glmnet does not include ordinal models. 
nonzero.glmnetSTIR.mask <- abs(glmnetSTIR.cc.results.mat[,1])>0.05  
as.matrix(glmnetSTIR.cc.results.mat[nonzero.glmnetSTIR.mask,],ncol=1)

# Naively remove negative coefficients, but would be better to modify shrinkage model.
pos.glmnetSTIR.mask <- glmnetSTIR.cc.results.mat[,1]>0.05  
as.matrix(glmnetSTIR.cc.results.mat[pos.glmnetSTIR.mask,],ncol=1)

# functional attribute detection stats
glmnetSTIR.cc.positives <- names(glmnetSTIR.cc.results.mat[nonzero.glmnetSTIR.mask,]) # p.adj<.05
glmnetSTIR.cc.detect.stats <- detectionStats(functional.case.control, glmnetSTIR.cc.positives)
cat(glmnetSTIR.cc.detect.stats$report)

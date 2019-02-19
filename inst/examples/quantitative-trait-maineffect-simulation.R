library(privateEC)
library(broom)
library(dplyr)

## npdr install
library(devtools)
install_github("insilico/npdr")
library(npdr)

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
# npdr utulity function
univariate.results <- uniReg(outcome="qtrait", dataset=qtrait.data, regression.type="lm")
univariate.results[1:10,]
# fdr significant.
univariate.05.fdr <- univariate.results[univariate.results[,"p.adj"]<.05,]
univariate.05.fdr
cat(detectionStats(functional.qtrait, rownames(univariate.05.fdr))$report)

##### Run npdr
npdr.qtrait.results <- npdr("qtrait", qtrait.data, regression.type="lm", attr.diff.type="numeric-abs",  
                            nbd.method="multisurf", nbd.metric = "manhattan", msurf.sd.frac=.5,
                            padj.method="bonferroni")
# attributes with npdr adjusted p-value less than .05 
npdr.qtrait.results[npdr.qtrait.results$pval.adj<.05,] # pval.adj, first column
# attributes with npdr raw/nominal p-value less than .05
#rownames(npdr.qtrait.results)[npdr.qtrait.results$pval.attr<.05] # pval.attr, second column

# functional attribute detection stats
npdr.qtrait.positives <- npdr.qtrait.results %>% filter(pval.adj<.05) %>% pull(att)
  #row.names(npdr.qtrait.results[npdr.qtrait.results$pval.adj<.05,]) # p.adj<.05
npdr.qtrait.detect.stats <- detectionStats(functional.qtrait, npdr.qtrait.positives)
cat(npdr.qtrait.detect.stats$report)

##### original (pseudo t-test) STIR
# original STIR not relevant for regression

##### CORElearn ReliefF with surf fixed k
# impression is that npdr ranks the attributes better than RReleifF.
# fixed k with theoretical surf value
library(CORElearn)
core.learn.qtrait <- CORElearn::attrEval("qtrait", data = qtrait.data,
                                      estimator = "RReliefFequalK",
                                      costMatrix = NULL,
                                      outputNumericSplits=FALSE,
                                      kNearestEqual = knnSURF(n.samples.qtrait,.5))
core.learn.qtrait.order <- order(core.learn.qtrait, decreasing = T)
t(t(core.learn.qtrait[core.learn.qtrait.order[1:20]]))

arbitrary_threshold = .005
t(t(core.learn.qtrait[core.learn.qtrait>arbitrary_threshold]))

# functional attribute detection stats
#core.learn.qtrait.detect <- detectionStats(functional.qtrait, 
#                                          names(core.learn.qtrait)[core.learn.qtrait.order[1:20]])
core.learn.qtrait.detect <- detectionStats(functional.qtrait, 
                                           names(core.learn.qtrait)[core.learn.qtrait>arbitrary_threshold])
cat(core.learn.qtrait.detect$report)

### Compare corelearn and npdr
corelearn.df <- data.frame(vars=names(core.learn.qtrait),rrelief=core.learn.qtrait)
npdr.beta.df <- data.frame(vars=npdr.qtrait.results$att,npdr.beta=(npdr.qtrait.results$beta.raw.att))

corelearn.cutoff <- arbitrary_threshold
npdr.pcutoff <- (npdr.qtrait.results$beta.attr[which(npdr.qtrait.results$pval.adj>.05)[1]-1])

library(ggplot2)
test.df <- merge(corelearn.df,npdr.beta.df)
functional <- factor(c(rep("Func",length(functional.qtrait)),rep("Non-Func",n.variables-length(functional.qtrait))))
ggplot(test.df, aes(x=rrelief,y=npdr.beta)) + geom_point(aes(colour = functional), size=4) +
  theme(text = element_text(size = 20)) +
  geom_vline(xintercept=corelearn.cutoff, linetype="dashed") +
  geom_hline(yintercept=npdr.pcutoff, linetype="dashed") +
  xlab("RRelief Scores") + ylab("NPDR Coefficients") 

##### Consensus Nested Cross Validation with ReliefF with surf fixed k
# selects features and learns regression model.

cncv.qtrait <- consensus_nestedCV(train.ds = qtrait.data, 
                                  validation.ds =  NULL, 
                                  label = "qtrait",
                                  method.model = "regression",
                                  is.simulated = TRUE,
                                  ncv_folds = c(10, 10),
                                  param.tune = FALSE,
                                  learning_method = "rf", 
                                  importance.algorithm = "RReliefFequalK",
                                  relief.k.method = "k_half_sigma",     # surf k
                                  num_tree = 500,
                                  verbose = F)

cat("\n Train R^2 [",cncv.qtrait$cv.acc,"]\n")
cat("\n Validation R^2 [",cncv.qtrait$Validation,"]\n")
cat("\n Selected Features \n [",cncv.qtrait$Features,"]\n")
cat("\n Elapsed Time [",cncv.qtrait$Elapsed,"]\n")
cat(detectionStats(functional.qtrait, cncv.qtrait$Features)$report)


##### Regular Nested Cross Validation with ReliefF with surf fixed k
# selects features and learns regression model.

rncv.qtrait <- regular_nestedCV(train.ds = qtrait.data, 
                                  validation.ds =  qtrait.3sets$validation, 
                                  label = "qtrait",
                                  method.model = "regression",
                                  is.simulated = TRUE,
                                  ncv_folds = c(5, 5),
                                  param.tune = FALSE,
                                  learning_method = "rf", 
                                  importance.algorithm = "RReliefFequalK",
                                  relief.k.method = "k_half_sigma",     # surf k
                                  num_tree = 500,
                                  verbose = F)

cat("\n Train R^2 [",rncv.qtrait$cv.acc,"]\n")
cat("\n Validation R^2 [",rncv.qtrait$Validation,"]\n")
cat("\n Selected Features \n [",rncv.qtrait$Features,"]\n")
cat("\n Elapsed Time [",rncv.qtrait$Elapsed,"]\n")
cat(detectionStats(functional.qtrait, rncv.qtrait$Features)$report)

##### GLMnet (penalized regression) comparison. 
# Impression for main effects is that TP is similar npdr, but npdr has higher FP

library(glmnet)
predictors.qtrait.mat <- qtrait.data[, - which(colnames(qtrait.data) == "qtrait")]
pheno.qtrait <- qtrait.data[, "qtrait"]

glmnet.qtrait.model<-cv.glmnet(as.matrix(predictors.qtrait.mat), pheno.qtrait, alpha=.1, type.measure="mse")
glmnet.qtrait.coeffs<-predict(glmnet.qtrait.model,type="coefficients")
#glmnet.cc.coeffs  # maybe 3 is most important, Excess kurtosis
model.qtrait.terms <- colnames(predictors.qtrait.mat)  # glmnet includes an intercept but we are going to ignore
nonzero.glmnet.qtrait.coeffs <- model.qtrait.terms[glmnet.qtrait.coeffs@i[which(glmnet.qtrait.coeffs@i!=0)]] # skip intercept if there, 0-based counting
nonzero.glmnet.qtrait.coeffs
cat(detectionStats(functional.qtrait, nonzero.glmnet.qtrait.coeffs)$report)

##### Run npdrNET, penalized npdr
npdrNET.qtrait.results <- npdr("qtrait", qtrait.data, regression.type="glmnet", attr.diff.type="numeric-abs",
                                 nbd.method="multisurf", nbd.metric = "manhattan", msurf.sd.frac=.5,
                                 glmnet.alpha=1, glmnet.family="gaussian",
                                 padj.method="bonferroni", verbose=T)
# attributes with npdr adjusted p-value less than .05 
npdrNET.qtrait.results.mat <- as.matrix(npdrNET.qtrait.results)
# .05 regression coefficient threshold is arbitrary
# not sure why glment did not force zeros
# Finds more interactions than regular glmnet, but not nearly as good as regular npdr
nonzero.npdrNET.qtrait.mask <- abs(npdrNET.qtrait.results.mat[,1])>0  
as.matrix(npdrNET.qtrait.results.mat[nonzero.npdrNET.qtrait.mask,],ncol=1)

# functional attribute detection stats
npdrNET.cc.positives <- names(npdrNET.cc.results.mat[nonzero.npdrNET.mask,]) # p.adj<.05
npdrNET.cc.detect.stats <- detectionStats(functional.case.control, npdrNET.cc.positives)
cat(npdrNET.cc.detect.stats$report)

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
# npdr utulity function
univariate.cc.results <- uniReg(outcome="class", dataset=case.control.data, regression.type="binomial")
#univariate.cc.results[1:10,]
# don't expect any less than .05 for interaction simulations
univariate.cc.results[univariate.cc.results[,"p.adj"]<.05,]

##### Run npdr

npdr.cc.results <- npdr("class", case.control.data, regression.type="binomial", attr.diff.type="numeric-abs",
                            nbd.method="multisurf", nbd.metric = "manhattan", msurf.sd.frac=.5, 
                            padj.method="bonferroni", verbose=T)
# attributes with npdr adjusted p-value less than .05 
npdr.cc.results[npdr.cc.results$pval.adj<.05,] # pval.adj, first column
# attributes with npdr raw/nominal p-value less than .05
#rownames(npdr.cc.results)[npdr.cc.results$pval.attr<.05] # pval.attr, second column

# functional attribute detection stats
#npdr.cc.positives <- row.names(npdr.cc.results[npdr.cc.results$pval.adj<.05,]) # p.adj<.05
npdr.cc.positives <- npdr.cc.results %>% filter(pval.adj<.05) %>% pull(att)
npdr.cc.detect.stats <- detectionStats(functional.case.control, npdr.cc.positives)
cat(npdr.cc.detect.stats$report)

### Compare univariate and npdr
univ.log10.df <- data.frame(vars=rownames(univariate.cc.results),univ.log10=-log10(univariate.cc.results[,"pval"]))
npdr.cc.log10.df <- data.frame(vars=npdr.cc.results$att,npdr.cc.log10=-log10(npdr.cc.results$pval.att))

univ.cc.pcutoff <- max(-log10(univariate.cc.results[,"pval"]))
#univ.pcutoff <- -log10(t_sorted_multisurf$t.pval.stir[which(t_sorted_multisurf$t.pval.adj.stir>.05)[1]-1])
npdr.cc.pcutoff <- -log10(npdr.cc.results$pval.att[which(npdr.cc.results$pval.adj>.05)[1]-1])

library(ggplot2)
test.cc.df <- merge(univ.log10.df,npdr.cc.log10.df)
functional <- factor(c(rep("Func",length(functional.case.control)),rep("Non-Func",n.variables-length(functional.case.control))))
ggplot(test.cc.df, aes(x=univ.log10,y=npdr.cc.log10)) + geom_point(aes(colour = functional), size=4) +
  theme(text = element_text(size = 20)) +
  #geom_vline(xintercept=univ.cc.pcutoff, linetype="dashed") +
  geom_hline(yintercept=npdr.cc.pcutoff, linetype="dashed") +
  xlab("Logistic Regression -log10(P)") + ylab("NPDR -log10(P)") 

##### original (pseudo t-test) STIR
# impression is that npdr gives same resutls as original t-STIR
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

### Compare STIR and npdr
stir.log10.df <- data.frame(vars=rownames(t_sorted_multisurf),stir.log10=-log10(t_sorted_multisurf$t.pval.stir))
npdr.log10.df <- data.frame(vars=npdr.cc.results$att,npdr.log10=-log10(npdr.cc.results$pval.att))

stir.pcutoff <- -log10(t_sorted_multisurf$t.pval.stir[which(t_sorted_multisurf$t.pval.adj.stir>.05)[1]-1])
npdr.pcutoff <- -log10(npdr.cc.results$pval.attr[which(npdr.cc.results$pval.adj>.05)[1]-1])

library(ggplot2)
test.df <- merge(stir.log10.df,npdr.log10.df)
functional <- factor(c(rep("Func",length(functional.case.control)),rep("Non-Func",n.variables-length(functional.case.control))))
ggplot(test.df, aes(x=stir.log10,y=npdr.log10)) + geom_point(aes(colour = functional), size=4) +
  theme(text = element_text(size = 20)) +
  geom_vline(xintercept=stir.pcutoff, linetype="dashed") +
  geom_hline(yintercept=npdr.pcutoff, linetype="dashed") +
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

### Compare corelearn and npdr
corelearn.cc.df <- data.frame(vars=names(core.learn.case.control),rrelief=core.learn.case.control)
npdr.cc.beta.df <- data.frame(vars=npdr.cc.results$att,npdr.beta=npdr.cc.results$beta.Z.att)

corelearn.cc.cutoff <- arbitrary.cc.threshold
npdr.cc.pcutoff <- (npdr.cc.results$beta.Z.att[which(npdr.cc.results$pval.adj>.05)[1]-1])


library(ggplot2)
test.cc.df <- merge(corelearn.cc.df,npdr.cc.beta.df)
functional <- factor(c(rep("Func",length(functional.case.control)),rep("Non-Func",n.variables-length(functional.case.control))))
ggplot(test.cc.df, aes(x=rrelief,y=npdr.beta)) + geom_point(aes(colour = functional), size=4) +
  theme(text = element_text(size = 20)) +
  geom_vline(xintercept=corelearn.cc.cutoff, linetype="dashed") +
  geom_hline(yintercept=npdr.cc.pcutoff, linetype="dashed") +
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
##### Run npdrNET, penalized npdr
npdrNET.cc.results <- npdr("class", case.control.data, regression.type="glmnet", attr.diff.type="numeric-abs",
                               nbd.method="multisurf", nbd.metric = "manhattan", msurf.sd.frac=.5,
                               glmnet.alpha=1, glmnet.family="binomial",
                               padj.method="bonferroni", verbose=T)
# attributes with npdr adjusted p-value less than .05 
npdrNET.cc.results.mat <- as.matrix(npdrNET.cc.results)
# .05 regression coefficient threshold is arbitrary
# not sure why glment did not force zeros
# Negative coefficients mean irrelevant attributes for Relief scores.
# However, glmnet does not include ordinal models. 
nonzero.npdrNET.mask <- abs(npdrNET.cc.results.mat[,1])>0.05  
as.matrix(npdrNET.cc.results.mat[nonzero.npdrNET.mask,],ncol=1)

# Naively remove negative coefficients, but would be better to modify shrinkage model.
pos.npdrNET.mask <- npdrNET.cc.results.mat[,1]>0.05  
as.matrix(npdrNET.cc.results.mat[pos.npdrNET.mask,],ncol=1)

# functional attribute detection stats
npdrNET.cc.positives <- names(npdrNET.cc.results.mat[nonzero.npdrNET.mask,]) # p.adj<.05
npdrNET.cc.detect.stats <- detectionStats(functional.case.control, npdrNET.cc.positives)
cat(npdrNET.cc.detect.stats$report)

## Random Forest
library(randomForest)
ranfor.cc.fit <- randomForest(as.factor(class) ~ ., data = case.control.data) 
rf.cc.importance <- importance(ranfor.cc.fit)  
rf.cc.sorted<-sort(rf.cc.importance, decreasing=T, index.return=T)
#rf.cc.sorted$ix
#rownames(rf.cc.importance)[rf.cc.sorted$ix]
#cbind(rownames(rf.cc.importance)[rf.cc.sorted$ix],rf.cc.importance[rf.cc.sorted$ix])[1:25,]

### Compare corelearn and npdr
rf.cc.df <- data.frame(vars=rownames(rf.cc.importance),MeanDecreaseGini=rf.cc.importance)
npdr.cc.beta.df <- data.frame(vars=npdr.cc.results$att,npdr.beta=npdr.cc.results$beta.Z.att)

rf.cc.cutoff <- 0
npdr.cc.pcutoff <- (npdr.cc.results$beta.Z.att[which(npdr.cc.results$pval.adj>.05)[1]-1])


library(ggplot2)
test.cc.df <- merge(rf.cc.df,npdr.cc.beta.df)
functional <- factor(c(rep("Func",length(functional.case.control)),rep("Non-Func",n.variables-length(functional.case.control))))
ggplot(test.cc.df, aes(x=MeanDecreaseGini,y=npdr.beta)) + geom_point(aes(colour = functional), size=4) +
  theme(text = element_text(size = 20)) +
  geom_vline(xintercept=rf.cc.cutoff, linetype="dashed") +
  geom_hline(yintercept=npdr.cc.pcutoff, linetype="dashed") +
  xlab("Random Forest Score") + ylab("NPDR Standardized Coefficient") 


## Testing out penalized neighbor idea

my.attrs <- case.control.data[,colnames(case.control.data)!="class"]
my.pheno <- as.numeric(as.character(case.control.data[,colnames(case.control.data)=="class"]))
neighbor.pairs.idx <- nearestNeighbors(my.attrs, 
                                       nb.method="relieff", nb.metric="manhattan", 
                                       sd.frac = .5, k=0,
                                       attr_removal_vec_from_dist_calc=NULL)

Ridx_vec <- neighbor.pairs.idx[,"Ri_idx"]
NNidx_vec <- neighbor.pairs.idx[,"NN_idx"]

attr.idx <- 1
my.attr <- my.attrs[,attr.idx] 

num.samp <- nrow(my.attrs)
knnSURF(num.samp,.5)
neighborhood.betas <- rep(0,num.samp)
neighborhood.pvals <- rep(0,num.samp)
for (Ridx in 1:num.samp){
  #Ridx <- 51
  Ri.attr.vals <- my.attr[Ridx]
  NN.attr.vals <- my.attr[NNidx_vec[Ridx_vec==Ridx]]
  attr.diff.vec <- npdrDiff(Ri.attr.vals, NN.attr.vals, diff.type="numeric-abs")

  Ri.pheno.vals <- my.pheno[Ridx]
  NN.pheno.vals <- my.pheno[NNidx_vec[Ridx_vec==Ridx]]
  pheno.diff.vec <- npdrDiff(Ri.pheno.vals, NN.pheno.vals, diff.type="match-mismatch")
  pheno.diff.vec <- factor(pheno.diff.vec, levels=c(0,1))
  mod <- glm(pheno.diff.vec ~ attr.diff.vec, family=binomial(link=logit))
  fit <- summary(mod)
  beta_a <- coef(fit)[2, 1]         # raw beta coefficient, slope (not standardized)
  beta_zscore_a <- coef(fit)[2, 3]  # standardized beta coefficient (col 3)
  ## use one-side p-value to test H1: beta>0 for case-control npdr scores
  pval_beta_a <- pt(beta_zscore_a, mod$df.residual, lower = FALSE)  # one-sided p-val
  neighborhood.betas[Ridx] <- beta_zscore_a
  neighborhood.pvals[Ridx] <- pval_beta_a
}
cbind(neighborhood.betas, neighborhood.pvals, my.pheno)
beta_zscore_ave <- mean(neighborhood.betas)
mean(neighborhood.pvals)
pt(beta_zscore_ave, knnSURF(num.samp,.5), lower = FALSE) 
pnorm(beta_zscore_ave, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)

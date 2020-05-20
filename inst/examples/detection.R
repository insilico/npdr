# devtools::install_github("insilico/npdr") # npdr install
#library(devtools)
#install_github("insilico/privateEC")
library(npdr)
library(ggplot2)
library(reshape2)

#==================== Simulate Data with npdr simulate2 ===================================#
##### simulate case-control interaction effect data 
n.samples <- 300     # 100 samples in train/holdout/test
n.variables <- 100   # 100 features
label <- "class" # tells simulator to do case/control and adds this colname
type <- "interactionErdos" # or mainEffect
#type <-"mainEffect"
bias <- .8           # .8 strong main effect, can also use for interactions 
prob.connect <- .25  # increase this probability for stronger interaction effect
pct.signals <- 0.1   # pct functional features
verbose <- FALSE
imbal <- .5   # .75 imbalance 

case.control.3sets <- createSimulation2(num.samples=n.samples,
                                        num.variables=n.variables,
                                        pct.imbalance=imbal, 
                                        pct.signals=pct.signals,
                                        main.bias=bias,
                                        interaction.bias=bias,
                                        hi.cor=0.8, #default
                                        lo.cor=0.2, #default
                                        prob.connected = prob.connect,
                                        out.degree = NULL,  # use for scalefree
                                        mix.type=NULL,
                                        save.file=NULL,
                                        label="class",
                                        sim.type=type,
                                        pct.mixed=0.5, # if mix.type not null
                                        pct.train=1/3,
                                        pct.holdout=1/3,
                                        pct.validation=1/3,
                                        plot.graph=F,
                                        verbose=verbose,
                                        use.Rcpp=F)

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
univ.cc.results <- uniReg(outcome="class", dataset=case.control.data, regression.type="binomial")
univ.cc.results[1:10,]
# don't expect any less than .05 for interaction simulations
univ.cc.positives<-univ.cc.results[univ.cc.results[,"p.adj"]<.05,]

##### Run npdr
npdr.cc.results <- npdr("class", case.control.data, regression.type="binomial", attr.diff.type="numeric-abs",
                        nbd.method="multisurf", nbd.metric = "manhattan", msurf.sd.frac=.5, k=0,
                        neighbor.sampling="none", dopar.nn = F,
                        padj.method="bonferroni", verbose=T)

npdr.cc.results[npdr.cc.results$pval.adj<.05,] # pval.adj, first column
npdr.cc.positives <- npdr.cc.results %>% filter(pval.adj<.05) %>% pull(att)

# sensitivity = TP/(TP+FN)
npdr.cc.detect.stats <- detectionStats(functional.case.control, npdr.cc.positives)
univ.cc.detect.stats <- detectionStats(functional.case.control, univ.cc.positives)
npdr.cc.detect.stats$TPR
univ.cc.detect.stats$TPR #nan

npdrDetected <- function(results.df,functional,top.pct){
  # results is the output of npdr, functional are known functional vars
  top.num <- floor(top.pct * nrow(results.df))
  top.vars <- results.df %>% top_n(-top.num, pval.att) %>% pull(att) 
  power <- detectionStats(functional,top.vars)$TP  # how many of top.pct are true
  ifelse(is.nan(power),0,power)/length(functional)
}

pcts <- seq(0,1,.05)
npdr.detected <- sapply(pcts,function(p){npdrDetected(npdr.cc.results,functional.case.control,p)})
plot(pcts,npdr.detected)

core.learn.cc <- CORElearn::attrEval("class", data = case.control.data,
                                     estimator = "ReliefFequalK",
                                     costMatrix = NULL,
                                     outputNumericSplits=FALSE,
                                     kNearestEqual = knnSURF(nrow(case.control.data),.5))
core.learn.cc.order <- order(core.learn.cc, decreasing = T)
corelearn.df <- data.frame(att=names(core.learn.cc)[core.learn.cc.order], rrelief=core.learn.cc[core.learn.cc.order])

reliefDetected <- function(results.df,functional,top.pct){
  # results is the output of npdr, functional are known functional vars
  top.num <- floor(top.pct * nrow(results.df))
  top.vars <- results.df %>% top_n(top.num, rrelief) %>% pull(att) 
  power <- detectionStats(functional,top.vars)$TP  # how many of top.pct are true
  ifelse(is.nan(power),0,power)/length(functional)
}

relief.detected <- sapply(pcts,function(p){reliefDetected(corelearn.df,functional.case.control,p)})
plot(pcts,relief.detected)

df <- data.frame(pcts=pcts, NPDR=npdr.detected, Relief=relief.detected)
melt.df <- melt(data = df, measure.vars = c("NPDR", "Relief"))
ggplot(melt.df, aes(x=pcts, y=value, group=variable)) +
  geom_line(aes(linetype=variable)) +
  geom_point(aes(shape=variable, color=variable), size=2) +
  scale_color_manual(values = c("#FC4E07", "#E7B800")) +
  xlab("Percent Selected") +
  ylab("Percent Correct") +
  ggtitle("Power to Detect Functional Variables")

# auDC: area under the detection curve
sum(npdr.detected)/length(npdr.detected)       # .83
sum(relief.detected)/length(relief.detected)   # .72


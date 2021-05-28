# test gwas added to createSimulation2()
library(npdr)
library(CORElearn)
library(randomForest)
library(reshape)

#main.bias <- 0.4      # this is probably a large effect
#pct.mixed <- 0.5
#nbias <- round(pct.mixed*num.variables)
#num.main <- round(pct.mixed*nbias) 
#num.int <- round((1 - pct.mixed)*nbias)
#sim.type <- "mainEffect"
#mix.type="main-interactionErdos"

sim.type="interactionErdos"
num.samples=100
num.variables=1000
pct.signals=0.1
interaction.bias=0.8
hi.cor=0.8
lo.cor=0.1
prob.connected=0.95 # for erdos
avg.maf=0.2

#notes: test imbalance, main effect, mixed effect, what interaction parameter range?

dataset <- createSimulation2(num.samples=num.samples,
                             num.variables=num.variables,
                             pct.imbalance=0.5,
                             pct.signals=pct.signals,
                             main.bias=main.bias,
                             interaction.bias=interaction.bias,
                             hi.cor=hi.cor,
                             lo.cor=lo.cor,
                             mix.type=mix.type,
                             label="class",
                             sim.type=sim.type,
                             pct.mixed=0.5,
                             pct.train=0.5,
                             pct.holdout=0.5,
                             pct.validation=0,
                             plot.graph=FALSE,
                             verbose=TRUE,
                             use.Rcpp=FALSE,
                             prob.connected=prob.connected, #erdos
                             out.degree=out.degree,      #scale-free
                             data.type="discrete",       #snps
                             avg.maf=avg.maf)            #snps
dats <- rbind(dataset$train, dataset$holdout, dataset$validation)
dats <- dats[order(dats[,ncol(dats)]),]

case.control.data <- dats[,-ncol(dats)]  # remove class

# calcluate average minor allele frequency
n.samp<-nrow(case.control.data)
n.attr<-ncol(case.control.data)
maf <- function(i){sum(case.control.data[,i])/(2*n.samp)}
cat("average minor allele frequency\n")
mean(sapply(1:n.attr,maf))

#notes: compare AM and GM

npdr.gwas.results <- npdr(dats[,ncol(dats)], case.control.data, regression.type="binomial", 
                        attr.diff.type="allele-sharing",   # AM
                        #attr.diff.type="match-mismatch",  # GM
                        nbd.method="relieff", 
                        #nbd.method="multisurf", 
                        nbd.metric = "manhattan", msurf.sd.frac=.5, k=0,
                        neighbor.sampling="none",
                        dopar.nn = F, padj.method="bonferroni", verbose=T)

functional.vars<-dataset$signal.names
npdr.gwas.positives <- npdr.gwas.results %>% filter(pval.adj<.05) %>% pull(att)
cat("Positives\n")
cat(npdr.gwas.positives)
cat("\n")
npdr.gwas.detect.stats <- detectionStats(functional.vars, npdr.gwas.positives)
cat(npdr.gwas.detect.stats$report)

# does e^beta have a meaning for npdr?
exp(npdr.gwas.results  %>% filter(pval.adj<.05) %>% pull(beta.raw.att))

#### Random Forest
ranfor.gwas.fit <- randomForest(as.factor(class) ~ ., data = dats) 
rf.gwas.importance <- importance(ranfor.gwas.fit)  
rf.gwas.sorted<-sort(rf.gwas.importance, decreasing=T, index.return=T)
rf.gwas.df <-data.frame(att=rownames(rf.gwas.importance)[rf.gwas.sorted$ix],rf.scores=rf.gwas.sorted$x)

##### Regular Relief
relief.gwas <- CORElearn::attrEval(as.factor(class) ~ ., data = dats, #"Class", data = sim2.df,
                                       estimator = "ReliefFequalK",
                                       costMatrix = NULL,
                                       outputNumericSplits=FALSE,
                                       kNearestEqual = knnSURF(nrow(dats),.5)) # fn from npdr

relief.gwas.order <- order(relief.gwas, decreasing = T)
relief.gwas.df <- data.frame(att=names(relief.gwas)[relief.gwas.order], rrelief=relief.gwas[relief.gwas.order])
#relief.gwas.df[1:10,]

##### Functions for creating detection rate vs percentage variables selected

reliefDetected <- function(results.df,functional,top.pct){
  # results.df is the results output of relief from corelearn 
  # functional are known functional var names
  # top.pct is percentile of top relief variables to compare with functional
  # results order is high to low score
  top.num <- floor(top.pct * nrow(results.df))
  top.vars <- results.df %>% top_n(top.num, rrelief) %>% pull(att) # rrelief column of results.df
  power <- detectionStats(functional,top.vars)$TP  # npdr:: fn, how many of top.pct are true
  ifelse(is.nan(power),0,power)/length(functional) # if nan, return 0, normalize by num of functional
}

rfDetected <- function(results.df,functional,top.pct){
  # results.df is the results output of random forest 
  # functional are known functional var names
  # top.pct is percentile of top random forest variables to compare with functional
  # results order is high to low score
  top.num <- floor(top.pct * nrow(results.df))
  top.vars <- results.df %>% top_n(top.num, rf.scores) %>% pull(att)  # rf.scores column for rf results
  power <- detectionStats(functional,top.vars)$TP  # npdr:: fn, how many of top.pct are true
  ifelse(is.nan(power),0,power)/length(functional) # if nan, return 0, normalize by num of functional
}

npdrDetected <- function(results.df,functional,top.pct){
  # results.df is the output of npdr
  # functional are known functional var names
  # top.pct is percentile of top npdr variables to compare with functional
  # results order is low P value to high P value
  top.num <- floor(top.pct * nrow(results.df))
  top.vars <- results.df %>% top_n(-top.num, pval.att) %>% pull(att)  # pval.att is npdr specific
  power <- detectionStats(functional,top.vars)$TP  # npdr:: fn, how many of top.pct are true
  ifelse(is.nan(power),0,power)/length(functional) # if nan, return 0, normalize by num of functional
}


##### compute detection of functional variables for the methods
pcts <- seq(0,1,.05)
rf.detected <- sapply(pcts,function(p){rfDetected(rf.gwas.df,functional.vars,p)})
relief.detected <- sapply(pcts,function(p){reliefDetected(relief.gwas.df,functional.vars,p)})
npdr.detected <- sapply(pcts,function(p){npdrDetected(npdr.gwas.results,functional.vars,p)})


# plot detection curves (DC) for three methods
df <- data.frame(pcts=pcts, NPDR=npdr.detected, Relief=relief.detected,RForest=rf.detected)
melt.df <- melt(data = df, measure.vars = c("NPDR", "Relief","RForest"))
ggplot(melt.df, aes(x=pcts, y=value, group=variable)) +
  geom_line(aes(linetype=variable)) +
  geom_point(aes(shape=variable, color=variable), size=2) +
  scale_color_manual(values = c("#FC4E07", "#E7B800", "#228B22")) +
  xlab("Percent Selected") +
  ylab("Percent Correct") +
  ggtitle("Power to Detect Functional Variables") + 
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw() 

# auDC: area under the detection curve for three methods
sum(rf.detected)/length(rf.detected)       
sum(npdr.detected)/length(npdr.detected)       
sum(relief.detected)/length(relief.detected)   

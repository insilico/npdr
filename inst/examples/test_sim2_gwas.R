# test gwas added to createSimulation2()
#library(npdr)

num.samples <- 1000
num.variables <- 100
main.bias <- 0.4      # this is probably a large effect
pct.mixed <- 0.5
nbias <- round(pct.mixed*num.variables)
num.main <- round(pct.mixed*nbias) 
num.int <- round((1 - pct.mixed)*nbias)

sim.type <- "interactionErdos"
#sim.type <- "mainEffect"
mix.type="main-interactionErdos"

#notes: test imbalance, main effect, mixed effect, what interaction parameter range?

dataset <- createSimulation2(num.samples=num.samples,
                             num.variables=num.variables,
                             pct.imbalance=0.5,
                             pct.signals=0.1,
                             main.bias=main.bias,
                             interaction.bias=.8,
                             hi.cor=0.8,
                             lo.cor=0.1,
                             mix.type=mix.type,
                             label="class",
                             sim.type=sim.type,
                             pct.mixed=0.5,
                             pct.train=0.5,
                             pct.holdout=0.5,
                             pct.validation=0,
                             plot.graph=F,
                             verbose=T,
                             use.Rcpp=F,
                             prob.connected=0.95,
                             out.degree=98,
                             data.type="discrete")
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
                        attr.diff.type="allele-sharing",
                        #attr.diff.type="match-mismatch",
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

npdr.gwas.results[1:10,]
exp(npdr.gwas.results  %>% filter(pval.adj<.05) %>% pull(beta.raw.att))

library(randomForest)
ranfor.gwas.fit <- randomForest(as.factor(class) ~ ., data = dats) 
rf.gwas.importance <- importance(ranfor.gwas.fit)  
rf.gwas.sorted<-sort(rf.gwas.importance, decreasing=T, index.return=T)
data.frame(vars=rownames(rf.gwas.importance)[rf.gwas.sorted$ix],scores=rf.gwas.sorted$x)[1:10,]

# generate replicates of auRC and auPRC from simulated data with varying degrees of class imbalance
#
# generates data for plotting script: make_auRC-auPRC_boxplots(imbalanced-hitmiss-nbds).R

library(npdr)

library(reshape2)
library(ggplot2)
library(PRROC)

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

# sim.type (options)
#
# "mainEffect": simple main effects
# "mainEffect_Erdos-Renyi": main effects with added correlation from Erdos-Renyi network
# "mainEffect_Scalefree": main effects with added correlation from Scale-free network
# "interactionErdos": interaction effects from Erdos-Renyi network
# "interactionScalefree": interaction effects from Scale-free network
# "mixed": main effects and interaction effects
#     mix.type (options)
#
#     "main-interactionErdos": main effects and interaction effects from Erdos-Renyi network
#     "main-interactionScalefree": main effects and interaction effects from Scale-free network

# data.type (options)
#
# "continuous": random normal data N(0,1) (e.g., gene expression data)
# "discrete": random binomial data B(n=2,prob) (e.g., GWAS data)

num.samples <- 100                         # number of samples
num.variables <- 100                       # number of variables
pct.imbalance <- 0.5                       # fraction of instances that are cases (class=1)
pct.signals <- 0.1                         # fraction of num.variables that are functional
main.bias <- 0.1                           # effect size parameter for main effects
interaction.bias <- 0.95                   # effect size parameter for interaction effects
mix.type <- "mainEffect-interactionErdos"  # mixed simulation type
sim.type <- "interactionErdos"             # simulation type
#sim.type <- "mainEffect"
data.type <- "discrete"

###################################################################################################
#
separate.hitmiss.nbds <- T # run with both T/F to generate all files for plots
#
###################################################################################################

num.iter <- 1 # generate num.iter replicates for each level of imbalance (THIS WILL TAKE AWHILE!!!)
imbalances <- c(0.1, 0.2, 0.3, 0.4, 0.5)
for(iter in 1:length(imbalances)){
  cat("Class Imbalance: ",imbalances[iter],"\n")
  
  pct.imbalance <- imbalances[iter]
  
  chosen.k.mat <- matrix(0,nrow=num.iter,ncol=num.variables)
  accu.vec <- numeric()
  auPRC.vec <- numeric()

  set.seed(1989)
  for(i in 1:num.iter){
    cat("Replicate Data Set: ",i,"\n")
    
    # simulated data
    dataset <- createSimulation2(num.samples=num.samples,
                                 num.variables=num.variables,
                                 pct.imbalance=pct.imbalance,
                                 pct.signals=pct.signals,
                                 main.bias=main.bias,
                                 interaction.bias=1,
                                 hi.cor=0.85,
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
                                 out.degree=(num.variables-2),
                                 data.type=data.type)
    dats <- rbind(dataset$train, dataset$holdout, dataset$validation)
    dats <- dats[order(dats[,ncol(dats)]),]
    
    # run vwak on data set to give beta and p-value matrices
    out <- vwak(dats=dats,                       
                k.grid=NULL,                                 # can provide custom grid of k's
                verbose=T,                       
                attr.diff.type="allele-sharing",             # attribute diff used by npdr
                separate.hitmiss.nbds=separate.hitmiss.nbds, # set = T for equal size hit/miss nbds
                label="class")
    
    betas <- out$beta.mat # num.variables x (num.samples - 1) matrix of betas
    pvals <- out$pval.mat # num.variables x (num.samples - 1) matrix of p-values 
    best.ks <- apply(betas,1,which.max) # best k's computed from max beta for each attribute
    best.betas <- numeric() # betas corresponding to best k's 
    best.pvals <- numeric() # p-values corresponding to best k's
    for(j in 1:nrow(betas)){
      best.betas[j] <- betas[j,as.numeric(best.ks[j])]
      best.pvals[j] <- pvals[j,as.numeric(best.ks[j])]
    }
    
    # data frame of attributes, betas, and p-values from vwak
    df.betas <- data.frame(att=row.names(betas),
                           betas=best.betas,
                           pval.att=best.pvals)
    df.betas <- df.betas[order(df.betas[,2],decreasing=T),] # sort data frame by decreasing beta
    
    functional.vars <- dataset$signal.names # functional variable names

    pcts <- seq(0,1,.05)
    npdr.detected <- sapply(pcts,function(p){npdrDetected(df.betas,functional.vars,p)})
    
    chosen.k.mat[i,] <- best.ks
    if(i==1){
      colnames(chosen.k.mat) <- names(best.ks)
    }
    
    accu.vec[i] <- sum(npdr.detected)/length(npdr.detected) # area under the recall curve
  
    idx.func <- which(c(as.character(df.betas[,"att"]) %in% functional.vars)==T)
    func.betas <- df.betas[idx.func,"betas"] # functional variable betas
    neg.betas <- df.betas[-idx.func,"betas"] # noise variable betas
    
    # precision-recall curve and area
    pr.npdr <- PRROC::pr.curve(scores.class0 = func.betas,
                               scores.class1 = neg.betas, 
                               curve = T)
    
    plot(pr.npdr) # plot precision-recall curve
    
    auPRC.vec[i] <- pr.npdr$auc.integral # area under the precision-recall curve
    
    
  }
  accu.df <- data.frame(auRC=accu.vec,auPRC=auPRC.vec)
  #setwd("C:/Users/bdawk/Documents/KNN_project_output") will need to change to desired directory
  
  balanced.hit.miss <- strsplit(as.character(separate.hitmiss.nbds), split="")[[1]][1]
  file <- paste("separate-hitmiss-nbds-",balanced.hit.miss,"_",sim.type,"_",data.type,"_imbalance-",iter,".csv",sep="")
  write.csv(accu.df,file,row.names=F)
  
  file <- paste("separate-hitmiss-nbds-",balanced.hit.miss,"_",sim.type,"_k-matrix_",data.type,"_imbalance-",iter,".csv",sep="")
  write.csv(chosen.k.mat,file,row.names=F)
  
}















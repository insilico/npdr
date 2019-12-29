# compute auRC and auPRC for NPDR, Relief, and Random Forest from 30 replicated data sets 
# Discrete (SNP) Data
library(npdr)
library(CORElearn)
library(randomForest)
library(reshape2)
library(ggplot2)
library(PRROC)

show.plots = T
save.files = F
run.pairs = F
num.iter <- 1  # just run one simulation
#num.iter <- 30 # 30 replicate simulations will take several minutes
if (save.files){
  cat("Results files for ",num.iter, " replicate simulation(s) will be saved in ", getwd(),".", sep="")
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

num.samples <- 100
num.variables <- 100
main.bias <- 0.5      
pct.mixed <- .5 
pct.imbalance <- .5  # percent of effects that are main effects, must also use sim.type = "mixed" 
mix.type <-  "main-interactionErdos"
#pct.imbalance <- 0.25 # 75% case - 25% ctrl
pct.signals <- 0.1 
nbias <- round(pct.signals*num.variables)
sim.type <- "interactionErdos"  # "mixed" for mixed main and interactions
data.type <- "discrete"  # or "continuous" for standard normal data

auRC.npdr.multisurf <- numeric()
auRC.npdr.fixedk <- numeric()
auRC.relief <- numeric()
auRC.RF <- numeric()

auPRC.vec.npdr.multisurf <- numeric()
auPRC.vec.npdr.fixedk <- numeric()
auPRC.vec.relief <- numeric()
auPRC.vec.RF <- numeric()
set.seed(1989)

for(iter in 1:num.iter){
  cat("Iteration: ",iter,"\n")
  
  dataset <- createSimulation2(num.samples=num.samples,
                               num.variables=num.variables,
                               pct.imbalance=pct.imbalance,
                               pct.signals=pct.signals,
                               main.bias=main.bias,
                               interaction.bias=1,
                               hi.cor=0.8,
                               lo.cor=0.1,
                               mix.type=mix.type,
                               label="class",
                               sim.type=sim.type,
                               pct.mixed=pct.mixed,
                               pct.train=0.5,
                               pct.holdout=0.5,
                               pct.validation=0,
                               plot.graph=F,
                               verbose=T,
                               use.Rcpp=F,
                               prob.connected=0.3,
                               out.degree=(num.variables-2),
                               data.type=data.type)
  dats <- rbind(dataset$train, dataset$holdout, dataset$validation)
  dats <- dats[order(dats[,ncol(dats)]),]
  
  if (run.pairs){  # skip this
    # pairwise interaction p-values or beta's if you want from inbixGAIN.R
    intPairs.mat <- getInteractionEffects("class", dats, 
                                      regressionFamily = "binomial", 
                                      numCovariates = 0,
                                      writeBetas = FALSE, 
                                      excludeMainEffects = FALSE, 
                                      interactOutput = "Pvals", 
                                      transformMethod = "", 
                                      numCores = 1, 
                                      verbose = T) 
    colnames(intPairs.mat) <- colnames(dats)[1:(ncol(dats)-1)]
    rownames(intPairs.mat) <- colnames(dats)[1:(ncol(dats)-1)]
    row <- 1
    for (var in rownames(intPairs.mat)){
      cat(var,": ", intPairs.mat[row,intPairs.mat[row,]>0 & intPairs.mat[row,]<.001],"\n")
      row <- row + 1
    }
  }
  
  # npdr - multisurf
  npdr.results1 <- npdr("class", dats, regression.type="binomial", 
                        attr.diff.type="allele-sharing",   
                        #nbd.method="relieff", 
                        nbd.method="multisurf", 
                        nbd.metric = "manhattan", msurf.sd.frac=.5, k=0,
                        neighbor.sampling="none", separate.hitmiss.nbds=F,
                        dopar.nn = T, dopar.reg=T, padj.method="bonferroni", verbose=T)
  df1 <- data.frame(att=npdr.results1$att,
                    beta=npdr.results1$beta.Z.att,
                    pval=npdr.results1$pval.adj)
  
  functional.vars <- dataset$signal.names
  npdr.positives1 <- npdr.results1 %>% filter(pval.adj<.05) %>% pull(att)
  
  df1 <- na.omit(df1)
  idx.func <- which(c(as.character(df1[,"att"]) %in% functional.vars))
  func.betas1 <- df1[idx.func,"beta"]
  neg.betas1 <- df1[-idx.func,"beta"]
  
  pr.npdr1 <- PRROC::pr.curve(scores.class0 = func.betas1,
                              scores.class1 = neg.betas1, 
                              curve = T)
  
  if (show.plots){plot(pr.npdr1)}
  
  npdr.detect.stats1 <- detectionStats(functional.vars, npdr.positives1)
  
  # npdr - fixed k
  npdr.results2 <- npdr("class", dats, regression.type="binomial", 
                        attr.diff.type="allele-sharing",   
                        nbd.method="relieff", 
                        #nbd.method="multisurf", 
                        nbd.metric = "manhattan", msurf.sd.frac=.5, k=0,
                        neighbor.sampling="none", separate.hitmiss.nbds=T,
                        dopar.nn = T, dopar.reg=T, padj.method="bonferroni", verbose=T)
  
  npdr.positives2 <- npdr.results2 %>% filter(pval.adj<.05) %>% pull(att)
  
  df2 <- data.frame(att=npdr.results2$att,
                    beta=npdr.results2$beta.Z.att,
                    pval=npdr.results2$pval.adj)
  df2 <- na.omit(df2)
  
  idx.func <- which(c(as.character(df2[,"att"]) %in% functional.vars)==TRUE)
  func.betas2 <- df2[idx.func,"beta"]
  neg.betas2 <- df2[-idx.func,"beta"]
  
  pr.npdr2 <- PRROC::pr.curve(scores.class0 = func.betas2,
                              scores.class1 = neg.betas2, 
                              curve = T)
  
  if (show.plots){plot(pr.npdr2)}
  
  npdr.detect.stats2 <- detectionStats(functional.vars, npdr.positives2)
  
  #### Random Forest
  ranfor.fit <- randomForest(as.factor(class) ~ ., data = dats) 
  rf.importance <- importance(ranfor.fit)  
  rf.sorted<-sort(rf.importance, decreasing=T, index.return=T)
  rf.df <-data.frame(att=rownames(rf.importance)[rf.sorted$ix],rf.scores=rf.sorted$x)
  
  idx.func <- which(c(as.character(rf.df$att) %in% functional.vars)==TRUE)
  func.scores.rf <- rf.df[idx.func,"rf.scores"]
  neg.scores.rf <- rf.df[-idx.func,"rf.scores"]
  
  pr.rf <- PRROC::pr.curve(scores.class0 = func.scores.rf,
                           scores.class1 = neg.scores.rf, 
                           curve = T)
  if (show.plots){plot(pr.rf)}
  
  ##### Regular Relief
  relief <- CORElearn::attrEval(as.factor(class) ~ ., data = dats, 
                                estimator = "ReliefFequalK",
                                costMatrix = NULL,
                                outputNumericSplits=FALSE,
                                kNearestEqual = floor(knnSURF(nrow(dats),.5)/2)) # fn from npdr
  
  relief.order <- order(relief, decreasing = T)
  relief.df <- data.frame(att=names(relief)[relief.order], rrelief=relief[relief.order])
  
  idx.func <- which(c(as.character(relief.df$att) %in% functional.vars)==TRUE)
  func.scores.relief <- relief.df[idx.func,"rrelief"]
  neg.scores.relief <- relief.df[-idx.func,"rrelief"]
  
  pr.relief <- PRROC::pr.curve(scores.class0 = func.scores.relief,
                               scores.class1 = neg.scores.relief, 
                               curve = T)
  if (show.plots){plot(pr.relief)}
  
  pcts <- seq(0,1,.05)
  rf.detected <- sapply(pcts,function(p){rfDetected(rf.df,functional.vars,p)})
  relief.detected <- sapply(pcts,function(p){reliefDetected(relief.df,functional.vars,p)})
  npdr.detected.multisurf <- sapply(pcts,function(p){npdrDetected(npdr.results1,functional.vars,p)})
  npdr.detected.fixedk <- sapply(pcts,function(p){npdrDetected(npdr.results2,functional.vars,p)})
  
  if (show.plots){
  # plot recall curves (RC) for several methods
  df <- data.frame(pcts=pcts, NPDR.MultiSURF=npdr.detected.multisurf, 
                   NPDR.Fixed.k=npdr.detected.fixedk,
                   Relief=relief.detected,
                   RForest=rf.detected)
  melt.df <- melt(data = df, measure.vars = c("NPDR.MultiSURF", "NPDR.Fixed.k", "Relief","RForest"))
  gg <- ggplot(melt.df, aes(x=pcts, y=value, group=variable)) +
    geom_line(aes(linetype=variable)) +
    geom_point(aes(shape=variable, color=variable), size=2) +
    scale_color_manual(values = c("#FC4E07", "brown", "#E7B800", "#228B22")) +
    xlab("Percent Selected") +
    ylab("Percent Correct") +
    ggtitle("Power to Detect Functional Variables") + 
    theme(plot.title = element_text(hjust = 0.5)) + theme_bw() 
  print(gg)
  
  # plot precision-recall curves (PRC) for several methods
  method.vec <- c(rep('NPDR.MultiSURF',length=nrow(pr.npdr1$curve)),
                  rep('NPDR.Fixed.k',length=nrow(pr.npdr2$curve)),
                  rep('Relief',length=nrow(pr.relief$curve)),
                  rep('RForest',length=nrow(pr.rf$curve)))
  df <- data.frame(recall=c(pr.npdr1$curve[,1],pr.npdr2$curve[,1],pr.relief$curve[,1],pr.rf$curve[,1]),
                   precision=c(pr.npdr1$curve[,2],pr.npdr2$curve[,2],pr.relief$curve[,2],pr.rf$curve[,2]),
                   method=method.vec)
  gg <- ggplot(df, aes(x=recall, y=precision, group=method)) +
    geom_line(aes(linetype=method)) +
    geom_point(aes(shape=method, color=method), size=2) +
    scale_color_manual(values= c("#FC4E07", "brown", "#E7B800", "#228B22")) +
    xlab("Recall") +
    ylab("Precision") +
    ggtitle("Precision-Recall Curves: Comparison of Methods") +
    theme(plot.title = element_text(hjust=0.5)) + theme_bw()
  print(gg)
  }
  
  # auRC: area under the recall curve for several methods
  auRC.RF[iter] <- sum(rf.detected)/length(rf.detected)       
  auRC.npdr.multisurf[iter] <- sum(npdr.detected.multisurf)/length(npdr.detected.multisurf)     
  auRC.npdr.fixedk[iter] <- sum(npdr.detected.fixedk)/length(npdr.detected.fixedk)
  auRC.relief[iter] <- sum(relief.detected)/length(relief.detected)
  
  # auPRC: area under the precision-recall curve for several methods
  auPRC.vec.npdr.multisurf[iter] <- pr.npdr1$auc.integral
  auPRC.vec.npdr.fixedk[iter] <- pr.npdr2$auc.integral
  auPRC.vec.relief[iter] <- pr.relief$auc.integral
  auPRC.vec.RF[iter] <- pr.rf$auc.integral
  
}

if (save.files){
# save results
df.save <- data.frame(cbind(iter=seq(1,30,by=1),
                            auRC.RForest=auRC.RF,
                            auRC.NPDR.MultiSURF=auRC.npdr.multisurf,
                            auRC.NPDR.Fixed.k=auRC.npdr.fixedk,
                            auRC.Relief=auRC.relief,
                            auPRC.RForest=auPRC.vec.RF,
                            auPRC.NPDR.MultiSURF=auPRC.vec.npdr.multisurf,
                            auPRC.NPDR.Fixed.k=auPRC.vec.npdr.fixedk,
                            auPRC.Relief=auPRC.vec.relief))
df.save <- apply(df.save,2,as.numeric)
#setwd("C:/Users/bdawk/Documents/KNN_project_output") will need to change to desired directory
file <- paste("auRC-auPRC_iterates_methods-comparison_",data.type,"-",sim.type,".csv",sep="")
write.csv(df.save,file,row.names=F)
}

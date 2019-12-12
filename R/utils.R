#=========================================================================#
#' checkPackages
#'
#' Utility to check for package installation
#'
#' @param pkg character vector of package names to install if not already.
#' @return null There is nothing returned.  
#' @examples
#' packages <- c("ggplot2", "CORElearn", "reshape2", "dplyr", "pROC", "plotROC")
#' checkPackages(packages)  
#' @export
checkPackages <- function(pkg){
  # check.packages function: install and load multiple R packages.
  # Check to see if packages are installed. Install them if they are not, 
  # then load them into the R session.
  # https://gist.github.com/smithdanielle/9913897
  
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, repos = "http://cran.us.r-project.org", dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#=========================================================================#
#' knnSURF
#'
#' Theoretical value for the number of expected neighbors for SURF or multiSURF
#'
#' @param m.samples number of samples in data.
#' @param sd.frac fraction of the standard deviation from the mean of all pairwise distances, dead-band. The default value used by the SURF and multiSURF algorithms is 1/2.  
#' @return knn Number of neighbors.  
#' @examples
#' k.surf <- knnSURF(200,.5)
#' @export
knnSURF <- function(m.samples,sd.frac=.5){
  # error function
  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
  # theoretical SURF knn formulat
  knn <- floor((m.samples-1)*(1-erf(sd.frac/sqrt(2)))/2)
  return(knn)
  }

#=========================================================================#
#' uniReg
#'
#' Univariate logistic or linear regression for a dataset.
#'
#' @param outcome string with name of class column.
#' @param dataset data matrix with predictor columns and outcome column.
#' @param regression.type "lm" or "binomial"
#' @param padj.method for p.adjust (\code{"fdr"}, \code{"bonferroni"}, ...) 
#' @param covars optional vector or matrix of covariate columns for correction. Or separate data matrix of covariates.
#' @return matrix of beta, p-value and adjusted p-value, sorted by p-value.  
#' @examples
#' lr.results <- uniReg(outcome="class", dataset=case.control.data, regression.type="binomial")
#  lr.results[lr.results[,"p.adj"]<.05] 
#' @export
uniReg <- function(outcome, dataset, regression.type="lm", padj.method="fdr", covars="none"){
  ## parse input
  if (length(outcome)==1){
    # e.g., outcome="qtrait" or outcome=101 (pheno col index) and data.set is data.frame including outcome variable
    pheno.vec <- dataset[,outcome] # get phenotype
    if (is.character(outcome)){ # example column name: outcome="qtrait"
      attr.mat <- dataset[ , !(names(dataset) %in% outcome)]  # drop the outcome/phenotype
    } else { # example column index: outcome=101
      attr.mat <- dataset[ , -outcome]  # drop the outcome/phenotype  
    }
   } else { # user specifies a separate phenotype vector
    pheno.vec <- outcome # assume users provides a separate outcome data vector
    attr.mat <- dataset # assumes data.set only contains attributes/predictors
  }
  ## set up model
  if (regression.type=="lm"){
    if (length(covars)>1){
      # 2 below means attribte stats (1 would be the intercept)
      model.func <- function(x) {as.numeric(summary(lm(pheno.vec ~ attr.mat[,x] + covars))$coeff[2,])}
    } else { # covar=="none"
      model.func <- function(x) {
        # if nrow(summary(lm(pheno.vec ~ attr.mat[,x]))$coeff) < 2 then there was an issue with the attribute
        as.numeric(summary(lm(pheno.vec ~ attr.mat[,x]))$coeff[2,])}  
    } 
    } else { # "binomial"
    if (length(covars)>1){
      #model.func <- function(x) {tidy(glm(pheno.vec ~ attr.mat[,x] + covars, family=binomial))[2,4:5]}
      model.func <- function(x) {summary(glm(pheno.vec ~ attr.mat[,x] + covars, family=binomial))$coeff[2,]}
    } else { # covar=="none"
      #model.func <- function(x) {tidy(glm(pheno.vec ~ attr.mat[,x], family=binomial))[2,4:5]}
      model.func <- function(x) {summary(glm(pheno.vec ~ attr.mat[,x], family=binomial))$coeff[2,]}
    } } # end else binomial
  #class.col <- which(colnames(dataset)==outcome)
  #predictor.cols <- which(colnames(dataset)!=outcome)
  num.attr <- ncol(attr.mat)
  if (is.null(num.attr)){ # if there is just one attribute
    attr.mat <- as.matrix(attr.mat)
    num.attr <- ncol(attr.mat) # num.attr <- 1
  }
  beta_pvals <- t(sapply(1:num.attr, model.func)) # stats for all predictors
  univariate.padj <- p.adjust(beta_pvals[,4], method=padj.method) # fdr
  univariate.padj <- as.numeric(format(univariate.padj, scientific = T, digits=5))
  betas <- as.numeric(format(beta_pvals[,1], scientific = F, digits=5))
  betas.Z.att <- as.numeric(format(beta_pvals[,3], scientific = F, digits=5))
  pvals <- as.numeric(format(beta_pvals[,4], scientific = T, digits=5))
  beta_pvals <- cbind(betas,betas.Z.att,pvals,univariate.padj) # adjusted p-val column
  row.names(beta_pvals)<- colnames(attr.mat) # add predictor names
  #row.names(beta_pvals)<- colnames(dataset)[-class.col] # add predictor names
  beta_pvals_sorted <- beta_pvals[order(as.numeric(beta_pvals[,3]), decreasing = F),] # sort by pval
  if (num.attr==1){ # special case of only 1 attribute
    names(beta_pvals_sorted) <- c("beta", "beta.Z.att", "pval", "p.adj")
  } else{  # multiple attributes, typical
    colnames(beta_pvals_sorted) <- c("beta", "beta.Z.att", "pval", "p.adj")
  }
  return(beta_pvals_sorted)
}

#=========================================================================#
#' detectionStats
#'
#' Given a vector functional (true) attribute names and a vector of positive
#' association attribute names, returns detection statistics like recall and precision.
#'
#' @param functional character vector of functional/true attribute names.
#' @param positives character vector of attribute names of positive associations (null hypothesis rejected or some threshold).
#' @return list with elements TP, FP, FN, TPR, FPR, precision, recall and summary message (string).  
#' @examples
#' functional <- case.control.3sets$signal.names  
#' positives <- row.names(npdr.cc.results.df[npdr.cc.results.df[,1]<.05,]) # p.adj<.05
#' npdr.cc.detect.stats <- detectionStats(functional.case.control, positives)
#' cat(npdr.cc.detect.stats$summary.msg)on(outcome="class", dataset=case.control.data)
#' @export
detectionStats <- function(functional, positives){
  TP <- sum(positives %in% functional)
  FP <- sum((positives %in% functional)==F)
  FN <- length(functional) - TP
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  num.positives <- length(positives)
  TPR <- TP/num.positives #rate, aka power or sensitivity
  FPR <- FP/num.positives #rate
  # summary message
  report <- paste(
    "Given ", length(functional)," functional (true) attributes.\n",
    "Given ", length(positives)," selected (positive) attributes.\n",
    "True Positives: ", TP," true out of ", length(positives)," positives. TP rate = ", TPR, ".\n",
    "False Positives: ", FP," false out of ", length(positives)," positives. FP rate = ", FPR, ".\n",
    "Precision: ", precision,".\n",
    "Recall: ", recall,".\n",
    sep="")
  return(list(TP=TP, FP=FP, FN=FN, TPR=TPR, FPR=FPR, 
              precision=precision, recall=recall, report=report))
}

#=========================================================================#
#' geneLowVarianceFilter
#'
#' Low variance mask and filtered data for gene expression matrix.
#'
#' @param dataMatrix data matrix with predictors only, sample x gene
#' @param pct percentile of low variance removed 
#' @return mask and filtered data  
#' @examples
#' pct <- 0.7 # higher value more strict
#' filter <- geneLowVarianceFilter(unfiltered.predictors.mat, pct)
#' filtered.data.df <- data.frame(filter$fdata, class = rnaseq.mdd.phenotype)
#' @export
geneLowVarianceFilter <- function(dataMatrix, percentile=0.5) {
  variances <- apply(as.matrix(dataMatrix), 2, var)
  threshold <- quantile(variances, c(percentile))
  # remove variable columns with lowest percentile variance
  mask <- apply(dataMatrix, 2, function(x) var(x) > threshold)
  fdata <- dataMatrix[, mask]
  # return the row mask and filtered data
  list(mask=mask, fdata=fdata)
}

#=========================================================================#
#' vwak
#' 
#' Variable-Wise Adaptive k
#' method for optimizing NPDR scores for each attribute as a function of k
#' Computes p x k beta and P value matrices for a data set with p attributes
#' 
#' @param dats m x (p+1) data set of m instances and p attributes with 1 binary outcome or m x [p(p - 1) + 1] with p(p-1) correlations and 1 outcome. Outcome is last column for standard m x (p + 1) and first column for m x [p(p - 1) + 1] (no good reason for the difference).
#' @param k.grid increasing sequence of k values used as looping index. Default is seq(1,(nrow(dats)-1),by=1).
#' @param verbose logical indicating whether to print progress with loop. Default is FALSE, but TRUE also does not give anything useful.
#' @param attr.diff.type character indicating the type of attribute diff to use. Default is 'numeric-abs' for standard continuous data. Use 'correlation-data' for rs-fMRI data.
#' @param corr.attr.names character indicating names of ROIs for attr.diff.type='correlation-data'. Default is NULL.
#' @param separate.hitmiss.nbds logical indicating whether to compute hit/miss neighborhoods separately. Default is FALSE.
#' @param label character indicating type of response. Default is "class" and should not change as of yet.
#' 
#' @return A list with:
#' \describe{
#'   \item{beta.mat}{p x k matrix of beta coefficients from NPDR}
#'   \item{pval.mat}{p x k matrix of P values corresponding to beta coefficients from NPDR}
#' }
#' 
#' @examples
#' # Main effect simulation for standard m x p data set
#' num.samples <- 100
#' num.variables <- 100
#' pct.imbalance <- 0.5
#' pct.signals <- 0.1
#' main.bias <- 0.5
#' sim.type <- "mainEffect"
#'
#' dataset <- createSimulation2(num.samples=num.samples,
#'                              num.variables=num.variables,
#'                              pct.imbalance=pct.imbalance,
#'                              pct.signals=pct.signals,
#'                              main.bias=main.bias,
#'                              label="class",
#'                              sim.type=sim.type,
#'                              pct.mixed=0.5,
#'                              pct.train=0.5,
#'                              pct.holdout=0.5,
#'                              pct.validation=0,
#'                              verbose=T,
#'                              data.type="continuous")
#' dats <- rbind(dataset$train, dataset$holdout, dataset$validation)
#' dats <- dats[order(dats[,ncol(dats)]),]
#'
#' # run Variable-Wise Adaptive k function
#' out <- vwak(dats=dats,
#'             k.grid=NULL,
#'             verbose=T,
#'             attr.diff.type="numeric-abs",
#'             label="class")
#' @export
#' 
vwak <- function(dats=NULL,
                 k.grid=NULL,
                 verbose=F,
                 attr.diff.type="numeric-abs",
                 corr.attr.names=NULL,
                 separate.hitmiss.nbds=FALSE,
                 label="class"){
  
  start_time <- Sys.time()
  
  m <- dim(dats)[1]
  if(attr.diff.type!="correlation-data"){
    p <- dim(dats[,-ncol(dats)])[2]
  }else{
    p <- dim(dats[,-1])[2]
  }
  
  if(attr.diff.type=="correlation-data"){   # corrdata
    mynum <- dim(dats[,-1])[2]               # corrdata
    for(i in seq(1,mynum-1,by=1)){          # corrdata
      mydiv <- i                            # corrdata
      if((mydiv*(mydiv - 1)) == mynum){     # corrdata
        my.dimension <- mydiv               # corrdata
        break                               # corrdata
      }                                     # corrdata
    }                                       # corrdata
    num.attr <- my.dimension                # corrdata
  }else{
    num.attr <- ncol(dats[,-ncol(dats)])
  }
  num.samp <- nrow(dats)
  
  if(attr.diff.type=="correlation-data"){    # corrdata
    attr.idx.list <- list()                  # corrdata
    for(i in 1:num.attr){                    # corrdata
      lo.idx <- (i - 1)*(num.attr-1) + 1     # corrdata
      hi.idx <- i*(num.attr-1)               # corrdata
      attr.idx.list[[i]] <- c(lo.idx:hi.idx) # corrdata
    }                                        # corrdata
  }
  
  if(is.null(k.grid)){
    #ks <- seq(1,floor((m-1)/2),by=1) # if computing hit/miss neighbors separately
    ks <- seq(1,(m-1),by=1) # if not separating hit/miss
  }else{
    ks <- k.grid
  }
  
  if(attr.diff.type=="correlation-data"){
    pheno.vec <- as.factor(dats[,1])
    dats <- as.matrix(dats[,-1])
  }
  
  avai.cors <- parallel::detectCores() - 2
  cl <- parallel::makeCluster(avai.cors)
  doParallel::registerDoParallel(cl)
  
  beta.mat <- matrix(0,nrow=num.attr,ncol=length(ks))
  pval.mat <- matrix(0,nrow=num.attr,ncol=length(ks))
  out.list <- foreach(k=ks, 
                      .export=c('npdr'), 
                      .packages=c("npdr")) %dopar% {
                        
                        #cat("k = ",k,"\n")
                        
                        if(attr.diff.type != "correlation-data"){
                          npdr.cc.results <- npdr::npdr(label, dats, regression.type="binomial", attr.diff.type="numeric-abs",
                                                        nbd.method="relieff", nbd.metric = "manhattan", msurf.sd.frac=.5, knn=k,
                                                        neighbor.sampling="none", separate.hitmiss.nbds=separate.hitmiss.nbds,
                                                        padj.method="bonferroni", verbose=verbose)
                        }else{
                          #pheno.vec <- dats[,1]
                          #dats <- dats[,-1]
                          npdr.cc.results <- npdr::npdr(pheno.vec, dats, regression.type="binomial", attr.diff.type="correlation-data",
                                                        nbd.method="relieff", nbd.metric = "manhattan", msurf.sd.frac=.5, knn=k,
                                                        neighbor.sampling="none", separate.hitmiss.nbds=separate.hitmiss.nbds,
                                                        padj.method="bonferroni", verbose=verbose,
                                                        corr.attr.names=corr.attr.names)
                          
                        }
                        
                        out.mat <- data.frame(cbind(att=npdr.cc.results$att,
                                                    beta=npdr.cc.results$beta.Z.att,
                                                    pval=npdr.cc.results$pval.adj))
                        out.mat <- out.mat[order(as.character(out.mat[,1])),]
                        
                        beta.mat <- as.numeric(as.character(out.mat[,"beta"]))
                        pval.mat <- as.numeric(as.character(out.mat[,"pval"]))
                        att.mat <- as.character(out.mat[,1])
                        
                        list(betas=beta.mat,
                             pvals=pval.mat,
                             atts=att.mat)
                        
                      }
  parallel::stopCluster(cl)
  
  avai.cors <- parallel::detectCores() - 2
  cl <- parallel::makeCluster(avai.cors)
  doParallel::registerDoParallel(cl)
  
  beta.mat <- NULL
  out.betas <- foreach(k=1:length(ks),.combine='cbind') %dopar% {
    
    out.list[[k]]$betas
    
  }
  parallel::stopCluster(cl)
  
  avai.cors <- parallel::detectCores() - 2
  cl <- parallel::makeCluster(avai.cors)
  doParallel::registerDoParallel(cl)
  
  pval.mat <- NULL
  out.pvals <- foreach(k=1:length(ks),.combine='cbind') %dopar% {
    
    out.list[[k]]$pvals
    
  }
  parallel::stopCluster(cl)
  
  beta.mat <- as.data.frame(out.betas)
  pval.mat <- as.data.frame(out.pvals)
  row.names(beta.mat) <- as.character(out.list[[1]]$atts)
  row.names(pval.mat) <- as.character(out.list[[1]]$atts)
  colnames(beta.mat) <- paste("k.",ks,sep="")
  colnames(pval.mat) <- paste("k.",ks,sep="")
  
  end_time <- Sys.time()
  elapsed <- end_time - start_time
  cat("Elapsed: ",elapsed,"\n")
  
  list(beta.mat=beta.mat,pval.mat=pval.mat)
  
}

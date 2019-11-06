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

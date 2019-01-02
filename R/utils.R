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
#' univariateRegression
#'
#' Univariate logistic or linear regression for a dataset.
#'
#' @param outcome string with name of class column.
#' @param dataset data matrix with predictor columns and outcome column.
#' @return matrix of beta, p-value and adjusted p-value, sorted by p-value.  
#' @examples
#' lr.results <- univariateRegression(outcome="class", dataset=case.control.data, regression.type="glm")
#  lr.results[lr.results[,"p.adj"]<.05] 
#' @export
univariateRegression <- function(outcome, dataset, regression.type="lm", covars="none"){
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
    if (length(covars>1)){
      model.func <- function(x) {as.numeric(tidy(lm(pheno.vec ~ attr.mat[,x] + covars))[2,4:5])}
    } else { # covar=="none"
      model.func <- function(x) {as.numeric(tidy(lm(pheno.vec ~ attr.mat[,x]))[2,4:5])}  
    } 
    } else { # "glm"
    if (length(covars>1)){
      #model.func <- function(x) {tidy(glm(pheno.vec ~ attr.mat[,x] + covars, family=binomial))[2,4:5]}
      model.func <- function(x) {glm(pheno.vec ~ attr.mat[,x] + covars, family=binomial)[2,4:5]}
    } else { # covar=="none"
      #model.func <- function(x) {tidy(glm(pheno.vec ~ attr.mat[,x], family=binomial))[2,4:5]}
      model.func <- function(x) {glm(pheno.vec ~ attr.mat[,x], family=binomial)[2,4:5]}
    }
  } # end else glm
  #class.col <- which(colnames(dataset)==outcome)
  #predictor.cols <- which(colnames(dataset)!=outcome)
  beta_pvals <- t(sapply(1:ncol(attr.mat), model.func)) # stats for all predictors
  univariate.padj <- p.adjust(beta_pvals[,2]) # fdr
  univariate.padj <- as.numeric(format(univariate.padj, scientific = T, digits=5))
  betas <- as.numeric(format(beta_pvals[,1], scientific = F, digits=5))
  pvals <- as.numeric(format(beta_pvals[,2], scientific = T, digits=5))
  beta_pvals <- cbind(betas,pvals,univariate.padj) # adjusted p-val column
  row.names(beta_pvals)<- colnames(attr.mat) # add predictor names
  #row.names(beta_pvals)<- colnames(dataset)[-class.col] # add predictor names
  beta_pvals_sorted <- beta_pvals[order(as.numeric(beta_pvals[,2]), decreasing = F),] # sort by pval
  colnames(beta_pvals_sorted) <- c("beta", "pval", "p.adj")
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
#' positives <- row.names(glm.stir.results.df[glm.stir.results.df[,1]<.05,]) # p.adj<.05
#' glm.stir.detect.stats <- detectionStats(functional.case.control, glm.stir.positives)
#' cat(glm.stir.detect.stats$summary.msg)on(outcome="class", dataset=case.control.data)
#' @export
detectionStats <- function(functional, positives){
  TP <- sum(positives %in% functional)
  FP <- sum((positives %in% functional)==F)
  FN <- length(functional) - TP
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  num.positives <- length(positives)
  TPR <- TP/num.positives #rate
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

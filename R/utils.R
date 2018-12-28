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
univariateRegression <- function(outcome, dataset, regression.type="lm"){
  if (regression.type=="lm"){
    model.func <- function(x) {tidy(lm(dataset[,outcome] ~ dataset[,x]))[2,4:5]}
  } else { # "glm"
    model.func <- function(x) {tidy(glm(dataset[,outcome] ~ dataset[,x], family=binomial))[2,4:5]}
  }
  class.col <- which(colnames(dataset)==outcome)
  predictor.cols <- which(colnames(dataset)!=outcome)
  beta_pvals <- t(sapply(predictor.cols, model.func)) # stats for all predictors
  row.names(beta_pvals)<- colnames(dataset)[-class.col] # add predictor names
  univariate.padj <- p.adjust(beta_pvals[,2]) # fdr
  beta_pvals <- cbind(beta_pvals,univariate.padj) # adjusted p-val column
  beta_pvals_sorted <- beta_pvals[order(as.numeric(beta_pvals[,2]), decreasing = F),] # sort by pval
  colnames(beta_pvals_sorted) <- c("beta", "pval", "p.adj")
  return(as.matrix(beta_pvals_sorted))
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
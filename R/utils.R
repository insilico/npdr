#=========================================================================#
#' check.packages
#'
#' Utility to check for package installation
#'
#' @param pkg character vector of package names to install if not already.
#' @return null There is nothing returned.  
#' @examples
#' packages <- c("ggplot2", "CORElearn", "reshape2", "dplyr", "pROC", "plotROC")
#' check.packages(packages)  
#' @export
check.packages <- function(pkg){
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
#' univariateLogisticRegression
#'
#' Univariate logistic regression for a data set.
#'
#' @param outcome string with name of class column.
#' @param dataset data matrix with predictor columns and outcome column.
#' @return matrix of beta, p-value and adjusted p-value, sorted by p-value.  
#' @examples
#' lr.results <- univariateLogisticRegression(outcome="class", dataset=case.control.data)
#  lr.results[lr.results[,"p.adj"]<.05] 
#' @export
univariateLogisticRegression <- function(outcome, dataset){
  glm.func2 <- function(x) {tidy(glm(dataset[,outcome] ~ dataset[,x], family=binomial))[2,4:5]}
  class.col <- which(colnames(dataset)==outcome)
  predictor.cols <- which(colnames(dataset)!=outcome)
  beta_pvals <- t(sapply(predictor.cols, glm.func2)) # stats for all predictors
  row.names(beta_pvals)<- colnames(dataset)[-class.col] # add predictor names
  univariate.padj <- p.adjust(beta_pvals[,2]) # fdr
  beta_pvals <- cbind(beta_pvals,univariate.padj) # adjusted p-val column
  beta_pvals_sorted <- beta_pvals[order(as.numeric(beta_pvals[,2]), decreasing = F),] # sort by pval
  colnames(beta_pvals_sorted) <- c("beta", "pval", "p.adj")
  return(as.matrix(beta_pvals_sorted))
}
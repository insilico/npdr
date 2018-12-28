#=========================================================================#
#' diffRegression
#'
#' Wrapper for lm and glm to run regression for a phenotype diff vector and one attribute diff vector. 
#' Organize regression statistics into a vector.
#'
#' @param pheno.diffs outcome variable as vector of diffs 
#' @param predictor.diffs one predictor varialbe as vector of diffs 
#' @param regression.type (\code{"lm"}) todo: \code{"glm"} and covariates
#' @return vector or regression stats to put into list for reSTIR and combine into matrix
#'
#' @examples
#'
#' @export
diffRegression <- function(pheno.diffs, predictor.diffs, regression.type="lm") {
  # 
  if (regression.type=="lm"){
  fit <- summary(lm(pheno.diffs ~ predictor.diffs))
  stats.vec <- c(
    fit$coefficients[2,4], # p-value for attribute beta, pval.a
    fit$coefficients[2,3], # beta_hat_a, standardize beta for attribute, Ba
    fit$r.squared,         # R^2 of fit, R.sqr
    fit$fstatistic[1],     # F-stat and next is its p-value, F.stat
    1 - pf(fit$fstatistic[1], fit$fstatistic[2], fit$fstatistic[3]), # Fstat.pval
    fit$coefficients[1,3], # beta_hat_0, intercept, B0
    fit$coefficients[1,4] # p-value for intercept, B0.pval
  ) 
  } else{ #regression.type=="glm"
    fit <- summary(glm(pheno.diffs ~ predictor.diffs, family=binomial(link=logit)))
    stats.vec <- c(
      fit$coefficients[2,4], # p-value for attribute beta, pval.a
      fit$coefficients[2,3], # beta_hat_a, standardize beta for attribute, Ba
      fit$coefficients[1,3], # beta_hat_0, intercept, B0
      fit$coefficients[1,4]  # p-value for intercept, B0.pval
    )
  }
  return(stats.vec)
}

#=========================================================================#
#' glmSTIR
#'
#' generalized linear model (GLM) based STatistical Inference Relief (STIR)
#' Computes Relief-based attribute statistical signficance for case/control or quantitative outcomes.
#' Allows for categorical (SNP) or numeric (expession) predictor data types. 
#' Allows for covariate correction.   
#'
#' @param outcome character name or length-m numeric outcome vector for linear regression, factor for logistic regression 
#' @param data.set m x p matrix of m instances and p attributes, May also include outcome vector but then outcome should be name. Include attr names as colnames. 
#' @param covar vector of names or indices of covariates in data.set for covariate correction. Or separate data matrix of covariates. 
#' @param regression.type (\code{"lm"} or \code{"glm"})
#' @param nbd.metric used in stirDistances for distance matrix between instances, default: \code{"manhattan"} (numeric). Used by nearestNeighbors().
#' @param nbd.method neighborhood method [\code{"multisurf"} or \code{"surf"} (no k) or \code{"relieff"} (specify k)]. Used by nearestNeighbors().
#' @param neighbor.pairs.idx nearest hit/miss matrix (2 columns), output from \code{find.neighbors}
#' @param attr.diff.type diff type for attributes (\code{"manhattan"} or \code{"euclidean"} for numeric). Phenotype diff same as attr.diff.type when lm regression. For glm, phenotype diff is match/mismatch. 
#' @param covar.diff.type diff type for covariates (\code{"manhattan"} or \code{"euclidean"} for numeric). 
#' @param k number of constant nearest hits/misses for \code{"relieff"} (fixed-k). Used by nearestNeighbors().
#' The default k=0 means use the expected SURF theoretical k with sd.frac (.5 by default) 
#' @param sd.frac multiplier of the standard deviation from the mean distances; subtracted from mean for SURF or multiSURF.
#' The multiSURF default is sd.frac=0.5: mean - sd/2. Used by nearestNeighbors(). 
#' @param fdr.method for p.adjust (\code{"fdr"}, \code{"bonferroni"}, ...) 
#' @return glmSTIR.stats.df: glmSTIR regression coefficients and p-values for each attribute
#'
#' @examples
#' # Data interface options.
#' # Specify name ("qtrait") of outcome and data.set, which is a data frame including the outcome column.
#' # ReliefF fixed-k neighborhood, uses surf theoretical default (with sd.frac=.5) if you do not specify k or let k=0
#' restir.results.df <- glmSTIR("qtrait", train.data, regression.type="lm", nbd.method="relieff", nbd.metric = "manhattan", attr.diff.type="manhattan", covar.diff.type="manhattan", sd.frac=0.5, fdr.method="bonferroni")
#'
#' # Specify column index (101) of outcome and data.set, which is a data frame including the outcome column.
#  # ReliefF fixed-k nbd, choose a k (k=10). Or choose sd.frac
#' restir.results.df <- glmSTIR(101, train.data, regression.type="lm", nbd.method="relieff", nbd.metric = "manhattan", attr.diff.type="manhattan", covar.diff.type="manhattan", k=10, fdr.method="bonferroni")
#'
#' # if outcome vector (pheno.vec) is separate from attribute matrix
#' # multisurf
#' restir.results.df <- glmSTIR(pheno.vec, predictors.mat, regression.type="lm", nbd.method="multisurf", nbd.metric = "manhattan", attr.diff.type="manhattan", covar.diff.type="manhattan", sd.frac=0.5, fdr.method="bonferroni")
#' restir.positives <- row.names(restir.results.df[restir.results.df[,1]<.05,]) # reSTIR p.adj<.05
#' @export
glmSTIR <- function(outcome, data.set, covars="none", regression.type="lm", nbd.method="multisurf", nbd.metric = "manhattan", 
                    attr.diff.type="numeric-abs", covar.diff.type="numeric-abs", 
                    k=0, sd.frac=0.5, fdr.method="bonferroni"){
  ##### parse the commandline 
  #if (is.character(outcome)){
  if (length(outcome)==1){
    # e.g., outcome="qtrait" or outcome=101 (pheno col index) and data.set is data.frame including outcome variable
    pheno.vec <- data.set[,outcome] # get phenotype
    if (is.character(outcome)){ # example column name: outcome="qtrait"
      attr.mat <- data.set[ , !(names(data.set) %in% outcome)]  # drop the outcome/phenotype
    } else { # example column index: outcome=101
      attr.mat <- data.set[ , -outcome]  # drop the outcome/phenotype  
    }
  } else { # user specifies a separate phenotype vector
    pheno.vec <- outcome # assume users provides a separate outcome data vector
    attr.mat <- data.set # assumes data.set only contains attributes/predictors
  }
  rm(data.set)  # cleanup memory
  
  num.attr <- ncol(attr.mat)
  num.samp <- nrow(attr.mat)
  
  ##### get Neighbors (no phenotype used)
  # nbd.method (relieff, multisurf...), nbd.metric (manhattan...), k (for relieff nbd, theoerical surf default) 
  # sd.frac used by surf/multisurf relieff for theoretical k
  neighbor.pairs.idx <- nearestNeighbors(attr.mat, nbd.method=nbd.method, nbd.metric = nbd.metric, 
                                         sd.vec = NULL, sd.frac = sd.frac, k=k)
  
  ##### run glmSTIR, each attribute is a list, then we do.call rbind to a matrix
  glmSTIR.stats.list <- vector("list",num.samp)
  for (attr.idx in seq(1, num.attr)){
    attr.vals <- attr.mat[, attr.idx]
    Ri.attr.vals <- attr.vals[neighbor.pairs.idx[,1]]
    NN.attr.vals <- attr.vals[neighbor.pairs.idx[,2]]
    attr.diff.vec <- stirDiff(Ri.attr.vals, NN.attr.vals, diff.type=attr.diff.type)
    # create pheno diff vector for linear regression (numeric)  
    if (regression.type=="lm"){
    Ri.pheno.vals <- pheno.vec[neighbor.pairs.idx[,1]]
    NN.pheno.vals <- pheno.vec[neighbor.pairs.idx[,2]]
    pheno.diff.vec <- stirDiff(Ri.pheno.vals, NN.pheno.vals, diff.type="numeric-abs")
    } else { #regression.type=="glm"
      # create pheno diff vector for logistic regression (match-mismatch or hit-miss)  
      Ri.pheno.vals <- pheno.vec[neighbor.pairs.idx[,1]]
      NN.pheno.vals <- pheno.vec[neighbor.pairs.idx[,2]]
      pheno.diff.vec <- stirDiff(Ri.pheno.vals, NN.pheno.vals, diff.type="match-mismatch")
      pheno.diff.vec <- as.factor(ifelse(pheno.diff.vec=="TRUE",1,0))
    }
    # utility function: RUN regression
    glmSTIR.stats.list[[attr.idx]] <- diffRegression(pheno.diff.vec, attr.diff.vec, regression.type=regression.type) 
  }
  # combine lists into matrix
  glmSTIR.stats.attr_ordered.mat <- do.call(rbind, glmSTIR.stats.list)
  
  # rownames
  if (!is.null(colnames(attr.mat))){
    # add attribute names to stats/results matrix if the data matrix contains them
    rownames(glmSTIR.stats.attr_ordered.mat) <- colnames(attr.mat)
  } else {
    message("If you have attribute names, add them to colnames of input data.")
  }
  
  # attribute p-values
  attr.pvals <- glmSTIR.stats.attr_ordered.mat[, 1]
  # order-index for sorted attribute-beta p-values
  attr.pvals.order.idx <- order(attr.pvals, decreasing = F)
  # adjust p-values using Benjamini-Hochberg (default)
  attr.pvals.adj <- p.adjust(attr.pvals[attr.pvals.order.idx], method=fdr.method)
  
  # order by attribute p-value
  glmSTIR.stats.pval_ordered.mat <- glmSTIR.stats.attr_ordered.mat[attr.pvals.order.idx, ]
  # prepend adjused attribute p-values to first column
  glmSTIR.stats.pval_ordered.mat <- cbind(attr.pvals.adj, glmSTIR.stats.pval_ordered.mat)
  if (regression.type=="lm"){# different stats colnames for lm and glm
    colnames(glmSTIR.stats.pval_ordered.mat) <- c("pval.adj", "pval.attr", "beta.attr", "R.sqr", 
                                                  "F.stat", "Fstat.pval", "beta.0", "pval.0")
  } else{ # "glm"
    colnames(glmSTIR.stats.pval_ordered.mat) <- c("pval.adj", "pval.attr", "beta.attr", "beta.0", "pval.0")
  }
  # dataframe it
  glmSTIR.stats.df <- data.frame(glmSTIR.stats.pval_ordered.mat)
  return(glmSTIR.stats.df)
}

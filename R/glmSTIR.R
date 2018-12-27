#=========================================================================#
#' stirDiff
#'
#' A diff is a function that computes the diffrence of values for an attribute between two instances.
#' It is used for attribute selection for attribute diffs and phenotype diffs.  
#' This function is vectorized: input a and b can be two vectors of values for one attribute.   
#'
#' @param a value of attribute for first instance
#' @param b value of attribute for second instance
#' @param type diff rule for the given attribute data type, such as numeric or categorical.
#' @return val diff or vector of diffs 
#' @examples
#' Example
#' @export
stirDiff <- function(a, b, diff.type = "numeric-abs", norm.fac = 1){
  # compute the difference between two vectors elementwise
  if (diff.type=="numeric-sqr"){ # numeric squared difference
    val <- abs(a - b)^2/norm.fac
  } else if (diff.type=="allele-sharing"){ # snps
    val <- abs(a-b)/2
  } else if (diff.type=="match-mismatch"){ # snps
    val <- as.character(a==b) # convert this to factor in reSTIR
  } else{ # numeric abs difference
    val <- abs(a - b)/norm.fac
  }
  return(val)
}

#=========================================================================#
#' stirDistances
#'
#' Create m x m distance matrix from m instances and p attributes using different metrics. Used by nearestNeighbors(). 
#' Note: Probably best to standardize data before manhattan and euclidean.
#'
#' @param attr.mat m x p matrix of m instances and p attributes 
#' @param metric for distance matrix between instances (default: \code{"manhattan"}, others include \code{"euclidean"}, 
#' versions scaled by max-min, \code{"relief-scaled-manhattan"} and \code{"relief-scaled-euclidean"}, and for GWAS \code{"allele-sharing-manhattan"}).
#' @return  distancesmat, matrix of m x m (instances x intances) pairwise distances.
#' @examples
#' dist.mat <- stirDistances(predictors.mat, metric = "manhattan")
#' @export
stirDistances <- function(attr.mat, metric="manhattan"){
  # Compute distance matrix between all samples (rows)
  # reSTIR default is numeric manhattan ("manhattan"), max-min scaling is not needed for stir
  if (metric == "hamming"){
    distance.mat <- hamming.binary(attr.mat)
  } else if (metric == "allele-sharing-manhattan"){
    # allele-sharing-manhattan, AM for SNPs
    attr.mat.scale <- attr.mat / 2
    distance.mat <- as.matrix(dist(attr.mat.scale, method = "manhattan"))
  } else if (metric == "relief-scaled-manhattan"){
    # value of metric, euclidean, manhattan or maximum
    maxminVec <- attr.range(attr.mat)
    minVec <- apply(attr.mat, 2, function(x) {min(x)})
    attr.mat.centered <- t(attr.mat) - minVec
    attr.mat.scale <- t(attr.mat.centered / maxminVec)
    distance.mat <- as.matrix(dist(attr.mat.scale, method = "manhattan"))
  } else if (metric == "relief-scaled-euclidean"){
    # value of metric, euclidean, manhattan or maximum
    maxminVec <- attr.range(attr.mat)
    minVec <- apply(attr.mat, 2, function(x) {min(x)})
    attr.mat.centered <- t(attr.mat) - minVec
    attr.mat.scale <- t(attr.mat.centered / maxminVec)
    distance.mat <- as.matrix(dist(attr.mat.scale, method = "euclidean"))
  } else if (metric=="euclidean"){
    distance.mat <- as.matrix(dist(attr.mat, method = "euclidean"))
  } else {
    distance.mat <- as.matrix(dist(attr.mat, method = "manhattan"))
  }
  distance.mat
}

#=========================================================================#
#' nearestNeighbors
#'
#' Find nearest neighbors of each instance using relief.method
#' Used for regression stir (reSTIR) (no hits or misses specified in function).
#' Also used in reSTIR for case/control, but hit/miss is used in reSTIR function. 
#'
#' @param attr.mat m x p matrix of m instances and p attributes 
#' @param nbd.metric used in stirDistances for distance matrix between instances, default: \code{"manhattan"} (numeric)
#' @param nbd.method neighborhood method [\code{"multisurf"} or \code{"surf"} (no k) or \code{"relieff"} (specify k)]
#' @param k number of constant nearest hits/misses for \code{"relieff"} (fixed k). 
#' The default k=0 means use the expected SURF theoretical k with sd.frac (.5 by default) for relieff nbd.
#' @param sd.frac multiplier of the standard deviation from the mean distances, subtracted from mean distance to create for SURF or multiSURF radius.
#' The multiSURF default "dead-band radius" is sd.frac=0.5: mean - sd/2 
#' @return  Ri_NN.idxmat, matrix of Ri's (first column) and their NN's (second column)
#'
#' @examples
#' # multisurf neighborhood with sigma/2 (sd.frac=0.5) "dead-band" boundary
#' neighbor.pairs.idx <- nearestNeighbors(predictors.mat, nbd.method="multisurf", nbd.metric = "manhattan", sd.frac = 0.5)
#' # reliefF (fixed-k) neighborhood using default k equal to theoretical surf expected value
#' # One can change the theoretical value by changing sd.frac (default 0.5)
#' neighbor.pairs.idx <- nearestNeighbors(predictors.mat, nbd.method="relieff", nbd.metric = "manhattan")
#' # reliefF (fixed-k) neighborhood with a user-specified k
#' neighbor.pairs.idx <- nearestNeighbors(predictors.mat, nbd.method="relieff", nbd.metric = "manhattan", k=10)
#'
#' @export
nearestNeighbors <- function(attr.mat, nbd.method="multisurf", nbd.metric = "manhattan", sd.vec = NULL, sd.frac = 0.5, k=0){
  # create a matrix with num.samp rows, two columns
  # first column is sample Ri, second is Ri's nearest neighbors
  
  dist.mat <- stirDistances(attr.mat, metric = nbd.metric)
  num.samp <- nrow(attr.mat)
  
  if (nbd.method == "relieff"){  
    if (k==0){ # if no k specified or value 0
    # replace k with the theoretical expected value for SURF (close to multiSURF)
    erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
    # theoretical surf k (sd.frac=.5) for regression problems (does not depend on a hit/miss group)
    k <- floor((num.samp-1)*(1-erf(sd.frac/sqrt(2)))/2)  # uses sd.frac
    }
    Ri_NN.idxmat <- matrix(0, nrow = num.samp * k, ncol = 2)
    colnames(Ri_NN.idxmat) <- c("Ri_idx","NN_idx")
    for (Ri in seq(1:num.samp)){ # for each sample Ri
      Ri.distances <- dist.mat[Ri,] # all distances to sample Ri
      Ri.nearest <- order(Ri.distances, decreasing = F) # closest to furthest
      ## bam_add
      Ri.nearest.idx <- Ri.nearest[2:(k+1)] # skip Ri self
      # stack matrix of neighbor indices
      row.start <- (Ri-1)*k + 1
      row.end <- row.start + k - 1
      Ri_NN.idxmat[row.start:row.end, 1] <- rep(Ri, k)     # col of Ri's
      Ri_NN.idxmat[row.start:row.end, 2] <- Ri.nearest.idx  # col of knn's of Ri's
    }
  } else {
    if (nbd.method == "surf"){
      num.pair <- num.samp * (num.samp-1) / 2 # number of paired distances
      radius.surf <- sum(dist.mat)/(2*num.pair) # const r = mean(all distances)
      sd.const <- sd(dist.mat[upper.tri(dist.mat)])  
      # bam: orignal surf does not subtract sd-frac but should for fair multisurf comparison
      Ri.radius <- rep(radius.surf - sd.frac*sd.const, num.samp) 
    }
    if (nbd.method == "multisurf"){
      if (is.null(sd.vec)) sd.vec <- sapply(1:num.samp, function(x) sd(dist.mat[-x, x]))
      Ri.radius <- colSums(dist.mat)/(num.samp - 1) - sd.frac*sd.vec # use adaptive radius
    }
    
    # put each Ri's nbd in a list then rbind them at the end with do.call(rbind, List)
    # initialize list:
    Ri.nearestPairs.list <- vector("list",num.samp)
    for (Ri in seq(1:num.samp)){ # for each sample Ri
      Ri.distances <- sort(dist.mat[Ri,], decreasing = F)
      Ri.nearest <- Ri.distances[Ri.distances < Ri.radius[Ri]] # within the threshold
      Ri.nearest <- Ri.nearest[-1] # skip Ri self
      Ri.nearest.idx <- match(names(Ri.nearest), row.names(attr.mat))
      if (length(Ri.nearest.idx) > 1){ # if neighborhood not empty
        # cbind automatically repeats Ri
        Ri.nearestPairs.list[[Ri]] <- cbind(Ri, Ri.nearest.idx) 
      } 
    } # end for, now stack lists into matrix, do.call rbind
    
    Ri_NN.idxmat <- do.call(rbind, Ri.nearestPairs.list)
    colnames(Ri_NN.idxmat) <- c("Ri_idx","NN_idx")
  }
  # matrix of Ri's (first column) and their NN's (second column)
  return(Ri_NN.idxmat)
}

#=========================================================================#
#' diffRegression
#'
#' Utility to run regression for a phenotype diff vector and one attribute diff vector. 
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
#' reSTIR
#'
#' regression-based STIR statistics for quantitative or case-control phenotypes
#' 
#'
#' @param outcome character name or length-m numeric outcome vector for linear regression, factor for logistic regression 
#' @param data.set m x p matrix of m instances and p attributes, May also include outcome vector but then outcome should be name. Include attr names as colnames. 
#' @param regression.type (\code{"lm"} or \code{"glm"})
#' @param nbd.metric used in stirDistances for distance matrix between instances, default: \code{"manhattan"} (numeric). Used by nearestNeighbors().
#' @param nbd.method neighborhood method [\code{"multisurf"} or \code{"surf"} (no k) or \code{"relieff"} (specify k)]. Used by nearestNeighbors().
#' @param neighbor.pairs.idx nearest hit/miss matrices, output from \code{find.neighbors}
#' @param attr.diff.type diff type for attributes (\code{"manhattan"} or \code{"euclidean"} for numeric). Phenotype diff same as attr.diff.type when lm regression. For glm, phenotype diff is match/mismatch. 
#' @param covar.diff.type diff type for covariates (\code{"manhattan"} or \code{"euclidean"} for numeric). 
#' @param k number of constant nearest hits/misses for \code{"relieff"} (fixed-k). Used by nearestNeighbors().
#' The default k=0 means use the expected SURF theoretical k with sd.frac (.5 by default) 
#' @param sd.frac multiplier of the standard deviation from the mean distances; subtracted from mean for SURF or multiSURF.
#' The multiSURF default is sd.frac=0.5: mean - sd/2. Used by nearestNeighbors(). 
#' @param fdr.method for p.adjust (\code{"fdr"}, \code{"bonferroni"}, ...) 
#' @return reSTIR.stats.df: reSTIR regression coefficients and p-values for each attribute
#'
#' @examples
#' # Specifies name ("qtrait") of outcome and data.set, which is a data frame including the outcome
#  # ReliefF fixed-k nbd, choose a k (k=10). Or choose sd.frac
#' restir.results.df <- reSTIR("qtrait", train.data, regression.type="lm", nbd.method="relieff", nbd.metric = "manhattan", attr.diff.type="manhattan", covar.diff.type="manhattan", k=10, fdr.method="bonferroni")

#' # Specifies index (101) of outcome and data.set, which is a data frame including the outcome
#' # ReliefF fixed-k neighborhood, surf theoretical default (sd.frac=.5) by not specifying k or let k=0
#' restir.results.df <- reSTIR(101, train.data, regression.type="lm", nbd.method="relieff", nbd.metric = "manhattan", attr.diff.type="manhattan", covar.diff.type="manhattan", sd.frac=0.5, fdr.method="bonferroni")

#' # 101 is column index of outcome variable
#' restir.results.df <- reSTIR(101, train.data, regression.type="lm", neighbor.pairs.idx, attr.diff.type="manhattan", covar.diff.type="manhattan", fdr.method="bonferroni")
#' row.names(restir.results.df[restir.results.df[,1]<.05,]) # reSTIR p.adj<.05
#'
#' # separate outcome vector and attribute matrix
#' # multisurf
#' restir.results.df <- reSTIR(pheno.vec, predictors.mat, regression.type="lm", nbd.method="multisurf", nbd.metric = "manhattan", attr.diff.type="manhattan", covar.diff.type="manhattan", sd.frac=0.5, fdr.method="bonferroni")
#' restir.positives <- row.names(restir.results.df[restir.results.df[,1]<.05,]) # reSTIR p.adj<.05
#' @export
reSTIR <- function(outcome, data.set, regression.type="lm", nbd.method="multisurf", nbd.metric = "manhattan", attr.diff.type="numeric-abs", covar.diff.type="numeric-abs", k=0, sd.frac=0.5, fdr.method="bonferroni"){

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
  neighbor.pairs.idx <- nearestNeighbors(attr.mat, nbd.method=nbd.method, nbd.metric = nbd.metric, sd.vec = NULL, sd.frac = sd.frac, k=k)
  
  ##### run reSTIR, each attribute is a list, then we do.call rbind to a matrix
  reSTIR.stats.list <- vector("list",num.samp)
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
    reSTIR.stats.list[[attr.idx]] <- diffRegression(pheno.diff.vec, attr.diff.vec, regression.type=regression.type) 
  }
  # combine lists into matrix
  reSTIR.stats.attr_ordered.mat <- do.call(rbind, reSTIR.stats.list)
  
  # rownames
  if (!is.null(colnames(attr.mat))){
    # add attribute names to stats/results matrix if the data matrix contains them
    rownames(reSTIR.stats.attr_ordered.mat) <- colnames(attr.mat)
  } else {
    message("If you have attribute names, add them to colnames of input data.")
  }
  
  # attribute p-values
  attr.pvals <- reSTIR.stats.attr_ordered.mat[, 1]
  # order-index for sorted attribute-beta p-values
  attr.pvals.order.idx <- order(attr.pvals, decreasing = F)
  # adjust p-values using Benjamini-Hochberg (default)
  attr.pvals.adj <- p.adjust(attr.pvals[attr.pvals.order.idx], method=fdr.method)
  
  # order by attribute p-value
  reSTIR.stats.pval_ordered.mat <- reSTIR.stats.attr_ordered.mat[attr.pvals.order.idx, ]
  # prepend adjused attribute p-values to first column
  reSTIR.stats.pval_ordered.mat <- cbind(attr.pvals.adj,reSTIR.stats.pval_ordered.mat)
  if (regression.type=="lm"){# different stats colnames for lm and glm
    colnames(reSTIR.stats.pval_ordered.mat) <- c("pval.adj", "pval.attr", "beta.attr", "R.sqr", "F.stat", "Fstat.pval", "beta.0", "pval.0")
  } else{ # "glm"
    colnames(reSTIR.stats.pval_ordered.mat) <- c("pval.adj", "pval.attr", "beta.attr", "beta.0", "pval.0")
  }
  # dataframe it
  reSTIR.stats.df <- data.frame(reSTIR.stats.pval_ordered.mat)
  return(reSTIR.stats.df)
}

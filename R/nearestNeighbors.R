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
  } else if (diff.type=="match-mismatch"){ 
    # used automatically for case-control pheno, optional genotype mismatch diff for snps
    val <- as.character(a==b) # convert this to factor in glmSTIR
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
#' @param sd.frac multiplier of the standard deviation from the mean distances, subtracted from mean distance to create for SURF or multiSURF radius. The multiSURF default "dead-band radius" is sd.frac=0.5: mean - sd/2 
#' @param k number of constant nearest hits/misses for \code{"relieff"} (fixed k). 
#' The default k=0 means use the expected SURF theoretical k with sd.frac (.5 by default) for relieff nbd.
#' @param rm.attr.from.dist attributes for removal (possible confounders) from the distance matrix calculation. 
#' 
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
#' 
#' testparams = get("testparams", envir=environment(foo)
nearestNeighbors <- function(attr.mat, 
                             nbd.method="multisurf", 
                             nbd.metric = "manhattan", 
                             sd.vec = NULL, sd.frac = 0.5, k=0,
                             rm.attr.from.dist=c()){
  # create a matrix with num.samp rows and two columns
  # first column is sample Ri, second is Ri's nearest neighbors
  
  if (length(rm.attr.from.dist)>1){ # remove attributes (possible confounders) from distance matrix calculation 
    attr.mat <- attr.mat[,-which(colnames(attr.mat) %in% rm.attr.from.dist)]
  }

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

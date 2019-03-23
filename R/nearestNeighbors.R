#=========================================================================#
#' npdrDiff
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
npdrDiff <- function(a, b, diff.type = "numeric-abs", norm.fac = 1){
  # compute the difference between two vectors elementwise
  if (diff.type=="numeric-sqr"){ # numeric squared difference
    val <- abs(a - b)^2/norm.fac
  } else if (diff.type=="allele-sharing"){ # snps
    val <- abs(a-b)/2
  } else if (diff.type=="match-mismatch"){ 
    # used automatically for case-control pheno, optional genotype mismatch diff for snps
    val <- ifelse(a==b,0,1)  # hit pairs are 0 and miss pairs are 1
    #val <- as.character(a==b) # convert this to factor in glmSTIR
  } else { # numeric abs difference
    val <- abs(a - b)/norm.fac
  }
  return(val)
}

#=========================================================================#
#' npdrDistances
#'
#' Create m x m distance matrix from m instances and p attributes using different metrics. Used by nearestNeighbors(). 
#' Note: Probably best to standardize data before manhattan and euclidean.
#'
#' @param attr.mat m x p matrix of m instances and p attributes 
#' @param metric for distance matrix between instances (default: \code{"manhattan"}, others include \code{"euclidean"}, 
#' versions scaled by max-min, \code{"relief-scaled-manhattan"} and \code{"relief-scaled-euclidean"}, and for GWAS \code{"allele-sharing-manhattan"}).
#' @return  distancesmat, matrix of m x m (instances x intances) pairwise distances.
#' @examples
#' dist.mat <- npdrDistances(predictors.mat, metric = "manhattan")
#' @export
npdrDistances <- function(attr.mat, metric="manhattan"){
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
    minVec <- apply(attr.mat, 2, min)
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
#' Used for npdr (no hits or misses specified in neighbor function).
#'
#' @param attr.mat m x p matrix of m instances and p attributes 
#' @param nb.metric used in npdrDistances for distance matrix between instances, default: \code{"manhattan"} (numeric)
#' @param nb.method neighborhood method [\code{"multisurf"} or \code{"surf"} (no k) or \code{"relieff"} (specify k)]
#' @param sd.frac multiplier of the standard deviation from the mean distances, subtracted from mean distance to create for SURF or multiSURF radius. The multiSURF default "dead-band radius" is sd.frac=0.5: mean - sd/2 
#' @param k number of constant nearest hits/misses for \code{"relieff"} (fixed k). 
#' The default k=0 means use the expected SURF theoretical k with sd.frac (.5 by default) for relieff nbd.
#' @param neighbor.sampling "none" or \code{"unique"} if you want to return only unique neighbor pairs
#' @param attr_removal_vec_from_dist_calc attributes for removal (possible confounders) from the distance matrix calculation. 
#' 
#' @return  Ri_NN.idxmat, matrix of Ri's (first column) and their NN's (second column)
#'
#' @examples
#' # multisurf neighborhood with sigma/2 (sd.frac=0.5) "dead-band" boundary
#' neighbor.pairs.idx <- nearestNeighbors(predictors.mat, nb.method="multisurf", nb.metric = "manhattan", sd.frac = 0.5)
#' # reliefF (fixed-k) neighborhood using default k equal to theoretical surf expected value
#' # One can change the theoretical value by changing sd.frac (default 0.5)
#' neighbor.pairs.idx <- nearestNeighbors(predictors.mat, nb.method="relieff", nb.metric = "manhattan")
#' # reliefF (fixed-k) neighborhood with a user-specified k
#' neighbor.pairs.idx <- nearestNeighbors(predictors.mat, nb.method="relieff", nb.metric = "manhattan", k=10)
#'
#' @export
nearestNeighbors <- function(attr.mat, 
                             nb.method = "multisurf", 
                             nb.metric = "manhattan", 
                             sd.vec = NULL, sd.frac = 0.5, k=0,
                             neighbor.sampling = "none",
                             attr_removal_vec_from_dist_calc=c()){
  # create a matrix with num.samp rows and two columns
  # first column is sample Ri, second is Ri's nearest neighbors
  
  if (!is.null(attr_removal_vec_from_dist_calc)){ 
    # remove attributes (possible confounders) from distance matrix calculation
    tryCatch(
      attr.mat <- attr.mat %>% data.frame() %>% 
        select(- attr_removal_vec_from_dist_calc), 
      error = function(c) 'The attribute to remove does not exist.'
    )
  }

  dist.mat <- attr.mat %>% as.matrix() %>% unname() %>%
    npdrDistances(metric = nb.metric) %>%
    as.data.frame()
  num.samp <- nrow(attr.mat)

  if (nb.method == "relieff"){  
    if (k==0){ # if no k specified or value 0
      # replace k with the theoretical expected value for SURF (close to multiSURF)
      erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
      # theoretical surf k (sd.frac=.5) for regression problems (does not depend on a hit/miss group)
      k <- floor((num.samp-1)*(1-erf(sd.frac/sqrt(2)))/2)  # uses sd.frac
    }
    
    Ri.nearestPairs.list <- vector("list", num.samp)
    for (Ri in colnames(dist.mat)){ # for each sample Ri
      Ri.int <- as.integer(Ri)
      Ri.nearest.idx <- dist.mat %>%
        dplyr::select(!!Ri) %>% # select the column Ri, hopefully reduce processing power
        rownames_to_column() %>% # push the neighbors from rownames to columns
        top_n(-(k+1), !!sym(Ri)) %>% # select the k closest neighbors, include self
        pull(rowname) %>% # get the neighbors
        as.integer() # convert from string (rownames - not factors) to integers
  
      if (!is.null(Ri.nearest.idx)){ # if neighborhood not empty
        # bind automatically repeated Ri, make sure to skip Ri self
        Ri.nearestPairs.list[[Ri.int]] <- data.frame(Ri_idx = Ri.int, NN_idx = Ri.nearest.idx[-1])
      }
    }
    
  } else {
    
    if (nb.method == "surf"){
      num.pair <- num.samp * (num.samp-1) / 2 # number of paired distances
      radius.surf <- sum(dist.mat)/(2*num.pair) # const r = mean(all distances)
      sd.const <- sd(dist.mat[upper.tri(dist.mat)])  
      # bam: orignal surf does not subtract sd-frac but should for fair multisurf comparison
      Ri.radius <- rep(radius.surf - sd.frac*sd.const, num.samp) 
    }
    if (nb.method == "multisurf"){
      if (is.null(sd.vec)) sd.vec <- sapply(1:num.samp, function(x) sd(dist.mat[-x, x]))
      Ri.radius <- colSums(dist.mat)/(num.samp - 1) - sd.frac*sd.vec # use adaptive radius
    }
    
    # put each Ri's nbd in a list then rbind them at the end with bind_rows()
    Ri.nearestPairs.list <- vector("list", num.samp) # initialize list
    
    for (Ri in colnames(dist.mat)){ # for each sample Ri
      Ri.int <- as.integer(Ri)
      Ri.nearest.idx <- dist.mat %>%
        dplyr::select(!!Ri) %>%
        rownames_to_column() %>% 
        filter((!!sym(Ri)) < Ri.radius[Ri]) %>%
        pull(rowname) %>%
        as.integer()
  
      if (!is.null(Ri.nearest.idx)){ # similar to relieff
        Ri.nearestPairs.list[[Ri.int]] <- data.frame(Ri_idx = Ri.int, NN_idx = Ri.nearest.idx[-1])
      }
    }
  }
  
  Ri_NN.idxmat <- dplyr::bind_rows(Ri.nearestPairs.list)
  
  if (neighbor.sampling=="unique"){
    # if you only want to return unique neighbors
    Ri_NN.idxmat <- uniqueNeighbors(Ri_NN.idxmat)
  }
  
  # matrix of Ri's (first column) and their NN's (second column)
  return(Ri_NN.idxmat)
}

#=========================================================================#
#' uniqueNeighbors
#'
#' Find pairs of unique nearest neighbors pairs from possible redundant pairs. 
#' Used as options (neighbor.sampling="unique") in nearestNeighbors and npdr functions. 
#'
#' @param neighbor.pairs.idx two columns of (possibly redundant) "i,j" pairs from nearestNeighbors function 
#' @return new neighborhood pair matrix of only unique pairs 
#'
#' @examples
#' unique.neighbor.pairs.idx <- uniqueNeighbors(neighbor.pairs.idx)  # unique neighbor pairs
#'
#' @export
uniqueNeighbors <- function(neighbor.pairs){
  # input: two columns of redundant "i,j" pairs
  # return: two columns of unique pairs from the redundant input
  # sort and make create redundant vector of "i,j" pairs
  # e.g., pairs 1  36 and 36  1 both become 1  36
  sorted_pairs <- data.frame(xmin = pmin(neighbor.pairs[, 1], neighbor.pairs[, 2]), 
                             xmax = pmax(neighbor.pairs[, 1], neighbor.pairs[, 2]))
  pair_str <- tidyr::unite(sorted_pairs, 'combined', xmin, xmax, sep = ',')
  unique.idx <- !duplicated(pair_str)
  return(neighbor.pairs[unique.idx,])
}

#=========================================================================#
#' knnVec
#'
#' Number of neighbors for each sample (vector) from a neighbor-pair matrix.
#'
#' @param neighbor.pairs.idx two columns of redundant "i,j" pairs from nearestNeighbors function 
#' @return  knn.vec vector number of nearest neighbors for each instance
#'
#' @examples
#' mean(knnVec(neighbor.pairs.idx))  # average number of neighbors
#'
#' @export
knnVec <- function(neighbor.pairs.mat){
  knn.vec <- data.frame(neighbor.pairs.mat) %>% dplyr::count(Ri_idx) %>% pull(n)
  return(knn.vec)
}

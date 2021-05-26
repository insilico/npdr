# =========================================================================#
#' npdrDiff
#'
#' A diff is a function that computes the difference of values for an attribute between two instances.
#' It is used for attribute selection for attribute diffs and phenotype diffs.
#' This function is vectorized: input a and b can be two vectors of values for one attribute.
#'
#' @param a value of attribute for first instance. Vector for correlation-data.
#' @param b value of attribute for second instance. Vector for correlation-data.
#' @param type diff rule for the given attribute data type, such as numeric, categorical or correlation-data vector.
#' @return val diff or vector of diffs
#' @examples
#' Example
#' @export
npdrDiff <- function(a, b, diff.type = "numeric-abs", norm.fac = 1) {
  # compute the difference between two vectors element-wise
  if (diff.type == "numeric-sqr") { # numeric squared difference
    val <- abs(a - b)^2 / norm.fac
  } else if (diff.type == "allele-sharing") { # snps
    val <- abs(a - b) / 2
  } else if (diff.type == "match-mismatch") {
    # used automatically for case-control pheno, optional genotype mismatch diff for snps
    val <- ifelse(a == b, 0, 1) # hit pairs are 0 and miss pairs are 1
    # val <- as.character(a==b) # convert this to factor in glmSTIR
  } else if (diff.type == "correlation-data") { # correlation data (e.g., fmri)      # corrdata
    val <- rowSums(abs(a - b) / norm.fac) # a and b are vectors in this case   # corrdata
  } else { # numeric abs difference
    val <- abs(a - b) / norm.fac
  }
  return(val)
}

# =========================================================================#
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
npdrDistances <- function(attr.mat, metric = "manhattan", fast.dist = FALSE) {
  if (fast.dist) {
    npdr.dist.fn <- wordspace::dist.matrix
  } else {
    npdr.dist.fn <- dist
  }

  # Compute distance matrix between all samples (rows)
  # default is numeric manhattan ("manhattan"), max-min scaling is only needed for relief
  if (metric == "hamming") {
    distance.mat <- hamming.binary(attr.mat)
  } else if (metric == "allele-sharing-manhattan") {
    # allele-sharing-manhattan, AM for SNPs
    attr.mat.scale <- attr.mat / 2
    distance.mat <- npdr.dist.fn(attr.mat.scale, method = "manhattan")
  } else if (metric == "relief-scaled-manhattan") {
    # value of metric, euclidean, manhattan or maximum
    maxminVec <- attr.range(attr.mat)
    minVec <- apply(attr.mat, 2, function(x) {
      min(x)
    })
    attr.mat.centered <- t(attr.mat) - minVec
    attr.mat.scale <- t(attr.mat.centered / maxminVec)
    distance.mat <- npdr.dist.fn(attr.mat.scale, method = "manhattan")
  } else if (metric == "relief-scaled-euclidean") {
    # value of metric, euclidean, manhattan or maximum
    maxminVec <- attr.range(attr.mat)
    minVec <- apply(attr.mat, 2, min)
    attr.mat.centered <- t(attr.mat) - minVec
    attr.mat.scale <- t(attr.mat.centered / maxminVec)
    distance.mat <- npdr.dist.fn(attr.mat.scale, method = "euclidean")
  } else if (metric == "euclidean") {
    distance.mat <- npdr.dist.fn(attr.mat, method = "euclidean")
  } else {
    distance.mat <- npdr.dist.fn(attr.mat, method = "manhattan")
  }
  as.matrix(distance.mat)
}

# =========================================================================#
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
#' @param att_to_remove attributes for removal (possible confounders) from the distance matrix calculation.
#' @param fast.dist whether or not distance is computed by faster algorithm in wordspace, default as F
#' @param dopar.nn whether or not neighborhood is computed in parallel, default as F
#' @return  Ri_NN.idxmat, matrix of Ri's (first column) and their NN's (second column)
#'
#' @examples
#' # multisurf neighborhood with sigma/2 (sd.frac=0.5) "dead-band" boundary
#' neighbor.pairs.idx <- nearestNeighbors(predictors.mat, nb.method = "multisurf", nb.metric = "manhattan", sd.frac = 0.5)
#' # reliefF (fixed-k) neighborhood using default k equal to theoretical surf expected value
#' # One can change the theoretical value by changing sd.frac (default 0.5)
#' neighbor.pairs.idx <- nearestNeighbors(predictors.mat, nb.method = "relieff", nb.metric = "manhattan")
#' # reliefF (fixed-k) neighborhood with a user-specified k
#' neighbor.pairs.idx <- nearestNeighbors(predictors.mat, nb.method = "relieff", nb.metric = "manhattan", k = 10)
#' @export
nearestNeighbors <- function(attr.mat,
                             nb.method = "multisurf",
                             nb.metric = "manhattan",
                             sd.vec = NULL, sd.frac = 0.5, k = 0,
                             neighbor.sampling = "none",
                             att_to_remove = c(), fast.dist = FALSE, dopar.nn = FALSE) {
  # create a matrix with num.samp rows and two columns
  # first column is sample Ri, second is Ri's nearest neighbors
  num.samp <- nrow(attr.mat)

  if (!is.null(att_to_remove)) {
    # remove attributes (possible confounders) from distance matrix calculation
    tryCatch(
      attr.mat <- attr.mat %>% data.frame() %>%
        select(-att_to_remove),
      error = function(c) "The attribute to remove does not exist."
    )
  }

  dist.mat <- attr.mat %>%
    as.matrix() %>%
    unname() %>%
    npdrDistances(metric = nb.metric, fast.dist = fast.dist) %>%
    as.data.frame()
  colnames(dist.mat) <- seq.int(num.samp)

  if (nb.method == "relieff") {
    if (k == 0) { # if no k specified or value 0
      # replace k with the theoretical expected value for SURF (close to multiSURF)
      erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
      # theoretical surf k (sd.frac=.5) for regression problems (does not depend on a hit/miss group)
      k <- floor((num.samp - 1) * (1 - erf(sd.frac / sqrt(2))) / 2) # uses sd.frac
    }

    if (dopar.nn == TRUE) {
      avai.cors <- parallel::detectCores() - 2
      cl <- parallel::makeCluster(avai.cors)
      doParallel::registerDoParallel(cl)
      Ri_NN.idxmat <- foreach::foreach(
        Ri.int = seq.int(num.samp), .combine = "rbind", .packages = c("dplyr", "tibble")
      ) %dopar% {
        Ri <- as.character(Ri.int)
        Ri.nearest.idx <- dist.mat %>%
          dplyr::select(!!Ri) %>% # select the column Ri, hopefully reduce processing power
          tibble::rownames_to_column() %>% # push the neighbors from rownames to a column named rowname
          top_n(-(k + 1), !!sym(Ri)) %>% # select the k closest neighbors, include self
          filter((!!sym(Ri)) > 0) %>% # top_n does not sort output, so make sure remove self
          pull(rowname) %>% # get the neighbors
          as.integer() # convert from string (rownames - not factors) to integers

        return(data.frame(Ri_idx = Ri.int, NN_idx = Ri.nearest.idx))
      }
      ## [1] 1.000000 1.414214 1.732051
      parallel::stopCluster(cl)
    } else {
      Ri.nearestPairs.list <- vector("list", num.samp)
      for (Ri in colnames(dist.mat)) { # for each sample Ri
        Ri.int <- as.integer(Ri)
        Ri.nearest.idx <- dist.mat %>%
          dplyr::select(!!Ri) %>%
          # select the column Ri, hopefully reduce processing power
          tibble::rownames_to_column() %>%
          # push the neighbors from rownames to columns
          dplyr::top_n(-(k + 1), !!sym(Ri)) %>%
          # select the k closest neighbors, include self
          dplyr::filter((!!sym(Ri)) > 0) %>%
          dplyr::pull(rowname) %>%
          # get the neighbors
          as.integer() # convert from string (rownames - not factors) to integers

        if (!is.null(Ri.nearest.idx)) { # if neighborhood not empty
          # bind automatically repeated Ri, make sure to skip Ri self
          Ri.nearestPairs.list[[Ri.int]] <- data.frame(Ri_idx = Ri.int, NN_idx = Ri.nearest.idx)
        }
      }
      Ri_NN.idxmat <- dplyr::bind_rows(Ri.nearestPairs.list)
    }
  } else {
    if (nb.method == "surf") {
      num.pair <- num.samp * (num.samp - 1) / 2 # number of paired distances
      radius.surf <- sum(dist.mat) / (2 * num.pair) # const r = mean(all distances)
      sd.const <- sd(dist.mat[upper.tri(dist.mat)])
      # bam: orignal surf does not subtract sd-frac but should for fair multisurf comparison
      Ri.radius <- rep(radius.surf - sd.frac * sd.const, num.samp)
      names(Ri.radius) <- as.character(1:num.samp)
    }
    if (nb.method == "multisurf") {
      if (is.null(sd.vec)) sd.vec <- sapply(1:num.samp, function(x) sd(dist.mat[-x, x]))
      Ri.radius <- colSums(dist.mat) / (num.samp - 1) - sd.frac * sd.vec # use adaptive radius
    }
    if (dopar.nn == TRUE) {
      avai.cors <- parallel::detectCores() - 2
      cl <- parallel::makeCluster(avai.cors)
      doParallel::registerDoParallel(cl)
      Ri_NN.idxmat <- foreach::foreach(
        Ri.int = seq.int(num.samp), .combine = "rbind", .packages = c("dplyr", "tibble")
      ) %dopar% {
        Ri <- as.character(Ri.int)
        Ri.nearest.idx <- dist.mat %>%
          dplyr::select(!!Ri) %>%
          # select the column Ri, hopefully reduce processing power
          tibble::rownames_to_column() %>%
          # push the neighbors from rownames to columns
          dplyr::filter(((!!sym(Ri)) < Ri.radius[Ri]) & ((!!sym(Ri)) > 0)) %>%
          dplyr::pull(rowname) %>%
          # get the neighbors
          as.integer() # convert from string (rownames - not factors) to integers

        return(data.frame(Ri_idx = Ri.int, NN_idx = Ri.nearest.idx))
      }
      parallel::stopCluster(cl)
    } else {
      # put each Ri's nbd in a list then rbind them at the end with bind_rows()
      Ri.nearestPairs.list <- vector("list", num.samp) # initialize list

      for (Ri in colnames(dist.mat)) { # for each sample Ri
        Ri.int <- as.integer(Ri)
        Ri.nearest.idx <- dist.mat %>%
          dplyr::select(!!Ri) %>%
          rownames_to_column() %>%
          filter(((!!sym(Ri)) < Ri.radius[Ri]) & ((!!sym(Ri)) > 0)) %>%
          pull(rowname) %>%
          as.integer()

        if (!is.null(Ri.nearest.idx)) { # similar to relieff
          Ri.nearestPairs.list[[Ri.int]] <- data.frame(Ri_idx = Ri.int, NN_idx = Ri.nearest.idx)
        }
      }
      Ri_NN.idxmat <- dplyr::bind_rows(Ri.nearestPairs.list)
    }
  }


  if (neighbor.sampling == "unique") {
    # if you only want to return unique neighbors
    Ri_NN.idxmat <- uniqueNeighbors(Ri_NN.idxmat)
  }

  # matrix of Ri's (first column) and their NN's (second column)
  return(Ri_NN.idxmat)
}

# =========================================================================#
#' nearestNeighborsSeparateHitMiss
#'
#' Find nearest neighbors of each instance using relief.method.
#' Treat the hit and miss distributions separately to circument potential hit bias.
#' ReliefF version makes hit/miss neighborhoods balanced. Surf and MultiSurf are still imbalanced.
#' Used for npdr (no hits or misses specified in neighbor function).
#'
#' @param attr.mat m x p matrix of m instances and p attributes
#' @param pheno.vec vector of class values for m instances
#' @param nb.metric used in npdrDistances for distance matrix between instances, default: \code{"manhattan"} (numeric)
#' @param nb.method neighborhood method [\code{"multisurf"} or \code{"surf"} (no k) or \code{"relieff"} (specify k)]
#' @param sd.frac multiplier of the standard deviation from the mean distances, subtracted from mean distance to create for SURF or multiSURF radius. The multiSURF default "dead-band radius" is sd.frac=0.5: mean - sd/2
#' @param k number of constant nearest hits/misses for \code{"relieff"} (fixed k).
#' The default k=0 means use the expected SURF theoretical k with sd.frac (.5 by default) for relieff nbd.
#' @param neighbor.sampling "none" or \code{"unique"} if you want to return only unique neighbor pairs
#' @param att_to_remove attributes for removal (possible confounders) from the distance matrix calculation.
#' @param fast.dist whether or not distance is computed by faster algorithm in wordspace, default as F
#' @param dopar.nn whether or not neighborhood is computed in parallel, default as F
#' @return  Ri_NN.idxmat, matrix of Ri's (first column) and their NN's (second column)
#'
#' @examples
#' # reliefF (fixed-k) neighborhood using default k equal to theoretical surf expected value
#' # One can change the theoretical value by changing sd.frac (default 0.5)
#' neighbor.pairs.idx <- nearestNeighborsSeparateHitMiss(cc.attrs, cc.pheno, # need attributes and pheno
#'   nb.method = "relieff", nb.metric = "manhattan",
#'   sd.frac = .5, k = 0
#' )
#' @export
nearestNeighborsSeparateHitMiss <- function(attr.mat, pheno.vec,
                                            nb.method = "relieff",
                                            nb.metric = "manhattan",
                                            sd.vec = NULL, sd.frac = 0.5, k = 0,
                                            neighbor.sampling = "none",
                                            att_to_remove = c(), fast.dist = FALSE, dopar.nn = FALSE) {
  # create a matrix with num.samp rows and two columns
  # first column is sample Ri, second is Ri's nearest neighbors
  num.samp <- nrow(attr.mat)
  pheno.vec <- as.numeric(as.character(pheno.vec))
  majority.pheno <- which.min(table(pheno.vec)) %>%
    names() %>%
    as.integer()
  majority.frac <- max(table(pheno.vec)) / length(pheno.vec)

  if (!is.null(att_to_remove)) {
    # remove attributes (possible confounders) from distance matrix calculation
    tryCatch(
      attr.mat <- attr.mat %>% data.frame() %>%
        select(-att_to_remove),
      error = function(c) "The attribute to remove does not exist."
    )
  }

  dist.mat <- attr.mat %>%
    as.matrix() %>%
    unname() %>%
    npdrDistances(metric = nb.metric, fast.dist = fast.dist) %>%
    as.data.frame()
  colnames(dist.mat) <- seq.int(num.samp)

  if (nb.method == "relieff") {
    if (k == 0) { # if no k specified or value 0
      # replace k with the theoretical expected value for SURF (close to multiSURF)
      erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
      # theoretical surf k (sd.frac=.5) for regression problems (does not depend on a hit/miss group)
      k <- floor((num.samp - 1) * (1 - erf(sd.frac / sqrt(2))) / 2) # uses sd.frac
    }

    if (dopar.nn == TRUE) {
      avai.cors <- parallel::detectCores() - 2
      cl <- parallel::makeCluster(avai.cors)
      doParallel::registerDoParallel(cl)
      Ri.nearestPairs.list <- foreach::foreach(
        Ri.int = seq.int(num.samp), .packages = c("dplyr", "tibble")
      ) %dopar% {
        # Ri <- as.character(Ri.int)
        # Ri.int <- as.integer(Ri)
        Ri.distances <- dist.mat[Ri.int, ] # all distances to sample Ri
        Ri.nearest <- order(Ri.distances, decreasing = F) # closest to farthest
        # consider distance distributions of hits and misses separately
        Ri.hits <- Ri.nearest[pheno.vec[Ri.int] == pheno.vec[Ri.nearest]]
        Ri.misses <- Ri.nearest[pheno.vec[Ri.int] != pheno.vec[Ri.nearest]]
        # make hit and miss neighborhoods the same size
        # depending on whether Ri is majority or minority class, the number of hits/misses changes
        if (pheno.vec[Ri.int] == majority.pheno) {
          Ri.nearest.idx <- Ri.hits[2:floor(majority.frac * k + 1)] # (2) skip Ri self
          # concatenate misses
          Ri.nearest.idx <- c(Ri.nearest.idx, Ri.misses[1:floor((1 - majority.frac) * k + 1)])
        } else {
          Ri.nearest.idx <- Ri.hits[2:floor((1 - majority.frac) * k + 1)] # (2) skip Ri self
          # concatenate misses
          Ri.nearest.idx <- c(Ri.nearest.idx, Ri.misses[1:floor(majority.frac * k + 1)])
        }

        if (!is.null(Ri.nearest.idx)) { # if neighborhood not empty
          # bind automatically repeated Ri, make sure to skip Ri self
          return(data.frame(Ri_idx = Ri.int, NN_idx = Ri.nearest.idx))
        }
      }
      ## [1] 1.000000 1.414214 1.732051
      parallel::stopCluster(cl)
    } else {
      Ri.nearestPairs.list <- vector("list", num.samp)
      for (Ri in colnames(dist.mat)) { # for each sample Ri
        Ri.int <- as.integer(Ri)
        Ri.distances <- dist.mat[Ri.int, ] # all distances to sample Ri
        Ri.nearest <- order(Ri.distances, decreasing = F) # closest to farthest
        # consider distance distributions of hits and misses separately
        Ri.hits <- Ri.nearest[pheno.vec[Ri.int] == pheno.vec[Ri.nearest]]
        Ri.misses <- Ri.nearest[pheno.vec[Ri.int] != pheno.vec[Ri.nearest]]
        # for misses, option to use farthest is not a good idea because it makes all variables appear
        # different between groups, even null variables
        # if (miss.ordering=="farthest"){ # choose misses that are farthest from Ri
        #  Ri.misses <- rev(Ri.misses)
        #    }
        #
        # make hit and miss neighborhoods the same size (balanced)
        # depending on whether Ri is majority or minority class, the number of hits/misses changes
        if (pheno.vec[Ri.int] == majority.pheno) {
          Ri.nearest.idx <- Ri.hits[2:floor(majority.frac * k + 1)] # (2) skip Ri self
          # concatenate misses
          Ri.nearest.idx <- c(Ri.nearest.idx, Ri.misses[1:floor((1 - majority.frac) * k + 1)])
        } else {
          Ri.nearest.idx <- Ri.hits[2:floor((1 - majority.frac) * k + 1)] # (2) skip Ri self
          # concatenate misses
          Ri.nearest.idx <- c(Ri.nearest.idx, Ri.misses[1:floor(majority.frac * k + 1)])
        }

        if (!is.null(Ri.nearest.idx)) { # if neighborhood not empty
          # bind automatically repeated Ri, make sure to skip Ri self
          Ri.nearestPairs.list[[Ri.int]] <- data.frame(Ri_idx = Ri.int, NN_idx = Ri.nearest.idx)
        }
      } # end for
      # Ri_NN.idxmat <- dplyr::bind_rows(Ri.nearestPairs.list)
    } # end else dopar.nn

    Ri_NN.idxmat <- dplyr::bind_rows(Ri.nearestPairs.list)
  } else { # surf or multisurf...


    # For treating hit/miss distance distributions separately, compute separate hit and miss radii
    # User might want to shrink alpha standard deviation fraction. Unlike relieff, the hit and miss
    # neighborhoods are not balanced.
    if (nb.method == "surf") { # compute surf radii

      hit.dist.rows <- vector("list", num.samp)
      for (i in seq(1, num.samp)) {
        hit.mask <- pheno.vec[i] == pheno.vec
        hit.mask <- hit.mask[-i] # remove self
        hit.dist.rows[[i]] <- dist.mat[i, hit.mask]
      }
      hit.dist.vec <- unlist(hit.dist.rows)
      # average of all hit neighbors
      Ri.hit.radii <- rep(mean(hit.dist.vec) - sd.frac * sd(hit.dist.vec), num.samp)
      names(Ri.hit.radii) <- as.character(1:num.samp)

      miss.dist.rows <- vector("list", num.samp)
      for (i in seq(1, num.samp)) {
        miss.mask <- pheno.vec[i] != pheno.vec
        miss.dist.rows[[i]] <- dist.mat[i, miss.mask]
      }
      miss.dist.vec <- unlist(miss.dist.rows)
      # average of all miss neighbors
      Ri.miss.radii <- rep(mean(miss.dist.vec) - sd.frac * sd(miss.dist.vec), num.samp)
      names(Ri.miss.radii) <- as.character(1:num.samp)
    } # end surf radius calc

    if (nb.method == "multisurf") { # compute multisurf radii

      Ri.hit.radii <- vector("numeric", num.samp)
      Ri.miss.radii <- vector("numeric", num.samp)
      for (i in seq(1, num.samp)) {
        # grab neighbors that are hits of Ri
        hit.mask <- pheno.vec[i] == pheno.vec
        hit.mask <- hit.mask[-i] # remove self
        hit.dist.row <- as.numeric(dist.mat[i, hit.mask])
        Ri.hit.radii[i] <- mean(hit.dist.row) - sd.frac * sd(hit.dist.row)

        # grab neighbors that are misses of Ri
        miss.mask <- pheno.vec[i] != pheno.vec
        miss.dist.row <- as.numeric(dist.mat[i, miss.mask])
        Ri.miss.radii[i] <- mean(miss.dist.row) - sd.frac * sd(miss.dist.row)
      }

      names(Ri.hit.radii) <- as.character(1:num.samp)
      names(Ri.miss.radii) <- as.character(1:num.samp)
    } # end multisurf radii calc

    if (dopar.nn == TRUE) {
      avai.cors <- parallel::detectCores() - 2
      cl <- parallel::makeCluster(avai.cors)
      doParallel::registerDoParallel(cl)
      Ri.nearestPairs.list <- foreach::foreach(
        Ri.int = seq.int(num.samp), .packages = c("dplyr", "tibble")
      ) %dopar% {
        Ri.distances <- dist.mat[Ri.int, ]
        Ri.nearest.hits <- which((pheno.vec[Ri.int] == pheno.vec) & (Ri.distances < Ri.hit.radii[Ri.int]) &
          (Ri.distances > 0)) # skip Ri self (dist=0)
        Ri.nearest.misses <- which((pheno.vec[Ri.int] != pheno.vec) & (Ri.distances < Ri.miss.radii[Ri.int]))
        # join hit and miss into one nbd
        Ri.nearest.idx <- c(Ri.nearest.hits, Ri.nearest.misses)

        if (!is.null(Ri.nearest.idx)) { # similar to relieff
          return(data.frame(Ri_idx = Ri.int, NN_idx = Ri.nearest.idx))
        }
      }
      parallel::stopCluster(cl)
    } else {
      # put each Ri's nbd in a list then rbind them at the end with bind_rows()
      Ri.nearestPairs.list <- vector("list", num.samp) # initialize list

      for (Ri in colnames(dist.mat)) { # for each sample Ri
        Ri.int <- as.integer(Ri)
        Ri.distances <- dist.mat[Ri.int, ]
        Ri.nearest.hits <- which((pheno.vec[Ri.int] == pheno.vec) & (Ri.distances < Ri.hit.radii[Ri.int]) &
          (Ri.distances > 0)) # skip Ri self (dist=0)
        Ri.nearest.misses <- which((pheno.vec[Ri.int] != pheno.vec) & (Ri.distances < Ri.miss.radii[Ri.int]))
        # join hit and miss into one nbd
        Ri.nearest.idx <- c(Ri.nearest.hits, Ri.nearest.misses)

        if (!is.null(Ri.nearest.idx)) { # similar to relieff
          Ri.nearestPairs.list[[Ri.int]] <- data.frame(Ri_idx = Ri.int, NN_idx = Ri.nearest.idx)
        }
      } # end for
    } # end else dopar.nn
    Ri_NN.idxmat <- dplyr::bind_rows(Ri.nearestPairs.list)
  }


  if (neighbor.sampling == "unique") {
    # if you only want to return unique neighbors
    Ri_NN.idxmat <- uniqueNeighbors(Ri_NN.idxmat)
  }

  # matrix of Ri's (first column) and their NN's (second column)
  return(Ri_NN.idxmat)
}

# =========================================================================#
#' uniqueNeighbors
#'
#' Find pairs of unique nearest neighbors pairs from possible redundant pairs.
#' Used as options (neighbor.sampling="unique") in nearestNeighbors and npdr functions.
#'
#' @param neighbor.pairs.idx two columns of (possibly redundant) "i,j" pairs from nearestNeighbors function
#' @return new neighborhood pair matrix of only unique pairs
#'
#' @examples
#' unique.neighbor.pairs.idx <- uniqueNeighbors(neighbor.pairs.idx) # unique neighbor pairs
#' @export
uniqueNeighbors <- function(neighbor.pairs) {
  # input: two columns of redundant "i,j" pairs
  # return: two columns of unique pairs from the redundant input
  # sort and make create redundant vector of "i,j" pairs
  # e.g., pairs 1  36 and 36  1 both become 1  36
  sorted_pairs <- data.frame(
    xmin = pmin(neighbor.pairs[, 1], neighbor.pairs[, 2]),
    xmax = pmax(neighbor.pairs[, 1], neighbor.pairs[, 2])
  )
  pair_str <- tidyr::unite(sorted_pairs, "combined", xmin, xmax, sep = ",")
  unique.idx <- !duplicated(pair_str)
  return(neighbor.pairs[unique.idx, ])
}

# =========================================================================#
#' knnVec
#'
#' Number of neighbors for each sample (vector) from a neighbor-pair matrix.
#'
#' @param neighbor.pairs.idx two columns of redundant "i,j" pairs from nearestNeighbors function
#' @return  knn.vec vector number of nearest neighbors for each instance
#'
#' @examples
#' mean(knnVec(neighbor.pairs.idx)) # average number of neighbors
#' @export
knnVec <- function(neighbor.pairs.mat) {
  knn.vec <- data.frame(neighbor.pairs.mat) %>%
    dplyr::count(Ri_idx) %>%
    pull(n)
  return(knn.vec)
}

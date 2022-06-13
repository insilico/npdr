# ==============================================================#
#' \code{npdrLearnerCV}
#'
#' Tune a hyperparmeter that maximizes the cross-validation accuracy of a k-nearest-neighbors classifier. You can tune k, but keep in mind that the resulting k might be underestimated because the training sample size is smaller than the original sample size. When other hyperparameters are optimized, k is fixed to the npdr theoretical value that adapts to the training size (todo: make more flexible with knn alpha). You can tune the number of ICA or PCA components as the components are used as the space for calculating nearest neighbors. todo: create function interface that allows user to create their own sapply_hyper_fn.
#' @param x (m+1) x p dataframe of m instances, 1 class column and p attributes
#' @param label column label for class \code{"class"}
#' @param tune_grid vector of hyperparameter values to test for best classification accuracy
#' @param dist_metric for distance matrix between instances
#' (default: \code{"manhattan"}, others include \code{"euclidean"},
#' and for GWAS \code{"allele-sharing-manhattan"}).
#' @param tune_type type of hyperparmater to optimize. default: \code{"knn"}, others include \code{"ica"} (number of ica components for ica space transformation, and \code{"pca"} (number of components for PCA transformation.
#' @param num_folds number of cross-validation folds for tuning
#' @return list containing best hyperparameter (best_param), its highest accuracy (best_acc), and a table of fold and parameter accuracies (cv_table) 
#' @examples
#' library(flexclust) # need for npdrLearner knn classifier
#' library(fastICA)   # need if tuning ica tansformation
#' cv.out <- npdrLearnerCV(x=dats, label="class", 
#'               tune_grid = seq(20,90,5),   # tuning knn
#'               dist_metric = "manhattan",
#'               tune_type = "knn",
#'               num_folds=5, verbose=T)
#' cv.out$best_param
#' plot(cv.out$cv_table$hyp,cv.out$cv_table$means,
#'         xlab="hyperparameter", ylab="accuracy", 
#'         main="CV hyperparameter tuning", type="l")
#' text(cv.out$best_param,cv.out$best_acc,paste("max.loc =",cv.out$best_param))
#' Or you can tune number of knns 
#' cv.out <- npdrLearnerCV(x=dats, label="class", 
#'                      tune_grid = seq(20,90,5),   # tuning knn
#'                        dist_metric = "manhattan",
#'                        tune_type = "knn",
#'                        num_folds=5, verbose=T)
#' @export
npdrLearnerCV <- function(x, label="class", 
                          tune_grid = seq(10,90,10), #knn
                          dist_metric = "manhattan",
                          tune_type = "knn",
                          num_folds = 5, verbose=F) 
{
  # Cross-Validtation Hyperparameter optimization 
  # by knn classification accuracy
  class_idx <- which(colnames(x)==label)
  folds <- caret::createFolds(x[,class_idx], k = num_folds) 
  ### switch the thing we want to tune (tune_type)
  sapply_hyper_fn = switch(tune_type,
                           knn = function(hyper.param,te.idx) {
                             # tune_type = "knn"
                             # ex. tune_grid = seq(10,90,10)
                             # values should be less than num samples-1
                             test_results <- npdrLearner(train.outcome="class", 
                                                         train.data=x[-te.idx,], 
                                                         test.outcome="class", 
                                                         test.data=x[te.idx,],
                                                         nbd.method = "relieff", 
                                                         nbd.metric = dist_metric, 
                                                         msurf.sd.frac = 0.5, 
                                                         knn=hyper.param) 
                             return(test_results$accuracy) },
                           ica = function(hyper.param,te.idx) { 
                             #   tune_type = "ica"
                             # ex.tune_grid = seq(10,80,10) # num of ICs
                             # values should be less than number of variables
                             # compute ICs for all samples for given n.comp
                             ICs <- fastICA(x[,-class_idx], 
                                            n.comp=hyper.param, fun="exp",
                                            method = "C", verbose=F,
                                            maxit=100000,tol=0.00000001)
                             IC.src <- data.frame(ICs$S)
                             IC.src$class <- x[,class_idx]  # add class back
                             m.samp <- nrow(IC.src)  # full sample size
                             # use training sample size for theoretical knn
                             # b/c knn's calculated in training
                             m.train <- m.samp-m.samp/num_folds # num train samp
                             k.train <- npdr::knnSURF(m.train - 1, 0.5)
                             test_results <- npdrLearner(train.outcome="class", 
                                                         train.data=IC.src[-te.idx,], 
                                                         test.outcome="class", 
                                                         test.data=IC.src[te.idx,],
                                                         nbd.method = "relieff", 
                                                         nbd.metric = dist_metric, 
                                                         msurf.sd.frac = 0.5, 
                                                         knn=k.train) 
                             return(test_results$accuracy) },
                           pca = function(hyper.param,te.idx) { 
                             #   tune_type = "pca"
                             # ex.tune_grid = seq(5,50,5) # num of PCs
                             # values should be less than number of variables
                             # compute ICs for all samples for given n.comp
                             PCs<-prcomp(x[,-class_idx])
                             topPCs <- data.frame(PCs$x[,1:hyper.param]) 
                             topPCs$class <- x[,class_idx]  # add class back
                             m.samp <- nrow(topPCs)  # full sample size
                             # use training sample size for theoretical knn
                             # b/c knn's calculated in training
                             m.train <- m.samp-m.samp/num_folds # num train samp
                             k.train <- npdr::knnSURF(m.train - 1, 0.5)
                             test_results <- npdrLearner(train.outcome="class", 
                                                         train.data=topPCs[-te.idx,], 
                                                         test.outcome="class", 
                                                         test.data=topPCs[te.idx,],
                                                         nbd.method = "relieff", 
                                                         nbd.metric = dist_metric, 
                                                         msurf.sd.frac = 0.5, 
                                                         knn=k.train) 
                             return(test_results$accuracy) }
  ) # end switch
  #### Begin cross validation
  ## Iterate over cv folds
  cv.results <- list()
  for (fold.id in seq(1,num_folds)){
    te.idx <- folds[[fold.id]]
    if (verbose){cat("fold", fold.id, "of",num_folds,"\n")}
    # add tune-alpha and tune-pca
    if(verbose){cat("\t inner loop over hyperparameters...\n")}
    # iterate over hyperparameter
    scores <- sapply(tune_grid,        # hyp loop var
                     function(hyp){sapply_hyper_fn(hyp,te.idx=te.idx)}
    )  # end sapply hyp loop over hyperparameters
    cv.results[[fold.id]] <- scores  # scores vector
  } # end for folds loop
  cv.results <- data.frame(cv.results)  # turn list to df
  cv.results$means <- rowMeans(as.matrix(cv.results))
  cv.results$hyp <- tune_grid
  colnames(cv.results) <- c(names(folds),"means","hyp")
  #### Select best performance
  best.idx <- which.max(cv.results$means)  # accuracy
  
  return(list(
    best_param = tune_grid[best.idx],
    best_acc = cv.results$means[best.idx],
    cv_table = cv.results
  ))
} # end npdrLearnerCV

# =========================================================================#
#' npdrLearner
#'
#' Uses npdr neighorhoods to learn a nearest neighbor classification or
#' regression model (latter not implemented but easy).
#' Finds the nearest neighbors of test instances to a training dataset.
#' Uses majority class of training neighbors for test prediction.
#' Allows adaptive Relief neighborhoods or specify k.
#' Regression would simply use the average value of the neighbor phenotypes.
#' Uses functions npdrDistances2 and nearestNeighbors2.
#'
#' @param train.outcome character name or length-m numeric outcome vector for
#' train data
#' @param train.data m x p matrix of m instances and p attributes of train data.
#' May also include outcome vector but then outcome should be a name.
#' Include attr names as colnames.
#' @param test.outcome character name or length-m numeric outcome vector of test data
#' @param test.data m x p matrix of m instances and p attributes for test data.
#' May also include outcome vector but then outcome should be a name.
#' Include attr names as colnames.
#' @param nbd.method neighborhood method: `multisurf` or `surf` (no k) or
#' `relieff` (specify k). Used by nearestNeighbors2().
#' @param nbd.metric used in npdrDistances2 for distance matrix between
#' instances, default: `manhattan` (numeric). Used by nearestNeighbors2().
#' @param knn number of constant nearest hits/misses for `relieff` (fixed-k).
#' Used by nearestNeighbors2().
#' The default knn=0 means use the expected SURF theoretical `k` with
#' `msurf.sd.frac` (0.5 by default)
#' @param msurf.sd.frac multiplier of the standard deviation from the mean
#' distances; subtracted from mean for SURF or multiSURF.
#' The multiSURF default is `msurf.sd.frac=0.5`: mean - sd/2.
#' Used by nearestNeighbors2().
#' @param dopar.nn whether or not neighborhood is computed in parallel,
#' default as FALSE.
#' @return  list: neighborhoods for each test instance, prediction for each
#' test instance, accuracy on test set
#'
#' @examples
#' test.results <- npdrLearner(
#'   train.outcome = "class", train.data = case.control.3sets$train,
#'   test.outcome = "class", test.data = case.control.3sets$validation,
#'   nbd.method = "relieff", nbd.metric = "manhattan",
#'   dopar.nn = FALSE, knn = 0
#' )
#' test.results$accuracy
#' @export
npdrLearner <- function(train.outcome, train.data, test.outcome, test.data,
                        nbd.method = "relieff", nbd.metric = "manhattan",
                        msurf.sd.frac = 0.5, dopar.nn = FALSE, knn = 0) {
  if (length(train.outcome) == 1) {
    # e.g., outcome="class" or outcome=101 (pheno col index) and dataset is data.frame including outcome variable
    train.pheno <- train.data[, train.outcome, drop = TRUE] %>%
      as.character() %>%
      as.numeric() # get phenotype
    train.data <- train.data %>% select(-train.outcome) # outcome = "qtrait" or 101
  } else { # user specifies a separate phenotype vector
    train.pheno <- train.outcome # assume users provides a separate outcome data vector
    # train.data <- train.data # assumes dataset only contains attributes/predictors
  }
  if (length(test.outcome) == 1) {
    # e.g., outcome="class" or outcome=101 (pheno col index) and dataset is data.frame including outcome variable
    test.pheno <- test.data[, test.outcome] %>%
      as.character() %>%
      as.numeric() # get phenotype
    test.data <- test.data %>% select(-test.outcome) # outcome = "qtrait" or 101
  } else { # user specifies a separate phenotype vector
    test.pheno <- test.outcome # assume users provides a separate outcome data vector
    # test.data <- test.data # assumes dataset only contains attributes/predictors
  }
  # get ids of nearest training instances for each test instance
  # number of elements in test.neighbors list equal to test.data sample size
  # each element has nearest neighbors from train.data
  test.neighbors <- nearestNeighbors2(train.data, test.data,
    nbd.method = nbd.method,
    nbd.metric = nbd.metric,
    sd.vec = NULL, sd.frac = msurf.sd.frac,
    k = knn, dopar.nn = dopar.nn
  )
  # predict test instances based on the majority class of nearest training neighbors
  test.predict <- sapply(
    1:length(test.neighbors),
    function(i) {
      table(train.pheno[test.neighbors[[i]]]) %>%
        which.max() %>%
        names() %>%
        as.numeric()
    }
  )
  test.acc <- sum(test.predict == test.pheno) / length(test.pheno)
  list(neighborhoods = test.neighbors, prediction = test.predict, accuracy = test.acc)
}


# =========================================================================#
#' npdrDistances2
#'
#' Create m1 x m2 distance matrix between two datasets,
#' (m1 instances and p attributes in dataset1 and
#' m2 instances and p attributes in dataset2).
#' Datasets should not include phenotype column.
#' Uses function dist2 from flexclust.
#' Used by nearestNeighbors2().
#'
#' @param attr.mat1 m1 x p matrix of m instances and p attributes
#' @param attr.mat2 m2 x p matrix of m instances and p attributes
#' @param metric for distance matrix between instances
#' (default: \code{"manhattan"}, others include \code{"euclidean"},
#' and for GWAS \code{"allele-sharing-manhattan"}).
#' @return matrix of m1 instances x m2 instances pairwise distances.
#' @examples
#' train_dat <- case.control.3sets$train
#' valid_dat <- case.control.3sets$validation
#'
#' if (require("flexclust")) {
#'   dist.mat <- npdrDistances2(
#'     train_dat[, names(train_dat) != "class"],
#'     valid_dat[, names(valid_dat) != "class"],
#'     metric = "manhattan"
#'   )
#' }
#' @export
npdrDistances2 <- function(attr.mat1, attr.mat2, metric = "manhattan") {
  check_installed("flexclust", reason = "for matrix distance computation `dist2()`")

  # first mat is rows and second is columns
  npdr.dist.fn <- flexclust::dist2
  # Compute distance matrix between all samples (rows) between test and training data
  # default is numeric manhattan ("manhattan"), max-min scaling is only needed for relief
  if (metric == "allele-sharing-manhattan") {
    # allele-sharing-manhattan, AM for SNPs
    distance.mat <- npdr.dist.fn(attr.mat1 / 2, attr.mat2 / 2, method = "manhattan")
  } else if (metric == "euclidean") {
    distance.mat <- npdr.dist.fn(attr.mat1, attr.mat2, method = "euclidean")
  } else {
    distance.mat <- npdr.dist.fn(attr.mat1, attr.mat2, method = "manhattan")
  }
  as.matrix(distance.mat)
}

# =========================================================================#
#' nearestNeighbors2
#'
#' Find nearest neighbors of each instance in attr.mat2 (test) to instances in
#' attr.mat1 (train) using relief neighborhood methods.
#' Used by npdrLearner, nearest neighbor classifier.
#' Input data should not include phenotype column.
#'
#' @param attr.mat1 m1 x p matrix of m instances and p attributes (training data)
#' @param attr.mat2 m2 x p matrix of m instances and p attributes (test data)
#' @param nbd.metric used in npdrDistances2 for distance matrix between
#' instances, default: `manhattan` (numeric)
#' @param nbd.method neighborhood method: `multisurf` or `surf` (no k) or
#' `relieff` (specify k)
#' @param sd.vec vector of standard deviations
#' @param sd.frac multiplier of the standard deviation from the mean distances,
#' subtracted from mean distance to create for SURF or multiSURF radius.
#' The multiSURF default "dead-band radius" is sd.frac=0.5: mean - sd/2
#' @param k number of constant nearest hits/misses for `relieff` (fixed k).
#' The default k=0 means use the expected SURF theoretical k with sd.frac
#' (0.5 by default) for relieff nbd.
#' @param dopar.nn whether or not neighborhood is computed in parallel, default as F
#' @return list of Ri's (data2 test instances) NN's in data1 (train instances)
#'
#' @examples
#' train_dat <- case.control.3sets$train
#' valid_dat <- case.control.3sets$validation
#' test.neighbors <- nearestNeighbors2(
#'   train_dat[, names(train_dat) != "class"],
#'   valid_dat[, names(valid_dat) != "class"], # no phenotype column
#'   nbd.method = "relieff",
#'   nbd.metric = "manhattan",
#'   sd.vec = NULL, sd.frac = 0.5,
#'   k = 0, # uses multisurf k estimate
#'   dopar.nn = FALSE
#' )
#' @export
nearestNeighbors2 <- function(attr.mat1, attr.mat2,
                              nbd.method = "multisurf",
                              nbd.metric = "manhattan",
                              sd.vec = NULL, sd.frac = 0.5, dopar.nn = FALSE,
                              k = 0) {
  if (dopar.nn) {
    check_installed("foreach", reason = "for fast parallel computing with `foreach()` and `%dopar%`")
    check_installed("doParallel", reason = "for `registerDoParallel()`")
    check_installed("parallel", reason = "for `makeCluster()`, `detectCores()`, and `stopCluster()`")
    `%dopar%` <- foreach::`%dopar%`
  }

  # create a matrix with num.samp rows and two columns
  # first column is sample Ri, second is Ri's nearest neighbors
  num.samp1 <- nrow(attr.mat1) # training set, rows of dist mat
  num.samp2 <- nrow(attr.mat2) # testing set, cols of dist mat

  dist.mat <- npdrDistances2(
    as.matrix(attr.mat1),
    as.matrix(attr.mat2),
    metric = nbd.metric
  )
  dist.df <- dist.mat %>% 
    as.data.frame() %>% 
    `colnames<-`(seq.int(num.samp2)) %>% 
    `rownames<-`(seq.int(num.samp1))

  if (nbd.method == "relieff") {
    if (k == 0) { # if no k specified or value 0
      # replace k with the theoretical expected value for SURF (close to multiSURF)
      erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
      # theoretical surf k (sd.frac=.5) for regression problems (does not depend on a hit/miss group)
      k <- floor((num.samp1 - 1) * (1 - erf(sd.frac / sqrt(2))) / 2) # uses sd.frac
    }

    if (dopar.nn) {
      avai.cors <- parallel::detectCores() - 2
      cl <- parallel::makeCluster(avai.cors)
      doParallel::registerDoParallel(cl)
      Ri.nearestNeighbors.list <- vector("list", num.samp2)
      Ri.nearestNeighbors.list <-
        foreach::foreach(Ri.int = seq.int(num.samp2), .packages = c("dplyr")) %dopar% {

          # foreach return, makes Ri.nearestNeighbors.list
          return(get_Ri_nearest(dist.df, as.character(Ri.int), nbd.method, k = k))
        }
      parallel::stopCluster(cl)
    } else { # relieff, no parallel
      Ri.nearestNeighbors.list <- vector("list", num.samp2)
      # look down each column of dist.df for neighbors of Ri
      for (Ri in colnames(dist.df)) { # for each instance/column Ri
        Ri.int <- as.integer(Ri)
        Ri.nearest.idx <- get_Ri_nearest(dist.df, as.character(Ri.int), nbd.method, k = k)

        if (!is.null(Ri.nearest.idx)) { # if neighborhood not empty
          Ri.nearestNeighbors.list[[Ri.int]] <- Ri.nearest.idx
        }
      } # end for
    }
    return(Ri.nearestNeighbors.list)
  }
  
  if (nbd.method == "surf") {
    radius.surf <- mean(dist.mat) # const r = mean(all distances)
    sd.const <- sd(dist.mat)
    # bam: orignal surf does not subtract sd-frac but should for fair multisurf comparison
    Ri.radius <- rep(radius.surf - sd.frac * sd.const, num.samp2)
    names(Ri.radius) <- as.character(1:num.samp2)
  }
  
  if (nbd.method == "multisurf") {
    sd.vec <- sd.vec %||% sapply(1:num.samp2, function(i) sd(dist.mat[, i]))
    Ri.radius <- colMeans(dist.mat) - sd.frac * sd.vec # use adaptive radius
    names(Ri.radius) <- as.character(1:num.samp2)
  }
  
  if (dopar.nn) {
    avai.cors <- parallel::detectCores() - 2
    cl <- parallel::makeCluster(avai.cors)
    doParallel::registerDoParallel(cl)
    Ri.nearestNeighbors.list <- vector("list", num.samp2)
    Ri.nearestNeighbors.list <-
      foreach::foreach(
        Ri.int = seq.int(num.samp2), .packages = c("dplyr")
      ) %dopar% {
        return(get_Ri_nearest(dist.df, as.character(Ri.int), nbd.method, Ri.radius = Ri.radius))
      }
    parallel::stopCluster(cl)
  } else {
    # put each Ri's nbd in a list then rbind them at the end with bind_rows()
    Ri.nearestNeighbors.list <- vector("list", num.samp2) # initialize list

    for (Ri in colnames(dist.df)) { # for each sample Ri
      Ri.int <- as.integer(Ri)
      Ri.nearest.idx <- get_Ri_nearest(dist.df, Ri, nbd.method, Ri.radius = Ri.radius)

      if (!is.null(Ri.nearest.idx)) { # similar to relieff
        Ri.nearestNeighbors.list[[Ri.int]] <- Ri.nearest.idx
      }
    }
  }

  # list of Ri's (data2 test instances) NN's in data1 (train instances)
  Ri.nearestNeighbors.list
}

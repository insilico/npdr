# TODO
# detectionStats
#' examples
#' functional <- case.control.3sets$signal.names
#' positives <- row.names(npdr.cc.results.df[npdr.cc.results.df[, 1] < .05, ]) # p.adj<.05
#' npdr.cc.detect.stats <- detectionStats(functional.case.control, positives)
#' cat(npdr.cc.detect.stats$summary.msg) # on(outcome="class", dataset=case.control.data)
#' 
#' #' npdrDetected
#'
#' Given a vector functional (true) attribute names, a vector of sorted attribute names, and percentile
#' threshold, returns true positive rate.
#'
#' @param results.df dataframe of sorted (low to high P value) attribute names from NPDR
#' @param functional character vector of functional/true attribute names
#' @param p percentile of top relief scores compared with the functional list
#' @return True positive rate: number of true postives divided by the number of functional
#' @examples
#' functional.vars <- dataset$signal.names
#' npdr.results1 <- npdr("class", dats,
#'   regression.type = "binomial",
#'   attr.diff.type = "allele-sharing", # nbd.method="relieff",
#'   nbd.method = "multisurf",
#'   nbd.metric = "manhattan", msurf.sd.frac = .5, k = 0,
#'   neighbor.sampling = "none", separate.hitmiss.nbds = F,
#'   dopar.nn = T, dopar.reg = T, padj.method = "bonferroni", verbose = T
#' )
#' npdrDetected(npdr.results1, functional.vars, p = .1)
#' 
#' #' geneLowVarianceFilter
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
#' 
#' # =========================================================================#
#' rfDetected
#'
#' Given a vector functional (true) attribute names, a vector of sorted attribute names, and percentile
#' threshold, returns true positive rate.
#'
#' @param results.df dataframe of sorted (high to low) attribute names from randomForest
#' @param functional character vector of functional/true attribute names
#' @param p percentile of top relief scores compared with the functional list
#' @return True positive rate: number of true postives divided by the number of functional
#' @examples
#' functional.vars <- dataset$signal.names
#' ranfor.fit <- randomForest(as.factor(class) ~ ., data = dats)
#' rf.importance <- importance(ranfor.fit)
#' rf.sorted <- sort(rf.importance, decreasing = T, index.return = T)
#' rf.df <- data.frame(att = rownames(rf.importance)[rf.sorted$ix], rf.scores = rf.sorted$x)
#' rfDetected(rf.df, functional.vars, p = .1)
#' 
#' # =========================================================================#
#' reliefDetected
#'
#' Given a vector functional (true) attribute names, a vector of sorted attribute names, and percentile
#' threshold, returns true positive rate.
#'
#' @param results.df dataframe of relief-sorted (high to low) attribute names from CORElearn
#' @param functional character vector of functional/true attribute names
#' @param p percentile of top relief scores compared with the functional list
#' @return True positive rate: number of true postives divided by the number of functional
#' @examples
#' functional.vars <- dataset$signal.names
#' relief <- CORElearn::attrEval(as.factor(class) ~ .,
#'   data = dats,
#'   estimator = "ReliefFequalK",
#'   costMatrix = NULL,
#'   outputNumericSplits = FALSE,
#'   kNearestEqual = floor(knnSURF(nrow(dats), .5) / 2)
#' ) # fn from npdr
#' relief.order <- order(relief, decreasing = T)
#' relief.df <- data.frame(att = names(relief)[relief.order], rrelief = relief[relief.order])
#' reliefDetected(relief.df, functional.vars, p = .1)
#' 
#' # =========================================================================#
#' uniReg
#'
#' Univariate logistic or linear regression for a dataset.
#'
#' @param outcome string with name of class column or outcome vector.
#' @param dataset data matrix with predictor columns and outcome column.
#' @param regression.type "lm" or "binomial"
#' @param padj.method for p.adjust (\code{"fdr"}, \code{"bonferroni"}, ...)
#' @param covars optional vector or matrix of covariate columns for correction. Or separate data matrix of covariates.
#' @return matrix of beta, p-value and adjusted p-value, sorted by p-value.
#' @examples
#' lr.results <- uniReg(outcome = "class", dataset = case.control.data, regression.type = "binomial")
#' #  lr.results[lr.results[,"p.adj"]<.05]
#' 
#' 
#' # =========================================================================#
#' knnVec
#'
#' Number of neighbors for each sample (vector) from a neighbor-pair matrix.
#'
#' @param neighbor.pairs.idx two columns of redundant "i,j" pairs from nearestNeighbors function
#' @return  knn.vec vector number of nearest neighbors for each instance
#'
#' @examples
#' mean(knnVec(neighbor.pairs.idx)) # average number of neighbors
#' # =========================================================================#
#' uniqueNeighbors
#'
#' Find pairs of unique nearest neighbors pairs from possible redundant pairs.
#' Used as options (neighbor.sampling="unique") in nearestNeighbors and npdr functions.
#'
#' @param neighbor.pairs.idx two columns of (possibly redundant) "i,j" pairs from nearestNeighbors function
#' @return new neighborhood pair matrix of only unique pairs
#'
#' @examples
#' unique.neighbor.pairs.idx <- uniqueNeighbors(neighbor.pairs.idx) 
#' 
#'  nearestNeighborsSeparateHitMiss
#'  
#' #' @examples
#' # reliefF (fixed-k) neighborhood using default k equal to theoretical surf expected value
#' # One can change the theoretical value by changing sd.frac (default 0.5)
#' neighbor.pairs.idx <- nearestNeighborsSeparateHitMiss(cc.attrs, cc.pheno, # need attributes and pheno
#'   nb.method = "relieff", nb.metric = "manhattan",
#'   sd.frac = .5, k = 0
#' )
#' 
#' #'# =========================================================================#
#' nearestNeighbors
#' @examples
#' # multisurf neighborhood with sigma/2 (sd.frac=0.5) "dead-band" boundary
#' neighbor.pairs.idx <- nearestNeighbors(predictors.mat, nb.method = "multisurf", nb.metric = "manhattan", sd.frac = 0.5)
#' # reliefF (fixed-k) neighborhood using default k equal to theoretical surf expected value
#' # One can change the theoretical value by changing sd.frac (default 0.5)
#' neighbor.pairs.idx <- nearestNeighbors(predictors.mat, nb.method = "relieff", nb.metric = "manhattan")
#' # reliefF (fixed-k) neighborhood with a user-specified k
#' neighbor.pairs.idx <- nearestNeighbors(predictors.mat, nb.method = "relieff", nb.metric = "manhattan", k = 10)
#' # =========================================================================#
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
#' 
#' 
#' add npdrLearner back in.
#' 
#' #' npdr
#'
#' @return npdr.stats.df: npdr fdr-corrected p-value for each attribute ($pval.adj [1]), raw p-value ($pval.attr [2]), and regression coefficient (beta.attr [3])
#'
#' @importFrom utils capture.output combn write.table
#' @importFrom foreach foreach `%dopar%`
#' @importFrom glmnet glmnet cv.glmnet
#' @import dplyr
#'
#' @examples
#' # Data interface options.
#' # Specify name ("qtrait") of outcome and dataset, which is a data frame including the outcome column.
#' # ReliefF fixed-k neighborhood, uses surf theoretical default (with msurf.sd.frac=.5) if you do not specify k or let k=0
#' npdr.results.df <- npdr("qtrait", case.control.3sets$train, regression.type = "lm", nbd.method = "relieff", nbd.metric = "manhattan", attr.diff.type = "manhattan", covar.diff.type = "manhattan", msurf.sd.frac = 0.5, padj.method = "bonferroni")
#'
#' # Specify column index (101) of outcome and dataset, which is a data frame including the outcome column.
#' #  # ReliefF fixed-k nbd, choose a k (knn=10). Or choose msurf.sd.frac
#' npdr.results.df <- npdr(101, case.control.3sets$train, regression.type = "lm", nbd.method = "relieff", nbd.metric = "manhattan", attr.diff.type = "manhattan", covar.diff.type = "manhattan", knn = 10, padj.method = "bonferroni")
#'
#' # if outcome vector (pheno.vec) is separate from attribute matrix
#' # multisurf
#' npdr.results.df <- npdr(pheno.vec, predictors.mat, regression.type = "lm", nbd.method = "multisurf", nbd.metric = "manhattan", attr.diff.type = "manhattan", covar.diff.type = "manhattan", msurf.sd.frac = 0.5, padj.method = "bonferroni")
#' # attributes with npdr adjusted p-value less than .05
#' npdr.positives <- row.names(npdr.results.df[npdr.results.df$pva.adj < .05, ]) # npdr p.adj<.05
#' @export
#' 
#' # ======================================================================================#
#' Split a data set for machine learning classification
#'
#' Return data.sets as a list of training set, holdout set and validation set
#' according to the predefined percentage of each partition
#' default is a 50-50 split into training and holdout, no testing set
#' code class/label/phenotypes as 1 and -1.
#' User can manage the simulation data to be dichotomious/quantitative using label (class/qtrait)
#'
#' @param all.data A data frame of n rows by d colums of data plus a label column
#' @param pct.train A numeric percentage of samples to use for traning
#' @param pct.holdout A numeric percentage of samples to use for holdout
#' @param pct.validation A numeric percentage of samples to use for testing
#' @param label A character vector of the data column name for the outcome label. class for classification
#' and qtrait for regression.
#' @return A list containing:
#' \describe{
#'   \item{train}{traing data set}
#'   \item{holdout}{holdout data set}
#'   \item{validation}{validation data set}
#' }
#' @examples
#' data("rsfMRIcorrMDD") # from privateEC
#' data.sets <- splitDataset(rsfMRIcorrMDD)
#' 
#' #' vwok
#' 
#' 
#' @examples
#' # Main effect simulation for standard m x p data set
#' num.samples <- 100
#' num.variables <- 100
#' pct.imbalance <- 0.5
#' pct.signals <- 0.1
#' main.bias <- 0.5
#' sim.type <- "mainEffect"
#'
#' dataset <- createSimulation2(num.samples=num.samples,
#'                              num.variables=num.variables,
#'                              pct.imbalance=pct.imbalance,
#'                              pct.signals=pct.signals,
#'                              main.bias=main.bias,
#'                              label="class",
#'                              sim.type=sim.type,
#'                              pct.mixed=0.5,
#'                              pct.train=0.5,
#'                              pct.holdout=0.5,
#'                              pct.validation=0,
#'                              verbose=T,
#'                              data.type="continuous")
#' dats <- rbind(dataset$train, dataset$holdout, dataset$validation)
#' dats <- dats[order(dats[,ncol(dats)]),]
#'
#' # run Variable-Wise Optimized k function
#' out <- vwok(dats=dats,
#'             k.grid=NULL,
#'             verbose=T,
#'             attr.diff.type="numeric-abs",
#'             label="class")
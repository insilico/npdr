# =========================================================================#
#' knnSURF
#'
#' Theoretical value for the number of expected neighbors for SURF or multiSURF
#'
#' @param m.samples number of samples in data.
#' @param sd.frac fraction of the standard deviation from the mean of all pairwise distances, dead-band. The default value used by the SURF and multiSURF algorithms is 1/2.
#' @return knn Number of neighbors.
#' @examples
#' k.surf <- knnSURF(200, .5)
#' @export
knnSURF <- function(m.samples, sd.frac = .5) {
  # error function
  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
  # theoretical SURF knn formulat
  knn <- floor((m.samples - 1) * (1 - erf(sd.frac / sqrt(2))) / 2)
  return(knn)
}

# =========================================================================#
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
#' @export
uniReg <- function(outcome, dataset, regression.type = "lm", padj.method = "fdr", covars = "none") {
  ## parse input
  if (length(outcome) == 1) {
    # e.g., outcome="qtrait" or outcome=101 (pheno col index) and data.set is data.frame including outcome variable
    pheno.vec <- dataset[, outcome] # get phenotype
    if (is.character(outcome)) { # example column name: outcome="qtrait"
      attr.mat <- dataset[, !(names(dataset) %in% outcome)] # drop the outcome/phenotype
    } else { # example column index: outcome=101
      attr.mat <- dataset[, -outcome] # drop the outcome/phenotype
    }
  } else { # user specifies a separate phenotype vector
    pheno.vec <- outcome # assume users provides a separate outcome data vector
    attr.mat <- dataset # assumes data.set only contains attributes/predictors
  }
  ## set up model
  if (regression.type == "lm") {
    if (length(covars) > 1) {
      # 2 below means attribte stats (1 would be the intercept)
      model.func <- function(x) {
        as.numeric(summary(lm(pheno.vec ~ attr.mat[, x] + covars))$coeff[2, ])
      }
    } else { # covar=="none"
      model.func <- function(x) {
        # if nrow(summary(lm(pheno.vec ~ attr.mat[,x]))$coeff) < 2 then there was an issue with the attribute
        fit <- summary(lm(pheno.vec ~ attr.mat[, x]))
        coeffs <- fit$coeff
        ## create summary stats
        if (nrow(coeffs) < 2) {
          # for example, a monomorphic SNP might result in attribute stats (row 2) not being created
          message("Regression failure. Continuing to next variable.\n")
          return(rep(NA, 4))
        } else {
          return(as.numeric(coeffs[2, ]))
        }
      } # end this model.func
    }
  } else { # "binomial" model
    if (length(covars) > 1) {
      # model.func <- function(x) {tidy(glm(pheno.vec ~ attr.mat[,x] + covars, family=binomial))[2,4:5]}
      model.func <- function(x) {
        fit <- summary(glm(pheno.vec ~ attr.mat[, x] + covars, family = binomial))
        coeffs <- fit$coeff
        ## create summary stats
        if (nrow(coeffs) < 2) {
          # for example, a monomorphic SNP might result in attribute stats (row 2) not being created
          message("Regression failure. Continuing to next variable.\n")
          return(rep(NA, 4))
        } else {
          return(as.numeric(coeffs[2, ]))
        }
      } # end this model.func
    } else { # covar=="none"
      # model.func <- function(x) {tidy(glm(pheno.vec ~ attr.mat[,x], family=binomial))[2,4:5]}
      model.func <- function(x) {
        fit <- summary(glm(pheno.vec ~ attr.mat[, x], family = binomial))
        coeffs <- fit$coeff
        ## create summary stats
        if (nrow(coeffs) < 2) {
          # for example, a monomorphic SNP might result in attribute stats (row 2) not being created
          message("Regression failure. Continuing to next variable.\n")
          return(rep(NA, 4))
        } else {
          return(as.numeric(coeffs[2, ]))
        }
      } # end this model.func
    }
  } # end else binomial
  # class.col <- which(colnames(dataset)==outcome)
  # predictor.cols <- which(colnames(dataset)!=outcome)
  num.attr <- ncol(attr.mat)
  if (is.null(num.attr)) { # if there is just one attribute
    attr.mat <- as.matrix(attr.mat)
    num.attr <- ncol(attr.mat) # num.attr <- 1
  }
  beta_pvals <- t(sapply(1:num.attr, model.func)) # stats for all predictors
  pvals_coerce <- as.numeric(unlist(beta_pvals[, 4])) # p.adjust - no dataframe input
  univariate.padj <- p.adjust(pvals_coerce, method = padj.method) # fdr
  univariate.padj <- as.numeric(format(univariate.padj, scientific = T, digits = 5))
  betas <- as.numeric(format(beta_pvals[, 1], scientific = F, digits = 5))
  betas.Z.att <- as.numeric(format(beta_pvals[, 3], scientific = F, digits = 5))
  pvals <- as.numeric(format(beta_pvals[, 4], scientific = T, digits = 5))
  beta_pvals <- cbind(betas, betas.Z.att, pvals, univariate.padj) # adjusted p-val column
  row.names(beta_pvals) <- colnames(attr.mat) # add predictor names
  # row.names(beta_pvals)<- colnames(dataset)[-class.col] # add predictor names
  beta_pvals_sorted <- beta_pvals[order(as.numeric(beta_pvals[, 3]), decreasing = F), ] # sort by pval
  if (num.attr == 1) { # special case of only 1 attribute
    names(beta_pvals_sorted) <- c("beta", "beta.Z.att", "pval", "p.adj")
  } else { # multiple attributes, typical
    colnames(beta_pvals_sorted) <- c("beta", "beta.Z.att", "pval", "p.adj")
  }
  return(beta_pvals_sorted)
}

# =========================================================================#
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
#' positives <- row.names(npdr.cc.results.df[npdr.cc.results.df[, 1] < .05, ]) # p.adj<.05
#' npdr.cc.detect.stats <- detectionStats(functional.case.control, positives)
#' cat(npdr.cc.detect.stats$summary.msg) # on(outcome="class", dataset=case.control.data)
#' @export
detectionStats <- function(functional, positives) {
  TP <- sum(positives %in% functional)
  FP <- sum((positives %in% functional) == F)
  FN <- length(functional) - TP
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  num.positives <- length(positives)
  TPR <- TP / num.positives # rate, aka power or sensitivity
  FPR <- FP / num.positives # rate
  # summary message
  report <- paste(
    "Given ", length(functional), " functional (true) attributes.\n",
    "Given ", length(positives), " selected (positive) attributes.\n",
    "True Positives: ", TP, " true out of ", length(positives), " positives. TP rate = ", TPR, ".\n",
    "False Positives: ", FP, " false out of ", length(positives), " positives. FP rate = ", FPR, ".\n",
    "Precision: ", precision, ".\n",
    "Recall: ", recall, ".\n",
    sep = ""
  )
  return(list(
    TP = TP, FP = FP, FN = FN, TPR = TPR, FPR = FPR,
    precision = precision, recall = recall, report = report
  ))
}

# =========================================================================#
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
#' @export
reliefDetected <- function(results.df, functional, top.pct) {
  top.num <- floor(top.pct * nrow(results.df))
  top.vars <- results.df %>%
    top_n(top.num, rrelief) %>%
    pull(att) # rrelief column of results.df
  power <- detectionStats(functional, top.vars)$TP # npdr:: fn, how many of top.pct are true
  ifelse(is.nan(power), 0, power) / length(functional) # if nan, return 0, normalize by num of functional
}

# =========================================================================#
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
#' @export
rfDetected <- function(results.df, functional, top.pct) {
  top.num <- floor(top.pct * nrow(results.df))
  top.vars <- results.df %>%
    top_n(top.num, rf.scores) %>%
    pull(att) # rf.scores column for rf results
  power <- detectionStats(functional, top.vars)$TP # npdr:: fn, how many of top.pct are true
  ifelse(is.nan(power), 0, power) / length(functional) # if nan, return 0, normalize by num of functional
}

# =========================================================================#
#' npdrDetected
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
#' @export
npdrDetected <- function(results.df, functional, top.pct) {
  top.num <- floor(top.pct * nrow(results.df))
  top.vars <- results.df %>%
    top_n(-top.num, pval.att) %>%
    pull(att) # pval.att is npdr specific
  power <- detectionStats(functional, top.vars)$TP # npdr:: fn, how many of top.pct are true
  ifelse(is.nan(power), 0, power) / length(functional) # if nan, return 0, normalize by num of functional
}

# =========================================================================#
#' geneLowVarianceFilter
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
#' @export
geneLowVarianceFilter <- function(dataMatrix, percentile = 0.5) {
  variances <- apply(as.matrix(dataMatrix), 2, var)
  threshold <- quantile(variances, c(percentile))
  # remove variable columns with lowest percentile variance
  mask <- apply(dataMatrix, 2, function(x) var(x) > threshold)
  fdata <- dataMatrix[, mask]
  # return the row mask and filtered data
  list(mask = mask, fdata = fdata)
}

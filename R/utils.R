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
#' knnSURF.balanced
#'
#' Theoretical value for the number of expected neighbors for SURF or multiSURF 
#' (fixed or adaptive radius) neighborhoods, but adjusted for imbalanced data. 
#' We use our theoretical formula with twice the size of the minority class 
#' as the input sample size.
#' @param class.vec vector of class labels to determine the minimum class size. 
#' @param sd.frac fraction of the standard deviation from the mean of all pairwise distances, dead-band. The default value used by the SURF and multiSURF algorithms is 1/2.
#' @return knn Theoretical number of neighbors.
#' @examples
#' k.surf.bal <- knnSURF.balanced(class.vec, .5)
#' @export
  knnSURF.balanced <- function(class.vec, sd.frac = .5) {
    class_table <- as.numeric(table(class.vec))
    table_len <- length(class_table)
    if (table_len > 2){    # if the class variable is numeric or 
       min.class.size <- length(class.vec)
    } else{
      min.class.size <- min(class_table)   
    }
  # theoretical SURF knn formula but using twice the minority class size
  knn <- npdr::knnSURF(2*min.class.size - 1, 0.5)
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
#' @export
#' @examples
#' out_univariate <- uniReg(
#'   outcome = "class",
#'   dataset = case.control.3sets$train,
#'   regression.type = "binomial"
#' )
#' head(out_univariate)
#'
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
#' functional <- c("var1", "var2", "var3")
#' detected <- c("var1", "var2", "var4")
#' detectionStats(functional, detected)
#' @export
detectionStats <- function(functional, positives) {
  TP <- sum(positives %in% functional)
  FP <- sum(!(positives %in% functional))
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
#' detected
#'
#' Given a vector functional (true) attribute names, a vector of
#' sorted attribute names, and percentile threshold, returns true positive rate.
#'
#' @param results.df dataframe of sorted attribute names (column `att`)
#' (e.g., feature important scores or low to high P value)
#' from either CORElearn, NPDR, or random forest.
#' @param functional character vector of functional/true attribute names
#' @param top.pct percentile of top relief scores compared with the functional list
#' @param sort_col column to sort importance score on.
#' e.g., `rrelief`, `rf.scores`, or `pval.att`.
#' @param get_min Boolean. Whether to reverse the sort, e.g., for p value.
#' Default to FALSE (i.e., sort by importance scores).
#' @return True positive rate: number of true postives divided by the number of functional
#' @export
#' @examples
#' out_npdr <- npdr("class", case.control.3sets$train)
#' functional_feats <- case.control.3sets$signal.names
#' detected(0.2, out_npdr, functional_feats, "pval.att", TRUE)
#'
#' \dontrun{
#' ranfor_fit <- randomForest(as.factor(class) ~ ., case.control.3sets$train)
#' rf_imp <- data.frame(importance(ranfor_fit)) %>%
#'   rownames2columns("att") 
#' detected(0.1, rf_imp, functional_feats, "MeanDecreaseGini")
#' }
#' 
detected <- function(top.pct, results.df, functional, sort_col, get_min = FALSE) {
  sort_col <- sym(sort_col)
  slice_fun <- if (get_min) slice_min else slice_max
  top.num <- floor(top.pct * nrow(results.df))
  top.vars <- results.df %>%
    data.frame() %>%
    slice_fun(!!sort_col, n = top.num) %>%
    pull(att)
  power <- detectionStats(functional, top.vars)$TP # how many of top.pct are true
  ifelse(is.nan(power), 0, power) / length(functional)
}

# =========================================================================#
#' geneLowVarianceFilter
#'
#' Low variance mask and filtered data for gene expression matrix.
#'
#' @param dataMatrix data matrix with predictors only, sample x gene
#' @param percentile percentile of low variance removed
#' @return mask and filtered data
geneLowVarianceFilter <- function(dataMatrix, percentile = 0.5) {
  variances <- apply(as.matrix(dataMatrix), 2, var)
  threshold <- quantile(variances, c(percentile))
  # remove variable columns with lowest percentile variance
  mask <- apply(dataMatrix, 2, function(x) var(x) > threshold)
  fdata <- dataMatrix[, mask]
  # return the row mask and filtered data
  list(mask = mask, fdata = fdata)
}

# =========================================================================#
#' Hamming distance for a binary matrix
#' https://johanndejong.wordpress.com/2015/10/02/faster-hamming-distance-in-r-2/
#'
#' @param X Original matrix.
#'
#' @return Distance matrix with the Hamming metric.
#' @export
hamming.binary <- function(X) {
  D <- t(1 - X) %*% X
  D + t(D)
}

#' =========================================================================#
#' Compute denominator of the diff formula
#' for each attribute x (column) in my.mat, max(x) - min(x)
#'
#' @param my.mat attribute matrix
#'
#' @return Numeric vectors representing the denominators in the diff formula
#'
attr.range <- function(my.mat) {
  apply(as.matrix(my.mat), 2, function(x) {
    max(x) - min(x)
  })
}

#' Convert rownames to column.
#' 
#' Adapted from https://github.com/tidyverse/tibble/blob/516449bbb0c76925b93b703b68d7979a53d7cdee/R/rownames.R.
#'
#' @param df Input dataframe.
#' @param var Name of new column.
#'
#' @return A data frame with a new column of row names.
rownames2columns <- function(df, var = "rowname"){
  df <- df %>% 
    mutate(!!var := rownames(df)) %>% 
    select(!!var, everything())
  rownames(df) <- NULL
  df
}

`%||%` <- function(a, b) if (is.null(a)) b else a

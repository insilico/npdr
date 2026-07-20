# ========================================================================= #
#' surv_npdr
#'
#' Survival NPDR: Nearest-Neighbor Projected-Distance Regression for
#' right-censored survival outcomes.
#'
#' Builds pairwise differences in a nearest-neighbor neighborhood and
#' regresses survival-time differences (\code{"lm"}) or binary pair-ordering
#' labels (\code{"binomial"}) on attribute differences. Two key choices
#' control the formulation:
#'
#' \strong{diff.type} — how differences are computed:
#' \itemize{
#'   \item \code{"absolute"} (symmetric): \eqn{|\Delta T| \sim |\Delta X|}.
#'     Symmetric in pair direction; best for detecting interaction effects.
#'   \item \code{"signed"} (asymmetric): signed \eqn{\Delta T \sim \Delta X}.
#'     Preserves direction; best for main-effect detection in AFT models.
#' }
#'
#' \strong{regression} — the regression model:
#' \itemize{
#'   \item \code{"lm"}: linear regression on pairwise \eqn{\Delta T}.
#'     Use \code{km.impute = TRUE} to replace censored times with their
#'     Kaplan-Meier conditional mean. Use \code{covariates} for adjustment.
#'   \item \code{"binomial"}: logistic regression on the binary pair-ordering
#'     label (neighbor outlives anchor = 1). Set \code{regularize = TRUE}
#'     for joint elastic-net selection via \pkg{glmnet}.
#' }
#'
#' @param outcome Named character vector specifying outcome column names:
#'   \code{c(time_var = "time", status_var = "status")}.
#' @param dataset Data.frame with predictor columns and the two outcome columns.
#' @param nbd.method Neighborhood method: \code{"multisurf"} for adaptive k
#'   (default) or \code{"relieff"} for fixed k (specify \code{knn}).
#' @param nbd.metric Distance metric for neighbor search. Default
#'   \code{"manhattan"}.
#' @param knn Number of nearest neighbors. \code{0} lets MultiSURF choose k
#'   adaptively. Default \code{0}.
#' @param msurf.sd.frac SD fraction for MultiSURF bandwidth pruning.
#'   Default \code{0.5}.
#' @param diff.type Difference formulation: \code{"absolute"} (symmetric,
#'   default) or \code{"signed"} (asymmetric).
#' @param regression Regression model applied to the pairwise design matrix:
#'   \code{"lm"} (default) or \code{"binomial"}.
#' @param km.impute Logical; replace censored survival times with their
#'   Kaplan-Meier conditional mean \eqn{\hat{E}[T \mid T > t_c]} before
#'   computing pairwise differences. Applies to both \code{"lm"} and
#'   \code{"binomial"} paths. Default \code{FALSE}.
#' @param nbd.filter Logical; drop pairs where the neighbor is censored at or
#'   before the anchor (ordering is ambiguous). Default \code{TRUE}.
#' @param covariates Optional character vector of column names to include as
#'   adjustment covariates, e.g. \code{c("age", "sex", "BMI")}. Signed
#'   pairwise differences of each covariate are added as regressors so that
#'   feature scores are adjusted. Only applies when \code{regression = "lm"}.
#'   Default \code{NULL}.
#' @param regularize Logical; fit a joint elastic-net logistic model via
#'   \pkg{glmnet} instead of per-feature GLMs. Only applies when
#'   \code{regression = "binomial"}. Default \code{FALSE}.
#' @param alpha Elastic-net mixing parameter: \code{1} = LASSO, \code{0} =
#'   ridge, values in \code{(0, 1)} = elastic net. Used when
#'   \code{regularize = TRUE}. Default \code{1}.
#'
#' @return Data.frame of features sorted by p-value (or \eqn{|\beta|} for
#'   regularized models). Columns: \code{feature}, \code{beta},
#'   \code{stat} (t or z), \code{pval}, \code{pval_adj}.
#'   Regularized binomial returns \code{feature} and \code{beta} only
#'   (non-zero coefficients).
#'
#' @examples
#' \dontrun{
#' # Main-effect detection — signed differences, adaptive neighborhood
#' res <- surv_npdr(my_data, diff.type = "signed")
#' head(res)
#'
#' # Interaction detection — absolute differences, fixed k = 20
#' res <- surv_npdr(my_data, diff.type = "absolute",
#'                  nbd.method = "relieff", knn = 20)
#'
#' # With KM imputation (works for both lm and binomial)
#' res <- surv_npdr(my_data, diff.type = "signed",
#'                  km.impute = TRUE, covariates = c("age", "sex"))
#'
#' # Binomial with KM imputation
#' res <- surv_npdr(my_data, regression = "binomial", km.impute = TRUE)
#'
#' # Logistic regression with LASSO regularization
#' res <- surv_npdr(my_data, regression = "binomial",
#'                  regularize = TRUE, alpha = 1)
#' }
#'
#' @importFrom dplyr mutate
#' @importFrom stats glm binomial lm na.omit p.adjust
#' @importFrom rlang .data
#' @export
surv_npdr <- function(outcome       = c(time_var = "time", status_var = "status"),
                      dataset,
                      nbd.method    = "multisurf",
                      nbd.metric    = "manhattan",
                      knn           = 0,
                      msurf.sd.frac = 0.5,
                      diff.type     = c("absolute", "signed"),
                      regression    = c("lm", "binomial"),
                      km.impute     = FALSE,
                      nbd.filter    = TRUE,
                      covariates    = NULL,
                      regularize    = FALSE,
                      alpha         = 1,
                      glmnet.lam    = "lambda.1se") {
  
  diff.type  <- match.arg(diff.type)
  regression <- match.arg(regression)
  
  time.var   <- as.character(outcome["time_var"])
  status.var <- as.character(outcome["status_var"])
  
  if (!is.data.frame(dataset)) dataset <- as.data.frame(dataset)
  
  if (!is.null(covariates)) {
    missing.covs <- setdiff(covariates, colnames(dataset))
    if (length(missing.covs) > 0)
      stop("covariates not found in dataset: ", paste(missing.covs, collapse = ", "))
    if (regression != "lm")
      stop("covariates only supported for regression = 'lm'")
  }
  
  ##### parse outcome columns; optionally impute censored survival times
  orig.time   <- dataset[[time.var]]
  orig.status <- dataset[[status.var]]
  
  if (km.impute)
    dataset[[time.var]] <- kmImputeTimes(orig.time, orig.status)
  
  ##### build attribute matrix (drop outcome and covariate columns)
  attr.mat   <- dataset[, !colnames(dataset) %in% c(time.var, status.var, covariates),
                        drop = FALSE]
  num.samp   <- nrow(attr.mat)
  num.attr   <- ncol(attr.mat)
  attr.names <- colnames(attr.mat)
  
  if (knn > num.samp - 1) {
    warning("knn too large; setting knn = num.samp - 1.")
    knn <- num.samp - 1
  }
  
  ##### find nearest neighbors
  neighbor.pairs.idx <- npdr::nearestNeighbors(
    attr.mat   = attr.mat,
    sd.frac    = msurf.sd.frac,
    k          = knn,
    nbd.method = nbd.method,
    nbd.metric = nbd.metric)
  
  # attach original survival times and status to each (anchor Ri, neighbor NN) pair
  neighbor.pairs.idx <- neighbor.pairs.idx |>
    dplyr::mutate(
      Ri.time   = orig.time[.data$Ri_idx],
      NN.time   = orig.time[.data$NN_idx],
      Ri.status = orig.status[.data$Ri_idx],
      NN.status = orig.status[.data$NN_idx])
  
  ##### censor net filter
  # ambiguous: neighbor is censored at or before anchor's time — true ordering unknown
  if (nbd.filter) {
    ambiguous.pairs <- which(neighbor.pairs.idx$NN.time <= neighbor.pairs.idx$Ri.time &
                               neighbor.pairs.idx$NN.status == 0)
    if (length(ambiguous.pairs) > 0)
      neighbor.pairs.idx <- neighbor.pairs.idx[-ambiguous.pairs, ]
  }
  
  ##### compute pairwise attribute differences
  # "signed"   (asymmetric): Ri_x - NN_x  — preserves direction; best for main effects
  # "absolute" (symmetric):  |Ri_x - NN_x| — direction-free; best for interactions
  signed.diffs <- (diff.type == "signed")
  
  attr.diff.mat <- matrix(0.0, nrow = nrow(neighbor.pairs.idx), ncol = num.attr,
                          dimnames = list(NULL, attr.names))
  for (i in seq.int(num.attr)) {
    x <- attr.mat[[i]]
    attr.diff.mat[, i] <- if (signed.diffs) {
      x[neighbor.pairs.idx$Ri_idx] - x[neighbor.pairs.idx$NN_idx]
    } else {
      abs(x[neighbor.pairs.idx$Ri_idx] - x[neighbor.pairs.idx$NN_idx])
    }
  }
  
  # ==========================================================================#
  # Path A — LM: regress phenotype diff (ΔT) on attribute diff (ΔX)
  # ==========================================================================#
  if (regression == "lm") {
    
    # survival time differences — uses imputed times when km.impute = TRUE
    Ri.time.vec    <- dataset[[time.var]][neighbor.pairs.idx$Ri_idx]
    NN.time.vec    <- dataset[[time.var]][neighbor.pairs.idx$NN_idx]
    pheno.diff.vec <- if (signed.diffs) Ri.time.vec - NN.time.vec else abs(Ri.time.vec - NN.time.vec)
    
    # optional covariate adjustment: signed ΔCovariate added as regressors to each model
    if (!is.null(covariates)) {
      covar.diff.mat <- do.call(cbind, lapply(covariates, function(cv) {
        v <- dataset[[cv]]
        v[neighbor.pairs.idx$Ri_idx] - v[neighbor.pairs.idx$NN_idx]
      }))
      colnames(covar.diff.mat) <- covariates
    } else {
      covar.diff.mat <- NULL
    }
    
    npdr.stats.list <- vector("list", num.attr)
    
    for (i in seq.int(num.attr)) {
      attr.diff.vec <- attr.diff.mat[, i]
      
      if (!is.null(covar.diff.mat)) {
        design.matrix.df <- data.frame(
          pheno.diff.vec = pheno.diff.vec,
          attr.diff.vec  = attr.diff.vec,
          covar.diff.mat)
        mod <- tryCatch(lm(pheno.diff.vec ~ ., data = design.matrix.df),
                        error = function(e) NULL)
      } else {
        design.matrix.df <- data.frame(
          pheno.diff.vec = pheno.diff.vec,
          attr.diff.vec  = attr.diff.vec)
        mod <- tryCatch(lm(pheno.diff.vec ~ attr.diff.vec, data = design.matrix.df),
                        error = function(e) NULL)
      }
      
      if (!is.null(mod)) {
        cs <- summary(mod)$coefficients
        if (nrow(cs) >= 2)
          npdr.stats.list[[i]] <- data.frame(
            feature = attr.names[i],
            beta    = cs[2, 1],
            stat    = cs[2, 3],
            pval    = cs[2, 4],
            stringsAsFactors = FALSE)
      }
    }
    
    npdr.stats.df <- do.call(rbind, Filter(Negate(is.null), npdr.stats.list))
    if (is.null(npdr.stats.df) || nrow(npdr.stats.df) == 0) return(data.frame())
    npdr.stats.df$pval_adj <- p.adjust(npdr.stats.df$pval, method = "bonferroni", n = num.attr)
    rownames(npdr.stats.df) <- NULL
    return(npdr.stats.df[order(npdr.stats.df$pval), ])
    
    # ==========================================================================#
    # Path B — Binomial: regress pair-ordering label on attribute diff (ΔX)
    # ==========================================================================#
  } else {
    
    # binary label: 1 if neighbor survives at least as long as anchor, else 0
    # uses imputed times when km.impute = TRUE, giving resolved orderings for censored pairs
    pair.times     <- dataset[[time.var]]
    pheno.diff.vec <- as.integer(pair.times[neighbor.pairs.idx$NN_idx] >= pair.times[neighbor.pairs.idx$Ri_idx])
    
    ##### regularized: fit a joint LASSO / elastic-net model via glmnet
    if (regularize) {
      mod <- tryCatch(
        glmnet::cv.glmnet(attr.diff.mat, pheno.diff.vec, alpha = alpha,
                          family = "binomial", type.measure = "class"),
        error = function(e) NULL)
      if (is.null(mod)) return(data.frame())
      coefs <- as.numeric(coef(mod, s = glmnet.lam))[-1]  # drop intercept
      npdr.stats.df <- data.frame(feature = attr.names, beta = coefs, stringsAsFactors = FALSE)
      # keep all features — non-zero first (by |β|), then zeros — so full ranking is preserved
      nonzero <- npdr.stats.df[npdr.stats.df$beta != 0, ]
      nonzero <- nonzero[order(abs(nonzero$beta), decreasing = TRUE), ]
      zeros   <- npdr.stats.df[npdr.stats.df$beta == 0, ]
      zeros   <- zeros[sample(nrow(zeros)), ]  # randomise — no ordering info for unselected features
      npdr.stats.df <- rbind(nonzero, zeros)
      rownames(npdr.stats.df) <- NULL
      return(npdr.stats.df)
    }
    
    ##### per-feature: one univariate binomial GLM per attribute
    attr.diff.df       <- as.data.frame(attr.diff.mat)
    attr.diff.df$class <- pheno.diff.vec
    npdr.stats.list    <- vector("list", num.attr)
    
    for (i in seq.int(num.attr)) {
      subset.data <- attr.diff.df[, c("class", attr.names[i]), drop = FALSE]
      subset.data <- na.omit(subset.data)
      mod <- tryCatch(
        glm(class ~ ., data = subset.data, family = binomial(link = "logit")),
        error = function(e) NULL)
      if (!is.null(mod)) {
        cs <- as.data.frame(summary(mod)$coefficients)
        if (nrow(cs) >= 2)
          npdr.stats.list[[i]] <- data.frame(
            feature = attr.names[i],
            beta    = cs[2, 1],
            stat    = cs[2, 3],
            pval    = cs[2, 4],
            stringsAsFactors = FALSE)
      }
    }
    
    npdr.stats.df <- do.call(rbind, Filter(Negate(is.null), npdr.stats.list))
    if (is.null(npdr.stats.df) || nrow(npdr.stats.df) == 0) return(data.frame())
    npdr.stats.df$pval_adj <- p.adjust(npdr.stats.df$pval, method = "bonferroni", n = num.attr)
    rownames(npdr.stats.df) <- NULL
    return(npdr.stats.df[order(npdr.stats.df$pval), ])
  }
}

# =========================================================================#
#' detectionStats
#'
#' Given a vector functional (true) attribute names and a vector of positive
#' association attribute names, returns detection statistics like recall and precision.
#'
#' @param functional character vector of functional/true attribute names.
#' @param positives character vector of attribute names of positive associations (null hypothesis rejected or some threshold).
#' @param num_all_vars total number of variables, needed for TN calculation for MCC
#' @return list with elements TP, FP, FN, TPR, FPR, precision, recall, F1, MCC, and summary message (string).
#' @examples
#' functional <- c("var1", "var2", "var3")
#' detected <- c("var1", "var2", "var4")
#' detectionStats(functional, detected)
#' @export
detectionStats <- function(functional, positives, num_all_vars) {
  TP <- sum(positives %in% functional)
  FP <- sum(!(positives %in% functional))
  FN <- length(functional) - TP
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  num.positives <- length(positives)
  TPR <- TP / num.positives # rate, aka power or sensitivity
  FPR <- FP / num.positives # rate
  TN <- num_all_vars - (FP + FN + TP)
  # F1 score
  # Calculation based on binary classification, 
  # like case/control or functional/non-functional feature.
  F1 <- 2 * ((precision * recall)/(precision + recall))
  # Matthews Correlation Coefficient (MCC)
  MCC <- (TP * TN - FP * FN)/sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  # summary message
  report <- paste(
    "Given ", length(functional), " functional (true) attributes.\n",
    "Given ", length(positives), " selected (positive) attributes.\n",
    "True Positives: ", TP, " true out of ", length(positives), " positives. TP rate = ", TPR, ".\n",
    "False Positives: ", FP, " false out of ", length(positives), " positives. FP rate = ", FPR, ".\n",
    "True Negatives: ", TN, ".\n",
    "Precision: ", precision, ".\n",
    "Recall: ", recall, ".\n",
    "F1: ", F1, ".\n",
    "MCC: ", MCC, ".\n",
    sep = ""
  )
  return(list(
    TP = TP, FP = FP, FN = FN, TPR = TPR, FPR = FPR,
    precision = precision, recall = recall, 
    F1=F1, MCC=MCC, 
    report = report
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
#' auPRC
#'
#' Area under the Precision-Recall Curve (AUPRC)
#' Wrapper for PRROC library
#' @examples
#' npdr.prc <- auPRC(npdr_results$beta.Z.att, npdr_results$att, 
#'                                            dataset$signal.names)
#' npdr.prc$auc
#' ggplot() +
#'  geom_line(data = npdr.prc$curve.df, aes(x = Recall, 
#'                                          y = Precision, color = "r")) 
#' curve: first column x-axis is recall = TP / (TP + FN)
#'        second column y-axis is precision = TP / (TP + FP)
#'        third column, not returned, is threshold
#'      : left to right, the threshold is increasing. So you start with
#'        all features selected/positive, and then all fnl variables are 
#'        positive (TP=1) and recall = 1/(1+0) (0 are declared false). 
#' xaxis: recall, fraction of true positives at a given threshold out of only
#'        the actual positives/functional
#' yaxis prec is fraction of positives that are acutally functional out of the  
#'       number threshold chose to be positive.
#'       Tends to decrease as all TPs are found and more FPs get added, and the
#'       final y value will eventually be num_functional/num_features.
#'
#' @param scores_all numeric vector, scores of all features, order is matched with varnames character vector of functional/true attribute names.
#' @param varnames_all character vector, all feature names with order matched with scores_all
#' @param varnames_fnl character vector of functional variable names
#' @return list of two elements: area under the PR curve (numeric) and dataframe of recall and precision values for threshold scan.
#' @export
auPRC <- function(scores_all,varnames_all, varnames_fnl){
  fnl_mask <- varnames_all %in% varnames_fnl
  
  scores_fnl <- scores_all[fnl_mask]
  scores_noise <- scores_all[!fnl_mask]
  
  pr.out <- PRROC::pr.curve(scores.class0 = as.numeric(scores_fnl),
                            scores.class1 = as.numeric(scores_noise), 
                            curve = T)
  curve.df <- data.frame(Recall=pr.out$curve[,1], Precision=pr.out$curve[,2])
  return(list(auc=pr.out$auc.integral, curve.df=curve.df))
}
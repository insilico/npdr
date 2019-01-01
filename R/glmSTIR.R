#=========================================================================#
#' diffRegression
#'
#' Wrapper for lm and glm to run regression for a phenotype diff vector, one attribute diff vector with optional covariate adjustment. Organizes regression statistics into a vector and then all attribute statistics combined in glmSTIR.
#'
#' @param design.matrix.df Desgin matrix with variables: pheno.diff.vec (outcome variable as vector of diffs), attr.diff.vec (one predictor varialbe as vector of diffs) and optional covariates (regressors of non-interest) vector diffs.   
#' @param regression.type (\code{"lm"}, \code{"glm"}) 
#' @return vector of regression stats to put into list for glmSTIR and combine into matrix
#'
#' @examples
#'
#' @export
#diffRegression <- function(pheno.diffs, predictor.diffs, regression.type="glm") {
diffRegression <- function(design.matrix.df, regression.type="glm") {
  # if there are no covariates then ~. model is pheno.diff.vec ~ attr.diff.vec
  # otherwise ~. model is pheno.diff.vec ~ attr.diff.vec + covariates
  if (regression.type=="lm"){
    fit <- summary(lm(pheno.diff.vec ~ ., data=design.matrix.df))
    stats.vec <- c(
      fit$coefficients[2,4], # p-value for attribute beta, pval.a
      fit$coefficients[2,3], # beta_hat_a, standardize beta for attribute, Ba
      #fit$fstatistic[1],     # F-stat and next is its p-value, F.stat
      #(1.0 - pf(fit$fstatistic[1], fit$fstatistic[2], fit$fstatistic[3])), # Fstat.pval
      fit$coefficients[1,3], # beta_hat_0, intercept, B0
      fit$coefficients[1,4], # p-value for intercept, B0.pval
      fit$r.squared         # R^2 of fit, R.sqr
  ) 
  } else{ #regression.type=="glm"
    fit <- summary(glm(pheno.diff.vec ~ ., family=binomial(link=logit), data=design.matrix.df))
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
#' glmSTIR
#'
#' generalized linear model (GLM) based STatistical Inference Relief (STIR)
#' Computes Relief-based attribute statistical signficance for case/control or quantitative outcomes.
#' Allows for categorical (SNP) or numeric (expession) predictor data types. 
#' Allows for covariate correction.   
#'
#' @param outcome character name or length-m numeric outcome vector for linear regression, factor for logistic regression 
#' @param data.set m x p matrix of m instances and p attributes, May also include outcome vector but then outcome should be name. Include attr names as colnames. 
#' @param regression.type (\code{"lm"} or \code{"glm"})
#' @param attr.diff.type diff type for attributes (\code{"numeric-abs"} or \code{"numeric-sqr"} for numeric, \code{"allele-sharing"} or \code{"match-mismatch"} for SNP). Phenotype diff uses same numeric diff as attr.diff.type when lm regression. For glm, phenotype diff is \code{"match-mismatch"}. 
#' @param nbd.method neighborhood method [\code{"multisurf"} or \code{"surf"} (no k) or \code{"relieff"} (specify k)]. Used by nearestNeighbors().
#' @param nbd.metric used in stirDistances for distance matrix between instances, default: \code{"manhattan"} (numeric). Used by nearestNeighbors().
#' @param k number of constant nearest hits/misses for \code{"relieff"} (fixed-k). Used by nearestNeighbors().
#' The default k=0 means use the expected SURF theoretical k with sd.frac (.5 by default) 
#' @param sd.frac multiplier of the standard deviation from the mean distances; subtracted from mean for SURF or multiSURF.
#' The multiSURF default is sd.frac=0.5: mean - sd/2. Used by nearestNeighbors(). 
#' @param covars optional vector or matrix of covariate columns for correction. Or separate data matrix of covariates.
#' @param covar.diff.type string (or string vector) specifying diff type(s) for covariate(s) (\code{"numeric-abs"} for numeric or \code{"match-mismatch"} for categorical). 
#' @param fdr.method for p.adjust (\code{"fdr"}, \code{"bonferroni"}, ...) 
#' @return glmSTIR.stats.df: glmSTIR fdr-corrected p-value for each attribute ($pval.adj [1]), raw p-value ($pval.attr [2]), and regression coefficient (beta.attr [3]) 
#'
#' @examples
#' # Data interface options.
#' # Specify name ("qtrait") of outcome and data.set, which is a data frame including the outcome column.
#' # ReliefF fixed-k neighborhood, uses surf theoretical default (with sd.frac=.5) if you do not specify k or let k=0
#' glmstir.results.df <- glmSTIR("qtrait", train.data, regression.type="lm", nbd.method="relieff", nbd.metric = "manhattan", attr.diff.type="manhattan", covar.diff.type="manhattan", sd.frac=0.5, fdr.method="bonferroni")
#'
#' # Specify column index (101) of outcome and data.set, which is a data frame including the outcome column.
#  # ReliefF fixed-k nbd, choose a k (k=10). Or choose sd.frac
#' glmstir.results.df <- glmSTIR(101, train.data, regression.type="lm", nbd.method="relieff", nbd.metric = "manhattan", attr.diff.type="manhattan", covar.diff.type="manhattan", k=10, fdr.method="bonferroni")
#'
#' # if outcome vector (pheno.vec) is separate from attribute matrix
#' # multisurf
#' glmstir.results.df <- glmSTIR(pheno.vec, predictors.mat, regression.type="lm", nbd.method="multisurf", nbd.metric = "manhattan", attr.diff.type="manhattan", covar.diff.type="manhattan", sd.frac=0.5, fdr.method="bonferroni")
#' # attributes with glmSTIR adjusted p-value less than .05 
#' glmstir.positives <- row.names(glmstir.results.df[glmstir.results.df$pva.adj<.05,]) # glmSTIR p.adj<.05
#' @export
glmSTIR <- function(outcome, data.set, regression.type="glm", attr.diff.type="numeric-abs",
                    nbd.method="multisurf", nbd.metric = "manhattan", k=0, sd.frac=0.5, 
                    covars="none", covar.diff.type="match-mismatch", 
                    fdr.method="bonferroni", verbose=FALSE){
  ##### parse the commandline 
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
  neighbor.pairs.idx <- nearestNeighbors(attr.mat, nbd.method=nbd.method, nbd.metric = nbd.metric, 
                                         sd.vec = NULL, sd.frac = sd.frac, k=k)
  num.neighbor.pairs <- nrow(neighbor.pairs.idx)
  if (verbose){
    cat(num.neighbor.pairs, " neighbor pairs. ", num.neighbor.pairs/num.samp, " average neighbors per instance\n")
    }
  ##### run glmSTIR, each attribute is a list, then we do.call rbind to a matrix
  glmSTIR.stats.list <- vector("list",num.samp) # initialize
  for (attr.idx in seq(1, num.attr)){
    attr.vals <- attr.mat[, attr.idx]
    Ri.attr.vals <- attr.vals[neighbor.pairs.idx[,1]]
    NN.attr.vals <- attr.vals[neighbor.pairs.idx[,2]]
    attr.diff.vec <- stirDiff(Ri.attr.vals, NN.attr.vals, diff.type=attr.diff.type)
    ### pheno diff vector
    # create pheno diff vector for linear regression (numeric)  
    if (regression.type=="lm"){
    Ri.pheno.vals <- pheno.vec[neighbor.pairs.idx[,1]]
    NN.pheno.vals <- pheno.vec[neighbor.pairs.idx[,2]]
    pheno.diff.vec <- stirDiff(Ri.pheno.vals, NN.pheno.vals, diff.type="numeric-abs")
    # model data.frame to go into lm
    design.matrix.df <- data.frame(attr.diff.vec=attr.diff.vec,pheno.diff.vec=pheno.diff.vec)
    } else { #regression.type=="glm"
      # create pheno diff vector for logistic regression (match-mismatch or hit-miss)  
      Ri.pheno.vals <- pheno.vec[neighbor.pairs.idx[,1]]
      NN.pheno.vals <- pheno.vec[neighbor.pairs.idx[,2]]
      pheno.diff.vec <- stirDiff(Ri.pheno.vals, NN.pheno.vals, diff.type="match-mismatch")
      pheno.diff.vec <- as.factor(ifelse(pheno.diff.vec=="TRUE",1,0))
      # model data.frame to go into glm
      design.matrix.df <- data.frame(attr.diff.vec=attr.diff.vec,pheno.diff.vec=pheno.diff.vec)
    }
    ### diff vector for each covariate
    # optional covariates to model
    if (length(covars)>1){ # if covars is a vector or matrix 
      # default value is covar="none" (no covariates) which has length 1
      covars <- as.matrix(covars)  # if covars is just one vector, make sure it's a 1-column matrix
      # covar.diff.type can be a vector of strings because each column of covars may be a different data type
      #covar.diff.list <- vector("list",length(covar.diff.type)) # initialize
      for (covar.col in (1:length(covar.diff.type))){
        covar.vals <- covars[, covar.col]
        Ri.covar.vals <- covar.vals[neighbor.pairs.idx[,1]]
        NN.covar.vals <- covar.vals[neighbor.pairs.idx[,2]]
        covar.diff.vec <- stirDiff(Ri.covar.vals, NN.covar.vals, diff.type=covar.diff.type[covar.col])
        #covar.diff.list[[covar.col]] <- covar.diff.vec
        # add covar diff vector to data.frame
        # these covars will be included in each attribute's model
        if (is.null(colnames(covars)[covar.col])){  # if covar vector has no column name, give it one
          covar.name <- paste("cov",covar.col,sep="") # cov1, etc.
        } else{
          covar.name <- colnames(covars)[covar.col] # else get the name from covars
        }
        design.matrix.df$temp <- covar.diff.vec  # add the diff covar to the design matrix data frame
        colnames(design.matrix.df)[2+covar.col] <- covar.name # change variable name
      }
      #covar.diff.mat <- do.call(cbind, covar.diff.list)  # all diff covariates as cols of a matrix
    }
    # utility function: RUN regression
    #                                                design.matrix.df = pheno.diff ~ attr.diff + option covar.diff
    glmSTIR.stats.list[[attr.idx]] <- diffRegression(design.matrix.df, regression.type=regression.type)
  }
  if (verbose){cat("Size of design matrices (phenotype + attribute + covariates, does not count intercept): ")
               cat(nrow(design.matrix.df)," diff-pairs by ", ncol(design.matrix.df)," variables.\n", sep="")
  }
  # combine lists into matrix
  glmSTIR.stats.attr_ordered.mat <- do.call(rbind, glmSTIR.stats.list)
  
  # rownames
  if (!is.null(colnames(attr.mat))){
    # add attribute names to stats/results matrix if the data matrix contains them
    rownames(glmSTIR.stats.attr_ordered.mat) <- colnames(attr.mat)
  } else {
    message("If you have attribute names, add them to colnames of input data.")
  }
  
  # attribute p-values
  attr.pvals <- glmSTIR.stats.attr_ordered.mat[, 1]
  # order-index for sorted attribute-beta p-values
  attr.pvals.order.idx <- order(attr.pvals, decreasing = F)
  # adjust p-values using Benjamini-Hochberg (default)
  attr.pvals.adj <- p.adjust(attr.pvals[attr.pvals.order.idx], method=fdr.method)
  
  # order by attribute p-value
  glmSTIR.stats.pval_ordered.mat <- glmSTIR.stats.attr_ordered.mat[attr.pvals.order.idx, ]
  # prepend adjused attribute p-values to first column
  glmSTIR.stats.pval_ordered.mat <- cbind(attr.pvals.adj, glmSTIR.stats.pval_ordered.mat)
  if (regression.type=="lm"){# stats colnames for lm
    colnames(glmSTIR.stats.pval_ordered.mat) <- c("pval.adj", "pval.attr", "beta.attr",  
                                                  "beta.0", "pval.0", "R.sqr")
  } else{ # stats columns for glm
    colnames(glmSTIR.stats.pval_ordered.mat) <- c("pval.adj", "pval.attr", "beta.attr", "beta.0", "pval.0")
  }
  # dataframe it
  glmSTIR.stats.df <- data.frame(glmSTIR.stats.pval_ordered.mat)
  return(glmSTIR.stats.df)
}

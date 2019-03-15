#=========================================================================#
#' diffRegression
#'
#' Wrapper for lm and glm-binomial to run regression for a phenotype diff vector, one attribute diff vector with optional covariate adjustment. Organizes regression statistics into a vector and then all attribute statistics combined in npdr.
#'
#' @param design.matrix.df Desgin matrix with variables: pheno.diff.vec (outcome variable as vector of diffs), attr.diff.vec (one predictor varialbe as vector of diffs) and optional covariates (regressors of non-interest) vector diffs.   
#' @param regression.type (\code{"lm"}, \code{"binomial"}) 
#' @return vector of regression stats to put into list for npdr and combine into matrix
#'
#' @examples
#'
#' @export
# regression of the neighbor diff vector for one attribute
diffRegression <- function(design.matrix.df, regression.type="binomial") {
  # if there are no covariates then ~. model is pheno.diff.vec ~ attr.diff.vec
  # otherwise ~. model is pheno.diff.vec ~ attr.diff.vec + covariates
  if (regression.type=="lm"){
    mod <- lm(pheno.diff.vec ~ ., data=design.matrix.df)
    fit <- summary(mod)
    beta_a <- coef(fit)[2, 1]         # raw beta coefficient, slope (not standardized)
    beta_zscore_a <- coef(fit)[2, 3]  # standardized beta coefficient (col 3)
    ## use one-side p-value to test H1: beta>0 for case-control and continuous outcome
    pval_beta_a <- pt(beta_zscore_a, mod$df.residual, lower = FALSE)  # one-sided p-val
    stats.vec <- c(
      #fit$coefficients[2,4], # p-value for attribute beta, pval.a
      #fit$coefficients[2,3], # beta_hat_a, standardize beta for attribute, Ba
      pval_beta_a,            # one-sided p-value for attribute beta
      beta_a,                 # beta_a for attribute a
      beta_zscore_a,          # standardized beta for attribut a
      fit$coefficients[1,1],  # beta_0, intercept, row 1 is inercept, col 1 is raw beta
      fit$coefficients[1,4],  # p-value for intercept, row 1 is intercept, col 4 is p-val 
      fit$r.squared           # R^2 of fit, R.sqr
    )
  } else{ #regression.type=="binomial"
    mod <- glm(pheno.diff.vec ~ ., family=binomial(link=logit), data=design.matrix.df)
    fit <- summary(mod)
    beta_a <- coef(fit)[2, 1]         # raw beta coefficient, slope (not standardized)
    beta_zscore_a <- coef(fit)[2, 3]  # standardized beta coefficient (col 3)
    ## use one-side p-value to test H1: beta>0 for case-control npdr scores
    pval_beta_a <- pt(beta_zscore_a, mod$df.residual, lower = FALSE)  # one-sided p-val
    stats.vec <- c(
      #fit$coefficients[2,4], # p-value for attribute beta, pval.a
      #fit$coefficients[2,3], # beta_hat_a, standardize beta for attribute, Ba
      pval_beta_a,            # one-sided p-value for attribute beta
      beta_a,                 # beta_a for attribute a
      beta_zscore_a,          # standardized beta for attribut a
      fit$coefficients[1,1],  # beta_0, intercept, row 1 is inercept, col 1 is raw beta
      fit$coefficients[1,4]   # p-value for intercept, row 1 is intercept, col 4 is p-val 
    )
  }
  return(stats.vec)
}

#=========================================================================#
#' npdr
#'
#' Nearest-Neighbor Projected-Distance Regression (npdr) 
#' generalized linear model (GLM) extension of STatistical Inference Relief (STIR)
#' Computes attribute statistical signficance with logistic for case/control and linear model for quantitative outcomes.
#' NPDR allows for categorical (SNP) or numeric (expession) predictor data types. 
#' NPDR allows for covariate correction.
#' Observations in the model are projected-distance differences between neighbors. 
#'
#' @param outcome character name or length-m numeric outcome vector for linear regression, factor for logistic regression 
#' @param dataset m x p matrix of m instances and p attributes, May also include outcome vector but then outcome should be name. Include attr names as colnames. 
#' @param regression.type (\code{"lm"} or \code{"binomial"})
#' @param attr.diff.type diff type for attributes (\code{"numeric-abs"} or \code{"numeric-sqr"} for numeric, \code{"allele-sharing"} or \code{"match-mismatch"} for SNP). Phenotype diff uses same numeric diff as attr.diff.type when lm regression. For glm-binomial, phenotype diff is \code{"match-mismatch"}. 
#' @param nbd.method neighborhood method [\code{"multisurf"} or \code{"surf"} (no k) or \code{"relieff"} (specify k)]. Used by nearestNeighbors().
#' @param nbd.metric used in npdrDistances for distance matrix between instances, default: \code{"manhattan"} (numeric). Used by nearestNeighbors().
#' @param knn number of constant nearest hits/misses for \code{"relieff"} (fixed-k). Used by nearestNeighbors().
#' The default knn=0 means use the expected SURF theoretical k with msurf.sd.frac (.5 by default) 
#' @param msurf.sd.frac multiplier of the standard deviation from the mean distances; subtracted from mean for SURF or multiSURF.
#' The multiSURF default is msurf.sd.frac=0.5: mean - sd/2. Used by nearestNeighbors(). 
#' @param covars optional vector or matrix of covariate columns for correction. Or separate data matrix of covariates.
#' @param covar.diff.type string (or string vector) specifying diff type(s) for covariate(s) (\code{"numeric-abs"} for numeric or \code{"match-mismatch"} for categorical).
#' @param glmnet.alpha penalty mixture for npdrNET: default alpha=1 (lasso, L1) alpha=0 (ridge, L2) 
#' @param glmnet.lower lower limit for coefficients for npdrNET: lower.limits=0 npdrNET default 
#' @param glment.family "binomial" for logistic regression, "gaussian" for regression
#' @param rm.attr.from.dist attributes for removal (possible confounders) from the distance matrix calculation. Argument for nearestNeighbors. None by default c().
#' @param neighbor.sampling "none" or \code{"unique"} if you want to use only unique neighbor pairs (used in nearestNeighbors)
#' @param padj.method for p.adjust (\code{"fdr"}, \code{"bonferroni"}, ...) 
#' @return npdr.stats.df: npdr fdr-corrected p-value for each attribute ($pval.adj [1]), raw p-value ($pval.attr [2]), and regression coefficient (beta.attr [3]) 
#'
#' @examples
#' # Data interface options.
#' # Specify name ("qtrait") of outcome and dataset, which is a data frame including the outcome column.
#' # ReliefF fixed-k neighborhood, uses surf theoretical default (with msurf.sd.frac=.5) if you do not specify k or let k=0
#' npdr.results.df <- npdr("qtrait", train.data, regression.type="lm", nbd.method="relieff", nbd.metric = "manhattan", attr.diff.type="manhattan", covar.diff.type="manhattan", msurf.sd.frac=0.5, padj.method="bonferroni")
#'
#' # Specify column index (101) of outcome and dataset, which is a data frame including the outcome column.
#  # ReliefF fixed-k nbd, choose a k (knn=10). Or choose msurf.sd.frac
#' npdr.results.df <- npdr(101, train.data, regression.type="lm", nbd.method="relieff", nbd.metric = "manhattan", attr.diff.type="manhattan", covar.diff.type="manhattan", knn=10, padj.method="bonferroni")
#'
#' # if outcome vector (pheno.vec) is separate from attribute matrix
#' # multisurf
#' npdr.results.df <- npdr(pheno.vec, predictors.mat, regression.type="lm", nbd.method="multisurf", nbd.metric = "manhattan", attr.diff.type="manhattan", covar.diff.type="manhattan", msurf.sd.frac=0.5, padj.method="bonferroni")
#' # attributes with npdr adjusted p-value less than .05 
#' npdr.positives <- row.names(npdr.results.df[npdr.results.df$pva.adj<.05,]) # npdr p.adj<.05
#' @export
npdr <- function(outcome, dataset, regression.type="binomial", attr.diff.type="numeric-abs",
                    nbd.method="multisurf", nbd.metric = "manhattan", knn=0, msurf.sd.frac=0.5, 
                    covars="none", covar.diff.type="match-mismatch",
                    glmnet.alpha=1, glmnet.lower=0, glmnet.family="binomial", 
                    rm.attr.from.dist=c(), neighbor.sampling="none",
                    padj.method="bonferroni", verbose=FALSE){
  ##### parse the commandline 
  if (length(outcome)==1){
    # e.g., outcome="qtrait" or outcome=101 (pheno col index) and dataset is data.frame including outcome variable
    pheno.vec <- dataset[,outcome] # get phenotype
    if (is.character(outcome)){ # example column name: outcome="qtrait"
      attr.mat <- dataset[ , !(names(dataset) %in% outcome)]  # drop the outcome/phenotype
    } else { # example column index: outcome=101
      attr.mat <- dataset[ , -outcome]  # drop the outcome/phenotype  
    }
  } else { # user specifies a separate phenotype vector
    pheno.vec <- outcome # assume users provides a separate outcome data vector
    attr.mat <- dataset # assumes dataset only contains attributes/predictors
  }
  rm(dataset)  # cleanup memory
  
  num.attr <- ncol(attr.mat)
  num.samp <- nrow(attr.mat)
  
  ##### get Neighbors (no phenotype used)
  # nbd.method (relieff, multisurf...), nbd.metric (manhattan...), k (for relieff nbd, theoerical surf default) 
  # msurf.sd.frac used by surf/multisurf relieff for theoretical k
  
  start_time <- Sys.time()                 
  neighbor.pairs.idx <- nearestNeighbors(attr.mat, nb.method=nbd.method, nb.metric=nbd.metric, 
                                         sd.frac = msurf.sd.frac, k=knn,
                                         attr_removal_vec_from_dist_calc=rm.attr.from.dist)
  num.neighbor.pairs <- nrow(neighbor.pairs.idx)
  if (neighbor.sampling=="unique"){
      # if you only want to return unique neighbors
      neighbor.pairs.idx <- uniqueNeighbors(neighbor.pairs.idx)
  }
  end_time <- Sys.time()
  num.neighbor.pairs <- nrow(neighbor.pairs.idx)
  if (verbose){
    cat("Neighborhood calculation time. "); difftime(end_time, start_time); cat("\n",sep="")
    cat(num.neighbor.pairs, "total neighbor pairs.", num.neighbor.pairs/num.samp, "average neighbors per instance.\n")
    erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
    # theoretical surf k (sd.frac=.5) for regression problems (does not depend on a hit/miss group)
    k.msurf.theory <- floor((num.samp-1)*(1-erf(msurf.sd.frac/sqrt(2)))/2)
    cat("Theoretical predicted multiSURF average neighbors: ", k.msurf.theory,".\n",sep="")
    if (neighbor.sampling=="unique"){
      # if you only want to return unique neighbors
      num.neighbor.pairs <- nrow(neighbor.pairs.idx)
      cat(num.neighbor.pairs, "unique neighbor pairs.", num.neighbor.pairs/num.samp, "average neighbors (unique) per instance.\n")
    }
  }
  ### pheno diff vector for glm-binomial or lm to use in each attribute's diff regression in for loop.
  # Not needed in loop.
  # create pheno diff vector for linear regression (numeric)  
  if (regression.type=="lm"){
    Ri.pheno.vals <- pheno.vec[neighbor.pairs.idx[,1]]
    NN.pheno.vals <- pheno.vec[neighbor.pairs.idx[,2]]
    pheno.diff.vec <- npdrDiff(Ri.pheno.vals, NN.pheno.vals, diff.type="numeric-abs")
  } else { #regression.type=="binomial"
    # create pheno diff vector for logistic regression (match-mismatch or hit-miss)  
    Ri.pheno.vals <- pheno.vec[neighbor.pairs.idx[,1]]
    NN.pheno.vals <- pheno.vec[neighbor.pairs.idx[,2]]
    pheno.diff.vec <- npdrDiff(Ri.pheno.vals, NN.pheno.vals, diff.type="match-mismatch")
    # the reference group is the hit group, so the logistic probability is prob of a pair being a miss
    pheno.diff.vec <- as.factor(pheno.diff.vec)
  }
  ##### run npdr, each attribute is a list, then we do.call rbind to a matrix
  npdr.stats.list <- vector("list",num.samp) # initialize
  for (attr.idx in seq(1, num.attr)){
    attr.vals <- attr.mat[, attr.idx]
    Ri.attr.vals <- attr.vals[neighbor.pairs.idx[,1]]
    NN.attr.vals <- attr.vals[neighbor.pairs.idx[,2]]
    attr.diff.vec <- npdrDiff(Ri.attr.vals, NN.attr.vals, diff.type=attr.diff.type)
    # model data.frame to go into lm or glm-binomial
    design.matrix.df <- data.frame(attr.diff.vec=attr.diff.vec,pheno.diff.vec=pheno.diff.vec)
    ### diff vector for each covariate
    # optional covariates to add to design.matrix.df model
    if (length(covars)>1){ # if covars is a vector or matrix
      if (regression.type=="glment"){
        message("penalized npdrNET does not currently support covariates.")
      }
      # default value is covar="none" (no covariates) which has length 1
      covars <- as.matrix(covars)  # if covars is just one vector, make sure it's a 1-column matrix
      # covar.diff.type can be a vector of strings because each column of covars may be a different data type
      #covar.diff.list <- vector("list",length(covar.diff.type)) # initialize
      for (covar.col in (1:length(covar.diff.type))){
        covar.vals <- covars[, covar.col]
        Ri.covar.vals <- covar.vals[neighbor.pairs.idx[,1]]
        NN.covar.vals <- covar.vals[neighbor.pairs.idx[,2]]
        covar.diff.vec <- npdrDiff(Ri.covar.vals, NN.covar.vals, diff.type=covar.diff.type[covar.col])
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
    npdr.stats.list[[attr.idx]] <- diffRegression(design.matrix.df, regression.type=regression.type)
  } # end of for loop, regression done for each attribute
  if (verbose){cat("Size of design matrices (phenotype + attribute + covariates, including intercept): ")
               cat(nrow(design.matrix.df)," diff-pairs by ", ncol(design.matrix.df)," columns.\n", sep="")
  }
  if (regression.type!="glmnet"){ # sort and format output if you did regular npdr
    # combine non-glmnet result lists into a matrix
    npdr.stats.attr_ordered.mat <- do.call(rbind, npdr.stats.list)
    # rownames
    if (!is.null(colnames(attr.mat))){
      # add attribute names to stats/results matrix if the data matrix contains them
      rownames(npdr.stats.attr_ordered.mat) <- colnames(attr.mat)
    } else {
      message("If you have attribute names, add them to colnames of input data.")
    }
    
    # attribute p-values
    # relies on first column [, 1] being the p-value for now
    # later first columns becomes att and adjusted p-value
    attr.pvals <- npdr.stats.attr_ordered.mat[, 1]
    # order-index for sorted attribute-beta p-values
    attr.pvals.order.idx <- order(attr.pvals, decreasing = F)
    # adjust p-values using Benjamini-Hochberg (default)
    attr.pvals.adj <- p.adjust(attr.pvals[attr.pvals.order.idx], method=padj.method)
    
    # order by attribute p-value
    npdr.stats.pval_ordered.mat <- npdr.stats.attr_ordered.mat[attr.pvals.order.idx, ]
    # prepend adjused attribute p-values to first column
    npdr.stats.pval_ordered.mat <- cbind(attr.pvals.adj, npdr.stats.pval_ordered.mat)
    # prepend attribute column (att)
    attributes.col <- as.character(rownames(npdr.stats.pval_ordered.mat))
    npdr.stats.pval_ordered.mat <- cbind(data.frame(att=attributes.col, stringsAsFactors=FALSE), 
                                         data.frame(npdr.stats.pval_ordered.mat, row.names=NULL))
    if (regression.type=="lm"){# stats colnames for lm
      colnames(npdr.stats.pval_ordered.mat) <- c("att", "pval.adj", "pval.att", "beta.raw.att", "beta.Z.att",  
                                                    "beta.0", "pval.0", "R.sqr")
    } else{ # stats columns for glm-binomial
     colnames(npdr.stats.pval_ordered.mat) <- c("att", "pval.adj", "pval.att", "beta.raw.att", "beta.Z.att", "beta.0", "pval.0")
    }
    # dataframe final output for regular npdr
    npdr.stats.df <- data.frame(npdr.stats.pval_ordered.mat)
    
  } else{ # Here we add an option npdrNET regression.type="glmnet"
          # Need to create a data matrix with each column as a vector of diffs for each attribute.
          # Need matrix because npdrNET operates on all attributes at once.
    attr.diff.mat <- matrix(0,nrow=nrow(neighbor.pairs.idx),ncol=num.attr)
    for (attr.idx in seq(1, num.attr)){
      attr.vals <- attr.mat[, attr.idx]
      Ri.attr.vals <- attr.vals[neighbor.pairs.idx[,1]]
      NN.attr.vals <- attr.vals[neighbor.pairs.idx[,2]]
      attr.diff.vec <- npdrDiff(Ri.attr.vals, NN.attr.vals, diff.type=attr.diff.type)
      attr.diff.mat[,attr.idx] <- attr.diff.vec
    }
    #
    Ri.pheno.vals <- pheno.vec[neighbor.pairs.idx[,1]]
    NN.pheno.vals <- pheno.vec[neighbor.pairs.idx[,2]]
    if (glmnet.family=="binomial"){
      pheno.diff.vec <- npdrDiff(Ri.pheno.vals, NN.pheno.vals, diff.type="match-mismatch")
      pheno.diff.vec <- as.factor(pheno.diff.vec)
      # Run glmnet on the diff attribute columns
      npdrNET.model<-cv.glmnet(attr.diff.mat, pheno.diff.vec,alpha=glmnet.alpha,family="binomial",
                               lower.limits=glmnet.lower, type.measure="class")
    } else{ # "gaussian"
      pheno.diff.vec <- npdrDiff(Ri.pheno.vals, NN.pheno.vals, diff.type="numeric-abs")
      # Run glmnet on the diff attribute columns
      npdrNET.model<-cv.glmnet(attr.diff.mat, pheno.diff.vec,alpha=glmnet.alpha,family="gaussian",
                               lower.limits=glmnet.lower, type.measure="mse")
    }
    npdrNET.coeffs<-as.matrix(predict(npdrNET.model,type="coefficients"))
    row.names(npdrNET.coeffs) <- c("intercept", colnames(attr.mat))  # add variable names to results
    glmnet.sorted<-as.matrix(npdrNET.coeffs[order(abs(npdrNET.coeffs),decreasing = T),],ncol=1) # sort
    npdr.stats.df<-data.frame(scores=glmnet.sorted)
  } # end glmnetSTIR option
  
  return(npdr.stats.df)
}

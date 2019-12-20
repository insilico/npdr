#=========================================================================#
#' diffRegression
#'
#' Wrapper for lm and glm-binomial to run regression for a phenotype diff vector, one attribute diff vector with optional covariate adjustment. Organizes regression statistics into a vector and then all attribute statistics combined in npdr.
#'
#' @param design.matrix.df Desgin matrix with variables: pheno.diff.vec (outcome variable as vector of diffs), attr.diff.vec (one predictor varialbe as vector of diffs) and optional covariates (regressors of non-interest) vector diffs.   
#' @param regression.type (\code{"lm"}, \code{"binomial"}) 
#' @param fast.reg logical, whether regression is run with speedlm or speedglm, default as F
#' @return vector of regression stats to put into list for npdr and combine into matrix
#'
#' @examples
#'
#' @export
# regression of the neighbor diff vector for one attribute
diffRegression <- function(design.matrix.df, regression.type = 'binomial', fast.reg = FALSE) {
  # if there are no covariates then ~. model is pheno.diff.vec ~ attr.diff.vec
  # otherwise ~. model is pheno.diff.vec ~ attr.diff.vec + covariates
  # design.matrix.df must have column named 'pheno.diff.vec'
  if (fast.reg == TRUE){
    if (regression.type == "lm"){
      mod <- speedlm(pheno.diff.vec ~ ., data = design.matrix.df)
    } else { # regression.type == "binomial"
      mod <- speedglm(pheno.diff.vec ~ ., data = design.matrix.df, family = binomial(link = logit))
    }
    res_df <- mod$df
  } else { # non-speedy version -- but why?
    if (regression.type == "lm"){
      mod <- lm(pheno.diff.vec ~ ., data = design.matrix.df)
    } else { # regression.type == "binomial"
      mod <- glm(pheno.diff.vec ~ ., family = binomial(link = logit), data = design.matrix.df)
    }
    res_df <- mod$df.residual
  }
  fit <- summary(mod)
  coeffs <- fit$coefficients
  
  ## create output NPDR summary stats (stats.vec)
  if (nrow(coeffs)<2){
    # for example, a monomorphic SNP might result in attribute stats (row 2) not being created
    beta_a <- NA
    beta_zscore_a <- NA
    pval.att <- NA
    message("Regression failure. Possible monomorphic SNP.\n")
  } else {
    beta_a <- coeffs[2,1]
    beta_zscore_a <- coeffs[2,3]  # standardized beta coefficient (col 3)
    ## use one-side p-value to test H1: beta>0 for case-control and continuous outcome
    pval.att <- pt(beta_zscore_a, res_df, lower = FALSE)  # one-sided p-val
  }
  
  stats.vec <- c(
    pval.att,            # one-sided p-value for attribute beta
    beta_a,                 # beta_a for attribute a
    beta_zscore_a,          # standardized beta for attribut a
    coeffs[1,1],  # beta_0, intercept, row 1 is inercept, col 1 is raw beta
    coeffs[1,4]   # p-value for intercept, row 1 is intercept, col 4 is p-val 
  )
  
  if (regression.type=="lm"){
    stats.vec <- c(stats.vec, fit$r.squared)} # add R^2 of fit, R.sqr for continuous outcomes
  
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
#' @param attr.diff.type diff type for attributes (\code{"numeric-abs"} or \code{"numeric-sqr"} for numeric, \code{"allele-sharing"} or \code{"match-mismatch"} for SNP). Phenotype diff uses same numeric diff as attr.diff.type when lm regression. For glm-binomial, phenotype diff is \code{"match-mismatch"} For correlation data (e.g., rs-fMRI), use \code{"correlation-data"}; diffs between two variables (e.g., ROIs) are taken across all their pairs of correlations and the attribute importances are given for the overall variable (e.g,. brain ROI), not individual pairs. 
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
#' @param use.glmnet logical, whether glmnet is employed
#' @param rm.attr.from.dist attributes for removal (possible confounders) from the distance matrix calculation. Argument for nearestNeighbors. None by default c()
#' @param neighbor.sampling "none" or \code{"unique"} if you want to use only unique neighbor pairs (used in nearestNeighbors)
#' @param separate.hitmiss.nbds for case/control data, find neighbors for same (hit) and opposite (miss) classes separately (TRUE) or find nearest neighborhoods before assigning hit/miss groups (FALSE). Uses nearestNeighborsSeparateHitMiss function
#' @param corr.attr.names character vector of p variable names that correspond to the variables used to create the p(p-1) correlation-data predictors. If not specified, integer (1...p) labels used.  
#' @param padj.method for p.adjust (\code{"fdr"}, \code{"bonferroni"}, ...) 
#' @param fast.reg logical, whether regression is run with speedlm or speedglm, default as F
#' @param dopar.nn logical, whether or not neighborhood is computed in parallel, default as F
#' @param dopar.reg logical, whether or not regression is run in parallel, default as F
#' 
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
#' 
npdr <- function(outcome, dataset, 
                 regression.type = "binomial", attr.diff.type = "numeric-abs",
                 nbd.method = "multisurf", nbd.metric = "manhattan", 
                 knn = 0, msurf.sd.frac = 0.5, 
                 covars = "none", covar.diff.type = "match-mismatch",
                 padj.method = "bonferroni", verbose = FALSE, 
                 use.glmnet = FALSE, glmnet.alpha = 1, glmnet.lower = 0, 
                 rm.attr.from.dist = c(), neighbor.sampling = "none",
                 separate.hitmiss.nbds = FALSE,
                 corr.attr.names=NULL,
                 fast.reg = FALSE, fast.dist = FALSE,
                 dopar.nn = FALSE, dopar.reg = FALSE){
  ##### parse the commandline 
  if (length(outcome)==1){
    # e.g., outcome="qtrait" or outcome=101 (pheno col index) and dataset is data.frame including outcome variable
    pheno.vec <- dataset[, outcome] # get phenotype
    attr.mat <- dataset %>% dplyr::select(- outcome) # outcome = "qtrait" or 101
  } else { # user specifies a separate phenotype vector
    pheno.vec <- outcome # assume users provides a separate outcome data vector
    attr.mat <- dataset # assumes dataset only contains attributes/predictors
  }
  rm(dataset)  # cleanup memory
  
  if(attr.diff.type=="correlation-data"){   # corrdata
    mynum <- dim(attr.mat)[2]               # corrdata
    for(i in seq(1,mynum-1,by=1)){          # corrdata
      mydiv <- i                            # corrdata
      if((mydiv*(mydiv - 1)) == mynum){     # corrdata
        my.dimension <- mydiv               # corrdata
        break                               # corrdata
      }                                     # corrdata
    }                                       # corrdata
    num.attr <- my.dimension                # corrdata
  }else{
    num.attr <- ncol(attr.mat)
  }
  num.samp <- nrow(attr.mat)
  
  # create a list of attribute indices for selecting columns in stretched matrix
  ###################################################################################################
  if(attr.diff.type=="correlation-data"){    # corrdata
    attr.idx.list <- list()                  # corrdata
    for(i in 1:num.attr){                    # corrdata
      lo.idx <- (i - 1)*(num.attr-1) + 1     # corrdata
      hi.idx <- i*(num.attr-1)               # corrdata
      attr.idx.list[[i]] <- c(lo.idx:hi.idx) # corrdata
    }                                        # corrdata
  }                                          # corrdata
  
  ##### get Neighbors (no phenotype used)
  # nbd.method (relieff, multisurf...), nbd.metric (manhattan...), k (for relieff nbd, theoerical surf default) 
  # msurf.sd.frac used by surf/multisurf relieff for theoretical k
  if (verbose){
    cat("Finding nearest neighbor pairs.\n")
  }
  start_time <- Sys.time() 
  if (separate.hitmiss.nbds){ # separate hit and miss neighborhoods
    neighbor.pairs.idx <- nearestNeighborsSeparateHitMiss(attr.mat, pheno.vec, 
                                         nb.method = nbd.method,
                                         nb.metric = nbd.metric, 
                                         sd.frac = msurf.sd.frac, k = knn,
                                         att_to_remove = rm.attr.from.dist,
                                         fast.dist = fast.dist,
                                         dopar.nn = dopar.nn)
  } else{ # allow neighborhoods to be imbalanced, often nearest hits are closer then misses, 
          # which could dilute the effect of misses
    neighbor.pairs.idx <- nearestNeighbors(attr.mat, nb.method = nbd.method, 
                                           nb.metric = nbd.metric, 
                                           sd.frac = msurf.sd.frac, k = knn,
                                           att_to_remove = rm.attr.from.dist,
                                           fast.dist = fast.dist,
                                           dopar.nn = dopar.nn)
  }
  num.neighbor.pairs <- nrow(neighbor.pairs.idx)
  k.ave.empirical <- mean(knnVec(neighbor.pairs.idx))
  if (neighbor.sampling == "unique"){
    if (verbose){
      cat("Extracting unique neighbors.\n")
    }
      # if you only want to return unique neighbors
      neighbor.pairs.idx <- uniqueNeighbors(neighbor.pairs.idx)
  }
  end_time <- Sys.time()
  if (verbose){
    cat("Neighborhood calculation:", capture.output(end_time - start_time), "\n")
    cat(num.neighbor.pairs, "total neighbor pairs (possible repeats).\n")
    erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
    # theoretical surf k (sd.frac=.5) for regression problems (does not depend on a hit/miss group)
    k.msurf.theory <- knnSURF(num.samp,msurf.sd.frac)
    cat("Theoretical (predicted) multiSURF average neighbors: ", k.msurf.theory,".\n",sep="")
    cat("Empirical (computed from neighborhood) average neighbors: ", k.ave.empirical,".\n",sep="")
    if (neighbor.sampling=="unique"){
      # if you only want to return unique neighbors
      num.neighbor.pairs <- nrow(neighbor.pairs.idx)
      cat(num.neighbor.pairs, "unique neighbor pairs.\n")
      cat("\nPerforming projected distance regression.\n")
    }
  }
  ### pheno diff vector for glm-binomial or lm to use in each attribute's diff regression in for loop.
  # Not needed in loop.
  # create pheno diff vector for linear regression (numeric)
  Ri.pheno.vals <- pheno.vec[neighbor.pairs.idx[,1]]
  NN.pheno.vals <- pheno.vec[neighbor.pairs.idx[,2]]
  if (regression.type == "lm"){
    pheno.diff.vec <- npdrDiff(Ri.pheno.vals, NN.pheno.vals, diff.type="numeric-abs")
  } else { #regression.type == "binomial"
    # create pheno diff vector for logistic regression (match-mismatch or hit-miss)  
    pheno.diff.vec <- npdrDiff(Ri.pheno.vals, NN.pheno.vals, diff.type="match-mismatch")
    # the reference group is the hit group, so the logistic probability is prob of a pair being a miss
    pheno.diff.vec <- as.factor(pheno.diff.vec)
  }
  
  # ----------------------------------------
  # run npdr, each attribute produces a list
  npdr.stats.list <- vector("list", num.attr) # initialize
  attr.diff.mat <- matrix(0, nrow = nrow(neighbor.pairs.idx), ncol = num.attr)
  # for npdrnet later, need matrix because npdrNET operates on all attributes at once
  if (length(covars) > 1){ # if covars is a vector or matrix
    if (use.glmnet == TRUE){
      message("penalized npdrNET does not currently support covariates.")
    }
    # default value is covar="none" (no covariates) which has length 1
    covars <- as.matrix(covars)  # if covars is just one vector, make sure it's a 1-column matrix
    num.covs <- length(covar.diff.type) # or ncol(covars)
    # covar.diff.type can be a vector of strings because each column of covars may be a different data type
    if (is.null(colnames(covars))){  # if covar vector has no column name, give it one
      covar.names <- paste0("cov", seq.int(num.covs)) # cov1, etc.
    } else {
      covar.names <- colnames(covars) # else get the name from covars
    }
    
    covar.diff.df <- data.frame(matrix(, nrow = nrow(neighbor.pairs.idx), ncol = 0))
    for (covar.col in (1:num.covs)){
      covar.name <- covar.names[covar.col]
      covar.vals <- covars[, covar.col]
      Ri.covar.vals <- covar.vals[neighbor.pairs.idx[,1]]
      NN.covar.vals <- covar.vals[neighbor.pairs.idx[,2]]
      covar.diff.vec <- npdrDiff(Ri.covar.vals, NN.covar.vals, 
                                 diff.type = covar.diff.type[covar.col])
      # add covar diff vector to data.frame
      # these covars will be included in each attribute's model

      covar.diff.df <- covar.diff.df %>% 
        dplyr::mutate(!!covar.name := covar.diff.vec)
    }
  }
  
  if (use.glmnet == FALSE){ # combine non-glmnet result lists into a matrix
    if (dopar.reg) { # perform regressions in parallel
      avai.cors <- parallel::detectCores() - 2
      cl <- parallel::makeCluster(avai.cors)
      doParallel::registerDoParallel(cl)
      
      npdr.stats.attr.mat <- foreach::foreach(
        attr.idx = seq.int(num.attr), .combine = 'rbind', .packages = c('dplyr')) %dopar% {
          if(attr.diff.type=="correlation-data"){                                                 # corrdata
            attr.vals <- attr.mat[, attr.idx.list[[attr.idx]]]                                    # corrdata
            Ri.attr.vals <- attr.vals[neighbor.pairs.idx[,1], ]                                   # corrdata
            NN.attr.vals <- attr.vals[neighbor.pairs.idx[,2], ]                                   # corrdata
            attr.diff.vec <- npdrDiff(Ri.attr.vals, NN.attr.vals, diff.type = attr.diff.type)     # corrdata
          }else{  
            attr.vals <- attr.mat[, attr.idx]
            Ri.attr.vals <- attr.vals[neighbor.pairs.idx[,1]]
            NN.attr.vals <- attr.vals[neighbor.pairs.idx[,2]]
            attr.diff.vec <- npdrDiff(Ri.attr.vals, NN.attr.vals, diff.type = attr.diff.type)
          }
        attr.diff.mat[, attr.idx] <- attr.diff.vec
        design.matrix.df <- data.frame(attr.diff.vec = attr.diff.vec,
                                       pheno.diff.vec = pheno.diff.vec)
        if (length(covars) > 1){ # if covars is a vector or matrix
          design.matrix.df <- data.frame(design.matrix.df, covar.diff.df)
        }
        # design.matrix.df = pheno.diff ~ attr.diff + option covar.diff
        return(diffRegression(design.matrix.df, regression.type = regression.type, fast.reg = fast.reg))
      } # end of foreach loop, regression done in parallel
      parallel::stopCluster(cl)
      
      
    } else { # non parallel version
      for (attr.idx in seq(1, num.attr)){
        if(attr.diff.type=="correlation-data"){                                                 # corrdata
          attr.vals <- attr.mat[, attr.idx.list[[attr.idx]]]                                    # corrdata
          Ri.attr.vals <- attr.vals[neighbor.pairs.idx[,1], ]                                   # corrdata
          NN.attr.vals <- attr.vals[neighbor.pairs.idx[,2], ]                                   # corrdata
          attr.diff.vec <- npdrDiff(Ri.attr.vals, NN.attr.vals, diff.type = attr.diff.type)     # corrdata
        } else{ 
          attr.vals <- attr.mat[, attr.idx]
          Ri.attr.vals <- attr.vals[neighbor.pairs.idx[,1]]
          NN.attr.vals <- attr.vals[neighbor.pairs.idx[,2]]
          attr.diff.vec <- npdrDiff(Ri.attr.vals, NN.attr.vals, diff.type = attr.diff.type)
        }
        attr.diff.mat[, attr.idx] <- attr.diff.vec
        # model data.frame to go into lm or glm-binomial
        design.matrix.df <- data.frame(attr.diff.vec = attr.diff.vec,
                                       pheno.diff.vec = pheno.diff.vec)
        ### diff vector for each covariate
        # optional covariates to add to design.matrix.df model
        if (length(covars) > 1){ # if covars is a vector or matrix
          design.matrix.df <- data.frame(design.matrix.df, covar.diff.df)
          # design.matrix.df$temp <- covar.diff.vec  # add the diff covar to the design matrix data frame
          # colnames(design.matrix.df)[2+covar.col] <- covar.name # change variable name
        }
        # design.matrix.df = pheno.diff ~ attr.diff + option covar.diff
        npdr.stats.list[[attr.idx]] <- diffRegression(design.matrix.df, regression.type = regression.type, fast.reg = fast.reg) 
      } # end of for loop, regression done for each attribute
      
      #npdr.stats.attr.mat <- bind_rows(npdr.stats.list)
      npdr.stats.attr.mat <- data.frame(do.call(rbind, npdr.stats.list))
      colnames(npdr.stats.attr.mat) <- c("pval.att","beta.raw.att", "beta.Z.att", "beta.0", "pval.0")
    }

    if(attr.diff.type=="correlation-data"){                                            # corrdata    
      if(is.null(corr.attr.names)){                                                    # corrdata
        att.names <- as.character(c(1:num.attr))                                       # corrdata
      }else{                                                                           # corrdata
        att.names <- as.character(corr.attr.names)                                     # corrdata
      }                                                                                # corrdata
      
      #### Create Results Data Frame for NPDR
      # attribute p-values
      # relies on first column [, 1] being the p-value for now
      # later, first columns become att and adjusted p-value
      attr.pvals <- npdr.stats.attr.mat[, 1]
      # order-index for sorted attribute-beta p-values
      attr.pvals.order.idx <- order(attr.pvals, decreasing = F)
      # adjust p-values using Benjamini-Hochberg (default)
      attr.pvals.adj <- p.adjust(attr.pvals[attr.pvals.order.idx], method=padj.method)
      # order by attribute p-value
      npdr.stats.pval_ordered.mat <- npdr.stats.attr.mat[attr.pvals.order.idx, ]
      # prepend adjused attribute p-values to first column
      npdr.stats.pval_ordered.mat <- cbind(attr.pvals.adj, npdr.stats.pval_ordered.mat)
      # prepend attribute column (att)
      att = colnames(attr.mat)
      npdr.stats.pval_ordered.mat <- cbind(data.frame(att=att, stringsAsFactors=FALSE), 
                                           data.frame(npdr.stats.pval_ordered.mat, row.names=NULL))
      colnames(npdr.stats.pval_ordered.mat) <- c("att", "pval.adj", "pval.att", "beta.raw.att", "beta.Z.att",
                                                   "beta.0", "pval.0")
      # dataframe final output for regular npdr
      npdr.stats.df <- data.frame(npdr.stats.pval_ordered.mat)
      
      #npdr.stats.df <- npdr.stats.attr.mat %>%                                         # corrdata
      #  mutate(att = att.names, # add an attribute column                              # corrdata
      #         pval.adj = p.adjust(pval.att, method = padj.method) # adjust p-values   # corrdata
      #  ) %>% arrange(pval.att) %>% # order by attribute p-value                       # corrdata
      #  dplyr::select(att, pval.adj, everything()) %>% # reorder columns               # corrdata
      #  as.data.frame() # convert tibbles to df -- can we remove this step?            # corrdata
      
    }else{ # non-correlation matrix predictors
      
      # npdr.stats.df <- as.data.frame(npdr.stats.attr.mat) %>% 
      #   mutate(att = colnames(attr.mat), # add an attribute column
      #          pval.adj = p.adjust(pval.att, method = padj.method) # adjust p-values
      #   ) %>% arrange(pval.att) %>% # order by attribute p-value 
      #   dplyr::select(att, pval.adj, everything()) %>% # reorder columns
      #   as.data.frame() # convert tibbles to df -- can we remove this step?
      
      #### Create Results Data Frame for NPDR
      # attribute p-values
      # relies on first column [, 1] being the p-value for now
      # later, first columns become att and adjusted p-value
      attr.pvals <- npdr.stats.attr.mat[, 1]
      # order-index for sorted attribute-beta p-values
      attr.pvals.order.idx <- order(attr.pvals, decreasing = F)
      # adjust p-values using Benjamini-Hochberg (default)
      attr.pvals.adj <- p.adjust(attr.pvals[attr.pvals.order.idx], method=padj.method)
      # order by attribute p-value
      npdr.stats.pval_ordered.mat <- npdr.stats.attr.mat[attr.pvals.order.idx, ]
      # prepend adjused attribute p-values to first column
      npdr.stats.pval_ordered.mat <- cbind(attr.pvals.adj, npdr.stats.pval_ordered.mat)
      # prepend attribute column (att)
      att = colnames(attr.mat)
      npdr.stats.pval_ordered.mat <- cbind(data.frame(att=att, stringsAsFactors=FALSE), 
                                           data.frame(npdr.stats.pval_ordered.mat, row.names=NULL))
      if (regression.type=="lm"){# stats colnames for lm
        colnames(npdr.stats.pval_ordered.mat) <- c("att", "pval.adj", "pval.att", "beta.raw.att", "beta.Z.att",  
                                                   "beta.0", "pval.0", "R.sqr")
      } else{ # stats columns for glm-binomial
        colnames(npdr.stats.pval_ordered.mat) <- c("att", "pval.adj", "pval.att", "beta.raw.att", "beta.Z.att",
                                                   "beta.0", "pval.0")
      }
      # dataframe final output for regular npdr
      npdr.stats.df <- data.frame(npdr.stats.pval_ordered.mat)
    } # end-else non-correlation-based predictors
    
  } else { # use.glmnet = TRUE, Option npdrNET
    # Run glmnet on the diff attribute columns
    # Need to create a data matrix with each column as a vector of diffs for each attribute.
    # Need matrix because npdrNET operates on all attributes at once.
    attr.diff.mat <- matrix(0,nrow=nrow(neighbor.pairs.idx),ncol=num.attr)
    for (attr.idx in seq(1, num.attr)){
      if(attr.diff.type=="correlation-data"){                                                 # corrdata
        attr.vals <- attr.mat[, attr.idx.list[[attr.idx]]]                                    # corrdata
        Ri.attr.vals <- attr.vals[neighbor.pairs.idx[,1], ]                                   # corrdata
        NN.attr.vals <- attr.vals[neighbor.pairs.idx[,2], ]                                   # corrdata
        attr.diff.vec <- npdrDiff(Ri.attr.vals, NN.attr.vals, diff.type = attr.diff.type)     # corrdata
      }else{      
        attr.vals <- attr.mat[, attr.idx]
        Ri.attr.vals <- attr.vals[neighbor.pairs.idx[,1]]
        NN.attr.vals <- attr.vals[neighbor.pairs.idx[,2]]
        attr.diff.vec <- npdrDiff(Ri.attr.vals, NN.attr.vals, diff.type=attr.diff.type)
      }
      attr.diff.mat[,attr.idx] <- attr.diff.vec
    } # end for
    #
    Ri.pheno.vals <- pheno.vec[neighbor.pairs.idx[,1]]
    NN.pheno.vals <- pheno.vec[neighbor.pairs.idx[,2]]
    if (glmnet.alpha != "cluster"){ # temporary trick to cluster attributes by diff vector
      if (regression.type=="binomial"){
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
      npdr.stats.df <- data.frame(scores = glmnet.sorted) 
      # %>%
        # tibble::rownames_to_column('att')
    } else{ # glmnet.alpha == "cluster", so don't do regression and return the attribute diff vectors
      # might not need the phenotype diff for clustering, but add anyway.
      if (regression.type=="binomial"){
        pheno.diff.vec <- npdrDiff(Ri.pheno.vals, NN.pheno.vals, diff.type="match-mismatch")
      } else{ # "gaussian"
        pheno.diff.vec <- npdrDiff(Ri.pheno.vals, NN.pheno.vals, diff.type="numeric-abs")
      }
      colnames(attr.diff.mat) <- colnames(attr.mat)
      # not actually stats, contains a pairs x attr diff matrix 
      npdr.stats.df <- data.frame(attr.diff.mat, pheno.diff=pheno.diff.vec)  
    } # end cluster option
  } # end glmnetNPDR option
  return(npdr.stats.df)
}

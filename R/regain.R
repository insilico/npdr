#'
#' \code{getMainEffects}
#' @param dats \code{labelledDataFrame.frame} with variables in columns and samples in rows.
#' @param depVarName \code{string} name of the classLabel variable.
#' @param indVarNames Vector of independent variable names, to name output matrices
#' @param reg.fam \code{string} glm regression family name.
#' @param nCovs \code{numeric}  of covariates included.
#' @param excludeMainEffects \code{logical} indicating whether to exclude main effect terms.
#' @return regainMatrix \code{matrix} of variable by variable regression coefficients.
#' @export
regain <- function(dats, depVarName, indVarNames, reg.fam="binomial", nCovs=0, excludeMainEffects=F){
  
  # get the main effect raw coefs, std errs, pvals, standard coefs, and convergence indicators
  mains <- getMainEffects(dats, 
                          indVarNames=indVarNames, 
                          depVarName=depVarName, 
                          reg.fam=reg.fam, 
                          nCovs=nCovs)
  
  # get the interaction effect raw coefs, std errs, pvals, standard coefs, and convergence indicators
  ints <- getInteractionEffects(dats, 
                                indVarNames=indVarNames,
                                depVarName=depVarName, 
                                reg.fam=reg.fam, 
                                nCovs=nCovs, 
                                excludeMainEffects=excludeMainEffects)
  
  tmp <- ints$conv.mat # matrix of convergence indicators for interaction models - 0's on diagonal
  diag(tmp) <- as.numeric(mains["converged",]) # load diagonal with main effect convergence indicators
  conv.all.mat <- tmp
  
  tmp <- ints$coef.mat # matrix of raw coefs for interaction models - 0's on diagonal
  diag(tmp) <- as.numeric(mains["coef",]) # load diagonal with main effect raw coefs
  coef.all.mat <- tmp
  
  tmp <- ints$stderr.mat # matrix of standard errors for interaction model coefs - 0's on diagonal
  diag(tmp) <- as.numeric(mains["se",]) # load diagonal with standard errors for main effect coefs
  stderr.all.mat <- tmp
  
  tmp <- ints$pval.mat # matrix of pvals for interaction model coefs - 0's on diagonal
  diag(tmp) <- as.numeric(mains["pval",]) # load diagonal with pvals for main effect coefs
  pval.all.mat <- tmp
  
  tmp <- ints$stdbeta.mat # matrix of standardized interaction model coefs - 0's on diagonal
  diag(tmp)  <- as.numeric(mains["stdbeta",]) # load diagonal with standardized main effect coefs
  stdbeta.all.mat <- tmp
  
  list(convergence=conv.all.mat, raw.coefs=coef.all.mat, stderrs=stderr.all.mat, pvals=pval.all.mat, stdbetas=stdbeta.all.mat)
  
}

#'
#' \code{getMainEffects}
#' Wrapper for fitMainEffectModel with sapply over independent variable names. 
#' Does not include covariates, so if you have covariates use 
#' indVarNames = colnames(dats)[1:(ncol(dats) - nCovs - 1)] 
#' if the columns are organized as: independent vars, covars, pheno
#' @param dats \code{data.frame} with variables in columns and samples in rows.
#' @param indVarNames vector of independent variable names to loop over.
#' @param depVarName \code{string} name of the classLabel variable.
#' @param regressionFamily \code{string} glm regression family name.
#' @param nCovs \code{numeric}  of included covariates.
#' @return mainEffectValues \code{vector} of main effect values.
getMainEffects <- function(dats, indVarNames=NULL, depVarName=NULL, 
                           reg.fam="binomial", nCovs=0){
  theMains <- sapply(indVarNames,
                     function(x) fitMainEffectModel(dats, variableName=x, 
                                        depVarName=depVarName, 
                                        regressionFamily=reg.fam, 
                                        numCovariates=nCovs),
                     simplify="matrix", USE.NAMES=T)
  
  return(theMains)
}

#'
#' \code{fitMainEffectModel}
#' Get the main effect of a variable using generalized linear regression .
#' @param labelledDataFrame \code{data.frame} with variables in columns and samples in rows.
#' @param variableName \code{string} name of the variable to consider.
#' @param depVarName \code{string} name of the classLabel variable.
#' @param regressionFamily \code{string} glm regression family name.
#' @param numCovariates \code{numeric}  of included covariates.
#' @return \code{data.frame} with variable, convergence status, beta coefficient,
#' p-value, standard error and standardized beta columns.
fitMainEffectModel <- function(labelledDataFrame, variableName, depVarName,
                               regressionFamily, numCovariates) {
  if (numCovariates > 0) {
    covarsStart <- ncol(labelledDataFrame) - numCovariates
    covarNames <- colnames(labelledDataFrame)[covarsStart:(ncol(labelledDataFrame) - 1)]
    covarModelParts <- NULL
    for (i in 1:length(covarNames)) {
      covarModelParts <- c(covarModelParts, paste("`", covarNames, "`", sep = ""))
    }
    formulaString <- paste(depVarName, " ~ `", variableName,
                           "` + ", paste(covarModelParts, collapse = "+"), sep = "")
    regressionFormula <- as.formula(formulaString)
  } else {
    regressionFormula <- as.formula(paste(depVarName, "~",
                                          paste("`", variableName, "`", sep = ""),
                                          sep = " "))
  }
  regressionModel <- glm2::glm2(regressionFormula, data = labelledDataFrame,
                                family = regressionFamily, na.action=na.omit)
  
  #interceptCoeff <- summary(regressionModel)$coefficients[1, "Estimate"]
  mainCoeff <- summary(regressionModel)$coefficients[2, "Estimate"]
  mainStdErr <- summary(regressionModel)$coefficients[2, "Std. Error"]
  if (regressionFamily == "binomial") {
    #interceptPval <- summary(regressionModel)$coefficients[1, "Pr(>|z|)"]
    mainPval <- summary(regressionModel)$coefficients[2, "Pr(>|z|)"]
    mainStat <- summary(regressionModel)$coefficients[2, "z value"]
  }
  else {
    #interceptPval <- summary(regressionModel)$coefficients[1, "Pr(>|t|)"]
    mainPval <- summary(regressionModel)$coefficients[2, "Pr(>|t|)"]
    mainStat <- summary(regressionModel)$coefficients[2, "t value"]
  }
  
  data.frame(converged = regressionModel$converged,
             coef = mainCoeff,
             se = mainStdErr,
             pval = mainPval,
             stdbeta = mainStat)
}

#'
#' \code{getInteractionEffects}
#'
#' @param dats \code{labelledDataFrame.frame} with variables in columns and samples in rows.
#' @param depVarName \code{string} name of the classLabel variable.
#' @param reg.fam \code{string} glm regression family name.
#' @param nCovs \code{numeric}  of covariates included.
#' @param excludeMainEffects \code{logical} indicating whether to exclude main effect terms.
#' @param indVarNames Vector of independent variable names, to name output matrices
#' @return \code{data.frame} with variable, convergence status, beta coefficient,
#' p-value, standard error and standardized beta columns.
# i.e. indVarNames = colnames(dats)[1:(ncol(dats) - nCovs - 1)] if the columns are organized as: independent vars, covars, pheno
#
getInteractionEffects <- function(dats, depVarName, reg.fam="binomial", nCovs=0, excludeMainEffects=F, indVarNames){
  
  numVars <- length(indVarNames)
  conv.mat <- matrix(0, nrow=numVars, ncol=numVars)
  coef.mat <- matrix(0, nrow=numVars, ncol=numVars)
  stderr.mat <- matrix(0, nrow=numVars, ncol=numVars)
  pval.mat <- matrix(0, nrow=numVars, ncol=numVars)
  stdbeta.mat <- matrix(0, nrow=numVars, ncol=numVars)
  for(i in 1:(numVars - 1)){
    
    idx1 <- i
    
    for(j in (i + 1):numVars){
      
      idx2 <- j
      
      getInts <- fitInteractionModel(dats, 
                                     variableIndices=c(idx1, idx2), 
                                     depVarName=depVarName, 
                                     regressionFamily=reg.fam, 
                                     numCovariates=nCovs, 
                                     excludeMainEffects=excludeMainEffects)
      
      conv.mat[i, j] <- getInts$converged
      coef.mat[i, j] <- getInts$coef
      stderr.mat[i, j] <- getInts$stderr
      pval.mat[i, j] <- getInts$pval
      stdbeta.mat[i, j] <- getInts$stdbeta
      
    }
    
  }
  
  conv.mat <- conv.mat + t(conv.mat)
  colnames(conv.mat) <- indVarNames
  row.names(conv.mat) <- indVarNames
  
  coef.mat <- coef.mat + t(coef.mat)
  colnames(coef.mat) <- indVarNames
  row.names(coef.mat) <- indVarNames
  
  stderr.mat <- stderr.mat + t(stderr.mat)
  colnames(stderr.mat) <- indVarNames
  row.names(stderr.mat) <- indVarNames
  
  pval.mat <- pval.mat + t(pval.mat)
  colnames(pval.mat) <- indVarNames
  row.names(pval.mat) <- indVarNames
  
  stdbeta.mat <- stdbeta.mat + t(stdbeta.mat)
  colnames(stdbeta.mat) <- indVarNames
  row.names(stdbeta.mat) <- indVarNames
  
  list(conv.mat=conv.mat, coef.mat=coef.mat, stderr.mat=stderr.mat, pval.mat=pval.mat, stdbeta.mat=stdbeta.mat)
}

#'
#' \code{fitInteractionModel}
#'
#' @param labelledDataFrame \code{labelledDataFrame.frame} with variables in columns and samples in rows.
#' @param variableIndices \code{vector} of column indices of variable pairs.
#' @param depVarName \code{string} name of the classLabel variable.
#' @param regressionFamily \code{string} glm regression family name.
#' @param numCovariates \code{numeric}  of covariates included.
#' @param excludeMainEffects \code{logical} indicating whether to exclude main effect terms.
#' @return \code{data.frame} with variable, convergence status, beta coefficient,
#' p-value, standard error and standardized beta columns.
fitInteractionModel <- function(labelledDataFrame, variableIndices, depVarName,
                                regressionFamily, numCovariates,
                                excludeMainEffects) {
  variable1Idx <- variableIndices[1]
  variable2Idx <- variableIndices[2]
  variableNames <- colnames(labelledDataFrame)[1:(ncol(labelledDataFrame) - 1)]
  variable1Name <- variableNames[variable1Idx]
  variable2Name <- variableNames[variable2Idx]
  if (excludeMainEffects) {
    interactionTerm <- paste("`", variable1Name, "`", ":", "`",
                             variable2Name, "`", sep = "")
  } else {
    interactionTerm <- paste("`", variable1Name, "`", "*", "`",
                             variable2Name, "`", sep = "")
  }
  if (numCovariates > 0) {
    covarsStart <- ncol(labelledDataFrame) - numCovariates
    covarNames <- colnames(labelledDataFrame)[covarsStart:(ncol(labelledDataFrame)-1)]
    covarsModelParts <- paste(covarNames, collapse = " + ")
    regressionFormula <- as.formula(paste(depVarName, "~",
                                          paste(interactionTerm, " + ",
                                                covarsModelParts, sep = ""),
                                          sep = " "))
  }
  else {
    regressionFormula <- as.formula(paste(depVarName, "~", interactionTerm,
                                          sep = " "))
  }
  
  # use glm2: Fits generalized linear models using the same model specification
  # as glm in the stats package, but with a modified default fitting method that
  # provides greater stability for models that may fail to converge using glm
  # bcw - 10/18/15
  regressionModel <- glm2::glm2(regressionFormula,
                                family = binomial(link = "logit"),
                                data = labelledDataFrame, na.action=na.omit)
  
  if (numCovariates > 0) {
    if (excludeMainEffects) {
      interactionTermIndex <- 1 + numCovariates + 1
    }
    else {
      interactionTermIndex <- 3 + numCovariates + 1
    }
  }
  else {
    if (excludeMainEffects) {
      interactionTermIndex <- 2
    }
    else {
      interactionTermIndex <- 4
    }
  }
  #print(summary(regressionModel))
  
  coeffs <- summary(regressionModel)$coefficients
  #print(nrow(coeffs))
  #print(interactionTermIndex)
  
  # catch NA problems with interaction term fit
  if (nrow(coeffs)<interactionTermIndex){
    interactionCoeff <- NA
    interactionStdErr <- NA
    interactionPval <- NA
    interactionStat <- NA
  } else {
    interactionCoeff <- coeffs[interactionTermIndex, "Estimate"]
    interactionStdErr <- coeffs[interactionTermIndex, "Std. Error"]
    if (regressionFamily == "binomial") {
      interactionPval <- coeffs[interactionTermIndex, "Pr(>|z|)"]
      interactionStat <- coeffs[interactionTermIndex, "z value"]
    }
    else {
      interactionPval <- summary(regressionModel)$coefficients[interactionTermIndex, "Pr(>|t|)"]
      interactionStat <- summary(regressionModel)$coefficients[interactionTermIndex, "t value"]
    }
  }  # end catching fit NA problems
  
  data.frame(converged = regressionModel$converged,
             coef = interactionCoeff,
             stderr = interactionStdErr,
             pval = interactionPval,
             stdbeta = interactionStat)
}

#' \code{EpistasisRank}
#' 
#' @family feature selection functions
#' @param G \code{matrix} genetic association interaction network.
#' @param Gamma_vec \code{numeric} gamma vector, either a constant value or a vector
#' @return rankTable \code{data.frame} with variable, EpistasisRank, 
#' diagonal and degree columns. Sorted.
#' @export
EpistasisRank <- function(G = NULL, Gamma_vec = 0.85){
  n <- nrow(G)
  geneNames <- colnames(G)
  Gdiag <- diag(G)
  Gtrace <- sum(Gdiag)
  #colsum <- colSums(G)
  diag(G) <- 0
  Gtrace <- Gtrace * n
  colsumG <- colSums(G)
  #rowSumG <- rowSums(G)
  rowsum_denom <- matrix(0, n, 1)
  for (i in 1:n) {
    localSum <- 0
    for (j in 1:n) {
      factor <- ifelse(G[i, j] != 0, 1, 0)
      localSum <- localSum + (factor * colsumG[j])
    }
    rowsum_denom[i] <- localSum
  }
  if (length(Gamma_vec) == 1) {
    gamma_vec <- rep(Gamma_vec, n)
  } else {
    gamma_vec <- Gamma_vec
  }
  gamma_matrix <- matrix(nrow = n, ncol = n, data = rep(gamma_vec, n))
  if (Gtrace) {
    b <- ((1 - gamma_vec) / n) + (Gdiag / Gtrace)
  }
  else {
    b <- ((1 - gamma_vec) / n)
  }
  D <- matrix(nrow = n, ncol = n, data = c(0))
  diag(D) <- 1 / colsumG
  I <- diag(n)
  temp <- I - gamma_matrix * G %*% D
  r <- solve(temp, b)
  ERank <- r / sum(r)
  rankTable <- data.frame(gene = geneNames, ER = ERank)
  rankTable[order(rankTable$ER, decreasing = TRUE), ]
}
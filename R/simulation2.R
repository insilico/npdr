# Functions in simulation2.R file in order:
# not exported: r_to_z_fn, stretch_mat
# exported: generate_structured_corrmat, splitDataset, createMainEffects
#           createSimulation2, createfMRIsimulation
#
#
# Rcpp::sourceCpp('R/arma_getEigenValues.cpp')
# Rcpp::cppFunction(depends="RcppArmadillo",
#                    'arma::vec getEigenValues(arma::mat M) {
#                    return arma::eig_sym(M);
# }')

# fisher r-to-z
r_to_z_fn <- function(x) {
  r.z <- 0.5 * log((1 + x) / (1 - x))
  r.z
}

# vectorize subject correlation matrix and preserve order
stretch_mat <- function(M) {
  mat <- numeric()
  for (k in 1:nrow(M)) {
    mat <- c(mat, M[k, -k])
  }
  return(mat)
}

# =========================================================================================#
#' generate_structured_corrmat
#'
#' parameters:
#'
#' @param g random graph
#' @param num.variables number of independent variables in data matrix
#' @param hi.cor upper baseline pairwise functional correlation in control group
#' @param lo.cor lower baseline pairwise functional correlation in control group
#' @param graph.type either Erdos-Renyi or Scale-Free graph
#' @param plot.graph logical indicating whether to plot graph or not
#' @param make.diff.cors logical indicating whether case correlation matrix for differential correlation is being created or not
#' @param nbias number of functional interaction variables
#' @param hi.cor.fixed high correlation for network connected/non-functional features
#' @param hi.cor.tmp high correlation proportion
#' @param lo.cor.tmp low correlation proportion
#' @param use.Rcpp logical, whether to use Rcpp or not
#' 
#' 
#' @return A list containing:
#' \describe{
#'   \item{cor.mat}{structured correlation matrix}
#'   \item{deg.vec}{degree vector corresponding to network adjacency}
#'   \item{A.mat}{adjacency matrix corresponding to network}
#'   \item{sig.vars}{indices of functional variables}
#' }
#'
#' @export

generate_structured_corrmat <- function(g = NULL,
                                        num.variables = 100,
                                        hi.cor.tmp = 0.8,
                                        lo.cor.tmp = 0.2,
                                        hi.cor.fixed = 0.8,
                                        graph.type = "Erdos-Renyi",
                                        plot.graph = FALSE,
                                        make.diff.cors = FALSE,
                                        nbias = 1, use.Rcpp = FALSE) {
  if (use.Rcpp) {
    Rcpp::cppFunction(
      depends = "RcppArmadillo",
      "arma::vec getEigenValues(arma::mat M) {
                      return arma::eig_sym(M);
  }"
    )
  }
  if (abs(as.integer(num.variables) - num.variables) > 1e-9) {
    stop("generate_structured_corrmat: num.variables should be a positive integer")
  }

  p <- num.variables

  kvec <- igraph::degree(g)

  # make plot of random graph
  if (plot.graph) {
    plot(g, vertex.size = 1, vertex.label.color = kvec)
  }

  # which variables are not connected in network
  idx.not.connected <- which(kvec == 0)

  # which variables are connected in network
  if (length(idx.not.connected) > 0) { # if any are not connected
    idx.connected <- seq(1, length(kvec), by = 1)[-idx.not.connected]
  } else { # if all are connected
    idx.connected <- seq(1, length(kvec), by = 1)
  }

  # functional variables are randomly defined from connected variables
  if (length(idx.connected) >= nbias) { # if number connected at least nbias
    diff.cor.vars <- sample(idx.connected, size = nbias, replace = FALSE)
  } else { # if too few are connected
    n.add.vars <- nbias - length(idx.connected)
    diff.cor.vars <- c(idx.connected, sample(idx.not.connected, size = n.add.vars, replace = FALSE)) # add some from unconnected
  }

  ## make make noisy covariance matrix form structure of adjacency matrix
  Adj <- as.matrix(igraph::get.adjacency(g))

  # if simply generating a correlation matrix without between-group differential correlation
  if (!make.diff.cors) {
    tmp <- Adj[upper.tri(Adj)] # upper triangle of adjacency
    tmp2 <- (hi.cor.tmp + rnorm(length(tmp), sd = .1)) * tmp # connected (high) correlations
    tmp3 <- -(tmp - 1) # additive inverse of adjacency upper triangle
    tmp4 <- (lo.cor.tmp + rnorm(length(tmp3), sd = .1)) * tmp3 # non-connected (low) correlations

    upper.mat <- matrix(0, nrow = dim(Adj)[1], ncol = dim(Adj)[1]) # initialize connected correlations matrix
    lower.mat <- matrix(0, nrow = dim(Adj)[1], ncol = dim(Adj)[1]) # initialize non-connected correlations matrix

    upper.mat[upper.tri(upper.mat)] <- tmp2 # load upper triangle of connected matrix
    lower.mat[upper.tri(lower.mat)] <- tmp4 # load upper triangle of non-connected matrix
    new.mat <- upper.mat + lower.mat # add connected/non-connected to get full matrix
    new.mat <- new.mat + t(new.mat) # make symmetric
    diag(new.mat) <- 1 # 1's on the diagonal

    # if generating correlation matrix and between-group differential correlation
  } else {

    # find connections of functional attributes, store connections in list,
    # and create binary matrix for functional and connections only (Subset of Adj)
    #
    adj.tmp1 <- Adj * 0 # initialize binary matrix
    connected.list <- list() # list for functional connections
    for (i in 1:length(diff.cor.vars)) { # for each functional,
      connected.list[[i]] <- which(abs(Adj[diff.cor.vars[i], ] - 1) < 1e-9) # find connections,
      adj.tmp1[diff.cor.vars[i], connected.list[[i]]] <- 1 # make connections 1's in binary matrix, and
      adj.tmp1[connected.list[[i]], diff.cor.vars[i]] <- 1 # make symmetric
    }

    adj.tmp2 <- Adj - adj.tmp1 # binary matrix for all non-functional connections (Subset of Adj)

    adj.tmp3 <- -(Adj - 1) # additive inverse of Adjacency for all non-connections

    tmp <- adj.tmp1[upper.tri(adj.tmp1)] # functional connections only (upper triangle)
    tmp2 <- (hi.cor.tmp + rnorm(length(tmp), sd = .1)) * tmp # functional connections correlations matrix
    tmp3 <- (lo.cor.tmp + rnorm(length(adj.tmp3[upper.tri(adj.tmp3)]), sd = .1)) * adj.tmp3[upper.tri(adj.tmp3)] # non-connections correlations matrix
    tmp4 <- (hi.cor.fixed + rnorm(length(adj.tmp2[upper.tri(adj.tmp2)]), sd = .1)) * adj.tmp2[upper.tri(adj.tmp2)] # non-functional connections correlations matrix

    mat1 <- matrix(0, nrow = dim(Adj)[1], ncol = dim(Adj)[1]) # initialize functional correlation matrix
    mat2 <- matrix(0, nrow = dim(Adj)[1], ncol = dim(Adj)[1]) # initialize non-connection correlation matrix
    mat3 <- matrix(0, nrow = dim(Adj)[1], ncol = dim(Adj)[1]) # initialize non-functional connection correlation matrix

    mat1[upper.tri(mat1)] <- tmp2 # load upper triangle of functional matrix
    mat2[upper.tri(mat2)] <- tmp3 # load upper triangle of non-connection matrix
    mat3[upper.tri(mat3)] <- tmp4 # load upper triangle of non-functional matrix

    new.mat <- mat1 + mat2 + mat3 # combine all correlations in single matrix
    new.mat <- new.mat + t(new.mat) # make symmetric
    diag(new.mat) <- 1 # 1's on diagonal
  }
  R <- new.mat

  # correct for negative eigenvalues to make matrix positive definite
  #
  if (use.Rcpp) { # compute eigenvalues and make diag matrix
    R.d <- diag(sort(c(getEigenValues(R)), decreasing = TRUE))
  } else {
    R.d <- diag(eigen(R)$values)
  }

  tmp <- diag(R.d) # vector of eigenvalues

  if (any(tmp < 0)) { # if any eigenvalues are negative
    R.V <- eigen(R)$vectors # compute eigenvectors,
    tmp[tmp < 0] <- 1e-7 # make negative into small positive,
    diag(R.d) <- tmp # replace in R.d,
    R.fix <- R.V %*% R.d %*% t(R.V) # compute new correlation matrix, and
    R <- R.fix # store in R
  }
  R <- as.matrix(R)

  # make 1's on diagonal of R
  #
  inv.diag <- 1 / diag(R) # multiplicative inverse of diag(R)
  mydiag <- diag(length(inv.diag)) # initialize diagonal matrix for inv.diag
  diag(mydiag) <- inv.diag # swap 1's for inv.diag
  mydiag <- sqrt(mydiag) # take sqrt of diag matrix
  R <- mydiag %*% R %*% mydiag # compute corrected correlation matrix with 1's on diagonal (Still Pos. Def.)

  # return correlation matrix, degree vector, adjacency matrix, and functional variables
  list(corrmat = R, deg.vec = kvec, A.mat = Adj, sig.vars = diff.cor.vars)
}
######################################################


# ======================================================================================#
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
#' @family simulation
#' @export
##########################################################################################
splitDataset <- function(all.data = NULL,
                         pct.train = 0.5,
                         pct.holdout = 0.5,
                         pct.validation = 0,
                         label = "class") {
  if (is.null(all.data)) {
    # stop or warning and return list of length 0?
    stop("No data passed")
  }
  if (1.0 - (pct.train + pct.holdout + pct.validation) > 0.001) {
    stop("Proportions of training, holdout and testing must to sum to 1")
  }
  if (!(label %in% colnames(all.data))) {
    stop("Class label is not in the column names or used more than once in data set column names")
  }
  if (label == "class") {
    if (!is.factor(all.data[, label])) {
      all.data[, label] <- factor(all.data[, label])
    }
    if (nlevels(all.data[, label]) != 2) {
      stop("Cannot split data set with more than or less than 2 factor levels in the class label")
    }
  }

  if (label == "class") {
    class.levels <- levels(all.data[, label])
    ind.case <- rownames(all.data)[all.data[, label] == class.levels[1]]
    ind.ctrl <- rownames(all.data)[all.data[, label] == class.levels[2]]

    n.case <- length(ind.case)
    n.ctrl <- length(ind.ctrl)

    n.validation.case <- floor(pct.validation * n.case)
    n.holdo.case <- floor(pct.holdout * n.case)
    n.train.case <- n.case - n.validation.case - n.holdo.case
    partition.case <- sample(c(
      rep(3, n.validation.case), rep(2, n.holdo.case),
      rep(1, n.train.case)
    ), n.case)

    n.validation.ctrl <- floor(pct.validation * n.ctrl)
    n.holdo.ctrl <- floor(pct.holdout * n.ctrl)
    n.train.ctrl <- n.ctrl - n.validation.ctrl - n.holdo.ctrl
    partition.ctrl <- sample(c(
      rep(3, n.validation.ctrl),
      rep(2, n.holdo.ctrl),
      rep(1, n.train.ctrl)
    ), n.ctrl)

    all.data <- data.frame(all.data)
    all.data[, label] <- factor(all.data[, label])
    levels(all.data[, label]) <- c(-1, 1)
    data.case <- all.data[ind.case, ]
    data.ctrl <- all.data[ind.ctrl, ]
    X_train <- rbind(data.case[partition.case == 1, ], data.ctrl[partition.ctrl == 1, ])
    X_holdo <- rbind(data.case[partition.case == 2, ], data.ctrl[partition.ctrl == 2, ])
    X_validation <- rbind(data.case[partition.case == 3, ], data.ctrl[partition.ctrl == 3, ])
  } else {
    num_sample <- length(all.data[, label])
    n.train <- floor(pct.train * num_sample)
    n.holdout <- floor(pct.holdout * num_sample)
    n.validation <- floor(pct.validation * num_sample)
    partition <- sample(c(
      rep(3, n.validation),
      rep(2, n.holdout),
      rep(1, n.train)
    ))
    X_train <- rbind(all.data[partition == 1, ])
    X_holdo <- rbind(all.data[partition == 2, ])
    X_validation <- rbind(all.data[partition == 3, ])
  }

  # if(nrow(X_validation) == 0) {
  #   data.sets <- list(train=X_train, holdout=X_holdo)
  # } else {
  data.sets <- list(train = X_train, holdout = X_holdo, validation = X_validation)
  # }
  #
  data.sets
}

# main effect with quantitative outcome using label ("class" (dichotomious) or "qtrait"(quantitative))
# main effect with imbalanced outcome using pct.imbalance (decreasing pct. will generate less control sample)
# parameters in effect n.e=num.variables, n.db=num.samples, sd.b=bias, p.b=pct.signals

#' Create a simulated data set with main effects
#'
#' \eqn{X = BS + \Gamma G + U}
#'
#' S = Biological group                                                                                                                   m
#' G = Batch
#' U = random error
#'
#' NOTE:  If you use conf=TRUE, then you must have exactly two surrogate
#' variables in the database this function only allows for confounding in
#' the database, not confounding in the new samples
#'
#' @param n.e number of variables
#' @param n.db sample size in database
#' @param n.ns sample size in newsample
#' @param label A character vector for the name of the outcome column. class for classification
#' and qtrait for regression
#' @param pct.imbalance A numeric percentage to indicate proportion of the imbalaced samples.
#' 0 means all controls and 1 mean all cases.
#' @param sv.db batches
#' @param sv.ns batches
#' @param sd.b a numeric for sd.b?
#' @param sd.gam a numeric for sd.gam?
#' @param sd.u a numeric for sd.u?
#' @param conf a flag for conf?
#' @param distr.db a numeric for distr.db?
#' @param p.b a numeric for p.b?
#' @param p.gam a numeric for p.gam?
#' @param p.ov a numeric for p.ov?
#' @return A list with:
#' \describe{
#'   \item{db}{database creation variables}
#'   \item{vars}{variables used in simulation}
#' }
#' @references
#' Leek,J.T. and Storey,J.D. (2007) Capturing heterogeneity in gene expression
#' studies by surrogate variable analysis. PLoS Genet., 3, 1724–1735
#' @family simulation
#' @export
##########################################################################################
createMainEffects <- function(n.e = 1000,
                              n.db = 70,
                              n.ns = 30,
                              label = "class",
                              pct.imbalance = 0.5,
                              sv.db = c("A", "B"),
                              sv.ns = c("A", "B"),
                              sd.b = 1,
                              sd.gam = 1,
                              sd.u = 1,
                              conf = FALSE,
                              distr.db = NA,
                              p.b = 0.3,
                              p.gam = 0.3,
                              p.ov = p.b / 2) {
  if (!any(label == c("class", "qtrait"))) {
    stop("CreateMainEffects: name of the label should be class or qtrait")
  }
  n <- n.db + n.ns
  # Create random error
  U <- matrix(nrow = n.e, ncol = n, rnorm(n.e * n, sd = sd.u))

  # Create index for database vs. new sample #
  ind <- as.factor(c(rep("db", n.db), rep("ns", n.ns)))

  if (label == "class") {
    # Create outcome, surrogate variables #
    # Use distr option to show % overlap of outcome, surrogate variables.
    # Note that .5 means no confounding between outcome, surrogate variables.

    # biological variable (fixed at 50% for each outcome)

    ##############################
    ##### imbalanced ability #####
    # ----------------------------
    S.db <- c(rep(0, round(pct.imbalance * n.db)), rep(1, round((1 - pct.imbalance) * n.db)))
    S.ns <- c(rep(0, round(pct.imbalance * n.ns)), rep(1, round((1 - pct.imbalance) * n.ns)))
    S <- c(S.db, S.ns)

    len0 <- sum(S.db == 0)
    len1 <- sum(S.db == 1)

    if (conf == FALSE) {
      # surrogate variable (no confounding in this function)
      n.sv.db <- length(sv.db)
      prop.db <- 1 / n.sv.db
      # create surrogate variables for outcome 0 in database
      x1 <- c()
      for (i in 1:n.sv.db) {
        x1 <- c(x1, rep(sv.db[i], floor(prop.db * len0)))
      }
      # If the rounding has caused a problem, randomly assign to fill out vector
      while (length(x1) != len0) {
        x1 <- c(x1, sample(sv.db, 1))
      }
      # surrogate variables for outcome 1 will NOT be the same
      x2 <- c()
      for (i in 1:n.sv.db) {
        x2 <- c(x2, rep(sv.db[i], floor(prop.db * len1)))
      }
      # If the rounding has caused a problem, randomly assign to fill out vector
      while (length(x2) != len1) {
        x2 <- c(x2, sample(sv.db, 1))
      }
    }

    if (conf) {
      x1 <- c(
        rep("A", round(distr.db * len0)),
        rep("B", len0 - round(distr.db * len0))
      )
      x2 <- c(
        rep("A", round((1 - distr.db) * len1)),
        rep("B", len1 - round((1 - distr.db) * len1))
      )
    }

    # create surrogate variables for outcome 0 in new samples
    n.sv.ns <- length(sv.ns)
    prop.ns <- 1 / n.sv.ns

    len0 <- sum(S.ns == 0)
    len1 <- sum(S.ns == 1)

    x3 <- c()
    for (i in 1:n.sv.ns) {
      x3 <- c(x3, rep(sv.ns[i], floor(prop.ns * len0)))
    }
    # If the rounding has caused a problem, randomly assign to fill out vector
    while (length(x3) != len0) {
      x3 <- c(x3, sample(sv.ns, 1))
    }

    # surrogate variables for outcome 1 will NOT be the same
    x4 <- c()
    for (i in 1:n.sv.ns) {
      x4 <- c(x4, rep(sv.ns[i], floor(prop.ns * len1)))
    }
    # If the rounding has caused a problem, randomly assign to fill out vector
    while (length(x4) != len1) {
      x4 <- c(x4, sample(sv.ns, 1))
    }

    G <- c(x1, x2, x3, x4)
    G <- t(stats::model.matrix(~ as.factor(G)))[-1, ]
    if (is.null(dim(G))) {
      G <- matrix(G, nrow = 1, ncol = n)
    }
    # Determine which probes are affected by what:
    # 30% for biological, 30% for surrogate, 10% overlap
    # First 30% of probes will be affected by biological signal
    ind.B <- rep(0, n.e)
    ind.B[1:round(p.b * n.e)] <- 1
    # Probes 20% thru 50% will be affected by surrogate variable
    ind.Gam <- rep(0, n.e)
    ind.Gam[round((p.b - p.ov) * n.e):round((p.b - p.ov + p.gam) * n.e)] <- 1

    # figure out dimensions for Gamma

    # create parameters for signal, noise
    B <- matrix(nrow = n.e, ncol = 1, rnorm(n.e, mean = 0, sd = sd.b) * ind.B)
    Gam <- matrix(
      nrow = n.e, ncol = dim(G)[1],
      rnorm(n.e * dim(G)[1], mean = 0, sd = sd.gam) * ind.Gam
    )

    # simulate the data
    sim.dat <- B %*% S + Gam %*% G + U
    sim.dat <- sim.dat + abs(min(sim.dat)) + 0.0001

    # simulate data without batch effects
    sim.dat.nobatch <- B %*% S + U
    sim.dat.nobatch <- sim.dat.nobatch + abs(min(sim.dat)) + 0.0001

    # divide parts into database, new samples
    db <- list()
    db$dat <- sim.dat[, ind == "db"]
    db$datnobatch <- sim.dat.nobatch[, ind == "db"]
    db$U <- U[, ind == "db"]
    db$B <- B
    db$S <- S[ind == "db"]
    db$Gam <- Gam
    db$G <- G[ind == "db"]
  } else if (label == "qtrait") {
    ##########################
    # ------- Comment --------
    # Since, we have quantitative outcome and not dichotomize outcome, we may not be able to include conf.
    # Also, simulation algorithm return a simulate data without batch effects, so we do not also need Gam.
    ##########################
    S.db <- rnorm(n.db, 0, 1)
    S.ns <- rnorm(n.ns, 0, 1)
    S <- c(S.db, S.ns)

    # Determine which probes are affected by what:
    # 30% for biological, 30% for surrogate, 10% overlap
    # First 30% of probes will be affected by biological signal
    ind.B <- rep(0, n.e)
    ind.B[1:round(p.b * n.e)] <- 1
    # Probes 20% thru 50% will be affected by surrogate variable
    ind.Gam <- rep(0, n.e)
    ind.Gam[round((p.b - p.ov) * n.e):round((p.b - p.ov + p.gam) * n.e)] <- 1

    # figure out dimensions for Gamma

    # create parameters for signal, noise
    B <- matrix(nrow = n.e, ncol = 1, rnorm(n.e, mean = 0, sd = sd.b) * ind.B)
    # Gam <- matrix(nrow = n.e, ncol = dim(G)[1],
    #               stats::rnorm(n.e * dim(G)[1], mean = 0, sd = sd.gam) * ind.Gam)
    #
    # # simulate the data
    # sim.dat <- B %*% S + Gam %*% G + U
    # sim.dat <- sim.dat + abs(min(sim.dat)) + 0.0001

    # simulate data without batch effects
    # since, there is no G and Gam matrix for quantitative outcome, we use without batch effects.
    sim.dat.nobatch <- B %*% S + U
    sim.dat.nobatch <- sim.dat.nobatch + abs(min(sim.dat.nobatch)) + 0.0001

    # divide parts into database, new samples
    db <- list()
    # db$dat <- sim.dat[, ind == "db"]
    db$datnobatch <- sim.dat.nobatch[, ind == "db"]
    db$U <- U[, ind == "db"]
    db$B <- B
    db$S <- S[ind == "db"]
    # db$Gam <- Gam
    # db$G <- G[ind == "db"]
  }

  vars <- list(
    n.e = n.e, n.db = n.db, n.ns = n.ns, sv.db = sv.db,
    sv.ns = sv.ns, sd.b = sd.b, sd.gam = sd.gam, sd.u = sd.u,
    conf = conf, distr.db = distr.db, p.b = p.b, p.gam = p.gam,
    p.ov = p.ov
  )

  list(db = db, vars = vars)
}
##########################################################################################

# =========================================================================================#
#' createSimulation2
#' update to createSimulation2() that allows for rs-fMRI data generation.
#'
#' Changes: only an additional input with conditional if/else statements
#'          wherever an igraph object is generated. Will allow for
#'          user-supplied graph object as well, this functionality requires
#'          an igraph object specifically for now.
#'
#' New Parameter:
#'
#'    graph.structure - (igraph) graph structure generated from igraph. If NULL,
#'                      then createSimulation2() typically generates one on its own.
#'                      However, this will be handled by createfMRIsimulation() and
#'                      the user will not be using createSimulation2() directly. This
#'                      input is used to circumvent automatic simulation of graphs
#'                      by the previous iteration of createSimulation2().
#'
#' @param num.samples number of samples
#' @param num.variables number of variables (features)
#' @param pct.imbalance fraction of num.samples that are cases
#' @param pct.signals fraction of num.variables that are functional
#' @param main.bias approximate effect size for main effect simulations
#' @param interaction.bias approximate effect size for interaction effects
#' @param hi.cor parameter to use for network-connected pairwise correlations
#' @param lo.cor parameter to use for network-non-connected pairwise correlations
#' @param label should just be "class" for binary response
#' @param sim.type a character that determines the type of simulation:
#' mainEffect/mainEffect_Erdos-Renyi/mainEffect_Scalefree/interactionErdos/interactionScalefree/mixed
#' @param pct.train fraction of num.samples used for training
#' @param pct.holdout fraction of num.samples used for holdout
#' @param pct.validation fraction of num.samples used for validation
#' @param save.file logical but not currently being used
#' @param mix.type character that determines the type of mixed effects simulation:
#' main-interactionErdos/main-interactionScalefree
#' @param pct.mixed percent of functional variables that are main effects (1 - pct_interaction). Use with sim.type="mixed" and specify mix.type.
#' @param verbose logical indicating whether to display time required to generate simulation
#' @param plot.graph logical indicating whether to plot networks
#' @param use.Rcpp if true use Rcpp to correct negative eigenvalues
#' @param prob.connected probability of drawing an edge between two arbitrary vertices in Erdos-Renyi graph
#' @param out.degree out-degree of vertices in Scale-free graph
#' @param data.type character indicating if data is from a "continuous" or "discrete" distribution
#' @param avg.maf numeric in (0,1) indicating the desired average MAF for GWAS main effect simulations
#' @param graph.structure (igraph) graph structure generated from igraph, NULL default
#' @return A list with:
#' \describe{
#'   \item{train}{traing data set}
#'   \item{holdout}{holdout data set}
#'   \item{validation}{validation data set}
#'   \item{label}{the class label/column name}
#'   \item{signal.names}{the variable names with simulated signals}
#'   \item{elapsed}{total elapsed time}
#'   \item{A.mat}{adjacency matrix for random network (for interaction/mixed/main plus correlation)}
#'   \item{main.vars}{indices of main effect variables (for mixed simulations)}
#'   \item{int.vars}{indices of interaction effect variables (for mixed simulations)}
#' }
#' @examples
#'
#' # main effects with correlation from Erdos-Renyi network (continuous)
#'
#' num.samples <- 100
#' num.variables <- 10
#'
#' dataset <- createSimulation2(
#'   num.samples = num.samples,
#'   num.variables = num.variables,
#'   pct.imbalance = 0.5,
#'   pct.signals = 0.2,
#'   main.bias = 0.4,
#'   interaction.bias = 1,
#'   hi.cor = 0.8,
#'   lo.cor = 0.2,
#'   mix.type = "main-interactionErdos",
#'   label = "class",
#'   sim.type = "mainEffect_Erdos-Renyi",
#'   pct.mixed = 0.5,
#'   pct.train = 0.5,
#'   pct.holdout = 0.5,
#'   pct.validation = 0,
#'   plot.graph = FALSE,
#'   verbose = TRUE,
#'   data.type = "continuous"
#' )
#'
#' # main effects with correlation from Erdos-Renyi network
#'
#' num.samples <- 100
#' num.variables <- 10
#'
#' dataset <- createSimulation2(
#'   num.samples = num.samples,
#'   num.variables = num.variables,
#'   pct.imbalance = 0.5,
#'   pct.signals = 0.2,
#'   main.bias = 0.4,
#'   interaction.bias = 1,
#'   hi.cor = 0.8,
#'   lo.cor = 0.2,
#'   mix.type = "main-interactionErdos",
#'   label = "class",
#'   sim.type = "mainEffect_Erdos-Renyi",
#'   pct.mixed = 0.5,
#'   pct.train = 0.5,
#'   pct.holdout = 0.5,
#'   pct.validation = 0,
#'   plot.graph = FALSE,
#'   verbose = TRUE,
#'   data.type = "discrete"
#' )
#' @export
createSimulation2 <- function(num.samples = 100,
                              num.variables = 100,
                              pct.imbalance = 0.5,
                              pct.signals = 0.1,
                              main.bias = 0.4,
                              interaction.bias = 0.4,
                              hi.cor = 0.8,
                              lo.cor = 0.2,
                              label = "class",
                              sim.type = "mainEffect",
                              pct.train = 0.5,
                              pct.holdout = 0.5,
                              pct.validation = 0,
                              save.file = NULL,
                              mix.type = NULL,
                              pct.mixed = 0.5,
                              verbose = FALSE,
                              plot.graph = FALSE, use.Rcpp = FALSE,
                              prob.connected = NULL,
                              out.degree = NULL,
                              data.type = "continuous",
                              avg.maf = 0.2,
                              graph.structure = NULL) {
  ptm <- proc.time() # start time

  if (use.Rcpp) {
    Rcpp::cppFunction(
      depends = "RcppArmadillo",
      "arma::vec getEigenValues(arma::mat M) {
                      return arma::eig_sym(M);
  }"
    )
  }

  nbias <- pct.signals * num.variables # number of functional attributes

  if (sim.type == "mainEffect") { # simple main effect simulation

    if (data.type == "continuous") {

      # new simulation:
      # sd.b sort of determines how large the signals are
      # p.b=0.1 makes 10% of the variables signal, bias <- 0.5
      my.sim.data <- createMainEffects(
        n.e = num.variables,
        n.db = num.samples,
        pct.imbalance = pct.imbalance,
        label = label,
        sd.b = main.bias,
        p.b = pct.signals
      )$db
      dataset <- cbind(t(my.sim.data$datnobatch), my.sim.data$S)

      # make numeric matrix into a data frame for splitting and subsequent ML algorithms
      dataset <- as.data.frame(dataset)

      signal.names <- paste("mainvar", 1:nbias, sep = "") # functional names
      background.names <- paste("var", 1:(num.variables - nbias), sep = "") # non-functional names
      var.names <- c(signal.names, background.names, label) # all variable names
      colnames(dataset) <- var.names # replace column names with variable names
    } else if (data.type == "discrete") {
      if (main.bias >= 1 || main.bias <= 0) {
        stop("Please choose main.bias in (0,1)")
      }

      if (avg.maf <= 0 || avg.maf >= 0.5) {
        stop("Please choose avg.maf in (0,0.5)")
      }

      m.case <- round((1 - pct.imbalance) * num.samples) # size of case group
      m.ctrl <- round(pct.imbalance * num.samples) # size of ctrl group

      lb <- avg.maf - 0.5 * main.bias # potentially the average maf for control group
      ub <- avg.maf + 0.5 * main.bias # potentially the average maf for the case group

      if (lb <= 0) { # check to see if lower maf is possible, if not just make it small positive number

        eps <- 0.01

        lb <- eps

        ub <- 2 * avg.maf - eps # chosen so that 0.5*(lb + ub) = avg.maf
      } else if (ub >= 1) { # check to see if upper maf is possible, if not just make it a number slightly less than 1

        eps <- 0.01

        ub <- 1 - eps

        lb <- 2 * avg.maf - (1 - eps) # chosen so that 0.5*(lb + ub) = avg.maf
      }

      prob.case <- ub # case group mean
      prob.ctrl <- lb # ctrl group mean
      prob.noise <- runif((num.variables - nbias), min = prob.ctrl, max = prob.case) # background mafs

      # case group functional attributes
      main.case <- matrix(0, nrow = m.case, ncol = nbias)
      for (j in 1:m.case) {
        tmp <- rbinom(n = nbias, size = 2, prob = prob.case)

        main.case[j, ] <- tmp
      }

      # ctrl group functional attributes
      main.ctrl <- matrix(0, nrow = m.ctrl, ncol = nbias)
      for (j in 1:m.ctrl) {
        tmp <- rbinom(n = nbias, size = 2, prob = prob.ctrl)

        main.ctrl[j, ] <- tmp
      }

      # main effect attributes
      main.effects <- rbind(main.case, main.ctrl)

      # background attributes
      background.data <- matrix(rbinom(n = (num.variables - nbias) * num.samples, size = 2, prob = prob.noise),
        nrow = num.samples
      )

      class.vec <- c(rep(1, length = m.case), rep(0, length = m.ctrl))
      dataset <- cbind(main.effects, background.data)
      dataset <- cbind(dataset, class.vec)

      # make numeric matrix into a data frame for splitting and subsequent ML algorithms
      dataset <- as.data.frame(dataset)

      signal.names <- paste("mainvar", 1:nbias, sep = "") # functional names
      background.names <- paste("var", 1:(num.variables - nbias), sep = "") # non-functional names
      var.names <- c(signal.names, background.names, label) # all variable names
      colnames(dataset) <- var.names # replace column names with variable names
    }
  } else if (sim.type == "mainEffect_Erdos-Renyi") { # main + Erdos-Renyi network

    if (data.type == "continuous") {

      # create main-effect simulation
      my.sim.data <- createMainEffects(
        n.e = num.variables,
        n.db = num.samples,
        pct.imbalance = pct.imbalance,
        label = label,
        sd.b = main.bias,
        p.b = pct.signals
      )$db
      dataset <- cbind(t(my.sim.data$datnobatch), my.sim.data$S)

      e <- 1 # fudge factor to the number of nodes to avoid giant component
      if (is.null(prob.connected)) {
        prob <- 1 / (num.variables + e) # probability of a node being connected to another node is less than 1/N to avoid giant component
      } else {
        prob <- prob.connected
      }

      if (is.null(graph.structure)) {
        g <- igraph::erdos.renyi.game(num.variables, prob) # Erdos-Renyi network
      } else {
        g <- graph.structure
      }

      # generate correlation matrix from g
      network.atts <- generate_structured_corrmat(
        g = g,
        num.variables = num.variables,
        hi.cor.tmp = hi.cor,
        lo.cor.tmp = lo.cor,
        hi.cor.fixed = hi.cor,
        graph.type = "Erdos-Renyi",
        plot.graph = plot.graph,
        nbias = nbias, use.Rcpp = use.Rcpp
      )

      R <- as.matrix(network.atts$corrmat) # correlation matrix

      A.mat <- network.atts$A.mat # adjacency from graph

      U <- t(chol(R)) # upper tri cholesky
      tmp <- t(U %*% t(dataset[, -ncol(dataset)])) # correlated data
      tmp <- cbind(tmp, dataset[, ncol(dataset)]) # combine with phenotype
      dataset <- tmp # main-effect data with correlation

      # make numeric matrix into a data frame for splitting and subsequent ML algorithms
      dataset <- as.data.frame(dataset)

      signal.names <- paste("mainvar", 1:nbias, sep = "") # functional variable names
      background.names <- paste("var", 1:(num.variables - nbias), sep = "") # non-functional variable names
      var.names <- c(signal.names, background.names, label) # all variable names
      colnames(dataset) <- var.names # replace column names with variable names
    } else if (data.type == "discrete") {
      if (main.bias >= 1) {
        stop("Please choose main.bias in (0,1)")
      }

      if (avg.maf <= 0 || avg.maf >= 0.5) {
        stop("Please choose avg.maf in (0,0.5)")
      }

      m.case <- round((1 - pct.imbalance) * num.samples) # size of case group
      m.ctrl <- round(pct.imbalance * num.samples) # size of ctrl group


      null.dats <- matrix(rnorm(num.samples * num.variables), nrow = num.samples, ncol = num.variables)

      e <- 1 # fudge factor to the number of nodes to avoid giant component
      if (is.null(prob.connected)) {
        prob <- 1 / (num.variables + e) # probability of a node being connected to another node is less than 1/N to avoid giant component
      } else {
        prob <- prob.connected
      }

      if (is.null(graph.structure)) {
        g <- igraph::erdos.renyi.game(num.variables, prob) # Erdos-Renyi network
      } else {
        g <- graph.structure
      }

      # generate correlation matrix from g
      network.atts <- generate_structured_corrmat(
        g = g,
        num.variables = num.variables,
        hi.cor.tmp = hi.cor,
        lo.cor.tmp = lo.cor,
        hi.cor.fixed = hi.cor,
        graph.type = "Erdos-Renyi",
        plot.graph = plot.graph,
        nbias = nbias, use.Rcpp = use.Rcpp
      )

      R <- as.matrix(network.atts$corrmat) # correlation matrix

      A.mat <- network.atts$A.mat # adjacency from graph

      U <- t(chol(R)) # upper tri cholesky
      tmp <- t(U %*% t(null.dats)) # correlated data
      tmp <- cbind(tmp, class = c(rep(1, length = 50), rep(0, length = 50))) # combine with phenotype
      dataset <- tmp # main-effect data with correlation

      trans.dats <- pnorm(dataset[, -ncol(dataset)])
      case.dats <- trans.dats[c(1:m.case), ]
      ctrl.dats <- trans.dats[c((m.case + 1):nrow(trans.dats)), ]

      lb <- avg.maf - 0.5 * main.bias # potentially the average maf for control group
      ub <- avg.maf + 0.5 * main.bias # potentially the average maf for the case group

      if (lb <= 0) { # check to see if lower maf is possible, if not just make it small positive number

        eps <- 0.01

        lb <- eps

        ub <- 2 * avg.maf - eps # chosen so that 0.5*(lb + ub) = avg.maf
      } else if (ub >= 1) { # check to see if upper maf is possible, if not just make it a number slightly less than 1

        eps <- 0.01

        ub <- 1 - eps

        lb <- 2 * avg.maf - (1 - eps) # chosen so that 0.5*(lb + ub) = avg.maf
      }

      prob.case <- ub # case group mean
      prob.ctrl <- lb # ctrl group mean
      prob.noise <- runif((num.variables - nbias), min = prob.ctrl, max = prob.case) # background mafs

      # functional and background features for cases
      main.case <- matrix(0, nrow = m.case, ncol = nbias)
      background.case <- matrix(0, nrow = m.case, ncol = length(c((nbias + 1):num.variables)))
      for (j in 1:m.case) {
        tmp <- qbinom(case.dats[j, c(1:nbias)], size = 2, prob = prob.case)

        main.case[j, ] <- tmp
        background.case[j, ] <- qbinom(case.dats[j, c((nbias + 1):num.variables)], size = 2, prob = prob.noise)
      }

      # functional and background features for ctrls
      main.ctrl <- matrix(0, nrow = m.ctrl, ncol = nbias)
      background.ctrl <- matrix(0, nrow = m.ctrl, ncol = length(c((nbias + 1):num.variables)))
      for (j in 1:m.ctrl) {
        tmp <- qbinom(ctrl.dats[j, c(1:nbias)], size = 2, prob = prob.ctrl)

        main.ctrl[j, ] <- tmp
        background.ctrl[j, ] <- qbinom(ctrl.dats[j, c((nbias + 1):num.variables)], size = 2, prob = prob.noise)
      }

      # main effect features
      main.effects <- rbind(main.case, main.ctrl)

      # background features
      background.data <- rbind(background.case, background.ctrl)

      class.vec <- c(rep(1, length = m.case), rep(0, length = m.ctrl))
      dataset <- cbind(main.effects, background.data)
      dataset <- cbind(dataset, class.vec)

      # make numeric matrix into a data frame for splitting and subsequent ML algorithms
      dataset <- as.data.frame(dataset)

      signal.names <- paste("mainvar", 1:nbias, sep = "") # functional names
      background.names <- paste("var", 1:(num.variables - nbias), sep = "") # non-functional names
      var.names <- c(signal.names, background.names, label) # all variable names
      colnames(dataset) <- var.names # replace column names with variable names
    }
  } else if (sim.type == "mainEffect_Scalefree") { # main + Scale-Free network

    if (data.type == "continuous") {

      # create main-effect simulation
      my.sim.data <- createMainEffects(
        n.e = num.variables,
        n.db = num.samples,
        pct.imbalance = pct.imbalance,
        label = label,
        sd.b = main.bias,
        p.b = pct.signals
      )$db
      dataset <- cbind(t(my.sim.data$datnobatch), my.sim.data$S)

      e <- 1 # fudge factor to the number of nodes to avoid giant component
      prob <- 1 / (num.variables + e) # probability of a node being connected to another node is less than 1/N to avoid giant component

      if (is.null(out.degree)) {
        g <- igraph::barabasi.game(num.variables, directed = FALSE) # scale-free network
      } else {
        g <- igraph::barabasi.game(num.variables, m = out.degree, directed = FALSE)
      }

      if (!is.null(graph.structure)) {
        g <- graph.structure
      }

      # generate correlation matrix from g
      network.atts <- generate_structured_corrmat(
        g = g,
        num.variables = num.variables,
        hi.cor.tmp = hi.cor,
        lo.cor.tmp = lo.cor,
        hi.cor.fixed = hi.cor,
        graph.type = "Scalefree",
        plot.graph = plot.graph,
        nbias = nbias, use.Rcpp = use.Rcpp
      )
      R <- as.matrix(network.atts$corrmat) # correlation matrix

      A.mat <- network.atts$A.mat # adjacency from graph

      U <- t(chol(R)) # upper tri cholesky
      tmp <- t(U %*% t(dataset[, -ncol(dataset)])) # correlated data
      tmp <- cbind(tmp, dataset[, ncol(dataset)]) # combine with phenotype
      dataset <- tmp # main-effect data with correlation

      # make numeric matrix into a data frame for splitting and subsequent ML algorithms
      dataset <- as.data.frame(dataset)

      signal.names <- paste("mainvar", 1:nbias, sep = "") # functional variable names
      background.names <- paste("var", 1:(num.variables - nbias), sep = "") # non-functional variable names
      var.names <- c(signal.names, background.names, label) # all variable names
      colnames(dataset) <- var.names # replace column names with variable names
    } else if (data.type == "discrete") {
      if (main.bias >= 1) {
        stop("Please choose main.bias in (0,1)")
      }

      if (avg.maf <= 0 || avg.maf >= 0.5) {
        stop("Please choose avg.maf in (0,0.5)")
      }

      m.case <- round((1 - pct.imbalance) * num.samples) # size of case group
      m.ctrl <- round(pct.imbalance * num.samples) # size of ctrl group


      null.dats <- matrix(rnorm(num.samples * num.variables), nrow = num.samples, ncol = num.variables)

      e <- 1 # fudge factor to the number of nodes to avoid giant component
      prob <- 1 / (num.variables + e) # probability of a node being connected to another node is less than 1/N to avoid giant component

      if (is.null(out.degree)) {
        g <- igraph::barabasi.game(num.variables, directed = FALSE) # scale-free network
      } else {
        g <- igraph::barabasi.game(num.variables, m = out.degree, directed = FALSE)
      }

      if (!is.null(graph.structure)) {
        g <- graph.structure
      }

      # generate correlation matrix from g
      network.atts <- generate_structured_corrmat(
        g = g,
        num.variables = num.variables,
        hi.cor.tmp = hi.cor,
        lo.cor.tmp = lo.cor,
        hi.cor.fixed = hi.cor,
        graph.type = "Scalefree",
        plot.graph = plot.graph,
        nbias = nbias, use.Rcpp = use.Rcpp
      )

      R <- as.matrix(network.atts$corrmat) # correlation matrix

      A.mat <- network.atts$A.mat # adjacency from graph

      U <- t(chol(R)) # upper tri cholesky
      tmp <- t(U %*% t(null.dats)) # correlated data
      tmp <- cbind(tmp, class = c(rep(1, length = 50), rep(0, length = 50))) # combine with phenotype
      dataset <- tmp # main-effect data with correlation

      trans.dats <- pnorm(dataset[, -ncol(dataset)])
      case.dats <- trans.dats[c(1:m.case), ]
      ctrl.dats <- trans.dats[c((m.case + 1):nrow(trans.dats)), ]

      lb <- avg.maf - 0.5 * main.bias # potentially the average maf for control group
      ub <- avg.maf + 0.5 * main.bias # potentially the average maf for the case group

      if (lb <= 0) { # check to see if lower maf is possible, if not just make it small positive number

        eps <- 0.01

        lb <- eps

        ub <- 2 * avg.maf - eps # chosen so that 0.5*(lb + ub) = avg.maf
      } else if (ub >= 1) { # check to see if upper maf is possible, if not just make it a number slightly less than 1

        eps <- 0.01

        ub <- 1 - eps

        lb <- 2 * avg.maf - (1 - eps) # chosen so that 0.5*(lb + ub) = avg.maf
      }

      prob.case <- ub # case group mean
      prob.ctrl <- lb # ctrl group mean
      prob.noise <- runif((num.variables - nbias), min = prob.ctrl, max = prob.case) # background mafs

      # functional and background features for cases
      main.case <- matrix(0, nrow = m.case, ncol = nbias)
      background.case <- matrix(0, nrow = m.case, ncol = length(c((nbias + 1):num.variables)))
      for (j in 1:m.case) {
        tmp <- qbinom(case.dats[j, c(1:nbias)], size = 2, prob = prob.case)

        main.case[j, ] <- tmp
        background.case[j, ] <- qbinom(case.dats[j, c((nbias + 1):num.variables)], size = 2, prob = prob.noise)
      }

      # functional and background features for ctrls
      main.ctrl <- matrix(0, nrow = m.ctrl, ncol = nbias)
      background.ctrl <- matrix(0, nrow = m.ctrl, ncol = length(c((nbias + 1):num.variables)))
      for (j in 1:m.ctrl) {
        tmp <- qbinom(ctrl.dats[j, c(1:nbias)], size = 2, prob = prob.ctrl)

        main.ctrl[j, ] <- tmp
        background.ctrl[j, ] <- qbinom(ctrl.dats[j, c((nbias + 1):num.variables)], size = 2, prob = prob.noise)
      }

      # main effect features
      main.effects <- rbind(main.case, main.ctrl)

      # background features
      background.data <- rbind(background.case, background.ctrl)

      class.vec <- c(rep(1, length = m.case), rep(0, length = m.ctrl))
      dataset <- cbind(main.effects, background.data)
      dataset <- cbind(dataset, class.vec)

      # make numeric matrix into a data frame for splitting and subsequent ML algorithms
      dataset <- as.data.frame(dataset)

      signal.names <- paste("mainvar", 1:nbias, sep = "") # functional names
      background.names <- paste("var", 1:(num.variables - nbias), sep = "") # non-functional names
      var.names <- c(signal.names, background.names, label) # all variable names
      colnames(dataset) <- var.names # replace column names with variable names
    }
  } else if (sim.type == "interactionErdos") { # simple interaction with Erdos-Renyi network

    if (data.type == "continuous") {
      m.case <- round((1 - pct.imbalance) * num.samples) # size of case group
      m.ctrl <- round(pct.imbalance * num.samples) # size of ctrl group

      # random cases data set
      X.case <- matrix(rnorm(m.case * num.variables), nrow = m.case, ncol = num.variables)

      # approximate effect size for pairwise correlations between functional attributes in case group
      case.hi.cor <- -hi.cor * interaction.bias + (1 - interaction.bias) * hi.cor

      e <- 1 # fudge factor to the number of nodes to avoid giant component
      if (is.null(prob.connected)) {
        prob <- 1 / (num.variables + e) # probability of a node being connected to another node is less than 1/N to avoid giant component
      } else {
        prob <- prob.connected
      }

      if (is.null(graph.structure)) {
        g <- igraph::erdos.renyi.game(num.variables, prob) # Erdos-Renyi network
      } else {
        g <- graph.structure
      }

      # generate correlation matrix from g
      network.atts <- generate_structured_corrmat(
        g = g,
        num.variables = num.variables,
        hi.cor.tmp = case.hi.cor,
        lo.cor.tmp = lo.cor,
        hi.cor.fixed = hi.cor,
        graph.type = "Erdos-Renyi",
        plot.graph = plot.graph,
        make.diff.cors = TRUE,
        nbias = nbias, use.Rcpp = use.Rcpp
      )

      R <- as.matrix(network.atts$corrmat) # correlation matrix for cases

      sig.vars <- network.atts$sig.vars # column indices of functional features

      A.mat <- network.atts$A.mat # adjacency from graph object

      U <- t(chol(R)) # upper tri cholesky
      tmp <- t(U %*% t(X.case)) # correlated case data
      X.case <- tmp

      # random controls data set
      X.ctrl <- matrix(rnorm(m.ctrl * num.variables), nrow = m.ctrl, ncol = num.variables)

      # for each functional variable, find connected variables
      sig.connected.list <- list() # list for storing connected variable indices
      for (i in 1:length(sig.vars)) { # for each functional feature
        sig.connected.list[[i]] <- which(abs(A.mat[sig.vars[i], ] - 1) < 1e-9) # find and store connected features
      }

      # create control group correlation matrix
      for (i in 1:length(sig.vars)) { # for each functional feature
        n.connected <- length(sig.connected.list[[i]]) # how many connections
        for (j in 1:n.connected) { # for each connection
          R[sig.vars[i], sig.connected.list[[i]][j]] <- min(hi.cor + rnorm(1, 0, .1), 1) # replace low correlation by high for ctrls
        }
      }

      # ensure that R is symmetric
      #
      tmp <- diag(R) # diag of R
      diag(R) <- 0 # replace diag of R with 0's
      R[lower.tri(R)] <- 0 # make lower triangle zeros

      R <- R + t(R) # make symmetric
      diag(R) <- tmp # original diag of R

      # correct for negative eigenvalues so R is positive definite
      #
      if (use.Rcpp) { # compute eigenvalues and make diag matrix
        R.d <- diag(sort(c(getEigenValues(R)), decreasing = TRUE))
      } else {
        R.d <- diag(eigen(R)$values)
      }
      tmp <- diag(R.d) # vector of eigenvalues

      if (any(tmp < 0)) { # if any eigenvalues are negative
        R.V <- eigen(R)$vectors # compute eigenvectors,
        tmp[tmp < 0] <- 1e-7 # make negative into small positive,
        diag(R.d) <- tmp # replace diag of R.d with positive eigenvalues,
        R.fix <- R.V %*% R.d %*% t(R.V) # compute new R, and
        R <- R.fix # store in R
      }
      R <- as.matrix(R)

      # make 1's on diag of R
      inv.diag <- 1 / diag(R) # multiplicative inverse of diag(R)
      mydiag <- diag(length(inv.diag)) # initialize diagonal matrix for inv.diag
      diag(mydiag) <- inv.diag # swap 1's for inv.diag
      mydiag <- sqrt(mydiag) # compute sqrt of inv.diag
      R <- mydiag %*% R %*% mydiag # compute corrected R with 1's on diagonal (Still Pos. Def.)

      U <- t(chol(R)) # upper tri cholesky
      tmp <- t(U %*% t(X.ctrl)) # correlated ctrl data
      X.ctrl <- tmp # store in X.ctrl

      X.all <- rbind(X.case, X.ctrl) # case/ctrl data
      X.all <- as.data.frame(X.all) # make data.frame

      dataset <- X.all

      signal.names <- paste("intvar", 1:nbias, sep = "") # interaction variable names
      background.names <- paste("var", 1:(num.variables - nbias), sep = "") # non-functional variable names
      var.names <- c(signal.names, background.names, label) # all variable names
      colnames(dataset)[sig.vars] <- signal.names # replace functional colnames with interaction variable names
      colnames(dataset)[-sig.vars] <- background.names # replace non-functional colnames with non-functional colnames

      dataset <- cbind(dataset, c(rep(1, length = m.case), rep(-1, length = m.ctrl))) # full data set with phenotype
      colnames(dataset)[ncol(dataset)] <- "class" # replace last colname with "class"
    } else if (data.type == "discrete") {
      if (avg.maf <= 0 || avg.maf >= 0.5) {
        stop("Please choose avg.maf in (0,0.5)")
      }

      m.case <- round((1 - pct.imbalance) * num.samples) # size of case group
      m.ctrl <- round(pct.imbalance * num.samples) # size of ctrl group

      # approximate effect size for pairwise correlations between functional attributes in case group
      case.hi.cor <- -hi.cor * interaction.bias + (1 - interaction.bias) * hi.cor

      e <- 1 # fudge factor to the number of nodes to avoid giant component
      if (is.null(prob.connected)) {
        prob <- 1 / (num.variables + e) # probability of a node being connected to another node is less than 1/N to avoid giant component
      } else {
        prob <- prob.connected
      }

      if (is.null(graph.structure)) {
        g <- igraph::erdos.renyi.game(num.variables, prob) # Erdos-Renyi network
      } else {
        g <- graph.structure
      }

      # generate correlation matrix from g
      network.atts <- generate_structured_corrmat(
        g = g,
        num.variables = num.variables,
        hi.cor.tmp = case.hi.cor,
        lo.cor.tmp = lo.cor,
        hi.cor.fixed = hi.cor,
        graph.type = "Erdos-Renyi",
        plot.graph = plot.graph,
        make.diff.cors = TRUE,
        nbias = nbias, use.Rcpp = use.Rcpp
      )

      R <- as.matrix(network.atts$corrmat) # correlation matrix
      R.case <- R

      null.dats <- matrix(rnorm(num.variables * m.case), nrow = m.case, ncol = num.variables)

      U <- t(chol(R.case)) # upper tri cholesky
      null.dats <- t(U %*% t(null.dats)) # correlated data

      sig.vars <- network.atts$sig.vars # column indices of functional features

      A.mat <- network.atts$A.mat # adjacency from graph object

      eps <- 0.1
      lb <- eps
      ub <- 2 * avg.maf - eps
      f.a <- runif(num.variables, min = lb, max = ub) # minor allele frequencies (success probabilities for bernoulli trials)

      trans.dats <- pnorm(null.dats)

      bin.data <- matrix(0, nrow = nrow(trans.dats), ncol = ncol(trans.dats))
      for (j in 1:ncol(trans.dats)) {
        tmp <- qbinom(trans.dats[, j], size = 2, prob = f.a[j])

        bin.data[, j] <- tmp
      }
      X.case <- bin.data

      # for each functional variable, find connected variables
      sig.connected.list <- list() # list for storing connected variable indices
      for (i in 1:length(sig.vars)) { # for each functional feature
        sig.connected.list[[i]] <- which(abs(A.mat[sig.vars[i], ] - 1) < 1e-9) # find and store connected features
      }

      # create control group correlation matrix
      for (i in 1:length(sig.vars)) { # for each functional feature
        n.connected <- length(sig.connected.list[[i]]) # how many connections
        for (j in 1:n.connected) { # for each connection
          R[sig.vars[i], sig.connected.list[[i]][j]] <- min(hi.cor + rnorm(1, 0, .1), 1) # replace low correlation by high for ctrls
        }
      }

      # ensure that R is symmetric
      #
      tmp <- diag(R) # diag of R
      diag(R) <- 0 # replace diag of R with 0's
      R[lower.tri(R)] <- 0 # make lower triangle zeros

      R <- R + t(R) # make symmetric
      diag(R) <- tmp # original diag of R

      # correct for negative eigenvalues so R is positive definite
      #
      if (use.Rcpp) { # compute eigenvalues and make diag matrix
        R.d <- diag(sort(c(getEigenValues(R)), decreasing = TRUE))
      } else {
        R.d <- diag(eigen(R)$values)
      }
      tmp <- diag(R.d) # vector of eigenvalues

      if (any(tmp < 0)) { # if any eigenvalues are negative
        R.V <- eigen(R)$vectors # compute eigenvectors,
        tmp[tmp < 0] <- 1e-7 # make negative into small positive,
        diag(R.d) <- tmp # replace diag of R.d with positive eigenvalues,
        R.fix <- R.V %*% R.d %*% t(R.V) # compute new R, and
        R <- R.fix # store in R
      }
      R <- as.matrix(R)

      # make 1's on diag of R
      inv.diag <- 1 / diag(R) # multiplicative inverse of diag(R)
      mydiag <- diag(length(inv.diag)) # initialize diagonal matrix for inv.diag
      diag(mydiag) <- inv.diag # swap 1's for inv.diag
      mydiag <- sqrt(mydiag) # compute sqrt of inv.diag
      R <- mydiag %*% R %*% mydiag # compute corrected R with 1's on diagonal (Still Pos. Def.)
      R.ctrl <- R

      null.dats <- matrix(rnorm(num.variables * m.ctrl), nrow = m.ctrl, ncol = num.variables)

      U <- t(chol(R.ctrl)) # upper tri cholesky
      null.dats <- t(U %*% t(null.dats)) # correlated data

      trans.dats <- pnorm(null.dats)

      bin.data <- matrix(0, nrow = nrow(trans.dats), ncol = ncol(trans.dats))
      for (j in 1:ncol(trans.dats)) {
        tmp <- qbinom(trans.dats[, j], size = 2, prob = f.a[j])

        bin.data[, j] <- tmp
      }
      X.ctrl <- bin.data

      X.all <- rbind(X.case, X.ctrl) # case/ctrl data
      X.all <- as.data.frame(X.all) # make data.frame

      dataset <- X.all

      label <- "class"
      signal.names <- paste("intvar", 1:nbias, sep = "") # interaction variable names
      background.names <- paste("var", 1:(num.variables - nbias), sep = "") # non-functional variable names
      var.names <- c(signal.names, background.names, label) # all variable names
      colnames(dataset)[sig.vars] <- signal.names # replace functional colnames with interaction variable names
      colnames(dataset)[-sig.vars] <- background.names # replace non-functional colnames with non-functional colnames

      dataset <- cbind(dataset, c(rep(1, length = m.case), rep(-1, length = m.ctrl))) # full data set with phenotype
      colnames(dataset)[ncol(dataset)] <- "class" # replace last colname with "class"
    }
  } else if (sim.type == "interactionScalefree") { # simple interaction with Scale-Free network

    if (data.type == "continuous") {
      m.case <- round((1 - pct.imbalance) * num.samples) # size of case group
      m.ctrl <- round(pct.imbalance * num.samples) # size of ctrl group

      # random cases data set
      X.case <- matrix(rnorm(m.case * num.variables), nrow = m.case, ncol = num.variables)

      # approximate effect size for pairwise correlations between functional attributes in case group
      case.hi.cor <- -hi.cor * interaction.bias + (1 - interaction.bias) * hi.cor

      e <- 1 # fudge factor to the number of nodes to avoid giant component
      prob <- 1 / (num.variables + e) # probability of a node being connected to another node is less than 1/N to avoid giant component

      if (is.null(out.degree)) {
        g <- igraph::barabasi.game(num.variables, directed = FALSE) # scale-free network
      } else {
        g <- igraph::barabasi.game(num.variables, m = out.degree, directed = FALSE)
      }

      if (!is.null(graph.structure)) {
        g <- graph.structure
      }

      # generate correlation matrix from g
      network.atts <- generate_structured_corrmat(
        g = g,
        num.variables = num.variables,
        hi.cor.tmp = case.hi.cor,
        lo.cor.tmp = lo.cor,
        hi.cor.fixed = hi.cor,
        graph.type = "Scalefree",
        plot.graph = plot.graph,
        make.diff.cors = TRUE,
        nbias = nbias, use.Rcpp = use.Rcpp
      )

      R <- as.matrix(network.atts$corrmat)

      sig.vars <- network.atts$sig.vars # functional variable indices

      A.mat <- network.atts$A.mat # adjacency from graph

      U <- t(chol(R)) # upper tri cholesky
      tmp <- t(U %*% t(X.case)) # correlated case data
      X.case <- tmp # store in X.case

      # random controls data set
      X.ctrl <- matrix(rnorm(m.ctrl * num.variables), nrow = m.ctrl, ncol = num.variables)

      # for each functional variable, find connected variables
      #
      sig.connected.list <- list() # list for functional connections
      for (i in 1:length(sig.vars)) { # for each functional variable
        sig.connected.list[[i]] <- which(abs(A.mat[sig.vars[i], ] - 1) < 1e-9) # find connections and store in list
      }

      # create control group correlation matrix
      for (i in 1:length(sig.vars)) { # for each functional variable
        n.connected <- length(sig.connected.list[[i]]) # how many connections
        for (j in 1:n.connected) { # for each connection
          R[sig.vars[i], sig.connected.list[[i]][j]] <- min(hi.cor + rnorm(1, 0, .1), 1) # replace low correlations with high for ctrls
        }
      }

      # make sure R is symmetric
      #
      tmp <- diag(R) # diag of R
      diag(R) <- 0 # make diag(R) 0
      R[lower.tri(R)] <- 0 # make lower triangle of R 0

      R <- R + t(R) # make symmetric
      diag(R) <- tmp # replace diag(R) with original diag

      # correct for negative eigenvalues to make R positive definite
      #
      if (use.Rcpp) { # compute eigenvalues and make diag matrix
        R.d <- diag(sort(c(getEigenValues(R)), decreasing = TRUE))
      } else {
        R.d <- diag(eigen(R)$values)
      }
      tmp <- diag(R.d) # vector of eigenvalues

      if (any(tmp < 0)) { # if any eigenvalues are negative
        R.V <- eigen(R)$vectors # compute eigenvectors
        tmp[tmp < 0] <- 1e-7 # make negative into small positive
        diag(R.d) <- tmp # swap negative for positive eigenvalues
        R.fix <- R.V %*% R.d %*% t(R.V) # compute corrected correlation matrix
        R <- R.fix # store in R
      }
      R <- as.matrix(R)

      # make 1's on diagonal
      #
      inv.diag <- 1 / diag(R) # multiplicative inverse of diag(R)
      mydiag <- diag(length(inv.diag)) # initialize diag matrix for inv.diag
      diag(mydiag) <- inv.diag # swap 1's for inv.diag
      mydiag <- sqrt(mydiag) # sqrt of inv.diag
      R <- mydiag %*% R %*% mydiag # compute corrected R with 1's on diagonal (Still Pos. Def.)

      U <- t(chol(R)) # upper tri cholesky
      tmp <- t(U %*% t(X.ctrl)) # correlated ctrl data
      X.ctrl <- tmp # store in X.ctrl

      X.all <- rbind(X.case, X.ctrl) # case/ctrl data
      X.all <- as.data.frame(X.all) # make data.frame

      dataset <- X.all

      signal.names <- paste("intvar", 1:nbias, sep = "") # interaction variable names
      background.names <- paste("var", 1:(num.variables - nbias), sep = "") # non-functional variable names
      var.names <- c(signal.names, background.names, label) # all variable names
      colnames(dataset)[sig.vars] <- signal.names # replace functional colnames with interaction variable names
      colnames(dataset)[-sig.vars] <- background.names # replace non-functional colnames with non-functional variable names

      dataset <- cbind(dataset, c(rep(1, length = m.case), rep(-1, length = m.ctrl))) # case/ctrl data with phenotype
      colnames(dataset)[ncol(dataset)] <- "class" # replace last colname with "class"
    } else if (data.type == "discrete") {
      if (avg.maf <= 0 || avg.maf >= 0.5) {
        stop("Please choose avg.maf in (0,0.5)")
      }

      m.case <- round((1 - pct.imbalance) * num.samples) # size of case group
      m.ctrl <- round(pct.imbalance * num.samples) # size of ctrl group

      # approximate effect size for pairwise correlations between functional attributes in case group
      case.hi.cor <- -hi.cor * interaction.bias + (1 - interaction.bias) * hi.cor

      e <- 1 # fudge factor to the number of nodes to avoid giant component
      prob <- 1 / (num.variables + e) # probability of a node being connected to another node is less than 1/N to avoid giant component

      if (is.null(out.degree)) {
        g <- igraph::barabasi.game(num.variables, directed = FALSE) # scale-free network
      } else {
        g <- igraph::barabasi.game(num.variables, m = out.degree, directed = FALSE)
      }

      if (!is.null(graph.structure)) {
        g <- graph.structure
      }

      # generate correlation matrix from g
      network.atts <- generate_structured_corrmat(
        g = g,
        num.variables = num.variables,
        hi.cor.tmp = case.hi.cor,
        lo.cor.tmp = lo.cor,
        hi.cor.fixed = hi.cor,
        graph.type = "Scalefree",
        plot.graph = plot.graph,
        make.diff.cors = TRUE,
        nbias = nbias, use.Rcpp = use.Rcpp
      )

      R <- as.matrix(network.atts$corrmat) # correlation matrix
      R.case <- R

      null.dats <- matrix(rnorm(num.variables * m.case), nrow = m.case, ncol = num.variables)

      U <- t(chol(R.case)) # upper tri cholesky
      null.dats <- t(U %*% t(null.dats)) # correlated data

      sig.vars <- network.atts$sig.vars # column indices of functional features

      A.mat <- network.atts$A.mat # adjacency from graph object

      eps <- 0.1
      lb <- eps
      ub <- 2 * avg.maf - eps
      f.a <- runif(num.variables, min = lb, max = ub) # minor allele frequencies (success probabilities for bernoulli trials)

      trans.dats <- pnorm(null.dats)

      bin.data <- matrix(0, nrow = nrow(trans.dats), ncol = ncol(trans.dats))
      for (j in 1:ncol(trans.dats)) {
        tmp <- qbinom(trans.dats[, j], size = 2, prob = f.a[j])

        bin.data[, j] <- tmp
      }
      X.case <- bin.data

      # for each functional variable, find connected variables
      sig.connected.list <- list() # list for storing connected variable indices
      for (i in 1:length(sig.vars)) { # for each functional feature
        sig.connected.list[[i]] <- which(abs(A.mat[sig.vars[i], ] - 1) < 1e-9) # find and store connected features
      }

      # create control group correlation matrix
      for (i in 1:length(sig.vars)) { # for each functional feature
        n.connected <- length(sig.connected.list[[i]]) # how many connections
        for (j in 1:n.connected) { # for each connection
          R[sig.vars[i], sig.connected.list[[i]][j]] <- min(hi.cor + rnorm(1, 0, .1), 1) # replace low correlation by high for ctrls
        }
      }

      # ensure that R is symmetric
      #
      tmp <- diag(R) # diag of R
      diag(R) <- 0 # replace diag of R with 0's
      R[lower.tri(R)] <- 0 # make lower triangle zeros

      R <- R + t(R) # make symmetric
      diag(R) <- tmp # original diag of R

      # correct for negative eigenvalues so R is positive definite
      #
      if (use.Rcpp) { # compute eigenvalues and make diag matrix
        R.d <- diag(sort(c(getEigenValues(R)), decreasing = TRUE))
      } else {
        R.d <- diag(eigen(R)$values)
      }
      tmp <- diag(R.d) # vector of eigenvalues

      if (any(tmp < 0)) { # if any eigenvalues are negative
        R.V <- eigen(R)$vectors # compute eigenvectors,
        tmp[tmp < 0] <- 1e-7 # make negative into small positive,
        diag(R.d) <- tmp # replace diag of R.d with positive eigenvalues,
        R.fix <- R.V %*% R.d %*% t(R.V) # compute new R, and
        R <- R.fix # store in R
      }
      R <- as.matrix(R)

      # make 1's on diag of R
      inv.diag <- 1 / diag(R) # multiplicative inverse of diag(R)
      mydiag <- diag(length(inv.diag)) # initialize diagonal matrix for inv.diag
      diag(mydiag) <- inv.diag # swap 1's for inv.diag
      mydiag <- sqrt(mydiag) # compute sqrt of inv.diag
      R <- mydiag %*% R %*% mydiag # compute corrected R with 1's on diagonal (Still Pos. Def.)
      R.ctrl <- R

      null.dats <- matrix(rnorm(num.variables * m.ctrl), nrow = m.ctrl, ncol = num.variables)

      U <- t(chol(R.ctrl)) # upper tri cholesky
      null.dats <- t(U %*% t(null.dats)) # correlated data

      trans.dats <- pnorm(null.dats)

      bin.data <- matrix(0, nrow = nrow(trans.dats), ncol = ncol(trans.dats))
      for (j in 1:ncol(trans.dats)) {
        tmp <- qbinom(trans.dats[, j], size = 2, prob = f.a[j])

        bin.data[, j] <- tmp
      }
      X.ctrl <- bin.data

      X.all <- rbind(X.case, X.ctrl) # case/ctrl data
      X.all <- as.data.frame(X.all) # make data.frame

      dataset <- X.all

      label <- "class"
      signal.names <- paste("intvar", 1:nbias, sep = "") # interaction variable names
      background.names <- paste("var", 1:(num.variables - nbias), sep = "") # non-functional variable names
      var.names <- c(signal.names, background.names, label) # all variable names
      colnames(dataset)[sig.vars] <- signal.names # replace functional colnames with interaction variable names
      colnames(dataset)[-sig.vars] <- background.names # replace non-functional colnames with non-functional colnames

      dataset <- cbind(dataset, c(rep(1, length = m.case), rep(-1, length = m.ctrl))) # full data set with phenotype
      colnames(dataset)[ncol(dataset)] <- "class" # replace last colname with "class"
    }
  } else if (sim.type == "mixed") { # main + interaction effect simulation

    if (mix.type == "main-interactionErdos") { # main + interaction with Erdos-Renyi network

      if (data.type == "continuous") {
        m.case <- round((1 - pct.imbalance) * num.samples) # size of case group
        m.ctrl <- round(pct.imbalance * num.samples) # size of ctrl group

        num.main <- round(pct.mixed * nbias) # number of main effect attributes
        num.int <- round((1 - pct.mixed) * nbias) # number of interaction effect attributes

        # make cases data set
        X.case <- matrix(rnorm(m.case * (num.variables - num.main)), nrow = m.case, ncol = (num.variables - num.main))

        # approximate effect size for pairwise correlations between functional attributes in case group
        case.hi.cor <- -hi.cor * interaction.bias + (1 - interaction.bias) * hi.cor

        e <- 1 # fudge factor to the number of nodes to avoid giant component
        if (is.null(prob.connected)) {
          prob <- 1 / ((num.variables - num.main) + e) # probability of a node being connected to another node is less than 1/N to avoid giant component
        } else {
          prob <- prob.connected
        }

        # generate random Erdos-Renyi network
        if (is.null(graph.structure)) {
          g <- igraph::erdos.renyi.game((num.variables - num.main), prob)
        } else {
          g <- graph.structure
        }

        # generate correlation matrix from g
        network.atts <- generate_structured_corrmat(
          g = g,
          num.variables = (num.variables - num.main),
          hi.cor.tmp = case.hi.cor,
          lo.cor.tmp = lo.cor,
          hi.cor.fixed = hi.cor,
          graph.type = "Erdos-Renyi",
          plot.graph = plot.graph,
          make.diff.cors = TRUE,
          nbias = num.int, use.Rcpp = use.Rcpp
        )

        R <- as.matrix(network.atts$corrmat) # case correlation matrix

        sig.vars <- network.atts$sig.vars # functional attribute indices

        A.mat <- network.atts$A.mat # adjacency from graph

        U <- t(chol(R)) # upper tri cholesky
        tmp <- t(U %*% t(X.case)) # correlated case data
        X.case <- tmp # store in X.case

        # make controls data set
        X.ctrl <- matrix(rnorm(m.ctrl * (num.variables - num.main)), nrow = m.ctrl, ncol = (num.variables - num.main))

        # find connections for functional attributes
        #
        sig.connected.list <- list() # list for functional attribute connections
        for (i in 1:length(sig.vars)) { # for each functional feature
          sig.connected.list[[i]] <- which(abs(A.mat[sig.vars[i], ] - 1) < 1e-9) # find all connections
        }

        # make high correlations for ctrls
        #
        for (i in 1:length(sig.vars)) { # for each functional variable
          n.connected <- length(sig.connected.list[[i]]) # how many connections
          for (j in 1:n.connected) { # for each connection
            R[sig.vars[i], sig.connected.list[[i]][j]] <- min(hi.cor + rnorm(1, 0, .1), 1) # replace with high correlation
          }
        }

        # ensure that R is symmetric
        #
        tmp <- diag(R) # diag of R
        diag(R) <- 0 # make diag(R) 0
        R[lower.tri(R)] <- 0 # make lower triangle of R 0

        R <- R + t(R) # make symmetric
        diag(R) <- tmp # swap for original diag

        # correct for negative eigenvalues so R is positive definite
        #
        if (use.Rcpp) { # compute eigenvalues and make diag matrix
          R.d <- diag(sort(c(getEigenValues(R)), decreasing = TRUE))
        } else {
          R.d <- diag(eigen(R)$values)
        }
        tmp <- diag(R.d) # vector of eigenvalues

        if (any(tmp < 0)) { # if any eigenvalues are negative
          R.V <- eigen(R)$vectors # compute eigenvectors,
          tmp[tmp < 0] <- 1e-7 # make negative into small positive,
          diag(R.d) <- tmp # swap diag of R.d with positive eigenvalues,
          R.fix <- R.V %*% R.d %*% t(R.V) # compute new correlation matrix, and
          R <- R.fix # store in R
        }
        R <- as.matrix(R)

        # make 1's on diagonal
        #
        inv.diag <- 1 / diag(R) # multiplicative inverse of diag(R)
        mydiag <- diag(length(inv.diag)) # initialize diag matrix for inv.diag
        diag(mydiag) <- inv.diag # swap 1's for inv.diag
        mydiag <- sqrt(mydiag) # compute sqrt of inv.diag
        R <- mydiag %*% R %*% mydiag # compute corrected R with 1's on diagonal (Still Pos. Def.)

        U <- t(chol(R)) # upper tri cholesky
        tmp <- t(U %*% t(X.ctrl)) # correlated ctrl data
        X.ctrl <- tmp # store in X.ctrl

        X.all <- rbind(X.case, X.ctrl) # case/ctrl data
        X.all <- as.data.frame(X.all) # make data.frame

        dataset <- X.all

        # create main-effect simulation for num.main attributes
        my.sim.data <- createMainEffects(
          n.e = num.main,
          n.db = num.samples,
          pct.imbalance = pct.imbalance,
          label = label,
          sd.b = main.bias,
          p.b = 1
        )$db

        # check dimensions to see if my.sim.data is a matrix or just a vector
        check.dim <- is.null(dim(my.sim.data$datnobatch))
        if (check.dim) {
          main.cols <- my.sim.data$datnobatch
        } else {
          main.cols <- t(my.sim.data$datnobatch)
        }
        dataset.tmp <- cbind(main.cols, my.sim.data$S)

        dataset <- cbind(dataset, dataset.tmp[, -ncol(dataset.tmp)])
        main.vars <- (dim(dataset)[2] - num.main + 1):(dim(dataset)[2])

        main.names <- paste("mainvar", 1:num.main, sep = "") # main effect variable names
        int.names <- paste("intvar", 1:num.int, sep = "") # interaction effect variable names

        signal.names <- c(main.names, int.names) # all functional variable names
        background.names <- paste("var", 1:(num.variables - nbias), sep = "") # all non-functional variable names
        var.names <- c(signal.names, background.names, label) # all variable names

        colnames(dataset)[sig.vars] <- int.names # replace interaction colnames with interaction variable names
        colnames(dataset)[main.vars] <- main.names # replace main colnames with main variable names
        colnames(dataset)[-c(sig.vars, main.vars)] <- background.names # replace non-functional colnames with non-functional variable names

        dataset <- cbind(dataset, c(rep(1, length = m.case), rep(-1, length = m.ctrl))) # full data set with phenotype
        colnames(dataset)[ncol(dataset)] <- "class" # replace last colname with "class"
      } else if (data.type == "discrete") {
        if (main.bias >= 1) {
          stop("Please choose main.bias in (0,1)")
        }

        if (avg.maf <= 0 || avg.maf >= 0.5) {
          stop("Please choose avg.maf in (0,0.5)")
        }

        m.case <- round((1 - pct.imbalance) * num.samples) # size of case group
        m.ctrl <- round(pct.imbalance * num.samples) # size of ctrl group

        num.main <- round(pct.mixed * nbias) # number of main effect attributes
        num.int <- round((1 - pct.mixed) * nbias) # number of interaction effect attributes

        # make cases data set
        X.case <- matrix(rnorm(m.case * (num.variables - num.main)), nrow = m.case, ncol = (num.variables - num.main))

        # approximate effect size for pairwise correlations between functional attributes in case group
        case.hi.cor <- -hi.cor * interaction.bias + (1 - interaction.bias) * hi.cor

        e <- 1 # fudge factor to the number of nodes to avoid giant component
        if (is.null(prob.connected)) {
          prob <- 1 / ((num.variables - num.main) + e) # probability of a node being connected to another node is less than 1/N to avoid giant component
        } else {
          prob <- prob.connected
        }

        # generate random Erdos-Renyi network
        g <- igraph::erdos.renyi.game((num.variables - num.main), prob)

        # generate correlation matrix from g
        network.atts <- generate_structured_corrmat(
          g = g,
          num.variables = (num.variables - num.main),
          hi.cor.tmp = case.hi.cor,
          lo.cor.tmp = lo.cor,
          hi.cor.fixed = hi.cor,
          graph.type = "Erdos-Renyi",
          plot.graph = plot.graph,
          make.diff.cors = TRUE,
          nbias = num.int, use.Rcpp = use.Rcpp
        )

        R <- as.matrix(network.atts$corrmat) # case correlation matrix

        sig.vars <- network.atts$sig.vars # functional attribute indices

        A.mat <- network.atts$A.mat # adjacency from graph

        U <- t(chol(R)) # upper tri cholesky

        null.dats <- matrix(rnorm((num.variables - num.main) * m.case), nrow = m.case, ncol = (num.variables - num.main))

        null.dats <- t(U %*% t(null.dats)) # correlated data

        A.mat <- network.atts$A.mat # adjacency from graph object

        eps <- 0.1
        lb <- eps
        ub <- 2 * avg.maf - eps
        f.a <- runif((num.variables - num.main), min = lb, max = ub) # minor allele frequencies (success probabilities for bernoulli trials)

        trans.dats <- pnorm(null.dats)

        bin.data <- matrix(0, nrow = nrow(trans.dats), ncol = ncol(trans.dats))
        for (j in 1:ncol(trans.dats)) {
          tmp <- qbinom(trans.dats[, j], size = 2, prob = f.a[j])

          bin.data[, j] <- tmp
        }
        X.case <- bin.data

        # for each functional variable, find connected variables
        sig.connected.list <- list() # list for storing connected variable indices
        for (i in 1:length(sig.vars)) { # for each functional feature
          sig.connected.list[[i]] <- which(abs(A.mat[sig.vars[i], ] - 1) < 1e-9) # find and store connected features
        }

        # create control group correlation matrix
        for (i in 1:length(sig.vars)) { # for each functional feature
          n.connected <- length(sig.connected.list[[i]]) # how many connections
          for (j in 1:n.connected) { # for each connection
            R[sig.vars[i], sig.connected.list[[i]][j]] <- min(hi.cor + rnorm(1, 0, .1), 1) # replace low correlation by high for ctrls
          }
        }

        # ensure that R is symmetric
        #
        tmp <- diag(R) # diag of R
        diag(R) <- 0 # replace diag of R with 0's
        R[lower.tri(R)] <- 0 # make lower triangle zeros

        R <- R + t(R) # make symmetric
        diag(R) <- tmp # original diag of R

        # correct for negative eigenvalues so R is positive definite
        #
        if (use.Rcpp) { # compute eigenvalues and make diag matrix
          R.d <- diag(sort(c(getEigenValues(R)), decreasing = TRUE))
        } else {
          R.d <- diag(eigen(R)$values)
        }
        tmp <- diag(R.d) # vector of eigenvalues

        if (any(tmp < 0)) { # if any eigenvalues are negative
          R.V <- eigen(R)$vectors # compute eigenvectors,
          tmp[tmp < 0] <- 1e-7 # make negative into small positive,
          diag(R.d) <- tmp # replace diag of R.d with positive eigenvalues,
          R.fix <- R.V %*% R.d %*% t(R.V) # compute new R, and
          R <- R.fix # store in R
        }
        R <- as.matrix(R)

        # make 1's on diag of R
        inv.diag <- 1 / diag(R) # multiplicative inverse of diag(R)
        mydiag <- diag(length(inv.diag)) # initialize diagonal matrix for inv.diag
        diag(mydiag) <- inv.diag # swap 1's for inv.diag
        mydiag <- sqrt(mydiag) # compute sqrt of inv.diag
        R <- mydiag %*% R %*% mydiag # compute corrected R with 1's on diagonal (Still Pos. Def.)
        R.ctrl <- R

        null.dats <- matrix(rnorm((num.variables - num.main) * m.ctrl), nrow = m.ctrl, ncol = (num.variables - num.main))

        U <- t(chol(R.ctrl)) # upper tri cholesky
        null.dats <- t(U %*% t(null.dats)) # correlated data

        trans.dats <- pnorm(null.dats)

        bin.data <- matrix(0, nrow = nrow(trans.dats), ncol = ncol(trans.dats))
        for (j in 1:ncol(trans.dats)) {
          tmp <- qbinom(trans.dats[, j], size = 2, prob = f.a[j])

          bin.data[, j] <- tmp
        }
        X.ctrl <- bin.data

        dataset <- rbind(X.case, X.ctrl)

        lb <- avg.maf - 0.5 * main.bias # potentially the average maf for control group
        ub <- avg.maf + 0.5 * main.bias # potentially the average maf for the case group

        if (lb <= 0) { # check to see if lower maf is possible, if not just make it small positive number

          eps <- 0.01

          lb <- eps

          ub <- 2 * avg.maf - eps # chosen so that 0.5*(lb + ub) = avg.maf
        } else if (ub >= 1) { # check to see if upper maf is possible, if not just make it a number slightly less than 1

          eps <- 0.01

          ub <- 1 - eps

          lb <- 2 * avg.maf - (1 - eps) # chosen so that 0.5*(lb + ub) = avg.maf
        }

        prob.case <- ub # case group mean
        prob.ctrl <- lb # ctrl group mean

        # main effects for cases
        main.case <- matrix(0, nrow = m.case, ncol = num.main)
        for (j in 1:m.case) {
          tmp <- rbinom(n = num.main, size = 2, prob = prob.case)

          main.case[j, ] <- tmp
        }

        # main effects for ctrls
        main.ctrl <- matrix(0, nrow = m.ctrl, ncol = num.main)
        for (j in 1:m.ctrl) {
          tmp <- rbinom(n = num.main, size = 2, prob = prob.ctrl)

          main.ctrl[j, ] <- tmp
        }

        # main effect features
        main.effects <- rbind(main.case, main.ctrl)

        # interactions + main effects
        dataset <- cbind(dataset, main.effects)
        colnames(dataset) <- paste("var", 1:num.variables, sep = "")

        main.vars <- (dim(dataset)[2] - num.main + 1):(dim(dataset)[2])

        main.names <- paste("mainvar", 1:num.main, sep = "") # main effect variable names
        int.names <- paste("intvar", 1:num.int, sep = "") # interaction effect variable names

        signal.names <- c(main.names, int.names) # all functional variable names
        background.names <- paste("var", 1:(num.variables - nbias), sep = "") # all non-functional variable names
        var.names <- c(signal.names, background.names, label) # all variable names

        colnames(dataset)[sig.vars] <- int.names # replace interaction colnames with interaction variable names
        colnames(dataset)[main.vars] <- main.names # replace main colnames with main variable names
        colnames(dataset)[-c(sig.vars, main.vars)] <- background.names # replace non-functional colnames with non-functional variable names

        dataset <- cbind(dataset, c(rep(1, length = m.case), rep(-1, length = m.ctrl))) # full data set with phenotype
        colnames(dataset)[ncol(dataset)] <- "class" # replace last colname with "class"
        dataset <- as.data.frame(dataset)
      }
    } else if (mix.type == "main-interactionScalefree") { # main + interaction with Scale-Free network

      if (data.type == "continuous") {
        m.case <- round((1 - pct.imbalance) * num.samples) # number of cases
        m.ctrl <- round(pct.imbalance * num.samples) # number of ctrls

        num.main <- round(pct.mixed * nbias) # number of main effect attributes
        num.int <- round((1 - pct.mixed) * nbias) # number of interaction effect attributes

        # make cases data set
        X.case <- matrix(rnorm(m.case * (num.variables - num.main)), nrow = m.case, ncol = (num.variables - num.main))

        # approximate effect size for pairwise correlations between functional attributes in case group
        case.hi.cor <- -hi.cor * interaction.bias + (1 - interaction.bias) * hi.cor

        e <- 1 # fudge factor to the number of nodes to avoid giant component
        prob <- 1 / ((num.variables - num.main) + e) # probability of a node being connected to another node is less than 1/N to avoid giant component

        # generate random Scale-Free network
        if (is.null(out.degree)) {
          g <- igraph::barabasi.game((num.variables - num.main), directed = FALSE)
        } else {
          g <- igraph::barabasi.game((num.variables - num.main), m = out.degree, directed = FALSE)
        }

        # generate case group correlation matrix
        network.atts <- generate_structured_corrmat(
          g = g,
          num.variables = (num.variables - num.main),
          hi.cor.tmp = case.hi.cor,
          lo.cor.tmp = lo.cor,
          hi.cor.fixed = hi.cor,
          graph.type = "Scalefree",
          plot.graph = plot.graph,
          make.diff.cors = TRUE,
          nbias = num.int, use.Rcpp = use.Rcpp
        )

        R <- as.matrix(network.atts$corrmat) # cases correlation matrix

        sig.vars <- network.atts$sig.vars # functional attribute indices

        A.mat <- network.atts$A.mat # adjacency from graph

        U <- t(chol(R)) # upper tri cholesky
        tmp <- t(U %*% t(X.case)) # correlated case data
        X.case <- tmp # store in X.case

        # make controls data set
        X.ctrl <- matrix(rnorm(m.ctrl * (num.variables - num.main)), nrow = m.ctrl, ncol = (num.variables - num.main))

        # find all functional connections
        #
        sig.connected.list <- list() # list for storing functional connections
        for (i in 1:length(sig.vars)) { # for each functional variable
          sig.connected.list[[i]] <- which(abs(A.mat[sig.vars[i], ] - 1) < 1e-9) # find all connections
        }

        # make high correlations matrix for ctrls
        #
        for (i in 1:length(sig.vars)) { # for each functional variable
          n.connected <- length(sig.connected.list[[i]]) # how many connections
          for (j in 1:n.connected) { # for each connection
            R[sig.vars[i], sig.connected.list[[i]][j]] <- min(hi.cor + rnorm(1, 0, .1), 1) # replace low correlation with high correlation
          }
        }

        # ensure that R is symmetric
        #
        tmp <- diag(R) # diag of R
        diag(R) <- 0 # make diag(R) 0
        R[lower.tri(R)] <- 0 # make lower triangle of R 0

        R <- R + t(R) # make symmetric
        diag(R) <- tmp # swap for original diag

        # correct for negative eigenvalues to make R positive definite
        #
        if (use.Rcpp) { # compute eigenvalues and make diag matrix
          R.d <- diag(sort(c(getEigenValues(R)), decreasing = TRUE))
        } else {
          R.d <- diag(eigen(R)$values)
        }
        tmp <- diag(R.d) # vector of eigenvalues

        if (any(tmp < 0)) { # if any eigenvalues are negative
          R.V <- eigen(R)$vectors # compute eigenvectores,
          tmp[tmp < 0] <- 1e-7 # make negative into small positive,
          diag(R.d) <- tmp # replace diag R.d with positive eigenvalues
          R.fix <- R.V %*% R.d %*% t(R.V) # compute new correlation matrix, and
          R <- R.fix # store in R
        }
        R <- as.matrix(R)

        # make 1's on diagonal
        #
        inv.diag <- 1 / diag(R) # multiplicative inverse of diag(R)
        mydiag <- diag(length(inv.diag)) # initialize diag matrix for inv.diag
        diag(mydiag) <- inv.diag # swap 1's for inv.diag
        mydiag <- sqrt(mydiag) # sqrt of inv.diag
        R <- mydiag %*% R %*% mydiag # compute corrected R with 1's on diagonal (Still Pos. Def.)

        U <- t(chol(R)) # upper tri cholesky
        tmp <- t(U %*% t(X.ctrl)) # correlated ctrl data
        X.ctrl <- tmp # store in X.ctrl

        X.all <- rbind(X.case, X.ctrl) # case/ctrl data
        X.all <- as.data.frame(X.all) # make data.frame

        dataset <- X.all

        # create main effect simulation with num.main features
        my.sim.data <- createMainEffects(
          n.e = num.main,
          n.db = num.samples,
          pct.imbalance = pct.imbalance,
          label = label,
          sd.b = main.bias,
          p.b = 1
        )$db

        # check dimensions to determine if my.dim.data is a matrix or vector
        check.dim <- is.null(dim(my.sim.data$datnobatch))
        if (check.dim) {
          main.cols <- my.sim.data$datnobatch
        } else {
          main.cols <- t(my.sim.data$datnobatch)
        }
        dataset.tmp <- cbind(main.cols, my.sim.data$S)

        dataset <- cbind(dataset, dataset.tmp[, -ncol(dataset.tmp)])
        main.vars <- (dim(dataset)[2] - num.main + 1):(dim(dataset)[2])

        main.names <- paste("mainvar", 1:num.main, sep = "") # main effect variable names
        int.names <- paste("intvar", 1:num.int, sep = "") # interaction effect variable names

        signal.names <- c(main.names, int.names) # all functional variable names
        background.names <- paste("var", 1:(num.variables - nbias), sep = "") # non-functional variable names
        var.names <- c(signal.names, background.names, label) # all variable names

        colnames(dataset)[sig.vars] <- int.names # replace interaction colnames with interaction variable names
        colnames(dataset)[main.vars] <- main.names # replace main colnames with main variable names
        colnames(dataset)[-c(sig.vars, main.vars)] <- background.names # replace non-functional colnames with non-functional variable names

        dataset <- cbind(dataset, c(rep(1, length = m.case), rep(-1, length = m.ctrl))) # full data set with phenotype
        colnames(dataset)[ncol(dataset)] <- "class" # replace last colname with "class"
      } else if (data.type == "discrete") {
        if (main.bias >= 1) {
          stop("Please choose main.bias in (0,1)")
        }

        if (avg.maf <= 0 || avg.maf >= 0.5) {
          stop("Please choose avg.maf in (0,0.5)")
        }

        m.case <- round((1 - pct.imbalance) * num.samples) # size of case group
        m.ctrl <- round(pct.imbalance * num.samples) # size of ctrl group

        num.main <- round(pct.mixed * nbias) # number of main effect attributes
        num.int <- round((1 - pct.mixed) * nbias) # number of interaction effect attributes

        # make cases data set
        X.case <- matrix(rnorm(m.case * (num.variables - num.main)), nrow = m.case, ncol = (num.variables - num.main))

        # approximate effect size for pairwise correlations between functional attributes in case group
        case.hi.cor <- -hi.cor * interaction.bias + (1 - interaction.bias) * hi.cor

        e <- 1 # fudge factor to the number of nodes to avoid giant component
        prob <- 1 / ((num.variables - num.main) + e) # probability of a node being connected to another node is less than 1/N to avoid giant component

        # generate random Scale-Free network
        if (is.null(out.degree)) {
          g <- igraph::barabasi.game((num.variables - num.main), directed = FALSE)
        } else {
          g <- igraph::barabasi.game((num.variables - num.main), m = out.degree, directed = FALSE)
        }

        # generate correlation matrix from g
        network.atts <- generate_structured_corrmat(
          g = g,
          num.variables = (num.variables - num.main),
          hi.cor.tmp = case.hi.cor,
          lo.cor.tmp = lo.cor,
          hi.cor.fixed = hi.cor,
          graph.type = "Scalefree",
          plot.graph = plot.graph,
          make.diff.cors = TRUE,
          nbias = num.int, use.Rcpp = use.Rcpp
        )

        R <- as.matrix(network.atts$corrmat) # case correlation matrix

        sig.vars <- network.atts$sig.vars # functional attribute indices

        A.mat <- network.atts$A.mat # adjacency from graph

        U <- t(chol(R)) # upper tri cholesky

        null.dats <- matrix(rnorm((num.variables - num.main) * m.case), nrow = m.case, ncol = (num.variables - num.main))

        null.dats <- t(U %*% t(null.dats)) # correlated data

        A.mat <- network.atts$A.mat # adjacency from graph object

        eps <- 0.1
        lb <- eps
        ub <- 2 * avg.maf - eps
        f.a <- runif((num.variables - num.main), min = lb, max = ub) # minor allele frequencies (success probabilities for bernoulli trials)

        trans.dats <- pnorm(null.dats)

        bin.data <- matrix(0, nrow = nrow(trans.dats), ncol = ncol(trans.dats))
        for (j in 1:ncol(trans.dats)) {
          tmp <- qbinom(trans.dats[, j], size = 2, prob = f.a[j])

          bin.data[, j] <- tmp
        }
        X.case <- bin.data

        # for each functional variable, find connected variables
        sig.connected.list <- list() # list for storing connected variable indices
        for (i in 1:length(sig.vars)) { # for each functional feature
          sig.connected.list[[i]] <- which(abs(A.mat[sig.vars[i], ] - 1) < 1e-9) # find and store connected features
        }

        # create control group correlation matrix
        for (i in 1:length(sig.vars)) { # for each functional feature
          n.connected <- length(sig.connected.list[[i]]) # how many connections
          for (j in 1:n.connected) { # for each connection
            R[sig.vars[i], sig.connected.list[[i]][j]] <- min(hi.cor + rnorm(1, 0, .1), 1) # replace low correlation by high for ctrls
          }
        }

        # ensure that R is symmetric
        #
        tmp <- diag(R) # diag of R
        diag(R) <- 0 # replace diag of R with 0's
        R[lower.tri(R)] <- 0 # make lower triangle zeros

        R <- R + t(R) # make symmetric
        diag(R) <- tmp # original diag of R

        # correct for negative eigenvalues so R is positive definite
        #
        if (use.Rcpp) { # compute eigenvalues and make diag matrix
          R.d <- diag(sort(c(getEigenValues(R)), decreasing = TRUE))
        } else {
          R.d <- diag(eigen(R)$values)
        }
        tmp <- diag(R.d) # vector of eigenvalues

        if (any(tmp < 0)) { # if any eigenvalues are negative
          R.V <- eigen(R)$vectors # compute eigenvectors,
          tmp[tmp < 0] <- 1e-7 # make negative into small positive,
          diag(R.d) <- tmp # replace diag of R.d with positive eigenvalues,
          R.fix <- R.V %*% R.d %*% t(R.V) # compute new R, and
          R <- R.fix # store in R
        }
        R <- as.matrix(R)

        # make 1's on diag of R
        inv.diag <- 1 / diag(R) # multiplicative inverse of diag(R)
        mydiag <- diag(length(inv.diag)) # initialize diagonal matrix for inv.diag
        diag(mydiag) <- inv.diag # swap 1's for inv.diag
        mydiag <- sqrt(mydiag) # compute sqrt of inv.diag
        R <- mydiag %*% R %*% mydiag # compute corrected R with 1's on diagonal (Still Pos. Def.)
        R.ctrl <- R

        null.dats <- matrix(rnorm((num.variables - num.main) * m.ctrl), nrow = m.ctrl, ncol = (num.variables - num.main))

        U <- t(chol(R.ctrl)) # upper tri cholesky
        null.dats <- t(U %*% t(null.dats)) # correlated data

        trans.dats <- pnorm(null.dats)

        bin.data <- matrix(0, nrow = nrow(trans.dats), ncol = ncol(trans.dats))
        for (j in 1:ncol(trans.dats)) {
          tmp <- qbinom(trans.dats[, j], size = 2, prob = f.a[j])

          bin.data[, j] <- tmp
        }
        X.ctrl <- bin.data

        dataset <- rbind(X.case, X.ctrl)

        lb <- avg.maf - 0.5 * main.bias # potentially the average maf for control group
        ub <- avg.maf + 0.5 * main.bias # potentially the average maf for the case group

        if (lb <= 0) { # check to see if lower maf is possible, if not just make it small positive number

          eps <- 0.01

          lb <- eps

          ub <- 2 * avg.maf - eps # chosen so that 0.5*(lb + ub) = avg.maf
        } else if (ub >= 1) { # check to see if upper maf is possible, if not just make it a number slightly less than 1

          eps <- 0.01

          ub <- 1 - eps

          lb <- 2 * avg.maf - (1 - eps) # chosen so that 0.5*(lb + ub) = avg.maf
        }

        prob.case <- ub # case group mean
        prob.ctrl <- lb # ctrl group mean

        # main effects for cases
        main.case <- matrix(0, nrow = m.case, ncol = num.main)
        for (j in 1:m.case) {
          tmp <- rbinom(n = num.main, size = 2, prob = prob.case)

          main.case[j, ] <- tmp
        }

        # main effects for ctrls
        main.ctrl <- matrix(0, nrow = m.ctrl, ncol = num.main)
        for (j in 1:m.ctrl) {
          tmp <- rbinom(n = num.main, size = 2, prob = prob.ctrl)

          main.ctrl[j, ] <- tmp
        }

        # main effect features
        main.effects <- rbind(main.case, main.ctrl)

        # interactions + main effects
        dataset <- cbind(dataset, main.effects)
        colnames(dataset) <- paste("var", 1:num.variables, sep = "")

        main.vars <- (dim(dataset)[2] - num.main + 1):(dim(dataset)[2])

        main.names <- paste("mainvar", 1:num.main, sep = "") # main effect variable names
        int.names <- paste("intvar", 1:num.int, sep = "") # interaction effect variable names

        signal.names <- c(main.names, int.names) # all functional variable names
        background.names <- paste("var", 1:(num.variables - nbias), sep = "") # all non-functional variable names
        var.names <- c(signal.names, background.names, label) # all variable names

        colnames(dataset)[sig.vars] <- int.names # replace interaction colnames with interaction variable names
        colnames(dataset)[main.vars] <- main.names # replace main colnames with main variable names
        colnames(dataset)[-c(sig.vars, main.vars)] <- background.names # replace non-functional colnames with non-functional variable names

        dataset <- cbind(dataset, c(rep(1, length = m.case), rep(-1, length = m.ctrl))) # full data set with phenotype
        colnames(dataset)[ncol(dataset)] <- "class" # replace last colname with "class"
        dataset <- as.data.frame(dataset)
      }
    }
  }

  # split data into train, holdout, and validation sets
  split.data <- splitDataset(
    all.data = dataset,
    pct.train = pct.train,
    pct.holdout = pct.holdout,
    pct.validation = pct.validation,
    label = label
  )

  if (!is.null(save.file)) {
    if (verbose) cat("saving to data/", save.file, ".Rdata\n", sep = "")
    save(split.data, pct.signals, main.bias, interaction.bias, sim.type, file = save.file)
  }

  elapsed <- (proc.time() - ptm)[3]
  if (verbose) cat("createSimulation elapsed time:", elapsed, "\n")

  # depending on sim.type, we have different potentially important outputs
  if (sim.type %in% c("interactionErdos", "interactionScalefree")) {

    # for interaction simulation only
    list(
      train = split.data$train,
      holdout = split.data$holdout,
      validation = split.data$validation,
      label = label,
      signal.names = signal.names,
      elapsed = elapsed,
      sig.vars = sig.vars,
      A.mat = A.mat
    )
  } else if (sim.type == "mixed") {

    # for mixed simulation (main+interaction)
    list(
      train = split.data$train,
      holdout = split.data$holdout,
      validation = split.data$validation,
      label = label,
      signal.names = signal.names,
      elapsed = elapsed,
      int.vars = sig.vars,
      main.vars = main.vars,
      A.mat = A.mat
    )
  } else if (sim.type == "mainEffect_Erdos-Renyi") {

    # for main effect + correlation simulation (Erdos-Renyi Network)
    list(
      train = split.data$train,
      holdout = split.data$holdout,
      validation = split.data$validation,
      label = label,
      signal.names = signal.names,
      elapsed = elapsed,
      A.mat = A.mat
    )
  } else if (sim.type == "mainEffect_Scalefree") {

    # for main effect + correlation simulation (Scale-Free Network)
    list(
      train = split.data$train,
      holdout = split.data$holdout,
      validation = split.data$validation,
      label = label,
      signal.names = signal.names,
      elapsed = elapsed,
      A.mat = A.mat
    )
  } else {

    # for simple main effect simulation
    list(
      train = split.data$train,
      holdout = split.data$holdout,
      validation = split.data$validation,
      label = label,
      signal.names = signal.names,
      elapsed = elapsed
    )
  }
}
#############################################################################

# =========================================================================================#
#' createfMRIsimulation
#'
#' function for generating NPDR-formatted rs-fMRI data set
#'
#' Parameters: Very similar to createSimulation2() in NPDR
#'
#'  New Parameters:
#'
#' num.times - (numeric) number of time points for each ROI (same for cases and controls)
#' graph.structure - (igraph like object) user-supplied graphical structure for all subjects.
#' Currently only igraph object is acceptable. We really just need the ability
#' to compute degree vec and adjacency matrix in generate_structured_corrmat().
#' sim.graph.structure - (logical) set to TRUE for random graph from igraph package. If FALSE, must provide graph structure as input.
#'
#' Value:
#'
#' (list) with the following elements:
#'
#' corr.attr.names - (character) ordered vec of ROI names from simulation. Each subject's correlation matrix
#' will have this exact order in its rows and columns.
#' dataset - (matrix) of dimension num.samples x num.variables*(num.variables - 1), where each row
#' represents a subject's pairwise correlations. The first (num.variables - 1) columns
#' are all pairwise correlations (excluding self-correlation) with the first ROI in
#' corr.attr.names vec, the next subsequent group of (num.variables - 1) columns are
#' pairwise correlations with the second ROI in corr.attr.names vec, ... etc. This order
#' is preserved in all individual subject correlation matrices.
#'
#' Details: Generates (num.times x num.variables) matrix of (num.variables) ROI time series for each of the num.samples subjects. Each ROI has (num.times) time points in the simulated fMRI scan. Functional ROIs can be created using (sim.type), which should be one of c("interactionErdos","interactionScalefree"). Effect size is controlled exactly the same as in createSimulation2(), using a combination of (interaction.bias), (hi.cor), and (prob.connected) (or out.degree). For example, setting interaction.bias = 1, hi.cor ~ 1, and prob.connected ~ 1 will yield maximal effect sizes for functional ROIs (or out.degree = num.variables - 1 for interactionScalefree). The order of ROIs in (dataset) columns is given by (corr.attr.names), which is a primary input to npdr() that allows us to assign importance to each ROI.
#'
#' Notes: Must be careful when using real data to make sure all subject correlation matrices have the same row/column order. The (corr.attr.names) parameter contains ROI names in the same order as row/column in subject correlation matrices. If subjects have differing row/column order, then ROI importance is meaningless.
#'
#' @param num.samples number of samples
#' @param num.variables number of variables (features)
#' @param pct.imbalance fraction of num.samples that are cases
#' @param pct.signals fraction of num.variables that are functional
#' @param main.bias approximate effect size for main effect simulations
#' @param interaction.bias approximate effect size for interaction effects
#' @param hi.cor parameter to use for network-connected pairwise correlations
#' @param lo.cor parameter to use for network-non-connected pairwise correlations
#' @param label should just be "class" for binary response
#' @param sim.type a character that determines the type of simulation:
#' mainEffect/mainEffect_Erdos-Renyi/mainEffect_Scalefree/interactionErdos/interactionScalefree/mixed
#' @param mix.type character that determines the type of mixed effects simulation:
#' main-interactionErdos/main-interactionScalefree
#' @param pct.mixed percent of functional variables that are main effects (1 - pct_interaction). Use with sim.type="mixed" and specify mix.type.
#' @param plot.graph logical indicating whether to plot networks
#' @param use.Rcpp if true use Rcpp to correct negative eigenvalues
#' @param prob.connected probability of drawing an edge between two arbitrary vertices in Erdos-Renyi graph
#' @param out.degree out-degree of vertices in Scale-free graph
#' @param num.times - (numeric) number of time points for each ROI (same for cases and controls)
#' @param graph.structure (igraph) graph structure generated from igraph, NULL default
#' @param data.type character indicating if data is from a "continuous" or "discrete" distribution
#' @param sim.graph.structure - (logical) set to TRUE for random graph from igraph package. If FALSE, must provide graph structure as input.
#' @return A list with:
#' \describe{
#'   \item{dataset}{(matrix) of dimension num.samples x num.variables*(num.variables - 1), where each row}
#'   \item{corr.attr.names}{(character) ordered vec of ROI names from simulation. Each subject's correlation matrix will have this exact order in its rows and columns.}
#' }
#' @export
createfMRIsimulation <- function(num.samples = 100,
                                 num.variables = 100,
                                 num.times = 100,
                                 pct.imbalance = 0.5,
                                 pct.signals = 0.1,
                                 main.bias = 0.4,
                                 interaction.bias = 0.4,
                                 hi.cor = 0.8,
                                 lo.cor = 0.2,
                                 label = "class",
                                 sim.type = "interactionErdos",
                                 mix.type = NULL,
                                 pct.mixed = 0.5,
                                 plot.graph = FALSE, use.Rcpp = FALSE,
                                 prob.connected = NULL,
                                 out.degree = NULL,
                                 sim.graph.structure = TRUE,
                                 data.type = "continuous",
                                 graph.structure = NULL) {
  if (sim.graph.structure) {
    nbias <- pct.signals * num.variables # number of functional attributes

    e <- 1 # fudge factor to the number of nodes to avoid giant component
    if (is.null(prob.connected)) {
      prob <- 1 / (num.variables + e) # probability of a node being connected to another node is less than 1/N to avoid giant component
    } else {
      prob <- prob.connected
    }

    if (is.null(graph.structure)) {
      if (sim.type != "mixed") {
        if (c(sim.type %in% c("interactionErdos", "mainEffect_Erdos-Renyi"))) {
          g <- igraph::erdos.renyi.game(num.variables, prob) # Erdos-Renyi network
        } else if (c(sim.type %in% c("interactionScalefree", "mainEffect_Scalefree"))) {
          if (is.null(out.degree)) {
            g <- igraph::barabasi.game(num.variables, directed = FALSE) # scale-free network
          } else {
            g <- igraph::barabasi.game(num.variables, m = out.degree, directed = FALSE)
          }
        }
      } else if (sim.type == "mixed") {
        num.main <- round(pct.mixed * nbias) # number of main effect attributes
        num.int <- round((1 - pct.mixed) * nbias) # number of interaction effect attributes

        if (mix.type == "main-interactionErdos") {
          e <- 1 # fudge factor to the number of nodes to avoid giant component
          if (is.null(prob.connected)) {
            prob <- 1 / ((num.variables - num.main) + e) # probability of a node being connected to another node is less than 1/N to avoid giant component
          } else {
            prob <- prob.connected
          }

          # generate random Erdos-Renyi network
          g <- igraph::erdos.renyi.game((num.variables - num.main), prob)
        } else if (mix.type == "main-interactionScalefree") {
          # generate random Scale-Free network
          if (is.null(out.degree)) {
            g <- igraph::barabasi.game((num.variables - num.main), directed = FALSE)
          } else {
            g <- igraph::barabasi.game((num.variables - num.main), m = out.degree, directed = FALSE)
          }
        }
      }
    } else {
      g <- graph.structure # user supplied graph structure
    }
  } else {
    g <- graph.structure
  }

  m.case <- round((1 - pct.imbalance) * num.samples) # size of case group
  m.ctrl <- round(pct.imbalance * num.samples) # size of ctrl group

  m <- num.samples # m.case + m.ctrl
  p <- num.variables # number of ROIs
  D.case <- matrix(0, nrow = m.case, ncol = p * (p - 1)) # case block of transformed data

  for (k in 1:m.case) {
    cat("Case Subject: ", k, "\n")

    data.mat <- createSimulation2(
      num.samples = num.times,
      num.variables = p,
      pct.imbalance = pct.imbalance,
      pct.signals = pct.signals,
      main.bias = main.bias,
      interaction.bias = interaction.bias,
      hi.cor = hi.cor,
      lo.cor = lo.cor,
      mix.type = mix.type,
      label = "class",
      sim.type = sim.type,
      pct.mixed = pct.mixed,
      pct.train = 0.5,
      pct.holdout = 0.5,
      pct.validation = 0,
      plot.graph = FALSE,
      verbose = TRUE,
      use.Rcpp = TRUE,
      prob.connected = prob.connected,
      out.degree = out.degree,
      data.type = data.type,
      graph.structure = g
    )
    dats <- rbind(data.mat$train, data.mat$holdout, data.mat$validation)
    dats <- dats[order(dats[, ncol(dats)]), ]

    case.dats <- dats[which(as.character(dats[, label]) == "1"), ]
    case.dats.tmp <- case.dats[, -ncol(case.dats)]
    case.dats.tmp <- case.dats.tmp[, sort(colnames(case.dats.tmp))]
    R <- cor(case.dats.tmp) # correlation matrix
    print(R[1:5, 1:5])

    zcorr <- apply(matrix(stretch_mat(R), ncol = 1), 1, r_to_z_fn)
    D.case[k, ] <- matrix(zcorr, ncol = (p * (p - 1)), nrow = 1)
  }

  D.ctrl <- matrix(0, nrow = m.ctrl, ncol = p * (p - 1)) # ctrl block of transformed data

  for (k in 1:m.ctrl) {
    cat("Control Subject: ", k, "\n")

    data.mat <- createSimulation2(
      num.samples = num.times,
      num.variables = p,
      pct.imbalance = pct.imbalance,
      pct.signals = pct.signals,
      main.bias = main.bias,
      interaction.bias = interaction.bias,
      hi.cor = hi.cor,
      lo.cor = lo.cor,
      mix.type = mix.type,
      label = "class",
      sim.type = sim.type,
      pct.mixed = pct.mixed,
      pct.train = 0.5,
      pct.holdout = 0.5,
      pct.validation = 0,
      plot.graph = FALSE,
      verbose = TRUE,
      use.Rcpp = TRUE,
      prob.connected = prob.connected,
      out.degree = out.degree,
      data.type = data.type
    )
    dats <- rbind(data.mat$train, data.mat$holdout, data.mat$validation)
    dats <- dats[order(dats[, ncol(dats)]), ]

    ctrl.dats <- dats[which(as.character(dats[, label]) == "-1"), ]
    ctrl.dats.tmp <- ctrl.dats[, -ncol(ctrl.dats)]
    ctrl.dats.tmp <- ctrl.dats.tmp[, sort(colnames(ctrl.dats.tmp))]
    R <- cor(ctrl.dats.tmp) # correlation matrix

    zcorr <- apply(matrix(stretch_mat(R), ncol = 1), 1, r_to_z_fn)
    D.ctrl[k, ] <- matrix(zcorr, ncol = (p * (p - 1)), nrow = 1)
  }

  D <- rbind(D.case, D.ctrl)
  case.ctrl <- c(rep(1, length = m.case), rep(-1, length = m.ctrl))
  D <- cbind(D, class = case.ctrl)
  corr.attr.names <- colnames(ctrl.dats.tmp)

  list(dataset = D, corr.attr.names = corr.attr.names)
}
##############################################################################

#Rcpp::sourceCpp('R/arma_getEigenValues.cpp')
#Rcpp::cppFunction(depends="RcppArmadillo",
#                    'arma::vec getEigenValues(arma::mat M) {
#                    return arma::eig_sym(M);
#}')


#=========================================================================================#
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
#' @return A list containing:
#' \describe{
#'   \item{cor.mat}{structured correlation matrix}
#'   \item{deg.vec}{degree vector corresponding to network adjacency}
#'   \item{A.mat}{adjacency matrix corresponding to network}
#'   \item{sig.vars}{indices of functional variables}
#' }
#' @examples 
#' 
#' @export

generate_structured_corrmat <- function(g=NULL,
                                        num.variables=100, 
                                        hi.cor.tmp=0.8, 
                                        lo.cor.tmp=0.2, 
                                        graph.type="Erdos-Renyi", 
                                        plot.graph=F,
                                        make.diff.cors=F,
                                        nbias=1, use.Rcpp=F){
  if(use.Rcpp){
    Rcpp::cppFunction(depends="RcppArmadillo",
                      'arma::vec getEigenValues(arma::mat M) {
                      return arma::eig_sym(M);
    }')
  }
  if(abs(as.integer(num.variables) - num.variables) > 1e-9){
    stop("generate_structured_corrmat: num.variables should be a positive integer")
  }
  
  p <- num.variables
  
  kvec <- degree(g)
  
  # make plot of random graph
  if(plot.graph==T){
    plot(g,vertex.size=1,vertex.label.color=kvec)
  }
  
  # which variables are not connected in network
  idx.not.connected <- which(kvec == 0)
  
  # which variables are connected in network
  if(length(idx.not.connected) > 0){                              # if any are not connected
    idx.connected <- seq(1,length(kvec),by=1)[-idx.not.connected]
  }else{                                                          # if all are connected
    idx.connected <- seq(1,length(kvec),by=1)
  }
  
  # functional variables are randomly defined from connected variables
  if(length(idx.connected) >= nbias){                             # if number connected at least nbias
    diff.cor.vars <- sample(idx.connected,size=nbias,replace=F)
  }else{                                                          # if too few are connected
    n.add.vars <- nbias - length(idx.connected)
    diff.cor.vars <- c(idx.connected, sample(idx.not.connected,size=n.add.vars,replace=F)) # add some from unconnected
  }
  
  ## make make noisy covariance matrix form structure of adjacency matrix
  Adj <- as.matrix(get.adjacency(g))
  
  # if simply generating a correlation matrix without between-group differential correlation
  if(make.diff.cors==F){
    
    tmp <- Adj[upper.tri(Adj)]                                 # upper triangle of adjacency
    tmp2 <- (hi.cor.tmp + rnorm(length(tmp),sd=.1))*tmp            # connected (high) correlations
    tmp3 <- -(tmp - 1)                                         # additive inverse of adjacency upper triangle
    tmp4 <- (lo.cor.tmp + rnorm(length(tmp3),sd=.1))*tmp3          # non-connected (low) correlations
    
    upper.mat <- matrix(0, nrow=dim(Adj)[1], ncol=dim(Adj)[1]) # initialize connected correlations matrix
    lower.mat <- matrix(0, nrow=dim(Adj)[1], ncol=dim(Adj)[1]) # initialize non-connected correlations matrix
    
    upper.mat[upper.tri(upper.mat)] <- tmp2                    # load upper triangle of connected matrix
    lower.mat[upper.tri(lower.mat)] <- tmp4                    # load upper triangle of non-connected matrix
    new.mat <- upper.mat + lower.mat                           # add connected/non-connected to get full matrix
    new.mat <- new.mat + t(new.mat)                            # make symmetric
    diag(new.mat) <- 1                                         # 1's on the diagonal
    
    # if generating correlation matrix and between-group differential correlation
  }else{
    
    # find connections of functional attributes, store connections in list,
    # and create binary matrix for functional and connections only (Subset of Adj)
    #
    adj.tmp1 <- Adj*0                                                      # initialize binary matrix
    connected.list <- list()                                               # list for functional connections
    for(i in 1:length(diff.cor.vars)){                                     # for each functional,                                    
      connected.list[[i]] <- which(abs(Adj[diff.cor.vars[i],] - 1) < 1e-9) # find connections,
      adj.tmp1[diff.cor.vars[i], connected.list[[i]]] <- 1                 # make connections 1's in binary matrix, and
      adj.tmp1[connected.list[[i]], diff.cor.vars[i]] <- 1                 # make symmetric
    }
    
    adj.tmp2 <- Adj - adj.tmp1 # binary matrix for all non-functional connections (Subset of Adj)                                    
    
    adj.tmp3 <- -(Adj - 1)     # additive inverse of Adjacency for all non-connections
    
    tmp <- adj.tmp1[upper.tri(adj.tmp1)] # functional connections only (upper triangle)
    tmp2 <- (hi.cor.tmp + rnorm(length(tmp),sd=.1))*tmp # functional connections correlations matrix  
    tmp3 <- (lo.cor.tmp + rnorm(length(adj.tmp3[upper.tri(adj.tmp3)]),sd=.1))*adj.tmp3[upper.tri(adj.tmp3)] # non-connections correlations matrix
    tmp4 <- (hi.cor.tmp + rnorm(length(adj.tmp2[upper.tri(adj.tmp2)]),sd=.1))*adj.tmp2[upper.tri(adj.tmp2)] # non-functional connections correlations matrix
    
    mat1 <- matrix(0, nrow=dim(Adj)[1], ncol=dim(Adj)[1]) # initialize functional correlation matrix
    mat2 <- matrix(0, nrow=dim(Adj)[1], ncol=dim(Adj)[1]) # initialize non-connection correlation matrix
    mat3 <- matrix(0, nrow=dim(Adj)[1], ncol=dim(Adj)[1]) # initialize non-functional connection correlation matrix
    
    mat1[upper.tri(mat1)] <- tmp2 # load upper triangle of functional matrix
    mat2[upper.tri(mat2)] <- tmp3 # load upper triangle of non-connection matrix
    mat3[upper.tri(mat3)] <- tmp4 # load upper triangle of non-functional matrix
    
    new.mat <- mat1 + mat2 + mat3   # combine all correlations in single matrix
    new.mat <- new.mat + t(new.mat) # make symmetric
    diag(new.mat) <- 1              # 1's on diagonal
    
  }
  R <- new.mat
  
  # correct for negative eigenvalues to make matrix positive definite
  #
  if(use.Rcpp){ # compute eigenvalues and make diag matrix
    R.d <- Matrix(diag(sort(c(getEigenValues(R)),decreasing=T)),sparse=T)
  }else{ 
    R.d <- Matrix(diag(eigen(R)$values),sparse=T)
  }
  
  tmp <- c(diag(R.d))                                      # vector of eigenvalues
  
  if (any(tmp<0)){                # if any eigenvalues are negative
    R.V <- Matrix(eigen(R)$vectors,sparse=T) # compute eigenvectors,
    tmp[tmp<0] <- 1e-7            # make negative into small positive,
    diag(R.d) <- tmp              # replace in R.d,
    R.fix <- R.V%*%R.d%*%t(R.V)   # compute new correlation matrix, and
    R <- R.fix                    # store in R
  }
  R <- as.matrix(R)
  
  # make 1's on diagonal of R
  #
  inv.diag <- 1/diag(R)            # multiplicative inverse of diag(R)
  mydiag <- diag(length(inv.diag)) # initialize diagonal matrix for inv.diag
  diag(mydiag) <- inv.diag         # swap 1's for inv.diag
  mydiag <- sqrt(mydiag)           # take sqrt of diag matrix
  R <- mydiag %*% R %*% mydiag     # compute corrected correlation matrix with 1's on diagonal (Still Pos. Def.)
  
  # return correlation matrix, degree vector, adjacency matrix, and functional variables
  list(corrmat = R, deg.vec = kvec, A.mat = Adj, sig.vars = diff.cor.vars)
}

#=========================================================================================#
#' createSimulation2
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
#' @param pct.mixed fraction of functional variables with interaction effects
#' @param verbose logical indicating whether to display time required to generate simulation
#' @param plot.graph logical indicating whether to plot networks
#' @param use.Rcpp if true use Rcpp to correct negative eigenvalues 
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
#' num.samples <- 100
#' num.variables <- 10
#'
#' dataset <- createSimulation2(num.samples=num.samples,
#'                             num.variables=num.variables,
#'                             pct.imbalance=0.5,
#'                             pct.signals=0.2,
#'                             main.bias=0.4,
#'                             interaction.bias=1,
#'                             hi.cor=10,
#'                             lo.cor=0.2,
#'                             mix.type="main-interactionErdos",
#'                             label="class",
#'                             sim.type="mainEffect_Erdos-Renyi",
#'                             pct.mixed=0.5,
#'                             pct.train=0.5,
#'                             pct.holdout=0.5,
#'                             pct.validation=0,
#'                             plot.graph=F,
#'                             verbose=T)
#' @export
createSimulation2 <- function(num.samples=100,
                              num.variables=100,
                              pct.imbalance = 0.5,
                              pct.signals=0.1,
                              main.bias=0.4,
                              interaction.bias=0.4,
                              hi.cor=0.8,
                              lo.cor=0.2,
                              label="class",
                              sim.type="mainEffect",
                              pct.train=0.5,
                              pct.holdout=0.5,
                              pct.validation=0,
                              save.file=NULL,
                              mix.type=NULL,
                              pct.mixed=0.5,
                              verbose=FALSE,
                              plot.graph=F, use.Rcpp=F){
  
  ptm <- proc.time() # start time
  
  if (use.Rcpp){
    Rcpp::cppFunction(depends="RcppArmadillo",
                      'arma::vec getEigenValues(arma::mat M) {
                      return arma::eig_sym(M);
    }')
  }
  
  nbias <- pct.signals * num.variables # number of functional attributes
  
  if(sim.type=="mainEffect"){ # simple main effect simulation
    
    # new simulation:
    # sd.b sort of determines how large the signals are
    # p.b=0.1 makes 10% of the variables signal, bias <- 0.5
    my.sim.data <- privateEC::createMainEffects(n.e=num.variables,                   
                                                n.db=num.samples,              
                                                pct.imbalance=pct.imbalance,
                                                label = label,
                                                sd.b=main.bias,
                                                p.b=pct.signals)$db
    dataset <- cbind(t(my.sim.data$datnobatch), my.sim.data$S)
    
    # make numeric matrix into a data frame for splitting and subsequent ML algorithms
    dataset <- as.data.frame(dataset)
    
    signal.names <- paste("mainvar", 1:nbias, sep = "")                   # functional names
    background.names <- paste("var", 1:(num.variables - nbias), sep = "") # non-functional names
    var.names <- c(signal.names, background.names, label)                 # all variable names
    colnames(dataset) <- var.names                                        # replace column names with variable names
    
  }else if(sim.type=="mainEffect_Erdos-Renyi"){ # main + Erdos-Renyi network
    
    # create main-effect simulation
    my.sim.data <- privateEC::createMainEffects(n.e=num.variables,                   
                                                n.db=num.samples,              
                                                pct.imbalance=pct.imbalance,
                                                label = label,
                                                sd.b=main.bias,
                                                p.b=pct.signals)$db
    dataset <- cbind(t(my.sim.data$datnobatch), my.sim.data$S)
    
    e <- 1    # fudge factor to the number of nodes to avoid giant component
    prob <- 1/(num.variables+e) # probability of a node being connected to another node is less than 1/N to avoid giant component
    
    g <- igraph::erdos.renyi.game(num.variables, prob) # Erdos-Renyi network
    
    # generate correlation matrix from g
    network.atts <- generate_structured_corrmat(g=g,
                                                num.variables=num.variables, 
                                                hi.cor.tmp=hi.cor, 
                                                lo.cor.tmp=lo.cor, 
                                                graph.type="Erdos-Renyi",
                                                plot.graph=plot.graph,
                                                nbias=nbias, use.Rcpp=use.Rcpp)
    
    R <- as.matrix(network.atts$corrmat) # correlation matrix
    
    A.mat <- network.atts$A.mat          # adjacency from graph
    
    U <- t(chol(R))                             # upper tri cholesky
    tmp <- t(U %*% t(dataset[,-ncol(dataset)])) # correlated data
    tmp <- cbind(tmp,dataset[,ncol(dataset)])   # combine with phenotype
    dataset <- tmp                              # main-effect data with correlation
    
    # make numeric matrix into a data frame for splitting and subsequent ML algorithms
    dataset <- as.data.frame(dataset)
    
    signal.names <- paste("mainvar", 1:nbias, sep = "")                   # functional variable names
    background.names <- paste("var", 1:(num.variables - nbias), sep = "") # non-functional variable names
    var.names <- c(signal.names, background.names, label)                 # all variable names
    colnames(dataset) <- var.names                                        # replace column names with variable names
    
  }else if(sim.type=="mainEffect_Scalefree"){ # main + Scale-Free network
    
    # create main-effect simulation
    my.sim.data <- privateEC::createMainEffects(n.e=num.variables,                   
                                                n.db=num.samples,              
                                                pct.imbalance=pct.imbalance,
                                                label = label,
                                                sd.b=main.bias,
                                                p.b=pct.signals)$db
    dataset <- cbind(t(my.sim.data$datnobatch), my.sim.data$S)
    
    e <- 1    # fudge factor to the number of nodes to avoid giant component
    prob <- 1/(num.variables+e) # probability of a node being connected to another node is less than 1/N to avoid giant component
    
    g <- igraph::barabasi.game(num.variables, directed=F) # scale-free network
    
    # generate correlation matrix from g
    network.atts <- generate_structured_corrmat(g=g,
                                                num.variables=num.variables, 
                                                hi.cor.tmp=hi.cor, 
                                                lo.cor.tmp=lo.cor, 
                                                graph.type="Scalefree",
                                                plot.graph=plot.graph,
                                                nbias=nbias, use.Rcpp=use.Rcpp)
    R <- as.matrix(network.atts$corrmat) # correlation matrix
    
    A.mat <- network.atts$A.mat          # adjacency from graph
    
    U <- t(chol(R))                             # upper tri cholesky
    tmp <- t(U %*% t(dataset[,-ncol(dataset)])) # correlated data
    tmp <- cbind(tmp,dataset[,ncol(dataset)])   # combine with phenotype
    dataset <- tmp                              # main-effect data with correlation
    
    # make numeric matrix into a data frame for splitting and subsequent ML algorithms
    dataset <- as.data.frame(dataset)
    
    signal.names <- paste("mainvar", 1:nbias, sep = "")                   # functional variable names
    background.names <- paste("var", 1:(num.variables - nbias), sep = "") # non-functional variable names  
    var.names <- c(signal.names, background.names, label)                 # all variable names
    colnames(dataset) <- var.names                                        # replace column names with variable names
    
  }else if(sim.type=="interactionErdos"){ # simple interaction with Erdos-Renyi network
    
    m.case <- round((1 - pct.imbalance)*num.samples) # size of case group
    m.ctrl <- round(pct.imbalance*num.samples)       # size of ctrl group
    
    # random cases data set
    X.case <- matrix(rnorm(m.case*num.variables), nrow=m.case, ncol=num.variables)
    
    # approximate effect size for pairwise correlations between functional attributes in case group
    case.hi.cor <- -hi.cor*interaction.bias + (1 - interaction.bias)*hi.cor
    
    e <- 1    # fudge factor to the number of nodes to avoid giant component
    prob <- 1/(num.variables+e) # probability of a node being connected to another node is less than 1/N to avoid giant component
    
    g <- igraph::erdos.renyi.game(num.variables, prob) # Erdos-Renyi network
    
    # generate correlation matrix from g
    network.atts <- generate_structured_corrmat(g=g,
                                                num.variables=num.variables, 
                                                hi.cor.tmp=case.hi.cor,
                                                lo.cor.tmp=lo.cor,
                                                graph.type="Erdos-Renyi",
                                                plot.graph=plot.graph,
                                                make.diff.cors=T,
                                                nbias=nbias, use.Rcpp=use.Rcpp)
    
    R <- as.matrix(network.atts$corrmat) # correlation matrix for cases
    
    sig.vars <- network.atts$sig.vars    # column indices of functional features
    
    A.mat <- network.atts$A.mat          # adjacency from graph object
    
    U <- t(chol(R))           # upper tri cholesky
    tmp <- t(U %*% t(X.case)) # correlated case data
    X.case <- tmp
    
    # random controls data set
    X.ctrl <- matrix(rnorm(m.ctrl*num.variables), nrow=m.ctrl, ncol=num.variables)
    
    # for each functional variable, find connected variables
    sig.connected.list <- list()  # list for storing connected variable indices
    for(i in 1:length(sig.vars)){                                           # for each functional feature
      sig.connected.list[[i]] <- which(abs(A.mat[sig.vars[i],] - 1) < 1e-9) # find and store connected features
    }
    
    # create control group correlation matrix
    for(i in 1:length(sig.vars)){                                               # for each functional feature
      n.connected <- length(sig.connected.list[[i]])                            # how many connections
      for(j in 1:n.connected){                                                   # for each connection
        R[sig.vars[i],sig.connected.list[[i]][j]] <- min(hi.cor+rnorm(1,0,.1),1) # replace low correlation by high for ctrls
      }
    }
    
    # ensure that R is symmetric
    #
    tmp <- diag(R)       # diag of R
    diag(R) <- 0         # replace diag of R with 0's
    R[lower.tri(R)] <- 0 # make lower triangle zeros
    
    R <- R + t(R)        # make symmetric
    diag(R) <- tmp       # original diag of R
    
    # correct for negative eigenvalues so R is positive definite
    #
    if(use.Rcpp){ # compute eigenvalues and make diag matrix
      R.d <- Matrix(diag(sort(c(getEigenValues(R)),decreasing=T)),sparse=T)
    }else{ 
      R.d <- Matrix(diag(eigen(R)$values),sparse=T)
    }
    tmp <- c(diag(R.d))                                # vector of eigenvalues
    
    if (any(tmp<0)){              # if any eigenvalues are negative
      R.V <- eigen(R)$vectors     # compute eigenvectors,
      tmp[tmp<0] <- 1e-7          # make negative into small positive,
      diag(R.d) <- tmp            # replace diag of R.d with positive eigenvalues,
      R.fix <- R.V%*%R.d%*%t(R.V) # compute new R, and
      R <- R.fix                  # store in R
    }
    R <- as.matrix(R)
    
    # make 1's on diag of R
    inv.diag <- 1/diag(R)            # multiplicative inverse of diag(R)
    mydiag <- diag(length(inv.diag)) # initialize diagonal matrix for inv.diag
    diag(mydiag) <- inv.diag         # swap 1's for inv.diag
    mydiag <- sqrt(mydiag)           # compute sqrt of inv.diag 
    R <- mydiag %*% R %*% mydiag     # compute corrected R with 1's on diagonal (Still Pos. Def.)
    
    U <- t(chol(R))           # upper tri cholesky
    tmp <- t(U %*% t(X.ctrl)) # correlated ctrl data
    X.ctrl <- tmp             # store in X.ctrl
    
    X.all <- rbind(X.case, X.ctrl) # case/ctrl data
    X.all <- as.data.frame(X.all)  # make data.frame
    
    dataset <- X.all
    
    signal.names <- paste("intvar", 1:nbias, sep = "")                    # interaction variable names
    background.names <- paste("var", 1:(num.variables - nbias), sep = "") # non-functional variable names
    var.names <- c(signal.names, background.names, label)                 # all variable names
    colnames(dataset)[sig.vars] <- signal.names                           # replace functional colnames with interaction variable names
    colnames(dataset)[-sig.vars] <- background.names                      # replace non-functional colnames with non-functional colnames
    
    dataset <- cbind(dataset, c(rep(1,length=m.case), rep(-1,length=m.ctrl))) # full data set with phenotype
    colnames(dataset)[ncol(dataset)] <- "class"                               # replace last colname with "class"
    
  }else if(sim.type=="interactionScalefree"){ # simple interaction with Scale-Free network
    
    m.case <- round((1 - pct.imbalance)*num.samples) # size of case group
    m.ctrl <- round(pct.imbalance*num.samples)       # size of ctrl group
    
    # random cases data set
    X.case <- matrix(rnorm(m.case*num.variables), nrow=m.case, ncol=num.variables)
    
    # approximate effect size for pairwise correlations between functional attributes in case group 
    case.hi.cor <- -hi.cor*interaction.bias + (1 - interaction.bias)*hi.cor
    
    e <- 1    # fudge factor to the number of nodes to avoid giant component
    prob <- 1/(num.variables+e) # probability of a node being connected to another node is less than 1/N to avoid giant component
    
    g <- igraph::barabasi.game(num.variables, directed=F) # scale-free network
    
    # generate correlation matrix from g
    network.atts <- generate_structured_corrmat(g=g,
                                                num.variables=num.variables, 
                                                hi.cor.tmp=case.hi.cor,
                                                lo.cor.tmp=lo.cor,
                                                graph.type="Erdos-Renyi",
                                                plot.graph=plot.graph,
                                                make.diff.cors=T,
                                                nbias=nbias, use.Rcpp=use.Rcpp)
    
    R <- as.matrix(network.atts$corrmat)
    
    sig.vars <- network.atts$sig.vars # functional variable indices
    
    A.mat <- network.atts$A.mat       # adjacency from graph
    
    U <- t(chol(R))           # upper tri cholesky
    tmp <- t(U %*% t(X.case)) # correlated case data
    X.case <- tmp             # store in X.case
    
    # random controls data set
    X.ctrl <- matrix(rnorm(m.ctrl*num.variables), nrow=m.ctrl, ncol=num.variables)
    
    # for each functional variable, find connected variables
    #
    sig.connected.list <- list() # list for functional connections
    for(i in 1:length(sig.vars)){                                           # for each functional variable
      sig.connected.list[[i]] <- which(abs(A.mat[sig.vars[i],] - 1) < 1e-9) # find connections and store in list
    }
    
    # create control group correlation matrix
    for(i in 1:length(sig.vars)){                                                # for each functional variable 
      n.connected <- length(sig.connected.list[[i]])                             # how many connections
      for(j in 1:n.connected){                                                    # for each connection
        R[sig.vars[i],sig.connected.list[[i]][j]] <- min(hi.cor+rnorm(1,0,.1),1)  # replace low correlations with high for ctrls
      }
    }
    
    # make sure R is symmetric
    #
    tmp <- diag(R)       # diag of R
    diag(R) <- 0         # make diag(R) 0
    R[lower.tri(R)] <- 0 # make lower triangle of R 0
    
    R <- R + t(R)        # make symmetric
    diag(R) <- tmp       # replace diag(R) with original diag
    
    # correct for negative eigenvalues to make R positive definite
    #
    if(use.Rcpp){ # compute eigenvalues and make diag matrix
      R.d <- Matrix(diag(sort(c(getEigenValues(R)),decreasing=T)),sparse=T)
    }else{ 
      R.d <- Matrix(diag(eigen(R)$values),sparse=T)
    }
    tmp <- c(diag(R.d))                                # vector of eigenvalues
    
    if (any(tmp<0)){              # if any eigenvalues are negative
      R.V <- eigen(R)$vectors     # compute eigenvectors
      tmp[tmp<0] <- 1e-7          # make negative into small positive
      diag(R.d) <- tmp            # swap negative for positive eigenvalues
      R.fix <- R.V%*%R.d%*%t(R.V) # compute corrected correlation matrix
      R <- R.fix                  # store in R
    }
    R <- as.matrix(R)
    
    # make 1's on diagonal
    # 
    inv.diag <- 1/diag(R)            # multiplicative inverse of diag(R)
    mydiag <- diag(length(inv.diag)) # initialize diag matrix for inv.diag
    diag(mydiag) <- inv.diag         # swap 1's for inv.diag
    mydiag <- sqrt(mydiag)           # sqrt of inv.diag
    R <- mydiag %*% R %*% mydiag     # compute corrected R with 1's on diagonal (Still Pos. Def.)
    
    U <- t(chol(R))           # upper tri cholesky
    tmp <- t(U %*% t(X.ctrl)) # correlated ctrl data
    X.ctrl <- tmp             # store in X.ctrl
    
    X.all <- rbind(X.case, X.ctrl) # case/ctrl data
    X.all <- as.data.frame(X.all)  # make data.frame
    
    dataset <- X.all
    
    signal.names <- paste("intvar", 1:nbias, sep = "")                    # interaction variable names
    background.names <- paste("var", 1:(num.variables - nbias), sep = "") # non-functional variable names
    var.names <- c(signal.names, background.names, label)                 # all variable names
    colnames(dataset)[sig.vars] <- signal.names                           # replace functional colnames with interaction variable names
    colnames(dataset)[-sig.vars] <- background.names                      # replace non-functional colnames with non-functional variable names
    
    dataset <- cbind(dataset, c(rep(1,length=m.case), rep(-1,length=m.ctrl))) # case/ctrl data with phenotype
    colnames(dataset)[ncol(dataset)] <- "class"                               # replace last colname with "class"
    
  }else if(sim.type=="mixed"){ # main + interaction effect simulation
    
    if(mix.type=="main-interactionErdos"){ # main + interaction with Erdos-Renyi network
      
      m.case <- round((1 - pct.imbalance)*num.samples) # size of case group
      m.ctrl <- round(pct.imbalance*num.samples)       # size of ctrl group
      
      num.main <- round(pct.mixed*nbias)      # number of main effect attributes
      num.int <- round((1 - pct.mixed)*nbias) # number of interaction effect attributes
      
      # make cases data set
      X.case <- matrix(rnorm(m.case*num.variables), nrow=m.case, ncol=num.variables)
      
      # approximate effect size for pairwise correlations between functional attributes in case group
      case.hi.cor <- -hi.cor*interaction.bias + (1 - interaction.bias)*hi.cor
      
      e <- 1    # fudge factor to the number of nodes to avoid giant component
      prob <- 1/(num.variables+e) # probability of a node being connected to another node is less than 1/N to avoid giant component
      
      # generate random Erdos-Renyi network
      g <- igraph::erdos.renyi.game(num.variables, prob)
      
      # generate correlation matrix from g
      network.atts <- generate_structured_corrmat(g=g,
                                                  num.variables=num.variables, 
                                                  hi.cor.tmp=case.hi.cor,
                                                  lo.cor.tmp=lo.cor,
                                                  graph.type="Erdos-Renyi",
                                                  plot.graph=plot.graph,
                                                  make.diff.cors=T,
                                                  nbias=num.int, use.Rcpp=use.Rcpp)
      
      R <- as.matrix(network.atts$corrmat) # case correlation matrix
      
      sig.vars <- network.atts$sig.vars    # functional attribute indices
      
      A.mat <- network.atts$A.mat          # adjacency from graph
      
      # determine which variables are not connected in graph
      degs <- degree(g)
      n.unconnected.vars <- length(which(degs==0))
      if(n.unconnected.vars > 0){
        unconnected.vars <- which(degs==0)
      }else{
        unconnected.vars <- NA
      }
      
      U <- t(chol(R))           # upper tri cholesky
      tmp <- t(U %*% t(X.case)) # correlated case data
      X.case <- tmp             # store in X.case
      
      # make controls data set
      X.ctrl <- matrix(rnorm(m.ctrl*num.variables), nrow=m.ctrl, ncol=num.variables)
      
      # find connections for functional attributes
      #
      sig.connected.list <- list()  # list for functional attribute connections
      for(i in 1:length(sig.vars)){                                           # for each functional feature
        sig.connected.list[[i]] <- which(abs(A.mat[sig.vars[i],] - 1) < 1e-9) # find all connections
      }
      
      # make high correlations for ctrls
      #
      for(i in 1:length(sig.vars)){                                                # for each functional variable
        n.connected <- length(sig.connected.list[[i]])                             # how many connections
        for(j in 1:n.connected){                                                    # for each connection
          R[sig.vars[i],sig.connected.list[[i]][j]] <- min(hi.cor+rnorm(1,0,.1),1)  # replace with high correlation
        }
      }
      
      # ensure that R is symmetric
      #
      tmp <- diag(R)       # diag of R
      diag(R) <- 0         # make diag(R) 0
      R[lower.tri(R)] <- 0 # make lower triangle of R 0
      
      R <- R + t(R)  # make symmetric
      diag(R) <- tmp # swap for original diag
      
      # correct for negative eigenvalues so R is positive definite
      #
      if(use.Rcpp){ # compute eigenvalues and make diag matrix
        R.d <- Matrix(diag(sort(c(getEigenValues(R)),decreasing=T)),sparse=T)
      }else{ 
        R.d <- Matrix(diag(eigen(R)$values),sparse=T)
      }
      tmp <- c(diag(R.d))                                # vector of eigenvalues
      
      if (any(tmp<0)){              # if any eigenvalues are negative
        R.V <- eigen(R)$vectors     # compute eigenvectors,
        tmp[tmp<0] <- 1e-7          # make negative into small positive,
        diag(R.d) <- tmp            # swap diag of R.d with positive eigenvalues,
        R.fix <- R.V%*%R.d%*%t(R.V) # compute new correlation matrix, and 
        R <- R.fix                  # store in R
      }
      R <- as.matrix(R)
      
      # make 1's on diagonal
      #
      inv.diag <- 1/diag(R)            # multiplicative inverse of diag(R)
      mydiag <- diag(length(inv.diag)) # initialize diag matrix for inv.diag
      diag(mydiag) <- inv.diag         # swap 1's for inv.diag
      mydiag <- sqrt(mydiag)           # compute sqrt of inv.diag
      R <- mydiag %*% R %*% mydiag     # compute corrected R with 1's on diagonal (Still Pos. Def.)
      
      U <- t(chol(R))           # upper tri cholesky
      tmp <- t(U %*% t(X.ctrl)) # correlated ctrl data
      X.ctrl <- tmp             # store in X.ctrl
      
      X.all <- rbind(X.case, X.ctrl) # case/ctrl data
      X.all <- as.data.frame(X.all)  # make data.frame
      
      dataset <- X.all
      
      # create main-effect simulation for num.main attributes
      my.sim.data <- privateEC::createMainEffects(n.e=num.main,                   
                                                  n.db=num.samples,              
                                                  pct.imbalance=pct.imbalance,
                                                  label = label,
                                                  sd.b=main.bias,
                                                  p.b=1)$db
      
      # check dimensions to see if my.sim.data is a matrix or just a vector
      check.dim <- is.null(dim(my.sim.data$datnobatch))
      if(check.dim){
        main.cols <- my.sim.data$datnobatch
      }else{
        main.cols <- t(my.sim.data$datnobatch)
      }
      dataset.tmp <- cbind(main.cols, my.sim.data$S)
      
      # determine indices of main effect variables in full data set
      nonsig.vars <- seq(1,num.variables,by=1)[-sig.vars]
      alternate.vars <- seq(1,num.variables,by=1)[-unique(c(unlist(sig.connected.list),sig.vars))]
      if(length(na.omit(unconnected.vars)) >= num.main){                         # if there are enough non-connected to choose from
        main.vars <- sample(c(na.omit(unconnected.vars)), size=num.main, replace=F) # sample randomly from non-connected
      }else{ # if there are too few non-connected to choose from
        main.vars <- sample(unique(c(c(na.omit(unconnected.vars)), alternate.vars)), size=num.main, replace=F) # sample randomly from unconnected and non-interaction variables
      }
      
      # insert main effects into data set
      main.idx <- sample(c(1:num.main),size=num.main,replace=F)
      for(i in 1:length(main.vars)){
        dataset[, main.vars[i]] <- dataset.tmp[,main.idx[i]]
      }
      
      main.names <- paste("mainvar", 1:num.main, sep="") # main effect variable names
      int.names <- paste("intvar", 1:num.int, sep="")    # interaction effect variable names
      
      signal.names <- c(main.names, int.names)           # all functional variable names
      background.names <- paste("var", 1:(num.variables - nbias), sep = "") # all non-functional variable names
      var.names <- c(signal.names, background.names, label) # all variable names
      
      colnames(dataset)[sig.vars] <- int.names   # replace interaction colnames with interaction variable names
      colnames(dataset)[main.vars] <- main.names # replace main colnames with main variable names
      colnames(dataset)[-c(sig.vars, main.vars)] <- background.names # replace non-functional colnames with non-functional variable names
      
      dataset <- cbind(dataset, c(rep(1,length=m.case), rep(-1,length=m.ctrl))) # full data set with phenotype
      colnames(dataset)[ncol(dataset)] <- "class" # replace last colname with "class"
      
    }else if(mix.type=="main-interactionScalefree"){ # main + interaction with Scale-Free network
      
      m.case <- round((1 - pct.imbalance)*num.samples) # number of cases
      m.ctrl <- round(pct.imbalance*num.samples)       # number of ctrls
      
      num.main <- round(pct.mixed*nbias)      # number of main effect attributes
      num.int <- round((1 - pct.mixed)*nbias) # number of interaction effect attributes
      
      # make cases data set
      X.case <- matrix(rnorm(m.case*num.variables), nrow=m.case, ncol=num.variables)
      
      # approximate effect size for pairwise correlations between functional attributes in case group
      case.hi.cor <- -hi.cor*interaction.bias + (1 - interaction.bias)*hi.cor
      
      e <- 1    # fudge factor to the number of nodes to avoid giant component
      prob <- 1/(num.variables+e) # probability of a node being connected to another node is less than 1/N to avoid giant component
      
      # generate random Scale-Free network
      g <- igraph::barabasi.game(num.variables, directed=F)
      
      # generate case group correlation matrix
      network.atts <- generate_structured_corrmat(g=g,
                                                  num.variables=num.variables, 
                                                  hi.cor.tmp=case.hi.cor,
                                                  lo.cor.tmp=lo.cor,
                                                  graph.type="Scalefree",
                                                  plot.graph=plot.graph,
                                                  make.diff.cors=T,
                                                  nbias=num.int, use.Rcpp=use.Rcpp)
      
      R <- as.matrix(network.atts$corrmat) # cases correlation matrix
      
      sig.vars <- network.atts$sig.vars    # functional attribute indices
      
      A.mat <- network.atts$A.mat          # adjacency from graph
      
      # find non-connected variables from graph
      degs <- degree(g)
      n.unconnected.vars <- length(which(degs==0))
      if(n.unconnected.vars > 0){
        unconnected.vars <- which(degs==0)
      }else{
        unconnected.vars <- NA
      }
      
      U <- t(chol(R))           # upper tri cholesky
      tmp <- t(U %*% t(X.case)) # correlated case data
      X.case <- tmp             # store in X.case
      
      # make controls data set
      X.ctrl <- matrix(rnorm(m.ctrl*num.variables), nrow=m.ctrl, ncol=num.variables)
      
      # find all functional connections
      # 
      sig.connected.list <- list()  # list for storing functional connections
      for(i in 1:length(sig.vars)){                                           # for each functional variable
        sig.connected.list[[i]] <- which(abs(A.mat[sig.vars[i],] - 1) < 1e-9) # find all connections
      }
      
      # make high correlations matrix for ctrls
      #
      for(i in 1:length(sig.vars)){                                                # for each functional variable
        n.connected <- length(sig.connected.list[[i]])                             # how many connections
        for(j in 1:n.connected){                                                    # for each connection
          R[sig.vars[i],sig.connected.list[[i]][j]] <- min(hi.cor+rnorm(1,0,.1),1)  # replace low correlation with high correlation
        }
      }
      
      # ensure that R is symmetric
      #
      tmp <- diag(R)       # diag of R
      diag(R) <- 0         # make diag(R) 0
      R[lower.tri(R)] <- 0 # make lower triangle of R 0
      
      R <- R + t(R)  # make symmetric
      diag(R) <- tmp # swap for original diag
      
      # correct for negative eigenvalues to make R positive definite
      #
      if(use.Rcpp){ # compute eigenvalues and make diag matrix
        R.d <- Matrix(diag(sort(c(getEigenValues(R)),decreasing=T)),sparse=T)
      }else{ 
        R.d <- Matrix(diag(eigen(R)$values),sparse=T)
      }
      tmp <- c(diag(R.d))                                # vector of eigenvalues
      
      if (any(tmp<0)){              # if any eigenvalues are negative
        R.V <- eigen(R)$vectors     # compute eigenvectores,
        tmp[tmp<0] <- 1e-7          # make negative into small positive,
        diag(R.d) <- tmp            # replace diag R.d with positive eigenvalues
        R.fix <- R.V%*%R.d%*%t(R.V) # compute new correlation matrix, and 
        R <- R.fix                  # store in R
      }
      R <- as.matrix(R)
      
      # make 1's on diagonal
      #
      inv.diag <- 1/diag(R)            # multiplicative inverse of diag(R)
      mydiag <- diag(length(inv.diag)) # initialize diag matrix for inv.diag
      diag(mydiag) <- inv.diag         # swap 1's for inv.diag
      mydiag <- sqrt(mydiag)           # sqrt of inv.diag
      R <- mydiag %*% R %*% mydiag     # compute corrected R with 1's on diagonal (Still Pos. Def.)
      
      U <- t(chol(R))                  # upper tri cholesky
      tmp <- t(U %*% t(X.ctrl))        # correlated ctrl data
      X.ctrl <- tmp                    # store in X.ctrl
      
      X.all <- rbind(X.case, X.ctrl)   # case/ctrl data
      X.all <- as.data.frame(X.all)    # make data.frame
      
      dataset <- X.all
      
      # create main effect simulation with num.main features
      my.sim.data <- privateEC::createMainEffects(n.e=num.main,                   
                                                  n.db=num.samples,              
                                                  pct.imbalance=pct.imbalance,
                                                  label = label,
                                                  sd.b=main.bias,
                                                  p.b=1)$db
      
      # check dimensions to determine if my.dim.data is a matrix or vector
      check.dim <- is.null(dim(my.sim.data$datnobatch))
      if(check.dim){
        main.cols <- my.sim.data$datnobatch
      }else{
        main.cols <- t(my.sim.data$datnobatch)
      }
      dataset.tmp <- cbind(main.cols, my.sim.data$S)
      
      # determine indices of main effect variables
      nonsig.vars <- seq(1,num.variables,by=1)[-sig.vars]
      alternate.vars <- seq(1,num.variables,by=1)[-unique(c(unlist(sig.connected.list),sig.vars))]
      if(length(na.omit(unconnected.vars)) >= num.main){                         # if there are enough non-connected variables
        main.vars <- sample(c(na.omit(unconnected.vars)), size=num.main, replace=F) # randomly sample from non-connected variables
      }else{ # if there are not enough non-connected variables
        main.vars <- sample(unique(c(c(na.omit(unconnected.vars)), alternate.vars)), size=num.main, replace=F) # randomly sample from non-connected and non-functional variables
      }
      
      # insert main effects into data set
      main.idx <- sample(c(1:num.main),size=num.main,replace=F)
      for(i in 1:length(main.vars)){
        dataset[, main.vars[i]] <- dataset.tmp[,main.idx[i]]
      }
      
      main.names <- paste("mainvar", 1:num.main, sep="") # main effect variable names
      int.names <- paste("intvar", 1:num.int, sep="")    # interaction effect variable names
      
      signal.names <- c(main.names, int.names) # all functional variable names
      background.names <- paste("var", 1:(num.variables - nbias), sep = "") # non-functional variable names
      var.names <- c(signal.names, background.names, label) # all variable names
      
      colnames(dataset)[sig.vars] <- int.names   # replace interaction colnames with interaction variable names
      colnames(dataset)[main.vars] <- main.names # replace main colnames with main variable names
      colnames(dataset)[-c(sig.vars, main.vars)] <- background.names # replace non-functional colnames with non-functional variable names
      
      dataset <- cbind(dataset, c(rep(1,length=m.case), rep(-1,length=m.ctrl))) # full data set with phenotype
      colnames(dataset)[ncol(dataset)] <- "class" # replace last colname with "class"
      
    }
    
  }
  
  # split data into train, holdout, and validation sets
  split.data <- privateEC::splitDataset(all.data = dataset,
                                        pct.train = pct.train,
                                        pct.holdout = pct.holdout,
                                        pct.validation = pct.validation,
                                        label = label)
  
  if (!is.null(save.file)) {
    if (verbose) cat("saving to data/", save.file, ".Rdata\n", sep = "")
    save(split.data, pct.signals, main.bias, interaction.bias,sim.type, file = save.file)
  }
  
  elapsed <- (proc.time() - ptm)[3]
  if (verbose) cat("createSimulation elapsed time:", elapsed, "\n")
  
  # depending on sim.type, we have different potentially important outputs
  if(sim.type %in% c("interactionErdos","interactionScalefree")){
    
    # for interaction simulation only
    list(train = split.data$train,
         holdout = split.data$holdout,
         validation = split.data$validation,
         label = label,
         signal.names = signal.names,
         elapsed = elapsed,
         sig.vars=sig.vars,
         A.mat=A.mat)
    
  }else if(sim.type == "mixed"){
    
    # for mixed simulation (main+interaction)
    list(train = split.data$train,
         holdout = split.data$holdout,
         validation = split.data$validation,
         label = label,
         signal.names = signal.names,
         elapsed = elapsed,
         int.vars=sig.vars,
         main.vars=main.vars,
         A.mat=A.mat)
    
  }else if(sim.type == "mainEffect_Erdos-Renyi"){
    
    # for main effect + correlation simulation (Erdos-Renyi Network)
    list(train = split.data$train,
         holdout = split.data$holdout,
         validation = split.data$validation,
         label = label,
         signal.names = signal.names,
         elapsed = elapsed,
         A.mat=A.mat)
    
  }else if(sim.type == "mainEffect_Scalefree"){
    
    # for main effect + correlation simulation (Scale-Free Network)
    list(train = split.data$train,
         holdout = split.data$holdout,
         validation = split.data$validation,
         label = label,
         signal.names = signal.names,
         elapsed = elapsed,
         A.mat=A.mat)
    
  }else{
    
    # for simple main effect simulation
    list(train = split.data$train,
         holdout = split.data$holdout,
         validation = split.data$validation,
         label = label,
         signal.names = signal.names,
         elapsed = elapsed)
  }
}
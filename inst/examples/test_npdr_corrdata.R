#library(devtools)
#install_github("insilico/npdr")
library(npdr)

p <- 100
attr.idx.list <- list()
for(i in 1:p){
  
  lo.idx <- (i - 1)*(p-1) + 1
  hi.idx <- i*(p-1)
  
  attr.idx.list[[i]] <- c(lo.idx:hi.idx)
  
}

m <- 100
p <- 100

myfisher <- function(rho){
  
  z <- 0.5*log((1+rho)/(1-rho))
  return(z)
}

stretch_mat <- function(M){
  
  mat <- numeric()
  for(k in 1:nrow(M)){
    mat <- c(mat,M[k,-k])
  }
  return(mat)
}

data.mat <- matrix(rep(0,length=m*(p*(p-1))),nrow=(p*(p-1)),ncol=m)
for(k in 1:m){
  #print(k)
  
  X <- matrix(runif(p^2,min=-1,max=1), ncol=p) 
  
  cov <- X %*% t(X) 
  
  mycorr <- cov2cor(cov) 
  
  zcorr <- apply(matrix(stretch_mat(mycorr),ncol=1),1,myfisher)
  
  data.mat[,k] <- matrix(zcorr,nrow=(p*(p-1)),ncol=1)
}

var.names <- c(paste("main",1:10,sep=""),
               paste("var",1:90,sep=""))

class.vec <- c(rep(1,length=round(m/2)), rep(0,length=round(m/2)))
case.control.data <- t(data.mat)
#case.control.data <- cbind(case.control.data, class.vec)
colnames(case.control.data) <- c((1:(ncol(case.control.data))))
row.names(case.control.data) <- c(paste("case",1:floor(m/2),sep=""),
                                  paste("ctrl",1:floor(m/2),sep=""))

# m=100 subjects
# p=100 attributes, p*(p-1)=9,900 correlation predictors
npdr.cc.results <- npdr(class.vec, case.control.data, regression.type="binomial", 
                        attr.diff.type="correlation-data",
                        nbd.method="relieff", nbd.metric = "manhattan", msurf.sd.frac=.5, k=0,
                        neighbor.sampling="none",
                        dopar.nn = F, padj.method="bonferroni", verbose=T,
                        corr.attr.names=var.names)

data.frame(cbind(att=npdr.cc.results$att,
                 beta=npdr.cc.results$beta.Z.att,
                 pval=npdr.cc.results$pval.adj))[1:10,]

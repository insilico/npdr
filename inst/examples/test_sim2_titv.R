# test gwas added to createSimulation2()
library(npdr)

num.samples <- 100
num.variables <- 100
main.bias <- 0.4
pct.mixed <- 0.5
nbias <- round(pct.mixed*num.variables)
num.main <- round(pct.mixed*nbias) 
num.int <- round((1 - pct.mixed)*nbias)

sim.type <- "interactionErdos"
mix.type="main-interactionErdos"

dataset <- createSimulation2(num.samples=num.samples,
                             num.variables=num.variables,
                             pct.imbalance=0.5,
                             pct.signals=0.1,
                             main.bias=main.bias,
                             interaction.bias=1,
                             hi.cor=0.99,
                             lo.cor=0.1,
                             mix.type=mix.type,
                             label="class",
                             sim.type=sim.type,
                             pct.mixed=0.5,
                             pct.train=0.5,
                             pct.holdout=0.5,
                             pct.validation=0,
                             plot.graph=F,
                             verbose=T,
                             use.Rcpp=F,
                             prob.connected=0.95,
                             out.degree=98,
                             data.type="discrete")
dats <- rbind(dataset$train, dataset$holdout, dataset$validation)
dats <- dats[order(dats[,ncol(dats)]),]

case.control.data <- dats[,-ncol(dats)]

npdr.cc.results <- npdr(dats[,ncol(dats)], case.control.data, regression.type="binomial", attr.diff.type="allele-sharing",
                        nbd.method="relieff", nbd.metric = "manhattan", msurf.sd.frac=.5, k=0,
                        neighbor.sampling="none",
                        dopar.nn = F, padj.method="bonferroni", verbose=T)

data.frame(cbind(att=npdr.cc.results$att,
                 beta=npdr.cc.results$beta.Z.att,
                 pval=npdr.cc.results$pval.adj,
                 eff.size=exp(npdr.cc.results$beta.raw.att)))[1:10,]

d.mat <- npdrDistances(case.control.data,metric="allele-sharing-manhattan")
case.d <- d.mat[c(51:100),c(51:100)]
ctrl.d <- d.mat[c(1:50),c(1:50)]
hist(case.d[upper.tri(case.d)],breaks=30,col='red',freq=F)
hist(ctrl.d[upper.tri(ctrl.d)],breaks=30,col='blue',add=T,freq=F)

hist(d.mat[upper.tri(d.mat)],breaks=30)

p <- dim(case.control.data)[2]
m <- dim(case.control.data)[1]

my.eta <- 2
g1 <- 1/(my.eta+1)
g0 <- runif(1,min=0.01,max=(my.eta/(my.eta+1) - 0.01))
g2 <- my.eta/(my.eta + 1) - g0

PuPy.samp <- sample(c(3,4,5),size=p,prob=c(g0,g1,g2),replace=T)
data.mat <- as.matrix(dats[,-ncol(dats)])

d.mat <- matrix(rep(0,length=m*m),nrow=m,ncol=m)
for(i in 1:m){
  
  for(j in i:m){
    diff <- abs(data.mat[i,]-data.mat[j,])
    diff.TiTv <- diff
    for(k in 1:p){
      
      if(diff[k]==1){
        if(PuPy.samp[k]==3){
          diff.TiTv[k] <- 0.25
        }else if(PuPy.samp[k]==5){
          diff.TiTv[k] <- 0.25
        }else{
          diff.TiTv[k] <- 0.5
        }
      }else if(diff[k]==2){
        if(PuPy.samp[k]==3){
          diff.TiTv[k] <- 0.75
        }else if(PuPy.samp[k]==5){
          diff.TiTv[k] <- 0.75
        }else{
          diff.TiTv[k] <- 1
        }
      }
    }
    d.mat[i,j] <- sum(diff.TiTv)
    
  }
}
tmp <- d.mat + t(d.mat)

if(sim.type=="mainEffect" || sim.type=="mainEffect_Erdos-Renyi" || sim.type=="mainEffect_Scalefree"){
  probs.cases <- rep((1 + main.bias)/2,length=round(num.variables*0.1))
  probs.ctrls <- rep((1 - main.bias)/2,length=round(num.variables*0.1))
  probs.null <- rep(mean(c(probs.cases,probs.ctrls)),length=round(num.variables*0.9))
  
  probs.cases <- c(probs.cases,probs.null)
  F.a.case <- ((1 - probs.cases)^3)*probs.cases + (probs.cases^3)*(1 - probs.cases)
  G.a.case <- ((1 - probs.cases)^2)*(probs.cases^2)
  mean.case <- (g0+g2+2*g1)*sum(F.a.case) + (1.5*(g0 + g2) + 2*g1)*sum(G.a.case)
  
  probs.ctrls <- c(probs.ctrls,probs.null)
  F.a.ctrl <- ((1 - probs.ctrls)^3)*probs.ctrls + (probs.ctrls^3)*(1 - probs.ctrls)
  G.a.ctrl <- ((1 - probs.ctrls)^2)*(probs.ctrls^2)
  mean.ctrl <- (g0+g2+2*g1)*sum(F.a.ctrl) + (1.5*(g0 + g2) + 2*g1)*sum(G.a.ctrl)
}else if(sim.type=="interactionErdos" || sim.type=="interactionScalefree"){
  probs.cases <- runif(num.variables,min=0.1,max=0.4)
  probs.ctrls <- probs.cases
  
  F.a.case <- ((1 - probs.cases)^3)*probs.cases + (probs.cases^3)*(1 - probs.cases)
  G.a.case <- ((1 - probs.cases)^2)*(probs.cases^2)
  mean.case <- (g0+g2+2*g1)*sum(F.a.case) + (1.5*(g0 + g2) + 2*g1)*sum(G.a.case)
  
  F.a.ctrl <- ((1 - probs.ctrls)^3)*probs.ctrls + (probs.ctrls^3)*(1 - probs.ctrls)
  G.a.ctrl <- ((1 - probs.ctrls)^2)*(probs.ctrls^2)
  mean.ctrl <- (g0+g2+2*g1)*sum(F.a.ctrl) + (1.5*(g0 + g2) + 2*g1)*sum(G.a.ctrl)
}else if(sim.type=="mixed"){
  probs.cases <- runif((num.variables - num.main),min=0.1,max=0.4)
  probs.ctrls <- probs.cases
  
  probs.cases.tmp <- rep((1 + main.bias)/2,length=num.main)
  probs.ctrls.tmp <- rep((1 - main.bias)/2,length=num.main)
  
  probs.cases <- c(probs.cases,probs.cases.tmp)
  F.a.case <- ((1 - probs.cases)^3)*probs.cases + (probs.cases^3)*(1 - probs.cases)
  G.a.case <- ((1 - probs.cases)^2)*(probs.cases^2)
  mean.case <- (g0+g2+2*g1)*sum(F.a.case) + (1.5*(g0 + g2) + 2*g1)*sum(G.a.case)
  
  probs.ctrls <- c(probs.ctrls,probs.ctrls.tmp)
  F.a.ctrl <- ((1 - probs.ctrls)^3)*probs.ctrls + (probs.ctrls^3)*(1 - probs.ctrls)
  G.a.ctrl <- ((1 - probs.ctrls)^2)*(probs.ctrls^2)
  mean.ctrl <- (g0+g2+2*g1)*sum(F.a.ctrl) + (1.5*(g0 + g2) + 2*g1)*sum(G.a.ctrl)
}

dist.vec <- tmp[upper.tri(tmp)]
d.mat.cross <- tmp[1:50,51:100]
dist.vec.cross <- d.mat.cross[upper.tri(d.mat.cross)]

par(mfrow=c(1,1),mar=c(4.5,4.1,1.7,0.8))
hist(dist.vec,breaks=30,freq=F,ylim=c(0,0.25),
     main="Histogram of TiTv Distances with Main Effect",
     xlab="TiTv Distance",ylab="Density",font.lab=2,cex.lab=1.5,cex.main=1.7)
abline(v=mean.case,lty=2,lwd=2,col='red')
abline(v=mean.ctrl,lty=2,lwd=2,col='blue')
abline(v=mean(dist.vec.cross),lty=2,lwd=2,col='green')
legend("topleft",c("Predicted Case-Case","Predicted Ctrl-Ctrl","Sample Case-Ctrl"),
       lty=2,lwd=2,col=c('red','blue','green'),bg='white',cex=1.3)
box()

##############################################################


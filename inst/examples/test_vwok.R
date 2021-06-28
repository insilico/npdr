#library(devtools)
#install_github(insilico\npdr)
library(npdr)

####################################################################################
# mixed simulation parameters
####################################################################################
main.bias <- 0.6
interaction.bias <- 1

pct.imbalance <- 0.5

pct.signals <- 0.1

sim.type <- "interactionErdos"

mix.type <- "main-interactionErdos"
#mix.type <- "main-interactionScalefree"

prob.connected <- 0.95
out.degree <- 90

hi.cor <- 0.9
lo.cor <- 0.1

verbose <- T

data.type <- "continuous"

m <- 100
p <- 100

pct.mixed <- 0.5
####################################################################################
# (end) mixed simulation parameters
####################################################################################

dataset <- createSimulation2(num.samples=m,
                             num.variables=p,
                             pct.imbalance=pct.imbalance,
                             pct.signals=pct.signals,
                             main.bias=main.bias,
                             interaction.bias=interaction.bias,
                             hi.cor=hi.cor,
                             lo.cor=lo.cor,
                             mix.type=mix.type,
                             label="class",
                             sim.type=sim.type,
                             pct.mixed=pct.mixed,
                             pct.train=0.5,
                             pct.holdout=0.5,
                             pct.validation=0,
                             plot.graph=FALSE,
                             verbose=TRUE,
                             use.Rcpp=F,
                             prob.connected=prob.connected,
                             out.degree=out.degree,
                             data.type=data.type)
dats <- rbind(dataset$train, dataset$holdout, dataset$validation)
dats <- dats[order(dats[,ncol(dats)]),]

vars <- colnames(dats)[-ncol(dats)]
int.vars <- grep(vars,pattern="intvar",value=T)
main.vars <- grep(vars,pattern="mainvar",value=T)

signal.names <- c(int.vars,main.vars)

# run variable-wise adaptive k algorithm
out <- vwok(dats=dats,
            #k.grid=seq(1,(m-1),by=1), # default
            verbose=F,
            signal.names=signal.names,
            label="class")

idx.func <- which((as.character(out$vwok.out[,"att"]) %in% signal.names))
func.betas <- out$vwok.out[idx.func,"betas"]
neg.betas <- out$vwok.out[-idx.func,"betas"]
pr.vwok <- PRROC::pr.curve(scores.class0 = func.betas,
                           scores.class1 = neg.betas, 
                           curve = T)

auPRC.vwok <- pr.vwok$auc.integral # area under the precision-recall curve (vwok)
auPRC.vwok

plot(pr.vwok$curve[,1],pr.vwok$curve[,2],type="l",lwd=2,col='red',xlab="Recall",ylab="Precision",font.lab=2,
     ylim=c(0,1))

# average k that optimizes NPDR betas
vwok.k.all <- round(mean(out$vwok.out[,"best.ks"])) # avg of all atts
idx <- grep(as.character(out$vwok.out$att),pattern="intvar")
vwok.k.func <- round(mean(out$vwok.out[idx,"best.ks"])) # avg of func atts

# predicted multisurf k
predk.vec <- knnSURF(m,sd.frac=0.5)

# global k -- best k for auPRC if you know the functional attributes
globk.vec <- as.numeric(as.character(out$best.auPRC.k[,"k"]))

# top 10 features
print(out$vwok.out[1:10,]) # vwok

# auPRC for global k (optimizes auPRC directly)
out$best.auPRC.k[,"auPRC"]

# auPRC for vwok (optimizes attr scores only)
auPRC.vwok

# global k -- best k for auPRC if you know the functional attributes
globk.vec
# avg vwok k (all)
vwok.k.all
# avg vwok k (func) -- average of k's that give hightest beta for each attribute
vwok.k.func 

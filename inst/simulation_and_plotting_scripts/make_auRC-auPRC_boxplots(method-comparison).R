# make auRC and auPRC method comparison plots
#setwd("C:/Users/bdawk/Documents/KNN_project_output") will need to change to desired directory

# sim.type (options)
#
# "mainEffect": simple main effects
# "mainEffect_Erdos-Renyi": main effects with added correlation from Erdos-Renyi network
# "mainEffect_Scalefree": main effects with added correlation from Scale-free network
# "interactionErdos": interaction effects from Erdos-Renyi network
# "interactionScalefree": interaction effects from Scale-free network
# "mixed": main effects and interaction effects
#     mix.type (options)
#
#     "main-interactionErdos": main effects and interaction effects from Erdos-Renyi network
#     "main-interactionScalefree": main effects and interaction effects from Scale-free network

# data.type (options)
#
# "continuous": random normal data N(0,1) (e.g., gene expression data)
# "discrete": random binomial data B(n=2,prob) (e.g., GWAS data)

# parameters for file name
sim.type <- "interactionErdos" # type of simulated effects
data.type <- "continuous"  # or "discrete" for GWAS
data.type <- "discrete"

file <- paste("auRC-auPRC_iterates_methods-comparison_",data.type,"-",sim.type,".csv",sep="")

df <- read.csv(file,header=T)

num.iter <- nrow(df) # number of replicate data sets in file

method.vec <- rep(c('RForest','NPDR.MultiSURF','NPDR.Fixed.k','Relief'), each=num.iter)
auRC.vec <- c(df[,"auRC.RForest"], df[,"auRC.NPDR.MultiSURF"], df[,"auRC.NPDR.Fixed.k"], df[,"auRC.Relief"])
df.plot <- data.frame(method=method.vec,auRC=auRC.vec)
df.plot[,"auRC"] <- as.numeric(as.character(df.plot[,"auRC"]))
df.plot[,"method"] <- as.factor(df.plot[,"method"])
head(df.plot)

# area under recall curve box plots
num.methods <- length(unique(df.plot[,"method"]))
p1 <- ggplot(df.plot, aes(x=method, y=auRC, color=method)) + 
  geom_boxplot(outlier.shape=NA,colour=rep('#636363',length=num.methods)) +
  geom_jitter(shape=16,position=position_jitter(0.2),alpha=0.4,size=4.5, 
              aes(color=method), inherit.aes=T)+
  scale_y_continuous(name="Area under Recall Curve", limits=c(0,1),breaks=seq(0,1,by=0.2))+
  scale_x_discrete(breaks=c("NPDR.Fixed.k","NPDR.MultiSURF","Relief","RForest"),
                   labels=c("NPDR (fixed-k)","NPDR (MS)","Relief","Random Forest"))+
  scale_color_manual(breaks=c("NPDR.MultiSURF","NPDR.Fixed.k","Relief","RForest"),
                     values = c("magenta", "brown", "#0072B2", "#228B22"),
                     labels=c("NPDR.MultiSURF","NPDR.Fixed.k","Relief","RForest"))+
  theme_bw() +
  theme(axis.text=element_text(size=14,face='bold'),
        axis.title=element_text(size=16,face='bold'),
        axis.title.x=element_blank(),
        legend.position="none")
p1

method.vec <- rep(c('RForest','NPDR.MultiSURF','NPDR.Fixed.k','Relief'), each=num.iter)
auPRC.vec <- c(df[,"auPRC.RForest"], df[,"auPRC.NPDR.MultiSURF"], df[,"auPRC.NPDR.Fixed.k"], df[,"auPRC.Relief"])
df.plot <- data.frame(method=method.vec,auPRC=auPRC.vec)
df.plot[,"auPRC"] <- as.numeric(as.character(df.plot[,"auPRC"]))
df.plot[,"method"] <- as.factor(df.plot[,"method"])
head(df.plot)

# area under precision-recall curve box plots
p2 <- ggplot(df.plot, aes(x=method, y=auPRC, color=method)) + 
  geom_boxplot(outlier.shape=NA,colour=rep('#636363',length=num.methods)) +
  geom_jitter(shape=16,position=position_jitter(0.2),alpha=0.4,size=4.5, 
              aes(color=method), inherit.aes=T)+
  scale_y_continuous(name="Area under Precision-Recall curve", limits=c(0,1),breaks=seq(0,1,by=0.2))+
  scale_x_discrete(breaks=c("NPDR.Fixed.k","NPDR.MultiSURF","Relief","RForest"),
                   labels=c("NPDR (fixed-k)","NPDR (MS)","Relief","Random Forest"))+
  scale_color_manual(breaks=c("NPDR.MultiSURF","NPDR.Fixed.k","Relief","RForest"),
                     values = c("magenta", "brown", "#0072B2", "#228B22"),
                     labels=c("NPDR.MultiSURF","NPDR.Fixed.k","Relief","RForest"))+
  theme_bw() +
  theme(axis.text=element_text(size=14,face='bold'),
        axis.title=element_text(size=16,face='bold'),
        axis.title.x=element_blank(),
        legend.position="none")
p2

# side-by-side plots
cowplot::plot_grid(p1,p2,labels='AUTO',ncol=2, label_size=19)

file <- paste("auRC-auPRC_iterates_methods-comparison_",data.type,"-",sim.type,"_boxplots.pdf",sep="")
ggsave(file,height=7,width=14.2)

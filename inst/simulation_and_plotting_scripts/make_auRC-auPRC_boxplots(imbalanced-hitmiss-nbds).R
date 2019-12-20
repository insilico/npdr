# create box plots comparing balanced and imbalanced hit/miss NPDR nbds
library(ggplot2)
library(cowplot)
theme_set(theme_bw())

#setwd("C:/Users/bdawk/Documents/KNN_project_output") will need to change to desired directory

sim.type <- "interactionErdos" # simulation type
data.type <- "discrete"        # data type: discrete ro continuous
num.imbalances <- 5            # number of imbalances considered 5 = c(0.1, 0.2, 0.3, 0.4, 0.5)
num.replicates <- 1            # number of replicate data sets from simulations (line 50 in generate_auRC-auPRC_replicates(imbalanced-data).R)

# create data frame for plotting
df.plot <- NULL
for(iter in 1:num.imbalances){
  
  file <- paste("separate-hitmiss-nbds-T_",sim.type,"_",data.type,"_imbalance-",iter,".csv",sep="")
  df.T <- read.csv(file, header=T)
  file <- paste("separate-hitmiss-nbds-F_",sim.type,"_",data.type,"_imbalance-",iter,".csv",sep="")
  df.F <- read.csv(file, header=T)
  
  df.plot <- rbind(df.plot, df.T, df.F)
  
}

method <- c(rep('Balanced',length=num.replicates),rep('Imbalanced',length=num.replicates))
methods <- rep(method,length=nrow(df.plot))

mixed.pct <- c(rep(0.1,length=length(method)),
               rep(0.2,length=length(method)),
               rep(0.3,length=length(method)),
               rep(0.4,length=length(method)),
               rep(0.5,length=length(method)))
df.plot <- cbind(df.plot,methods=methods,pct=mixed.pct)

df.plot[,"auRC"] <- as.numeric(as.character(df.plot[,"auRC"]))
df.plot[,"methods"] <- as.factor(df.plot[,"methods"])
df.plot[,"pct"] <- as.factor(df.plot[,"pct"])
head(df.plot)

# area under recall curve replicate box plots

# title of plot based on simulation type
df.plot$title <- paste(sim.type,": Balanced vs Imbalanced NPDR",sep="")

ggplot(df.plot, aes(x=pct, y=auRC, fill=methods)) + 
  geom_boxplot(outlier.shape=NA,colour=rep('#636363',length=10)) +
  scale_fill_manual(values = rep("white",length=10))+
  geom_jitter(shape=16,position=position_jitterdodge(),alpha=0.4,size=3.5, 
              aes(x=pct,y=auRC,color=methods, group=methods), inherit.aes=F)+
  scale_color_manual(breaks=c('Balanced','Imbalanced'),
                     values = c("#228B22","#0072B2"),
                     labels=c('Balanced','Imbalanced'))+
  scale_y_continuous(name="Area under the Recall curve", limits=c(0,1),breaks=seq(0,1,by=0.2))+
  facet_grid(. ~ title)+
  theme(axis.text=element_text(size=14,face='bold'),
        axis.title=element_text(size=19,face='bold'),
        axis.title.x = element_text(size=19,face='bold'),
        strip.text=element_text(size=16,face='bold'),
        legend.title=element_text(size=12,face='bold'),
        legend.title.align=0.5)
ggsave(paste("auRC_comparison_",sim.type,"(balanced_and_imbalanced-NPDR-nbds).pdf",sep=""), height=6, width=7.2)


# area under precision-recall curve replicate box plots

# title of plot based on simulation type
df.plot$title <- paste(sim.type,": Balanced vs Imbalanced NPDR",sep="")

ggplot(df.plot, aes(x=pct, y=auPRC, fill=methods)) + 
  geom_boxplot(outlier.shape=NA,colour=rep('#636363',length=10)) +
  scale_fill_manual(values = rep("white",length=10))+
  geom_jitter(shape=16,position=position_jitterdodge(),alpha=0.4,size=3.5, 
              aes(x=pct,y=auPRC,color=methods, group=methods), inherit.aes=F)+
  scale_color_manual(breaks=c('Balanced','Imbalanced'),
                     values = c("#228B22","#0072B2"),
                     labels=c('Balanced','Imbalanced'))+
  scale_y_continuous(name="Area under Precision-Recall curve", limits=c(0,1),breaks=seq(0,1,by=0.2))+
  facet_grid(. ~ title)+
  theme(axis.text=element_text(size=14,face='bold'),
        axis.title=element_text(size=19,face='bold'),
        axis.title.x = element_text(size=19,face='bold'),
        strip.text=element_text(size=16,face='bold'),
        legend.title=element_text(size=12,face='bold'),
        legend.title.align=0.5)
ggsave(paste("auRC_comparison_",sim.type,"(balanced_and_imbalanced-NPDR-nbds).pdf",sep=""), height=6, width=7.2)

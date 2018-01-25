setwd("~/Data/Nurse_Larva/")

library(plyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(edgeR)
library(topGO)
library(scales)
library(ggmosaic)
library(questionr)
library(VennDiagram)

PS_palette = c("#deebf7","#9ecae1","#4292c6","#2171b5","#084594","darkblue")
tissue_palette = c("#9e9ac8","#6a51a3","cornflowerblue","royalblue4")

f.est <- read.csv("~/Data/Nurse_Larva/MKtestConstraintOneAlpha.csv")
colnames(f.est) = c("Gene","f")
ext <- read.csv("~/Downloads/msx123_Supp (1)/MpharAnn.csv") #load in MBE results

makePlot <- function(nurse){
  df <- read.csv(paste("JAN23",nurse,"GenieTabConn.csv",sep=""))
  df$tissue = gsub("_.*","",df$Gene)
  df$Gene = gsub(".*_","",df$Gene)
  df = df[,-c(1)]
  return(df)
}

codes <- c("CH","CG")
plotsG <- lapply(codes,makePlot)

allD <- ldply(plotsG,data.frame)

allDf <- merge(allD,f.est,by = "Gene",all.x = TRUE)
allDf$reg_diff = allDf$reg_between - allDf$reg_within
allDf$targ_diff = allDf$targ_between - allDf$targ_within
allDf$tissueC = as.factor(apply(allDf[,c('tissue','code')],1,paste,collapse="_"))
levels(allDf$tissueC) = c("larva -> nurse abdomen","larva -> nurse head","nurse abdomen","nurse head")

f2b <- ggplot(allDf,aes (x = tissueC, y = reg_within/max(allDf$reg_within),fill = tissueC))+
  geom_boxplot(notch = TRUE)+
  ylab("within-tissue regulatory strength")+
  xlab("tissue")+
  theme_bw()+theme4+
  scale_fill_manual(values = tissue_palette,name = "tissue")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

f2c <- ggplot(allDf,aes (x = tissueC, y = reg_between/max(allDf$reg_within),fill = tissueC))+
  geom_boxplot(notch = TRUE)+
  ylab("social regulatory strength")+
  xlab("tissue")+
  theme_bw()+theme4+
  scale_fill_manual(values = tissue_palette,name = "tissue")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

f2d <- ggplot(allDf, aes(x = reg_within/max(allDf$reg_within), y = reg_between/max(allDf$reg_within), color = tissueC))+
  geom_point(alpha = 0.5)+
  geom_smooth(method = "lm",se = FALSE)+
  theme_bw()+theme4+
  ylab("social regulatory strength")+
  xlab("within-tissue regulatory strength")+
  theme(legend.position = c(0.8,0.8))+
  scale_color_manual(values = tissue_palette,name = "tissue")

png("~/Writing/Figures/NurseLarva/corApproach/fig2.png",height = 3000,width = 4000,res =300)
grid.arrange(f2b,f2c,f2d,nrow=1)
dev.off()

allDN = droplevels(allDf[!grepl("larva",allDf$tissueC),])
levels(allDN$tissueC) = c("nurse abdomen","nurse head")

nurse_pallete=c("red","blue")

f3a <- ggplot(allDN,aes( x = reg_within/max(allDN$reg_within), y = 1 - log(f), color = tissueC))+
  geom_point(alpha = 0.5)+
  geom_smooth(method = "lm",se=FALSE)+
  theme_bw()+
  ylab("constraint")+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  theme4+
  xlab("within-tissue regulatory strength")+
  scale_color_manual(values = nurse_pallete,name = "tissue")+
  theme(legend.position = "none")

f3b <- ggplot(allDN,aes( x = reg_between/max(allDN$reg_between), y = 1 - log(f), color = tissueC))+
  geom_point(alpha = 0.5)+
  geom_smooth(method = "lm",se=FALSE)+
  theme_bw()+
  ylab("constraint")+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  theme4+
  xlab("social regulatory strength")+
  scale_color_manual(values = nurse_pallete,name = "tissue")+
  theme(legend.position = "none")



lmStat <- function(tissue){
  a1 <- summary(glm(1 - log(f) ~ reg_within, data = droplevels(allDf[allDf$tissueC == tissue,])))$coefficients[2,3:4]
  a2 <- summary(glm(1 - log(f) ~ reg_between, data = droplevels(allDf[allDf$tissueC == tissue,])))$coefficients[2,3:4]
  a3 <- summary(glm(1 - log(f) ~ targ_within, data = droplevels(allDf[allDf$tissueC == tissue,])))$coefficients[2,3:4]
  a4 <- summary(glm(1 - log(f) ~ targ_between, data = droplevels(allDf[allDf$tissueC == tissue,])))$coefficients[2,3:4]
  a5 <- summary(glm(1 - log(f) ~ reg_diff, data = droplevels(allDf[allDf$tissueC == tissue,])))$coefficients[2,3:4]
  a6 <- summary(glm(1 - log(f) ~ targ_diff, data = droplevels(allDf[allDf$tissueC == tissue,])))$coefficients[2,3:4]
  return(rbind(a1,a2,a3,a4,a5,a6))
}

lapply(levels(allDf$tissueC),lmStat)

all <- merge(allDN,ext,by="Gene",all.x=TRUE)
allP = all[!is.na(all$Raw.PS),]
allP$PS2 = factor(allP$PS2,levels = c("cellular",'eukaryote','bilaterian','insect','hymenopteran_ant'))
allP$PS2[allP$Raw.PS=="Neoptera"]="insect"

bootCI <- function(v,boots){
  v = v[!is.na(v)]
  bV = unlist(lapply(seq(1,boots),function(x){
    mean(v[sample(x = seq(1,length(v)),size = length(v),replace=TRUE)])
  }))
  c1 = quantile(bV, 0.025)
  c2 = quantile(bV, 0.975)
  m = mean(v)
  return(data.frame(mean = m, c1 = c1, c2 = c2))
}

fM = lapply(levels(allP$PS2), function(x){
  cbind(bootCI(1-allP$f[allP$PS2 == x],1000),bootCI(allP$reg_diff[allP$PS2 == x],1000))
})

fM = ldply(fM,data.frame)
colnames(fM)[c(1,4)] = c("f","reg_diff")
fM$phylostrata = factor(levels(allP$PS2),levels = levels(allP$PS2))

f3c <- ggplot(fM,aes(x=reg_diff,y=f,color = phylostrata))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymin = c1, ymax = c2))+
  geom_errorbarh(aes(xmin = c1.1, xmax = c2.1))+
  theme_bw()+
  scale_color_manual(values = PS_palette)+
  ylab("constraint")+
  xlab("social minus within-tissue regulatory strength")+theme4

png("~/Writing/Figures/NurseLarva/corApproach/fig3.png",height = 4000,width = 4000,res =300)
grid.arrange(f3a,f3b,f3c,nrow = 2)
dev.off()

lm <- glm(1 - log(f) ~ tissueC + reg_between + reg_within + PS2, data = allP[!grepl("larva",allP$tissueC),])
av <- aov(lm)
summary(av)
TukeyHSD(av,"PS2")



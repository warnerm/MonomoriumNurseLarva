cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
library(scales)
library(plyr)
library(reshape2)
library(gridExtra)
library(grid)
library(pBrackets)
library(ggplot2)

PS_palette = c("white","#deebf7","#9ecae1","#4292c6","#2171b5","#084594")
Soc_palette = c("#dadaeb","#9e9ac8","#6a51a3")

theme_all = theme(axis.text=element_text(size=13),
                  axis.title=element_text(size=17),
                  legend.text=element_text(size=13),
                  legend.title=element_text(size=17),
                  plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
                  axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)),
                  axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)),
                  legend.position = c(.05, .98),
                  legend.justification = c("left", "top"),
                  legend.box.just = "left",
                  legend.margin = margin(6, 6, 6, 6),
                  legend.background=element_rect(fill=FALSE))

bootMean <- function(data,boots){
  res = c()
  for (i in 1:boots){
    d = sample(data,length(data),replace=TRUE)
    res = c(res,mean(d,na.rm=TRUE))
  }
  return(res)
}

compareNgene <- function(nGene){
  df <- read.csv(paste("~/Data/Nurse_Larva/workerQR_1000000_",nGene,"_connStrengths.csv",sep=""))
  dfR <- read.csv(paste("~/Data/Nurse_Larva/random_1000000_",nGene,"_connStrengths.csv",sep=""))
  df = df[,-c(1:2)]
  dfR = dfR[,-c(1,2)]
  colnames(df)[2:10]=colnames(dfR)[2:10]=c("L_L","H_L","G_L","L_H","H_H","G_H","L_G","H_G","G_G")
  colnames(dfR)[2:10]=paste(colnames(dfR)[2:10],"random",sep="")
  
  df[,2:10]=df[,2:10]/max(df[,2:10])
  dfR[,2:10]=dfR[,2:10]/max(dfR[,2:10])
  dfB = merge(df,dfR,by="gene",all=TRUE)
  dfM = melt(dfB,id.vars="gene")
  dfM$cond = "focal"
  dfM$cond[grepl("random",dfM$variable)]="random"
  dfM$variable=gsub("random","",dfM$variable)
  dfM = droplevels(dfM)
  dfM$variable=factor(dfM$variable,levels=c("L_L","H_H","G_G","H_G","G_H","L_H","L_G","H_L","G_L"))
  levels(dfM$variable)=c("larva-larva","head-head","abd-abd",
                         "head-abd","abd-head","larva-head",
                         "larva-abd","head-larva","abd-larva")
  
  
  p <- ggplot(dfM,aes(x=variable,y=value,fill=cond))+
    geom_boxplot(notch=TRUE)+theme_bw()+theme(legend.position=c(0.8,0.95))+
    scale_fill_manual(name="nurse type",values=c("dodgerblue2","firebrick3"))+
    ylab("relative connection strength")+
    geom_vline(xintercept=3.5,linetype=2)+
    geom_vline(xintercept=5.5,linetype=2)+
    annotate("text",x=4.5,y=0.9,lineheight=0.9,vjust=0.75,hjust=0.5,label='atop(bold("between\n tissue"))',size=6,parse=TRUE)+
    annotate("text",x=7.5,y=0.9,label='atop(bold("social"))',size=6,parse=TRUE)+
    annotate("text",x=2,y=0.9,lineheight=0.9,vjust=0.75,label='atop(bold("within\ntissue"))',size=6,parse=TRUE)+
    xlab("connection type")+theme(legend.position="top",legend.justification = c(0.5,0))+
    theme(axis.title.x=element_text(margin=margin(15,0,0,0)))+
    theme(axis.title.y=element_text(margin=margin(0,15,0,0)))+
    theme(axis.text.x=element_text(angle=-45,vjust=1,hjust=0))+
    ggtitle(nGene)
  
    df$SIL_WH=df[,5] - df[,2]
    df$SIL_WG=df[,8] - df[,2]
    df$SIWH=df[,3]- df[,6]
    df$SIWG=df[,4] - df[,10]
    df$SIO = ((df$SIL_WH+df$SIL_WG)/2+df$SIWH+df$SIWG)/3
    df$SIL = (df[,5] + df[,8])/2 - df[,2]
    
    f.est <- read.csv("~/Data/MKtestConstraintOneAlpha.csv")
    colnames(f.est) = c("gene","f.est")
    dfAll <- merge(df,f.est,by="gene")
    
    return(list(p,dfAll))
}

tesCor <- function(dfAll){
  print(cor.test(dfAll$f.est,dfAll$SIL_WH))
  print(cor.test(dfAll$f.est,dfAll$SIL_WG))
  print(cor.test(dfAll$f.est,dfAll$SIWG))
  print(cor.test(dfAll$f.est,dfAll$SIWH))
}

tesCorConn <- function(dfAll){
  print(cor.test(dfAll$f.est,dfAll$L_L))
  print(cor.test(dfAll$f.est,dfAll$L_H))
  print(cor.test(dfAll$f.est,dfAll$L_G))
  print(cor.test(dfAll$f.est,dfAll$H_L))
  print(cor.test(dfAll$f.est,dfAll$H_G))
  print(cor.test(dfAll$f.est,dfAll$H_H))
  print(cor.test(dfAll$f.est,dfAll$G_L))
  print(cor.test(dfAll$f.est,dfAll$G_G))
  print(cor.test(dfAll$f.est,dfAll$G_H))
}

g500 = compareNgene(500)
g1000 = compareNgene(1000)
g2000 = compareNgene(2000)
g5000 = compareNgene(5000)
tesCor(g500[[2]])
tesCorConn(g500[[2]])

lm <- glm(log(f.est )~ L_H + L_G + H_H + H_G + H_L + G_L + G_H + G_G,data=g1000[[2]])
drop1(lm,.~.,test="Chi") 

lm <- glm(log(f.est) ~ SIWH + SIWG + SIL_WH + SIL_WG, data = g1000[[2]])

df <- read.csv("~/Data/Nurse_Larva/workerQR_1000000_500_connStrengths.csv")
dfR <- read.csv("~/Data/Nurse_Larva/random_1000000_500_connStrengths.csv")
df = df[,-c(1:2)]
dfR = dfR[,-c(1,2)]
colnames(df)[2:10]=colnames(dfR)[2:10]=c("L_L","H_L","G_L","L_H","H_H","G_H","L_G","H_G","G_G")
colnames(dfR)[2:10]=paste(colnames(dfR)[2:10],"random",sep="")

df[,2:10]=df[,2:10]/max(df[,2:10])
dfR[,2:10]=dfR[,2:10]/max(dfR[,2:10])
dfB = merge(df,dfR,by="gene",all=TRUE)
dfM = melt(dfB,id.vars="gene")
dfM$cond = "focal"
dfM$cond[grepl("random",dfM$variable)]="random"
dfM$variable=gsub("random","",dfM$variable)
dfM = droplevels(dfM)
dfM$variable=factor(dfM$variable,levels=c("L_L","H_H","G_G","H_G","G_H","L_H","L_G","H_L","G_L"))
levels(dfM$variable)=c("larva-larva","head-head","abd-abd",
                          "head-abd","abd-head","larva-head",
                          "larva-abd","head-larva","abd-larva")


p <- ggplot(dfM,aes(x=variable,y=value,fill=cond))+
  geom_boxplot(notch=TRUE)+theme_bw()+theme_all+theme(legend.position=c(0.8,0.95))+
  scale_fill_manual(name="nurse type",values=c("dodgerblue2","firebrick3"))+
  ylab("relative connection strength")+
  geom_vline(xintercept=3.5,linetype=2)+
  geom_vline(xintercept=5.5,linetype=2)+
  annotate("text",x=4.5,y=0.9,lineheight=0.9,vjust=0.75,hjust=0.5,label='atop(bold("between\n tissue"))',size=6,parse=TRUE)+
  annotate("text",x=7.5,y=0.9,label='atop(bold("social"))',size=6,parse=TRUE)+
  annotate("text",x=2,y=0.9,lineheight=0.9,vjust=0.75,label='atop(bold("within\ntissue"))',size=6,parse=TRUE)+
  xlab("connection type")+theme(legend.position="top",legend.justification = c(0.5,0))+
  theme(axis.title.x=element_text(margin=margin(15,0,0,0)))+
  theme(axis.title.y=element_text(margin=margin(0,15,0,0)))+
  theme(axis.text.x=element_text(angle=-45,vjust=1,hjust=0))

png("~/Writing/Figures/NurseLarva/Fig1b.png")
p
dev.off()

for (lev in levels(dfM$variable)){
  print(lev)
  print(wilcox.test(dfM$value[dfM$variable==lev&dfM$cond=="focal"],dfM$value[dfM$variable==lev&dfM$cond=="random"],paired=FALSE,alternative="greater"))
}

#################
##Continuing with only focal nurse data; using QR and QL together
#################
df <- read.csv("~/Data/Nurse_Larva/workerQR_1000000_5000_connStrengths.csv")
df = df[,-c(1:2)]
colnames(df)[2:10]=c("L_L","WH_L","WG_L","L_WH","WH_WH","WG_WH","L_WG","WH_WG","WG_WG")
df$SIL_WH=df[,5] - df[,2]
df$SIL_WG=df[,8] - df[,2]
df$SIWH=df[,3]- df[,6]
df$SIWG=df[,4] - df[,10]
df$SIO = ((df$SIL_WH+df$SIL_WG)/2+df$SIWH+df$SIWG)/3
df$SIL = (df[,5] + df[,8])/2 - df[,2]

##################
###Look at correlation with constraint
##################
f.est <- read.csv("~/Data/MKtestConstraintOneAlpha.csv")
colnames(f.est) = c("gene","f.est")
dfAll <- merge(df,f.est,by="gene")

p1 <- ggplot(dfAll,aes(x=SIWH,y=1-f.est))+
  geom_point(size=0.8)+
  theme_bw()+ylab("")+xlab("")+
  scale_x_continuous(breaks=c(-0.06,-0.03,0,0.03),limits=c(-0.08,0.04))+
  ggtitle("nurse head")+
  geom_smooth(method="lm",se=FALSE,color="red")+
  annotate("text",x=-0.045,y=0.22,label="r = - 0.096\np = 0.018",size=5)

p2 <- ggplot(dfAll,aes(x=SIWG,y=1-f.est))+
  geom_point(size=0.8)+
  theme_bw()+ylab("")+xlab("")+
  scale_x_continuous(breaks=c(-0.06,-0.03,0,0.03),limits=c(-0.08,0.04))+
  geom_smooth(method="lm",se=FALSE,color="red")+
  ggtitle("nurse abd")+
  annotate("text",x=-0.045,y=0.22,label="r = - 0.102\np = 0.011",size=5)

p4 <- ggplot(dfAll,aes(x=SIL_WH,y=1-f.est))+
  geom_point(size=0.8)+
  theme_bw()+
  scale_x_continuous(breaks=c(-0.09,-0.06,-0.03,0),limits=c(-0.1,0.02))+
  ylab("")+xlab("")+
  geom_smooth(method="lm",se=FALSE,color="red")+
  ggtitle("larva to head")+
  annotate("text",x=-0.06,y=0.22,label="r = - 0.055\np = 0.175",size=5)

p3 <- ggplot(dfAll,aes(x=SIL_WG,y=1-f.est))+
  geom_point(size=0.8)+
  theme_bw()+
  scale_x_continuous(breaks=c(-0.09,-0.06,-0.03,0),limits=c(-0.1,0.02))+ylab("")+xlab("")+
  geom_smooth(method="lm",se=FALSE,color="red")+
  ggtitle("larva to abd")+
  annotate("text",x=-0.06,y=0.22,label="r = - 0.113\np = 0.005",size=5)
theme_f2 = theme(plot.title=element_text(size=20,hjust=0.5),
                 plot.margin=margin(0.6,0,0,0),
                 axis.text=element_text(size=12),
                 axis.title=element_text(size=15),
                 legend.position = "none")
png("~/Writing/Figures/NurseLarva/Fig2.png",width=3000,height=1500,res=300)
grid.arrange(p1+theme_f2,
             p2+theme_f2,
             p3+theme_f2,
             p4+theme_f2,
             bottom = textGrob("sociality index", gp=gpar(fontsize=23,font=8),hjust=0.3,vjust=0),
             left = textGrob("constraint", rot = 90, vjust = 1,gp=gpar(fontsize=23,font=8)),nrow=1)
dev.off()

m = melt(dfAll[,c(1,11:14,17)],id.vars=c("gene","f.est"))

m$quantile="0-25"
for (i in 1:nrow(m)){
  if (m$value[i] > quantile(m$value[m$variable==m$variable[i]],0.25)){
    if (m$value[i] > quantile(m$value[m$variable==m$variable[i]],0.5)){
      if (m$value[i] > quantile(m$value[m$variable==m$variable[i]],0.75)){
        m$quantile[i]="75-100"
      } else {
        m$quantile[i]="50-75"
      }
    } else {
      m$quantile[i]="25-50"
    }
  }
}

ggplot(m,aes(x=variable,y=1-log(f.est),fill=quantile))+
  geom_boxplot(notch=TRUE)

pc <- prcomp()

###########
##Adding phylostrata
ext <- read.csv("~/Downloads/msx123_Supp (1)/MpharAnn.csv")
dfAll = merge(ext,df,by.x="Gene",by.y="gene")


dfAll$PS2 = factor(dfAll$PS2, levels = c("cellular","eukaryote","bilaterian","insect","hymenopteran_ant"))
dfAll = dfAll[!is.na(dfAll$PS2),]
levels(dfAll$PS2)[5]="hymenopteran"
#levels(dfAll$PS2)[4]="hymenopteran"
dfAll$PS2 = droplevels(dfAll$PS2)
dfAll$PS2=as.character(dfAll$PS2)
dfAll$PS2[dfAll$Raw.PS=="Monomorium pharaonis"] = "M. pharaonis"
dfAll$PS2 = factor(dfAll$PS2,levels = c("cellular","eukaryote","bilaterian","insect","hymenopteran","M. pharaonis"))

lm <- glm(SIWH ~ PS2, data=dfAll)
lm <- glm(SIWG ~ PS2, data=dfAll)

summary(glht(lm, mcp(PS2="Tukey"))) 

d2 <- melt(dfAll[,c(1,27,37:40)],id.vars=c("Gene","PS2"))
d2$tissue="larva-head"
d2$tissue[grepl("SIL_WG",d2$variable)]="larva-abd"
d2$tissue[grepl("SIWG",d2$variable)]="abdomen"
d2$tissue[grepl("SIWH",d2$variable)]="head"

rowN = c()
for (i in 1:6){
  rowN = c(rowN,paste(levels(dfAll$PS2)[i]," (",sum(dfAll$PS2==levels(dfAll$PS2)[i]),")",sep=""))
}

p1 <- ggplot(d2,aes(x=tissue,y=value,fill=PS2))+geom_boxplot(notch=TRUE)+
  theme_bw()+ylab("sociality index")+xlab("tissue")+ylim(-0.1,0.1)+
  scale_fill_manual(values=PS_palette,name = "phylostrata",labels=rowN)+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=25),
        legend.text=element_text(size=11),
        legend.title=element_text(size=14),
        legend.position = "right")

png("~/Writing/Figures/NurseLarva/Fig3a.png",width=2500,height=2000,res=300)
grid.arrange(p1 + theme(axis.text=element_text(size=17),
                        axis.title=element_text(size=22),
                        axis.title.x=element_text(margin=margin(13,0,0,0)),
                        legend.text=element_text(size=15),
                        legend.title=element_text(size=20),
                        legend.position = c(.05, .98),
                        legend.justification = c("left", "top"),
                        legend.box.just = "left",
                        legend.margin = margin(6, 6, 6, 6),
                        legend.background=element_rect(fill=FALSE),
                        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
                        axis.title.y=element_text(margin=margin(0,10,0,0))))
dev.off()

##Add f back in so dataframe has sociality, f and PS
dfAll <- merge(dfAll,f.est,by.x="Gene",by.y="gene")
d4 <- ddply(dfAll,~PS2,summarise,
            meanSI = mean(SIO),
            c1B = quantile(bootMean(SIO,1000),0.025),
            c2B = quantile(bootMean(SIO,1000),0.975),
            meanF = mean(f.est),
            c1W = quantile(bootMean(f.est,1000),0.025),
            c2W = quantile(bootMean(f.est,1000),0.975))

p2 <- ggplot(d4[!is.na(d4$PS2),],aes(x=meanSI,y=1-meanF,fill=PS2))+
  geom_point(aes(fill=PS2),pch=21,color="black",size=7)+
  geom_errorbarh(aes(xmin=c1B,xmax=c2B))+
  geom_errorbar(aes(ymin=1-c1W,ymax=1-c2W))+ylim(0.8,0.95)+
  theme_bw()+xlab("average sociality index")+ylab("constraint")+
  scale_fill_manual(values=PS_palette,name="phylostrata")+
  #scale_color_manual(values=PS_palette,name="phylostrata")+
  theme(legend.position="none")+
  annotate("text",x=-0.0325,y=0.93,label="eukaryote",size=6)+
  annotate("text",x=-0.034,y=0.879,label="bilaterian",size=6)+
  annotate("text",x=-0.0295,y=0.908,label="cellular",size=6)+
  annotate("text",x=-0.029,y=0.845,label="hymenopteran",size=6)

png("~/Writing/Figures/NurseLarva/Fig3b.png",width=2000,height=2000,res=300)
grid.arrange(p2+theme(plot.margin=unit(c(0.5,1,0.5,0.5),"cm"),
                      axis.text=element_text(size=17),
                      axis.title=element_text(size=22),
                      axis.title.x=element_text(margin=margin(10,0,0,0)),
                      axis.title.y=element_text(margin=margin(0,10,0,0))))
dev.off()

png("~/Writing/Figures/NurseLarva/Fig3.png",width=4000,height=2000,res=300)
grid.arrange(p1 + theme(axis.text=element_text(size=15),
                        axis.title=element_text(size=20),
                        axis.title.x=element_text(margin=margin(13,0,0,0)),
                        legend.text=element_text(size=13),
                        legend.title=element_text(size=17),
                        legend.position = c(.05, .98),
                        legend.justification = c("left", "top"),
                        legend.box.just = "left",
                        legend.margin = margin(6, 6, 6, 6),
                        legend.background=element_rect(fill=FALSE),
                        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
                        axis.title.y=element_text(margin=margin(0,10,0,0)))+
               annotate("text",x=3,y=0.08,label="A",size=11),
             p2+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
                      axis.text=element_text(size=15),
                      axis.title=element_text(size=20),
                      axis.title.x=element_text(margin=margin(10,0,0,0)),
                      axis.title.y=element_text(margin=margin(0,10,0,0)))+
               annotate("text",x=-0.04,y=0.95 - 0.02*(3/5),label="B",size=11),nrow=1)
dev.off()



p1 <- ggplot(dfAll,aes(x=SIWH,y=BSnIPRE.est))+
  geom_point(size=0.8)+
  theme_bw()+ylab("")+xlab("")+
  scale_x_continuous(breaks=c(-0.06,-0.03,0,0.03),limits=c(-0.08,0.04))+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),
        legend.position = "none")+
  ggtitle("a. nurse head")+
  geom_smooth(method="lm",se=FALSE,color="red")+
  annotate("text",x=-0.045,y=0.22,label="r = - 0.096\np = 0.018",size=5)

p2 <- ggplot(dfAll,aes(x=SIWG,y=BSnIPRE.est))+
  geom_point(size=0.8)+
  theme_bw()+ylab("")+xlab("")+
  scale_x_continuous(breaks=c(-0.06,-0.03,0,0.03),limits=c(-0.08,0.04))+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),legend.position = "none")+
  geom_smooth(method="lm",se=FALSE,color="red")+
  ggtitle("b. nurse abdomen")+
  annotate("text",x=-0.045,y=0.22,label="r = - 0.102\np = 0.011",size=5)

p3 <- ggplot(dfAll,aes(x=SIL_WH,y=BSnIPRE.est))+
  geom_point(size=0.8)+
  theme_bw()+
  scale_x_continuous(breaks=c(-0.09,-0.06,-0.03,0))+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),
        legend.position = "none")+ylab("")+xlab("")+
  geom_smooth(method="lm",se=FALSE,color="red")+
  ggtitle("c. larva")+
  annotate("text",x=-0.06,y=0.22,label="r = - 0.091\np = 0.025",size=5)

p4 <- ggplot(dfAll,aes(x=SIL_WG,y=BSnIPRE.est))+
  geom_point(size=0.8)+
  theme_bw()+
  scale_x_continuous(breaks=c(-0.09,-0.06,-0.03,0))+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),
        legend.position = "none")+ylab("")+xlab("")+
  geom_smooth(method="lm",se=FALSE,color="red")+
  ggtitle("c. larva")+
  annotate("text",x=-0.06,y=0.22,label="r = - 0.091\np = 0.025",size=5)

png("~/Writing/Figures/NurseLarva/FigS_alpha.png",width=3000,height=1500,res=300)
grid.arrange(p1+theme(plot.title=element_text(size=20,hjust=0.5)),
             p2+theme(plot.title=element_text(size=20,hjust=0.5)),
             p3+theme(plot.title=element_text(size=20,hjust=0.5)),
             p4+theme(plot.title=element_text(size=20,hjust=0.5)),
             bottom = textGrob("sociality index", gp=gpar(fontsize=23,font=8),hjust=0.3,vjust=-0.5),
             left = textGrob("constraint", rot = 90, vjust = 1.5,gp=gpar(fontsize=23,font=8)),nrow=1)
dev.off()



m = melt(dfAll[,c(1,19,37:40)],id.vars=c("Gene","BSnIPRE.est"))

m$quantile="0-25"
for (i in 1:nrow(m)){
  if (m$value[i] > quantile(m$value[m$variable==m$variable[i]],0.25)){
    if (m$value[i] > quantile(m$value[m$variable==m$variable[i]],0.5)){
      if (m$value[i] > quantile(m$value[m$variable==m$variable[i]],0.75)){
        m$quantile[i]="75-100"
      } else {
        m$quantile[i]="50-75"
      }
    } else {
      m$quantile[i]="25-50"
    }
  }
}

ggplot(m,aes(x=variable,y=BSnIPRE.est,fill=quantile))+
  geom_boxplot(notch=TRUE)



m = melt(dfAll[,c(1,19,27,37:40,43)],id.vars=c("Gene","BSnIPRE.est","f.est","PS2"))

m$quantile="0-25"
for (i in 1:nrow(m)){
  if (m$value[i] > quantile(m$value[m$variable==m$variable[i]],0.25)){
    if (m$value[i] > quantile(m$value[m$variable==m$variable[i]],0.5)){
      if (m$value[i] > quantile(m$value[m$variable==m$variable[i]],0.75)){
        m$quantile[i]="75-100"
      } else {
        m$quantile[i]="50-75"
      }
    } else {
      m$quantile[i]="25-50"
    }
  }
}

ggplot(m[m$variable=="SIWG",],
       aes(x=BSnIPRE.est,y=f.est,color=quantile,shape=PS2))+
  geom_point()

d2 <- melt(df[,c(1,3:11)],id.vars="gene")
d2$tissue="larva"
d2$tissue[grepl("^WG",d2$variable)]="nurse abdomen"
d2$tissue[grepl("^WH",d2$variable)]="nurse head"
d2$connection.type = "within tissue"
d2$connection.type[grepl("between",d2$variable)]="social"
d2$connection.type[grepl("WH.WG",d2$variable)|grepl("WG.WH",d2$variable)]="between tissue"
d2$connection.type=factor(d2$connection.type,levels=c("within tissue","social","between tissue"))

p1 <- ggplot(d2,aes(x=connection.type,y=value,fill=tissue))+
  geom_boxplot(notch=TRUE)+
  theme_bw()+xlab("connection type")+ylab("connection strength")+
  scale_fill_manual(values=Soc_palette,name="tissue")+
  theme(axis.text=element_text(size=21),
        axis.title=element_text(size=25),
        axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)),
        axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)),
        legend.text=element_text(size=20),
        legend.title=element_text(size=23),
        legend.position = c(.05, .98),
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.margin = margin(6, 6, 6, 6),
        legend.background=element_rect(fill=FALSE))+
  ylim(0,0.25)+
  guides(fill=guide_legend(
    keywidth=0.2,
    keyheight=0.30,
    default.unit="inch")
  )
png("~/Writing/Figures/NurseLarva/Fig1.png",width=3000,height=3000,res=300)
p1
dev.off()

p2 <- ggplot(df,aes(x=WHwithin,y=WHbetween))+geom_point(size=0.7)+
  theme_bw()+xlab("")+geom_smooth(method="lm",se=FALSE,color="red")+
  ylab("")+geom_abline(slope=1,linetype="dashed")+xlim(0,0.25)+ylim(0,0.25)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))+
  ggtitle("A. nurse head")+theme(plot.title=element_text(hjust=0.5))+
  annotate("text",x=0.225,y=0.025,label="n.s.",size=5)

p3 <- ggplot(df,aes(x=WGwithin,y=WGbetween))+geom_point(size=0.7)+
  theme_bw()+xlab("")+geom_smooth(method="lm",se=FALSE,color="red")+
  ylab("")+geom_abline(slope=1,linetype="dashed")+xlim(0,0.25)+ylim(0,0.25)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))+
  ggtitle("B. nurse gaster")+theme(plot.title=element_text(hjust=0.5))+
  annotate("text",x=0.225,y=0.025,label="***",size=5)


p4 <- ggplot(df,aes(x=Lwithin,y=Lbetween))+geom_point(size=0.7)+
  theme_bw()+xlab("")+geom_smooth(method="lm",se=FALSE,color="red")+
  ylab("")+geom_abline(slope=1,linetype="dashed")+xlim(0,0.25)+ylim(0,0.25)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))+
  ggtitle("C. larva")+theme(plot.title=element_text(hjust=0.5))+
  annotate("text",x=0.225,y=0.025,label="***",size=5)

###Add histograms as insets
p1s <- ggplot(df,aes(x=SIL))+geom_histogram(binwidth=0.005)+
  xlim(-0.15,0.15)+ylim(0,100)+ggtitle("sociality index")+
  xlab("")+ylab("")+
  theme(plot.title = element_text(hjust = 0.5,size=10,vjust=0.5),
        plot.margin=unit(c(0,0,-0.5,-0.5),"cm"))
p2s <- ggplot(df,aes(x=SIWH))+geom_histogram(binwidth=0.005)+
  xlim(-0.15,0.15)+ylim(0,100)+ggtitle("sociality index")+
  xlab("")+ylab("")+
  theme(plot.title = element_text(hjust = 0.5,size=10,vjust=0.5),
        plot.margin=unit(c(0,0,-0.5,-0.5),"cm"))

p3s <- ggplot(df,aes(x=SIWG))+geom_histogram(binwidth=0.005)+
  xlim(-0.15,0.15)+ylim(0,100)+ggtitle("sociality index")+
  xlab("")+ylab("")+
  theme(plot.title = element_text(hjust = 0.5,size=10,vjust=0.5),
        plot.margin=unit(c(0,0,-0.5,-0.5),"cm"))
png("~/Writing/Figures/NurseLarva/FigS1.png",width=3000,height=1000,res=300)
grid.arrange(p2+theme_all+theme(axis.text=element_text(size=11),
                                plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))+
               annotation_custom(grob=ggplotGrob(p1s),xmin=0.11,xmax=0.25,ymin=0.11,ymax=0.25),
             p3+theme_all+theme(axis.text=element_text(size=11),
                                plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))+
               annotation_custom(grob=ggplotGrob(p2s),xmin=0.11,xmax=0.25,ymin=0.11,ymax=0.25),
             p4+theme_all+theme(axis.text=element_text(size=11),
                                plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"))+
               annotation_custom(grob=ggplotGrob(p3s),xmin=0.11,xmax=0.25,ymin=0.11,ymax=0.25),
             bottom = textGrob("within tissue connection strength", gp=gpar(fontsize=15,font=8),hjust=0.4,vjust=-1.5),
             left = textGrob("social connection strength", rot = 90, vjust = 2.5,hjust=0.4,gp=gpar(fontsize=15,font=8)),nrow=1)
dev.off()

p1 <- ggplot(df,aes(y=WGwithin,x=WG.WH))+
  geom_point(size=0.8)+xlim(0.025,0.12)+ylim(0.025,0.15)+
  ylab("within tissue")+xlab("between tissue")+
  geom_smooth(method="lm",se=FALSE,color="red")+
  annotate("text",x=0.03,y=0.14,label="C",size=10)+
  annotate("text",x=0.1,y=0.12,label="***",size=8)+
  annotate("text",x=0.1,y=0.03,label="nurse gaster",size=8)

p2 <- ggplot(df,aes(x=WG.WH,y=WGbetween))+
  geom_point(size=0.8)+xlim(0.025,0.12)+ylim(0.025,0.15)+
  ylab("social")+xlab("between tissue")+
  geom_smooth(method="lm",se=FALSE,color="red")+
  annotate("text",x=0.03,y=0.14,label="A",size=10)+
  annotate("text",x=0.1,y=0.12,label="***",size=8)+
  annotate("text",x=0.1,y=0.03,label="nurse gaster",size=8)


p3 <- ggplot(df,aes(y=WHwithin,x=WH.WG))+
  geom_point(size=0.8)+xlim(0.025,0.12)+ylim(0.025,0.15)+
  ylab("within tissue")+xlab("between tissue")+
  geom_smooth(method="lm",se=FALSE,color="red")+
  annotate("text",x=0.03,y=0.14,label="D",size=10)+
  annotate("text",x=0.1,y=0.12,label="***",size=8)+
  annotate("text",x=0.1,y=0.03,label="nurse head",size=8)


p4 <- ggplot(df,aes(x=WH.WG,y=WHbetween))+
  geom_point(size=0.8)+xlim(0.025,0.12)+ylim(0.025,0.15)+
  ylab("social")+xlab("between tissue")+
  geom_smooth(method="lm",se=FALSE,color="red")+
  annotate("text",x=0.03,y=0.14,label="B",size=10)+
  annotate("text",x=0.1,y=0.12,label="***",size=8)+
  annotate("text",x=0.1,y=0.03,label="nurse head",size=8)


png("~/Writing/Figures/NurseLarva/FigS2.png",width=3000,height=3000,res=300)
grid.arrange(p2+theme_all,
             p4+theme_all,
             p1+theme_all,
             p3+theme_all,nrow=2,ncol=2)
dev.off()


p1 <- ggplot(df,aes(x=SIWH,y=SIWG))+
  geom_point(size=0.8)+xlim(-0.1,0.05)+ylim(-0.15,0.05)+
  ylab("nurse gaster")+xlab("nurse head")+
  geom_smooth(method="lm",se=FALSE,color="red")+
  annotate("text",x=0.03,y=-0.13,label="***",size=7)+
  annotate("text",x=-0.08,y=0.04,label="A",size=9)


p2 <- ggplot(df,aes(x=SIWH,y=SIL))+
  geom_point(size=0.8)+xlim(-0.1,0.05)+ylim(-0.15,0.05)+
  ylab("larva")+xlab("nurse head")+
  geom_smooth(method="lm",se=FALSE,color="red")+
  annotate("text",x=0.03,y=-0.13,label="n.s.",size=7)+
  annotate("text",x=-0.08,y=0.04,label="B",size=9)

p3 <- ggplot(df,aes(x=SIWG,y=SIL))+
  geom_point(size=0.8)+xlim(-0.1,0.05)+ylim(-0.15,0.05)+
  ylab("larva")+xlab("nurse gaster")+
  geom_smooth(method="lm",se=FALSE,color="red")+
  annotate("text",x=0.03,y=-0.13,label="**",size=7)+
  annotate("text",x=-0.08,y=0.04,label="C",size=9)

png("~/Writing/Figures/NurseLarva/FigS3.png",width=3000,height=1000,res=300)
grid.arrange(p1+theme_all,
             p2+theme_all,
             p3+theme_all,nrow=1)
dev.off()

###########
##Adding phylostrata
ext <- read.csv("~/Downloads/msx123_Supp (1)/MpharAnn.csv")
dfAll = merge(ext,df,by.x="Gene",by.y="gene")

dfAll$PS2 = factor(dfAll$PS2, levels = c("cellular","eukaryote","bilaterian","insect","hymenopteran_ant"))
dfAll = dfAll[!is.na(dfAll$PS2),]
levels(dfAll$PS2)[5]="hymenopteran"

lm <- glm(SIWH ~ PS2, data=dfAll)
lm <- glm(SIWG ~ PS2, data=dfAll)

summary(glht(lm, mcp(PS2="Tukey"))) 

d2 <- melt(dfAll[,c(1,27,37:39)],id.vars=c("Gene","PS2"))
d2$tissue="larva"
d2$tissue[grepl("WG",d2$variable)]="nurse gaster"
d2$tissue[grepl("WH",d2$variable)]="nurse head"

rowN = c()
for (i in 1:5){
  rowN = c(rowN,paste(levels(dfAll$PS2)[i]," (",sum(dfAll$PS2==levels(dfAll$PS2)[i]),")",sep=""))
}

p1 <- ggplot(d2,aes(x=tissue,y=value,fill=PS2))+geom_boxplot(notch=TRUE)+
  theme_bw()+ylab("sociality index")+xlab("tissue")+ylim(-0.15,0.1)+
  scale_fill_manual(values=PS_palette,name = "phylostrata",labels=rowN)+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=25),
        legend.text=element_text(size=11),
        legend.title=element_text(size=14),
        legend.position = "right")

png("~/Writing/Figures/NurseLarva/Fig3a.png",width=2500,height=2000,res=300)
grid.arrange(p1 + theme(axis.text=element_text(size=17),
                        axis.title=element_text(size=22),
                        axis.title.x=element_text(margin=margin(13,0,0,0)),
                        legend.text=element_text(size=15),
                        legend.title=element_text(size=20),
                        legend.position = c(.05, .98),
                        legend.justification = c("left", "top"),
                        legend.box.just = "left",
                        legend.margin = margin(6, 6, 6, 6),
                        legend.background=element_rect(fill=FALSE),
                        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
                        axis.title.y=element_text(margin=margin(0,10,0,0))))
dev.off()

##Add f back in so dataframe has sociality, f and PS
dfAll <- merge(dfAll,f.est,by.x="Gene",by.y="gene")

d4 <- ddply(dfAll,~PS2,summarise,
            meanSI = mean(SI_Overall),
            c1B = quantile(bootMean(SI_Overall,1000),0.025),
            c2B = quantile(bootMean(SI_Overall,1000),0.975),
            meanF = mean(f.est),
            c1W = quantile(bootMean(f.est,1000),0.025),
            c2W = quantile(bootMean(f.est,1000),0.975))

p2 <- ggplot(d4[!is.na(d4$PS2),],aes(x=meanSI,y=1-meanF,fill=PS2))+
  geom_point(aes(fill=PS2),pch=21,color="black",size=7)+
  geom_errorbarh(aes(xmin=c1B,xmax=c2B))+
  geom_errorbar(aes(ymin=1-c1W,ymax=1-c2W))+ylim(0.8,0.95)+
  theme_bw()+xlab("average sociality index")+ylab("constraint")+
  scale_fill_manual(values=PS_palette,name="phylostrata")+
  #scale_color_manual(values=PS_palette,name="phylostrata")+
  theme(legend.position="none")+
  annotate("text",x=-0.0545,y=0.925,label="eukaryote",size=6)+
  annotate("text",x=-0.0578,y=0.879,label="bilaterian",size=6)+
  annotate("text",x=-0.049,y=0.908,label="cellular",size=6)+
  annotate("text",x=-0.046,y=0.87,label="insect",size=6)+
  annotate("text",x=-0.039,y=0.845,label="hymenopteran",size=6)

png("~/Writing/Figures/NurseLarva/Fig3b.png",width=2000,height=2000,res=300)
grid.arrange(p2+theme(plot.margin=unit(c(0.5,1,0.5,0.5),"cm"),
         axis.text=element_text(size=17),
         axis.title=element_text(size=22),
         axis.title.x=element_text(margin=margin(10,0,0,0)),
         axis.title.y=element_text(margin=margin(0,10,0,0))))
dev.off()

png("~/Writing/Figures/NurseLarva/Fig3.png",width=4000,height=2000,res=300)
grid.arrange(p1 + theme(axis.text=element_text(size=15),
                        axis.title=element_text(size=20),
                        axis.title.x=element_text(margin=margin(13,0,0,0)),
                        legend.text=element_text(size=13),
                        legend.title=element_text(size=17),
                        legend.position = c(.05, .98),
                        legend.justification = c("left", "top"),
                        legend.box.just = "left",
                        legend.margin = margin(6, 6, 6, 6),
                        legend.background=element_rect(fill=FALSE),
                        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
                        axis.title.y=element_text(margin=margin(0,10,0,0)))+
               annotate("text",x=3,y=0.08,label="A",size=11),
             p2+theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
                      axis.text=element_text(size=15),
                      axis.title=element_text(size=20),
                      axis.title.x=element_text(margin=margin(10,0,0,0)),
                      axis.title.y=element_text(margin=margin(0,10,0,0)))+
               annotate("text",x=-0.04,y=0.95 - 0.02*(3/5),label="B",size=11),nrow=1)
dev.off()

####
##Fig 4: relationship between f and d
dfAll = merge(ext,df,by.x="Gene",by.y="gene")
dfAll$PS2 = factor(dfAll$PS2, levels = c("cellular","eukaryote","bilaterian","insect","hymenopteran_ant"))
dfAll = dfAll[!is.na(dfAll$PS2),]
levels(dfAll$PS2)[5]="hymenopteran"
dfAll <- merge(dfAll,f.est,by.x="Gene",by.y="gene")
corMat = pMat = matrix(nrow=6,ncol=10)
cols = c(29:34,37:40)

for (j in 1:10){
  for (i in 1:5){
    d2 = dfAll[dfAll$PS2==levels(dfAll$PS2)[i],]
    res= cor.test((1-d2$f.est),d2[,cols[j]])
    corMat[i,j] = res$estimate
    pMat[i,j] = res$p.value
  }
  res = cor.test((1-dfAll$f.est),dfAll[,cols[j]])
  corMat[6,j] = res$estimate
  pMat[6,j] = res$p.value
}

pMat = pMat[,c(5,1,3,6,2,4,7,9,8,10)]
corMat = corMat[,c(5,1,3,6,2,4,7,9,8,10)]
colN = c("larva","nurse gaster","nurse head",
         "larva","nurse gaster","nurse head",
         "larva","nurse gaster","nurse head","overall")
rownames(pMat) = rownames(corMat) = c(levels(dfAll$PS2),"overall")

##Quick heatmap
pMat = t(pMat)
corMat = t(corMat)
corM = melt(corMat)
corP = melt(pMat)
corM$p.value = corP$value
corM$logp = -log10(corM$p.value)
corM$logp[corM$value<0]=-corM$logp[corM$value<0]
corM$text = apply(round(corM[,c(3,4)],3),1,paste,collapse="\n")
corM$text[corM$p.value < 0.001]=gsub("\n0","\n<0.001",corM$text[corM$p.value < 0.001])

rowN = c()
for (i in 1:5){
  rowN = c(rowN,paste(colnames(pMat)[i]," (",sum(dfAll$PS2==colnames(pMat)[i]),")",sep=""))
}
rowN = c(rowN,paste("overall (",nrow(dfAll),")",sep=""))

p <- ggplot(corM,aes(x=Var2,y=factor(Var1)))+geom_tile(aes(fill=logp))+
  geom_text(aes(label=text))+
  scale_fill_gradient2(low="blue",high="red",mid="white")+
  xlab("")+ylab("")+
  guides(fill= guide_colorbar(barheight=2,
                              barwidth=30,nbin=40,
                              title=NULL,ticks=FALSE,label=FALSE))+
  theme_bw()+
  geom_hline(yintercept=3.5,color="black",size=2.5)+
  geom_hline(yintercept=6.5,color="black",size=2.5)+
  geom_vline(xintercept=5.5,color="black",size=1.5,linetype="dashed")+
  scale_y_discrete(labels=colN,expand=c(0,0))+
  scale_x_discrete(labels=rowN,expand=c(0,0))+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.position = "top",
        legend.text=element_text(size=15),
        legend.title=element_text(size=18),
        axis.text.x=element_text(angle = 45,hjust=1))

xax = 80

png("~/Writing/Figures/NurseLarva/FigS4.png",width=3000,height=3000,res=300)
p+theme(plot.margin=unit(c(2,0.5,0.5,2),"cm"),axis.text.x=element_text(margin=margin(10,0,0,0)))
grid.brackets(xax,295,xax,120,lwd=2,col="red")
grid.brackets(xax,435,xax,305,lwd=2,col="red")
grid.brackets(xax,575,xax,445,lwd=2,col="red")
grid.text("sociality index",x=0.035,y=0.75,rot=90,gp=gpar(fontsize=18,fontface="italic"))
grid.text("within tissue",x=0.035,y=0.52,rot=90,gp=gpar(fontsize=18,fontface="italic"))
grid.text("social",x=0.035,y=0.3,rot=90,gp=gpar(fontsize=18,fontface="italic"))
grid.text("correlation with constraint",x=0.6,y=0.97,gp=gpar(fontsize=18,fontface="bold"))
grid.text("negative",x=0.355,y=0.93,gp=gpar(fontsize=18))
grid.text("positive",x=0.865,y=0.93,gp=gpar(fontsize=18))
grid.text("phyostrata (number of genes)",x=0.55,y=0.025,gp=gpar(fontsize=20))
dev.off()

######
##Fig 5: add in queen expression
load("~/Dropbox/monomorium nurses/data.processed/cleandata.RData")
fpkm = log(fpkm + sqrt(fpkm ^ 2 + 1)) #hyperbolic sine transformation to normalize gene expression data

meanQH = rowSums(fpkm[,grepl("QH",colnames(fpkm))])
meanQG = rowSums(fpkm[,grepl("QG",colnames(fpkm))])

expr = data.frame(QH=meanQH,QG=meanQG,Gene=names(meanQH))

dfAll = merge(ext,df,by.x="Gene",by.y="gene")
dfAll$PS2 = factor(dfAll$PS2, levels = c("cellular","eukaryote","bilaterian","insect","hymenopteran_ant"))
dfAll = dfAll[!is.na(dfAll$PS2),]
levels(dfAll$PS2)[5]="hymenopteran"
d = merge(dfAll,expr,by="Gene")

p1 <- ggplot(d,aes(x=QH,y=SIWH,color=PS2))+
  geom_point(aes(fill=PS2),pch=21,color="black",size=1.5)+
  theme_bw()+
  annotate("text",x=40,y=0.05,label="***",size=7)+
  ylab("nurse head sociality index")+xlab("queen head normalized expression")+
  scale_fill_manual(values=PS_palette,name="phylostrata")+
  ylim(-0.12,0.06)+
  xlim(0,40)+
  geom_smooth(method="lm",se=FALSE,color="red")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        legend.position = "none",
        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))+
  annotate("text",x=5,y=0.04,label="A",size=11)

p2 <- ggplot(d,aes(x=QG,y=SIWG,color=PS2))+
  geom_point(aes(fill=PS2),pch=21,color="black",size=1.5)+
  theme_bw()+
  ylim(-0.12,0.06)+
  xlim(0,40)+
  geom_smooth(method="lm",se=FALSE,color="red")+
  annotate("text",x=40,y=0.05,label="***",size=7)+
  ylab("nurse gaster sociality index")+xlab("queen gaster normalized expression")+
  scale_fill_manual(values=PS_palette,name="phylostrata")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        legend.title=element_text(size=20),
        legend.position = c(0.83, 0.22),
        legend.background=element_rect(fill=FALSE),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))+
  annotate("text",x=5,y=0.04,label="B",size=11)

png("~/Writing/Figures/NurseLarva/FigS5.png",width=2500,height=3000,res=300)
grid.arrange(p1+
               annotate("text",x=16,y=0.06,label="myosin regulatory light chain 2",size=4)+
               geom_segment(aes(x=21,xend=22.9,y=0.055,yend=0.051),
                            arrow=arrow(length=unit(0.3,"cm")),
                            color="black")+
               annotate("text",x=35,y=0.035,label="cathD",size=4)+
               geom_segment(aes(x=36,xend=38.5,y=0.03,yend=0.02),
                            arrow=arrow(length=unit(0.3,"cm")),
                            color="black")+
               annotate("text",x=30,y=0.047,label="troponin I",size=4)+
               geom_segment(aes(x=29,xend=27.6,y=0.042,yend=0.02),
                            arrow=arrow(length=unit(0.3,"cm")),
                            color="black"),
             p2+
               annotate("text",x=7.6,y=-0.1,label="myosin regulatory light chain 2",size=4)+
               geom_segment(aes(x=11.5,xend=13.8,y=-0.09,yend=-0.053),
                            arrow=arrow(length=unit(0.3,"cm")),
                            color="black"),
             ncol=1)
dev.off()


########
##conceptual figure
library(rsvg)
image <- rsvg("~/Downloads/smaller larva with shadow.svg")
mult=dim(image)[1]/dim(image)[2]

#yT is the center of the cluster, layD is the number of points on each side of the grid,
#nGene is number of genes/cluster, rad is proportional to the size of the circle, in terms of proportion of grid
randCluster <- function(yT,layD,nGene,rad){
  rad = layD*rad
  ret = c()
  angle.inc = 2*pi/nGene
  for (i in 1:nGene){
    yadj = floor(sin(i*angle.inc)*rad+0.5)
    xadj = floor(cos(i*angle.inc)*rad+0.5)
    adj = yT+xadj+yadj*layD
    ret = c(ret,adj)
  }
  return(ret)
}


library(igraph)
nGene = 5
nGT = nGene*3
Soc_palette = c("#dadaeb","#9e9ac8","#6a51a3")


links <- data.frame(tissue1=c(rep(rep(c("larv","wh","wg"),each=nGene),nGT),"not","not"),
                    tissue2=c(rep(c("larv","wh","wg"),each=nGene*nGT),"not","not"),
                    ID1=paste("gene",c(rep(1:nGT,each=nGT),nGT+1,nGT+2)),
                    ID2=paste("gene",c(rep(seq(1:nGT),nGT),nGT+1,nGT+2)),weight=c(rep(rnorm(nGT*nGT,0,1)),5,5),
                    edgeCol = "gold")
links$edgeCol=as.character(links$edgeCol)
links$edgeCol[links$tissue1==links$tissue2]="#6a51a3"
links$edgeCol[links$tissue1!=links$tissue2]="#dadaeb"
links$edgeCol[(links$tissue1=="wh"&links$tissue2=="wg")|(links$tissue1=="wg"&links$tissue2=="wh")]="#9e9ac8"
links$weight[links$tissue1==links$tissue2]=rnorm(sum(links$tissue1==links$tissue2),0.9,1)


##Trim number of connections
keep1 = links[sample(rownames(links[links$edgeCol=="#6a51a3",]),24),]
keep2 = links[sample(rownames(links[links$edgeCol=="#dadaeb",]),5),]
keep3 = links[sample(rownames(links[links$edgeCol=="#9e9ac8",]),5),]
keep = rbind(keep1,keep2,keep3)
keep = keep[keep$ID1!=keep$ID2,]
missing = unique(links$ID1)[!unique(links$ID1) %in% c(as.character(keep$ID1),as.character(keep$ID2))]
for (gene in missing){
  keep = rbind(keep,links[sample(rownames(links[links$ID1==gene,]),1),])
}

keep$lty = "solid"
keep$lty[keep$edgeCol=="#dadaeb"]="dotted"
keep$lty[keep$edgeCol=="#9e9ac8"]="dashed"
keep$curved= FALSE
keep$curved[keep$edgeCol!="#6a51a3"]=0.2
keep = keep[,c(3:8)]

nodes <- paste("gene",seq(1:(nGT+2)))
net <- graph_from_data_frame(d=keep, vertices=nodes, directed=T) 

E(net)$width = rep(3,nrow(keep))
V(net)$color = rep("black",nGT+2)
E(net)$edge.color = rep("black",nrow(keep))

net <- simplify(net, remove.multiple = F, remove.loops = T) 


######
##Layout on lattice
######

layD = 100
lay <- layout_on_grid(make_lattice( c(layD,layD)))
lcenter = floor(layD/2)*layD+floor(layD/6)+4
whcenter = floor(3*layD/4)*layD+floor(3*layD/6)
wgcenter = floor(3*layD/5)*layD+layD+floor(5*layD/6)+1

larv = randCluster(lcenter,layD,nGene,0.1)
wh = randCluster(whcenter,layD,nGene,0.1)
wg = randCluster(wgcenter,layD,nGene,0.1)

layout =lay[c(larv,wh,wg,floor(layD/2),layD*layD-floor(layD/2)),]
V(net)$size=rep(4,nGT+2)
V(net)$color=c(rep(c("red","darkblue","lightblue"),each=nGene),"black","black")
V(net)$shape="square"

plot.new()
a=heightDetails(grid.text("within-tissue\n",x=0.15,y=0.88,gp=gpar(fontsize=20),just="left"))
lineh = as.numeric(substr(a,1,6))/8.3 ##8.3 is number of inches the figure is tall

png("~/Writing/Figures/NurseLarva/conceptFig.png",width=2500,height=2500,res=300)
par(mar=c(0,0,0,0))
plot(net,edge.color=E(net)$edge.color,edge.arrow.size=0.8,vertex.label=NA,
     layout=layout,edge.curved=E(net)$curved)
grid.raster(image,x=.45,y=0.45,width=0.8/mult,height=0.8)
grid.curve(0.13,0.28,0.04,0.43,curvature=0.2,gp=gpar(lwd=1.5,col="darkgrey"))
grid.curve(0.2,0.32,0.26,0.48,curvature=-0.1,gp=gpar(lwd=1,col="darkgrey"))
grid.curve(0.45,0.31,0.37,0.67,curvature=0.2,gp=gpar(lwd=1.5,col="darkgrey"))
grid.curve(0.53,0.33,0.58,0.7,curvature=-0.2,gp=gpar(lwd=1,col="darkgrey"))
grid.curve(0.75,0.39,0.7,0.55,curvature=0.2,gp=gpar(lwd=1.5,col="darkgrey"))
grid.curve(0.85,0.37,0.9,0.5,curvature=-0.2,gp=gpar(lwd=1,col="darkgrey"))
grid.rect(x=0.45,y=0.05,width=0.1,height=0.05,gp=gpar(col="white"))
grid.rect(x=0.45,y=0.95,width=0.1,height=0.05,gp=gpar(col="white"))
grid.text("within-tissue\nbetween-tissue\nsocial",x=0.15,y=0.86,gp=gpar(fontsize=20),just="left")
grid.lines(x=unit(c(0.04,0.13),"npc"),y=unit(c(0.86+lineh,0.86+lineh),"npc"),gp=gpar(lwd=1.5,lty="solid"))
grid.lines(x=unit(c(0.04,0.13),"npc"),y=unit(c(0.86,0.86),"npc"),gp=gpar(lwd=1.5,lty="dashed"))
grid.lines(x=unit(c(0.04,0.13),"npc"),y=unit(c(0.86-lineh,0.86-lineh),"npc"),gp=gpar(lwd=1.5,lty="dotted"))
grid.text(expression(underline("connection type")),x=0.15,y=0.96,just="left",gp=gpar(fontsize=22,fontface="bold"))
dev.off()


#################
##GO term analysis using drosophila orthologs
#################
library(topGO)
library(grid)
library(gridExtra)

#Use to get Dmel mappings if we want them
getDmel <- function(){
  m <- inverseList(annFUN.org("BP",mapping="org.Dm.eg.db",ID="entrez"))
  
  x <- org.Dm.egFLYBASEPROT
  # Get the entrez gene IDs that are mapped to a Flybase prot ID
  mapped_genes <- mappedkeys(x)
  # Convert to a list
  xx <- as.list(x[mapped_genes])
  new <- list()
  for (i in 1:length(m)){
    entrez = names(m)[i]
    flyBase = xx[entrez]
    if (is.na(names(flyBase))){
      m[[i]]=NULL
      next;
    }
    names(m)[i]=flyBase[[1]][1]
    if (length(flyBase[[1]]) > 1){
      for (j in 2:length(flyBase)){
        newL = list(m[i][[1]])
        names(newL) = flyBase[[1]][j]
        new = c(new,newL)
      }
    }
  }
  
  m = c(m,new)
  return(m)
}

geneList = dfO$SIWG
names(geneList)=dfO$gene_Dmel_FLYBASE

#Not necessary for GSEA, but necessary for identifying significant genes
select <- function(score){
  return(score > quantile(score,0.8))
}

#Previously derived GO annotations
go <- read.csv("~/Writing/Data/NurseSpecialization_transcriptomicData/GOannotation.csv")
new <- list()
for (gene in unique(go$gene)){
  d = go[go$gene %in% gene,]
  new[[gene]]=as.character(d$GO)
}

#Gene set enrichment analysis. Input is a named list of values
getGSEA <- function(stat){
  GOdata <- new("topGOdata",
                description="Simple session",ontology="BP",
                allGenes=stat,geneSel=select,
                nodeSize = 10,
                annot=annFUN.gene2GO,gene2GO=new)
  
  
  resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks",scoreOrder="decreasing")
  allRes <- GenTable(GOdata,KS=resultKS,topNodes=100)
  return(allRes[allRes$Significant>allRes$Expected,])
}

d = list(df$L_WH,df$L_WG,df$WH_L,df$WG_L,df$SIL_WH,df$SIL_WG,df$SIWH,df$SIWG)
res = list()
for (i in 1:8){
  stat = d[[i]]
  names(stat)=df$gene
  res[[i]]=getGSEA(stat)
}

names = c("larva-head","larva-abdomen","head-larva","abdomen-larva",
          "larva-head","larva-abdomen","head","abdomen")

tabGO <- function(res){
  return(res$Term[res$KS < 0.05])
}

GOres <- lapply(res,tabGO)
max = 10 #Only keep 10 GO terms
res = sapply(GOres,'[',1:max) #keep 10 GO terms
res <- as.data.frame(res)
colnames(res) = names

tt3 <- ttheme_minimal(
  core=list(
    fg_params=list(fontsize=8)),
  colhead=list(fg_params=list(fontface="bold",fontsize=14)))

#Make table grob out of results
resT <- tableGrob(res[,5:8],theme=tt3,rows=NULL)
grid.arrange(resT)
g <- gtable_add_grob(resT,
                     grobs = segmentsGrob( # line across the bottom
                       x0 = unit(0,"npc"),
                       y0 = unit(0,"npc"),
                       x1 = unit(1,"npc"),
                       y1 = unit(0,"npc"),
                       gp = gpar(lwd = 3.0)),
                     t = 1, b = 1, l = 1, r = 4)
separators <- replicate(ncol(g) - 1,
                        segmentsGrob(x1 = unit(0, "npc"), gp=gpar(lty=2)),
                        simplify=FALSE)
## add vertical lines on the left side of columns (after 2nd)
g <- gtable::gtable_add_grob(g, grobs = separators,
                             t = 2, b = nrow(g), l = seq_len(ncol(g)-1)+1)
p <- grid.arrange(g)
ggsave(p,file="~/Writing/Figures/NurseLarva/TableS1.png",width=10,height=4,dpi=300)







######
##Mock phylogeny image
######
library(ggdendro)
d <- matrix(runif(100,0,1),ncol=10)
hc <- hclust(dist(d))
dhc <- as.dendrogram(hc)
ddata <- dendro_data(dhc, type = "rectangle")
seg = ddata$segments
for (i in 1:nrow(seg)){
  if (seg$xend[i]!=seg$x[i]){
    if (seg$xend[i] > seg$x[i]){
      seg$xend[i]= seg$xend[i] + 0.065
    } else {
      seg$xend[i]= seg$xend[i] - 0.065
    }
  }
}
ggplot(seg) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend),color="black",size=2) + 
  coord_flip()+
  scale_y_reverse()+
  xlab("")+ylab("")+
  theme(panel.background = element_rect(fill="transparent",color=NA),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_rect(fill="transparent",color=NA),
        axis.text=element_blank(),
        axis.ticks=element_blank())

ggsave("~/Writing/Figures/NurseLarva/dummyPhylo.png",bg="transparent")




library(jsonlite)
d <- fromJSON("~/genome_assembly/mphar/histogram_contig.json")
f <- read.table("~/genome_assembly/scafs.txt",sep=',')
f$Length = f$V2
for (i in 1:nrow(f)){
  f$cum[i]=sum(f$Length[1:i])
}

n50 = min(f$Length[f$cum < 0.5*f$cum[nrow(f)]])
n90 = min(f$length[f$cum > 0.1*sum(f$nbase)])
l50 = sum(f$N[f$cum > 0.5*sum(f$nbase)])

d <- fromJSON("~/genome_assembly/mphar/histogram_scaffold.json")

f = data.frame(N = d$vals)
f$length=seq(from=d$min,to=d$max,by=d$binsize)
f$nbase = f$length*f$N
for (i in 1:nrow(f)){
  f$cum[i]=colSums(f[1:i,])[3]
}

n50 = min(f$length[f$cum > 0.5*sum(f$nbase)])
n90 = min(f$length[f$cum > 0.9*sum(f$nbase)])
l50 = sum(f$N[f$cum > 0.5*sum(f$nbase)])

load("~/Dropbox/monomorium nurses/data.processed/cleandata.RData")
write.csv(fpkm,file="~/Data/Nurse_Larva/fpkm.csv")

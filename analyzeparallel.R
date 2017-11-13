cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
library(scales)
library(plyr)
library(reshape2)
library(gridExtra)
library(grid)
library(pBrackets)
library(ggplot2)

PS_palette = c("#deebf7","#9ecae1","#4292c6","#2171b5","#084594")
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

df <- read.csv("~/Data/new/CompiledConnections.csv")

df <- read.csv("~/Data/new/TopExprWorkerNetSocialityDF.csv")
load("~/Dropbox/monomorium nurses/data.processed/cleandata.RData")
exprH = rowSums(fpkm[,grepl("C.*_WH",colnames(fpkm))])/sum(grepl("C.*_WH",colnames(fpkm)))
exprG = rowSums(fpkm[,grepl("C.*_WG",colnames(fpkm))])/sum(grepl("C.*_WG",colnames(fpkm)))
exprL = rowSums(fpkm[,grepl("W_L",colnames(fpkm))])/sum(grepl("W_L",colnames(fpkm)))
exprQH = rowSums(fpkm[,grepl("_QH",colnames(fpkm))])/sum(grepl("_QH",colnames(fpkm)))
exprQG = rowSums(fpkm[,grepl("_QG",colnames(fpkm))])/sum(grepl("_QG",colnames(fpkm)))
expr = data.frame(exprH=exprH,exprG=exprG,exprL=exprL,
                  exprQH=exprQH,exprQG=exprQG,gene=names(exprH))

df = df[,c(1,2,4:11)]
df$SIL = df$Lbetween - df$Lwithin
df$SIWH = df$WHbetween - df$WHwithin
df$SIWG = df$WGbetween - df$WGwithin
df$SI_Overall = (df$SIL + df$SIWH + df$SIWG)/3
df$SIWH_worker = df$WH.WG - df$WHwithin
df$SIWG_worker = df$WG.WH - df$WGwithin

df <- merge(df,expr,by="gene")
dfAll = df

d2 <- melt(df[,c(1,3:11)],id.vars="gene")
d2$tissue="larva"
d2$tissue[grepl("^WG",d2$variable)]="nurse gaster"
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

########
##Adding MKtest constraint
f.est <- read.csv("~/Data/MKtestConstraintOneAlpha.csv")
colnames(f.est) = c("gene","f.est")
dfAll <- merge(df,f.est,by="gene")

p1 <- ggplot(dfAll,aes(x=SIWH,y=1-f.est))+
  geom_point(size=0.8)+
  theme_bw()+ylab("")+xlab("")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),
        legend.position = "none")+
  ggtitle("nurse head")+
  geom_smooth(method="lm",se=FALSE,color="red")+
  annotate("text",x=-0.065,y=0.25,label="r = - 0.129\np = 0.001",size=5)

p2 <- ggplot(dfAll,aes(x=SIWG,y=1-f.est))+
  geom_point(size=0.8)+
  theme_bw()+ylab("")+xlab("")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),legend.position = "none")+
  geom_smooth(method="lm",se=FALSE,color="red")+
  ggtitle("nurse gaster")+
  annotate("text",x=-0.08,y=0.25,label="r = - 0.131\np = 0.001",size=5)

p3 <- ggplot(dfAll,aes(x=SIL,y=1-f.est))+
  geom_point(size=0.8)+
  theme_bw()+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),
        legend.position = "none")+ylab("")+xlab("")+
  geom_smooth(method="lm",se=FALSE,color="red")+
  ggtitle("larva")+
  annotate("text",x=-0.11,y=0.25,label="r = - 0.086\np = 0.034",size=5)
png("~/Writing/Figures/NurseLarva/Fig2.png",width=3000,height=1500,res=300)
grid.arrange(p1+theme(plot.title=element_text(size=20,hjust=0.5)),
             p2+theme(plot.title=element_text(size=20,hjust=0.5)),
             p3+theme(plot.title=element_text(size=20,hjust=0.5)),
             bottom = textGrob("sociality index", gp=gpar(fontsize=23,font=8),hjust=0.3,vjust=-0.8),
             left = textGrob("constraint", rot = 90, vjust = 1.5,gp=gpar(fontsize=23,font=8)),nrow=1)
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

##import annotation information
annNew <- read.csv("~/Dropbox/monomorium nurses/Amel_Mphar/ThreeWayOGGMap.csv",header=TRUE)
annMp <- annNew[,c(2,5)]
colnames(annMp) = c("Gene","uniprot drosophila")
annMp <- annMp[!duplicated(annMp$Gene),]
dA <- merge(d,annMp,by="Gene")
dATop = dA[dA$QH > 18 &dA$SIWH > 0,]

lm <- glm(log(f.est) ~ PS2 + SIWH + exprQG + exprH + exprG,
          data = dfAll)

lm <- glm(WGwithin~ PS2,data=dfAll)

drop1(lm,.~.,test="Chi") 
summary(glht(lm, mcp(PS2="Tukey"))) 


d2 <- melt(dfAll[,c(1,30:35,27,48)],id.vars=c("Gene","PS2","midF"))
d2$tissue="Larva"
d2$tissue[grepl("WG",d2$variable)]="WorkerGaster"
d2$tissue[grepl("WH",d2$variable)]="WorkerHead"

d2$connection.type = "within"
d2$connection.type[grepl("between",d2$variable)]="between"
d3 <- dcast(d2, PS2+tissue+Gene+midF ~ connection.type)

d4 <- ddply(d3,~PS2 + tissue,summarise,
            meanB = mean(between),
            c1B = quantile(bootMean(between,1000),0.025),
            c2B = quantile(bootMean(between,1000),0.975),
            meanF = mean(within),
            c1W = quantile(bootMean(within,1000),0.025),
            c2W = quantile(bootMean(within,1000),0.975))

ggplot(d4,aes(x=meanB,y=meanF,color=PS2,shape=tissue))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=c1W,ymax=c2W))+
  geom_errorbarh(aes(xmin=c1B,xmax=c2B))+
  xlab("between")+ylab("within")

########
##
corWH = cor(t(fpkm[,grepl("WH",colnames(fpkm))]))
connWH = data.frame(gene = rownames(corWH),connWH = rowSums(abs(corWH)))
corWG = cor(t(fpkm[,grepl("WG",colnames(fpkm))]))
connWG = data.frame(gene = rownames(corWG),connWG = rowSums(abs(corWG)))
corL = cor(t(fpkm[,grepl("W_L",colnames(fpkm))]))
connL = data.frame(gene = rownames(corL),connL = rowSums(abs(corL)))

df2 = merge(df,connWH,by="gene")
df2 = merge(df2,connWG,by="gene")
df2 = merge(df2,connL,by="gene")

#########
#GO enrichment of the 100 genes with the highest connectivity in each category, plus the highest social indices in each category
#Note that the "universe" is defined as the 1000 genes with the highest expression
load("~/Dropbox/workspace/GODB.RData")
universe = universe[universe %in% dfAll$Gene]

GOenrich <- function(column){
  df = dfAll[dfAll$Gene %in% universe,]
  genes = df$Gene[df[,column] > quantile(df[,column],0.8)]
  return(GOstatUniverse(genes,universe))
}

GOdat = list()
for (i in 30:43){
  column = colnames(dfAll)[i]
  GOdat[[column]] = GOenrich(column)
}

df <- data.frame(ConnType=names(GOdat),FirstTerm=NA,SecondTerm=NA,ThirdTerm=NA)
for (i in 1:length(GOdat)){
  for (j in 1:3){
    df[i,(j+1)] = GOdat[[names(GOdat)[i]]]$Term[j]
  }
}

write.csv(df,"~/Data/GObyConn.csv")
df <- data.frame(ConnType=names(GOdat),Gene1=NA,Gene2=NA,Gene3=NA,
                 Gene4=NA,Gene5=NA,Gene6=NA,Gene7=NA,Gene8=NA,Gene9=NA,Gene10=NA)
df2 = df
for (i in 1:length(GOdat)){
  d = dfAll[order(dfAll[,names(GOdat[i])],decreasing=TRUE),]
  for (j in 1:10){
    df[i,(j+1)] = as.character(d$SwissProt[j])
    df2[i,(j+1)] = as.character(d$UniProt[j])
  }
}

write.csv(df,file="~/Data/SwissProtbyConn.csv")
write.csv(df2,file="~/Data/UniProtbyConn.csv")

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


library(rsvg)
image <- rsvg("~/Downloads/smaller larva with shadow.svg")
mult=dim(image)[1]/dim(image)[2]

png("~/Writing/Figures/NurseLarva/sepNurseLarv.png",width=2500,height=2500,res=300)
grid.raster(image,x=.45,y=0.45,width=0.8/mult,height=0.8)
dev.off()

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

setwd("~/Writing/Data/NurseSpecialization_transcriptomicData/")
ogg <- read.csv("ThreeWayOGGMap.csv") #Import 3-way OGG map
df <- read.csv("Dmel_secreted.csv") #Derived from http://www.flyrnai.org/tools/glad/web/
df = df[,c(1,2)]
colnames(df)[2] = "Flybase"
key <- read.table("DmelKey.txt") #Generated from Drosophila melanogaster gff file
oggK = merge(ogg,key,by.x="gene_Dmel_FLYBASE",by.y="V1")
oggS = merge(oggK,df,by.x="V2",by.y="Flybase") #only secreted orthologs

dfAll$secreted = "no"
dfAll$secreted[dfAll$gene %in% oggS$gene_Mphar] = 'yes'

mean(dfAll$SI_Overall[dfAll$secreted=="no"])
mean(dfAll$SI_Overall[dfAll$secreted=="yes"])

d = melt(dfAll[,c(1,11:14,22)],id.vars=c("gene","secreted"))

ggplot(d,aes(x=variable,y=value,fill=secreted))+
  geom_boxplot(notch=TRUE)  

wilcox.test(dfAll$SIWG[dfAll$secreted=='no'],dfAll$SIWG[dfAll$secreted=='yes'])

####################
##Trying a purely correlational approach
####################
setInput <- function(codes,names,fpkm){
  AllExpr = list()
  for (i in 1:5){
    nSamp = getNsamp(codes,i)
    expr = list()
    for (j in 1:length(codes)){
      f = factors[grepl(codes[j],factors$sample.id)&factors$stage==i,]
      samps = sample(f$sample.id,nSamp,replace=FALSE) ##Want same number of samples per sample type
      expr[[j]] = fpkm[,samps]
      expr[[j]]$tissue=names[j]
      expr[[j]]$gene = rownames(expr[[j]])
      expr[[j]]$tissue_gene = with(expr[[j]],paste(tissue,gene,sep="_"))
      colnames(expr[[j]]) = c(paste("Stage",i,"_",seq(1,nSamp),sep=""),"tissue","gene","tissue_gene")
    }
    AllExpr[[i]] = ldply(expr,data.frame)
    rownames(AllExpr[[i]]) = AllExpr[[i]]$tissue_gene
    AllExpr[[i]] = AllExpr[[i]][,c(1:nSamp)]
  }
  input <- do.call(cbind,AllExpr)
  return(input)
}

getNsamp <- function(codes,stage){
  nSamp = c()
  for (code in codes){
    f = factors[grepl(code,factors$sample.id)&factors$stage==stage,]
    nSamp=c(nSamp,nrow(f))
  }
  return(min(nSamp))
}

deriveCodes <- function(samp){
  if (samp=="worker"){
    codes = c("W.*_L","C.*WH","C.*WG")
    names = c("WorkLarv","WorkNurseH","WorkNurseG")
  } else if (samp=="random"){
    codes = c("QW.*_L","R.*WH","R.*WG")
    names = c("WorkLarvQR","RandNurseH","RandNurseG")
  } else if (samp=="forager"){
    codes = c("W.*_L","F.*WH","F.*WG")
    names = c("WorkeLarv","ForagerH","ForagerG")
  } else if (samp=="workerQR"){
    codes = c("QW.*_L","QCH","QCG")
    names = c("WorkLarvQR","WorkNurseHQR","WorkNurseGQR")
  }
  return(list(codes,names))
}

sortData <- function(N,fpkm,codes){
  fpkm <- fpkm[,grepl(codes[1],colnames(fpkm))|grepl(codes[2],colnames(fpkm))|grepl(codes[3],colnames(fpkm))]
  rowS = rowSums(fpkm)
  keep = rowS[order(rowS,decreasing=TRUE)]
  keep = names(keep)[1:N]
  fpkm = fpkm[keep,]
  fpkm = log(fpkm + sqrt(fpkm ^ 2 + 1)) #hyperbolic sine transformation to normalize gene expression data
  
  return(fpkm)
}

getConns <- function(Ngenes,sample){
  load("~/Dropbox/monomorium nurses/data.processed/cleandata.RData")
  n <- deriveCodes(sample) #returns names, codes
  fpkm <- sortData(Ngenes,fpkm,n[[1]])
  d <- setInput(n[[1]],n[[2]],fpkm) #makes meta-sample expression dataframe (genes labeled by individual)
  cMat <- cor(t(d)) #calculate gene-gene correlation
  nGene = nrow(fpkm)
  
  res <- matrix(nrow=1,ncol=4)
  for (i in 1:nrow(cMat)){
    for (j in 1:3){
      if (i/nGene > j){ #Will iterate over the three blocks (matrix is symetric)
        next;
      } else {
        gene <- gsub(".*_","",rownames(cMat)[i])
        samp1 <- gsub("_.*","",rownames(cMat)[i])
        samp2 <- n[[2]][j] #n[[2]] holds the names list
        cAbs = mean(abs(cMat[i,(nGene*(j-1)+1):(nGene*j)])) #Subset so calculates over one tissue
        res = rbind(res,c(gene,samp1,samp2,cAbs))
      }
    }
  }
  
  res = as.data.frame(res[-c(1),])
  res$samp = as.factor(do.call(paste,c(res[,c(2,3)],list(sep="_"))))
  res$V4 = as.numeric(as.character(res$V4))
  res$V1 = as.character(res$V1)
  rh = dcast(res[,c(1,4,5)],V1~samp,value.var="V4")
  colnames(rh)[1] = "gene"
  return(rh)
}

rand = getConns(1000,"random")
work = getConns(1000,"workerQR")
rand$SIWH=rand$WorkLarvQR_RandNurseH-rand$RandNurseH_RandNurseH
rand$SIWG=rand$WorkLarvQR_RandNurseG-rand$RandNurseG_RandNurseG
rand$SIL=(rand$WorkLarvQR_RandNurseH+rand$WorkLarvQR_RandNurseG)/2-rand$WorkLarvQR_WorkLarvQR
work$SIWH=work$WorkLarvQR_WorkNurseHQR-work$WorkNurseHQR_WorkNurseHQR
work$SIWG=work$WorkLarvQR_WorkNurseGQR-work$WorkNurseGQR_WorkNurseGQR
work$SIL=(work$WorkLarvQR_WorkNurseHQR+work$WorkLarvQR_WorkNurseGQR)/2-work$WorkLarvQR_WorkLarvQR
wilcox.test(rand$SIL,work$SIL,paired=TRUE,alternative="less")
wilcox.test(rand$SIWH,work$SIWH,paired=TRUE,alternative="less")
wilcox.test(rand$SIWG,work$SIWG,paired=TRUE,alternative="less")
wilcox.test(rand$WorkLarvQR_WorkLarvQR,work$WorkLarvQR_WorkLarvQR,paired=TRUE,alternative="less")


rh$SIWH = -rh$RandNurseH_RandNurseH + rh$WorkLarvQR_RandNurseH
rh$SIWG = -rh$RandNurseG_RandNurseG + rh$WorkLarvQR_RandNurseG
rh$SIL = (rh$WorkLarvQR_RandNurseG+rh$WorkLarvQR_RandNurseH)/2 - rh$WorkLarvQR_WorkLarvQR
f.est <- read.csv("~/Data/MKtestConstraintOneAlpha.csv")
colnames(f.est) = c("gene","f.est")
dfAll <- merge(f.est,rh,by="gene")
dfAll$SI_Overall = (dfAll$SIWH+dfAll$SIWG+dfAll$SIL)/3




ext <- read.csv("~/Downloads/msx123_Supp (1)/MpharAnn.csv")
dfAll = merge(ext,rh,by.x="Gene",by.y="gene")

dfAll$PS2 = factor(dfAll$PS2, levels = c("cellular","eukaryote","bilaterian","insect","hymenopteran_ant"))
dfAll = dfAll[!is.na(dfAll$PS2),]
levels(dfAll$PS2)[5]="hymenopteran"

ggplot(dfAll,aes(x=PS2,y=SI_Overall))+geom_boxplot(notch=TRUE)


lm <- glm(SIWH ~ PS2, data=dfAll)
lm <- glm(SI_Overall ~ PS2, data=dfAll)
summary(glht(lm, mcp(PS2="Tukey"))) 


ggplot(res,aes(x=samp,y=V4))+
  geom_boxplot()

res$


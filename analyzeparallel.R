cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
library(scales)
library(plyr)
library(reshape2)
library(gridExtra)
library(grid)
library(pBrackets)

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

df <- read.csv("~/Data/TopExprWorkerNetNetSocialityDF.csv")

df$SIL = df$Lbetween - df$Lwithin
df$SIWH = df$WHbetween - df$WHwithin
df$SIWG = df$WGbetween - df$WGwithin
df$SI_Overall = (df$SIL + df$SIWH + df$SIWG)/3
df$SIWH_worker = df$WH.WG - df$WHwithin
df$SIWG_worker = df$WG.WH - df$WGwithin

d2 <- melt(df[,c(2:10)],id.vars="gene")
d2$tissue="larva"
d2$tissue[grepl("^WG",d2$variable)]="nurse gaster"
d2$tissue[grepl("^WH",d2$variable)]="nurse head"
d2$connection.type = "within tissue"
d2$connection.type[grepl("between",d2$variable)]="social"
d2$connection.type[grepl("WH.WG",d2$variable)|grepl("WG.WH",d2$variable)]="between tissue"

p1 <- ggplot(d2,aes(x=tissue,y=value,fill=connection.type))+
  geom_boxplot(notch=TRUE)+
  theme_bw()+xlab("tissue")+ylab("connection strength")+
  scale_fill_manual(values=Soc_palette,name="connection type")+
  theme(axis.text=element_text(size=17),
        axis.title=element_text(size=23),
        axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)),
        axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)),
        legend.text=element_text(size=14),
        legend.title=element_text(size=17),
        legend.position = c(.05, .98),
                  legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.margin = margin(6, 6, 6, 6),
        legend.background=element_rect(fill=FALSE))+
  ylim(0,0.25)

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
  ggtitle("A. nurse head")+
  geom_smooth(method="lm",se=FALSE,color="red")+
  annotate("text",x=-0.065,y=0.25,label="r = 0.110\np = 0.007",size=5)

p2 <- ggplot(dfAll,aes(x=SIWG,y=1-f.est))+
  geom_point(size=0.8)+
  theme_bw()+ylab("")+xlab("")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),legend.position = "none")+
  geom_smooth(method="lm",se=FALSE,color="red")+
  ggtitle("B. nurse gaster")+
  annotate("text",x=-0.08,y=0.25,label="r = 0.170\np < 0.001",size=5)

p3 <- ggplot(dfAll,aes(x=SIL,y=1-f.est))+
  geom_point(size=0.8)+
  theme_bw()+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),
        legend.position = "none")+ylab("")+xlab("")+
  geom_smooth(method="lm",se=FALSE,color="red")+
  ggtitle("C. larva")+
  annotate("text",x=-0.11,y=0.25,label="r = 0.123\np = 0.002",size=5)
png("~/Writing/Figures/NurseLarva/Fig2.png",width=3000,height=1500,res=300)
grid.arrange(p1+theme(plot.title=element_text(size=20,hjust=0.5)),
             p2+theme(plot.title=element_text(size=20,hjust=0.5)),
             p3+theme(plot.title=element_text(size=20,hjust=0.5)),
             bottom = textGrob("sociality index", gp=gpar(fontsize=23,font=8),hjust=0.3,vjust=-0.8),
             left = textGrob("constraint", rot = 90, vjust = 1.5,gp=gpar(fontsize=23,font=8)),nrow=1)
dev.off()

###########
##Adding phylostrata
ext <- read.csv("~/Downloads/msx123_Supp (1)/External Database S1.csv")
dfAll = merge(ext,df,by.x="Gene",by.y="gene")

dfAll$PS2 = factor(dfAll$PS2, levels = c("cellular","eukaryote","bilaterian","insect","hymenopteran_ant"))
dfAll = dfAll[!is.na(dfAll$PS2),]
levels(dfAll$PS2)[5]="hymenopteran"

lm <- glm(SIL ~ PS2, data=dfAll)
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
  annotate("text",x=-0.055,y=0.925,label="eukaryote",size=5)+
  annotate("text",x=-0.057,y=0.88,label="bilaterian",size=5)+
  annotate("text",x=-0.049,y=0.908,label="cellular",size=5)+
  annotate("text",x=-0.046,y=0.87,label="insect",size=5)+
  annotate("text",x=-0.040,y=0.845,label="hymenopteran",size=5)



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
    res= cor.test(d2$f.est,d2[,cols[j]])
    corMat[i,j] = res$estimate
    pMat[i,j] = res$p.value
  }
  res = cor.test(dfAll$f.est,dfAll[,cols[j]])
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
corM$text = apply(round(corM[,c(3,4)],3),1,paste,collapse="\n")
corM$text[corM$p.value < 0.001]=gsub("\n0","\n<0.001",corM$text[corM$p.value < 0.001])

rowN = c()
for (i in 1:5){
  rowN = c(rowN,paste(colnames(pMat)[i]," (",sum(dfAll$PS2==colnames(pMat)[i]),")",sep=""))
}
rowN = c(rowN,paste("overall (",nrow(dfAll),")",sep=""))

p <- ggplot(corM,aes(x=Var2,y=factor(Var1)))+geom_tile(aes(fill=value))+
  geom_text(aes(label=text))+
  scale_fill_gradient2(name="pearson r",low=muted("blue"),high=muted("red"),mid="white",
                       limits=c(-0.5,0.5))+
  xlab("phylostrata")+ylab("connectivity measurement")+
  guides(fill= guide_colorbar(barheight=2,
                              barwidth=30,nbin=40,
                              title.vjust=0.75))+
  theme_bw()+
  geom_hline(yintercept=3.5,color="black")+
  geom_hline(yintercept=6.5,color="black")+
  scale_y_discrete(labels=colN,expand=c(0,0))+
  scale_x_discrete(labels=rowN,expand=c(0,0))+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.position = "top",
        legend.text=element_text(size=15),
        legend.title=element_text(size=18),
        axis.text.x=element_text(angle = 45,hjust=1))

png("~/Writing/Figures/NurseLarva/FigS4.png",width=3000,height=3000,res=300)
p
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
  ylab("nurse head sociality index")+xlab("queen head normalized expression")+
  scale_fill_manual(values=PS_palette,name="phylostrata")+
  ylim(-0.12,0.06)+
  xlim(0,40)+
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
                                color="red")+
               annotate("text",x=35,y=0.035,label="cathD",size=4)+
               geom_segment(aes(x=36,xend=38.5,y=0.03,yend=0.02),
                            arrow=arrow(length=unit(0.3,"cm")),
                            color="red")+
               annotate("text",x=30,y=0.047,label="troponin I",size=4)+
               geom_segment(aes(x=29,xend=27.6,y=0.042,yend=0.02),
                            arrow=arrow(length=unit(0.3,"cm")),
                            color="red"),
             p2+
               annotate("text",x=7.6,y=-0.1,label="myosin regulatory light chain 2",size=4)+
               geom_segment(aes(x=11.5,xend=13.8,y=-0.09,yend=-0.053),
                            arrow=arrow(length=unit(0.3,"cm")),
                            color="red"),
             ncol=1)
dev.off()

##import annotation information
annNew <- read.csv("~/Dropbox/monomorium nurses/Amel_Mphar/ThreeWayOGGMap.csv",header=TRUE)
annMp <- annNew[,c(2,5)]
colnames(annMp) = c("Gene","uniprot drosophila")
annMp <- annMp[!duplicated(annMp$Gene),]
dA <- merge(d,annMp,by="Gene")
dATop = dA[dA$QH > 18 &dA$SIWH > 0,]

lm <- glm(log(midF) ~ PS2 + QH + QG + 
            SIWH + SIWG + SIL,
          data = d)

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




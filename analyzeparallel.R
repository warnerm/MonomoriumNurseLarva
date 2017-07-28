cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
library(scales)
library(plyr)
library(reshape2)
library(gridExtra)

ext <- read.csv("~/Downloads/msx123_Supp (1)/External Database S1.csv")

snipre <- read.csv("~/Dropbox/monomorium nurses/data/bayesian_results.csv")
sn = snipre[,c(1,10,11,15,18,21,22,25)]

dnds <- read.table("~/Dropbox/monomorium nurses/data/dnds.txt")
colnames(dnds) = c("gene","dnds")
ann=read.table("~/Dropbox/monomorium nurses/data/monomorium.annotation.txt",sep="\t",head=T,quote="\"",stringsAsFactors = FALSE)   # note, problem caused by "'" in file
ann$gene=dnds$gene
a1 = ann[,c("TransLength","gene")]
ann$HSPEvalueTR=as.numeric(ann$HSPEvalueTR)

#order list based on TR Evalue
ann=ann[order(ann$HSPEvalueTR),]

#remove duplicates so only isoform with highest HSPEvalueTR remaining; 12648 gene_id now
ann=ann[!duplicated(ann$gene),]   

dnds <- merge(dnds,a1,by="gene")
dnds <- dnds[order(dnds$TransLength,decreasing=TRUE),]
dnds <- dnds[!duplicated(dnds$gene),]

f.est <- read.csv("~/Data/MKtestConstraintOneAlpha.csv")
colnames(f.est) = c("gene","f.est")
ext = merge(ext,f.est,by.x="Gene",by.y="gene")

glmConn <- function(name){
  df <- read.csv(paste("~/Data/",name,".csv",sep=""))
  
  df$SIL = df$Lbetween - df$Lwithin
  df$SIWH = df$WHbetween - df$WHwithin
  df$SIWG = df$WGbetween - df$WGwithin
  df$SI_Overall = (df$SIL + df$SIWH + df$SIWG)/3
  
  dfAll = merge(ext,df,by.x="Gene",by.y="gene")
  dfAll = merge(dfAll,sn,by.x="Gene",by.y="gene")
  
  AllLm <- SixLm(dfAll)
  
  AllDf <- parseLm(AllLm)
  
  temp = list()
  for (i in 1:3){
    temp[[i]] =  rbind(AllDf[[i*2-1]],AllDf[[i*2]])
    colnames(temp[[i]]) = paste(colnames(temp[[i]]),gsub("\\..*","",names(AllDf)[i*2]))
  }
  
  df <- do.call(cbind,temp)
  return(list(df,AllLm))
}

parseLm <- function(AllLm){
  AllDf <- list()
  for (i in 1:length(AllLm)){
    df <- summary(AllLm[[i]])$coefficients
    AllDf[[i]] = df[2:nrow(df),c(3,4)]
  }
  names(AllDf) = names(AllLm)
  return(AllDf)
}

SixLm <- function(dfAll){
  lm1 <- glm(BSnIPRE.est.x ~ Lbetween+Lwithin+WHbetween+WHwithin+WGbetween+WGwithin,data=dfAll)
  lm2 <- glm(BSnIPRE.est.x ~ SIL+SIWH+SIWG,data=dfAll)
  lm3 <- glm(log(f.est) ~ Lbetween+Lwithin+WHbetween+WHwithin+WGbetween+WGwithin,data=dfAll)
  lm4 <- glm(log(f.est) ~ SIL+SIWH+SIWG,data=dfAll)
  lm5 <- glm(log(BSnIPRE.f) ~ Lbetween+Lwithin+WHbetween+WHwithin+WGbetween+WGwithin,data=dfAll)
  lm6 <- glm(log(BSnIPRE.f) ~ SIL+SIWH+SIWG,data=dfAll)

  
  AllLm <- list(lm1,lm2,lm3,lm4,lm5,lm6)
  names(AllLm) = c("est.conn","est.soc","MKf.conn","MKf.soc","BSnIPREf.conn","BSnIPREf.soc")
  return(AllLm)
}

resWexpr <- glmConn("TopExprWorkerNetNetSocialityDF")
resSexpr <- glmConn("TopExprSexualNetNetSocialityDF")
resRexpr <- glmConn("TopExprRandomNetNetSocialityDF")
resWconn <- glmConn("TopConnWorkerNetNetSocialityDF")
resSconn <- glmConn("TopConnSexualNetNetSocialityDF")
resRconn <- glmConn("TopConnRandomNetNetSocialityDF")


bootMean <- function(data,boots){
  res = c()
  for (i in 1:boots){
    d = sample(data,length(data),replace=TRUE)
    res = c(res,mean(d,na.rm=TRUE))
  }
  return(res)
}

df <- read.csv("~/Data/TopExprWorkerNetNetSocialityDF.csv")

d2 <- melt(dfAll[,c(1,30:35)],id.vars="Gene")
d2$tissue="Larva"
d2$tissue[grepl("WG",d2$variable)]="WorkerGaster"
d2$tissue[grepl("WH",d2$variable)]="WorkerHead"
d2$connection.type = "within tissue"
d2$connection.type[grepl("between",d2$variable)]="between tissue"

png("~/Writing/Figures/NurseLarva/Fig1a.png",width=3000,height=3000,res=300)
p1 <- ggplot(d2,aes(x=tissue,y=value,fill=connection.type))+geom_boxplot()+
  theme_bw()+xlab("tissue")+ylab("connection strength")+
  scale_fill_grey(start = 0.5, end = 0.8,name="connection type")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),
        legend.position = "top")
  
dev.off()
maxConn = max(c(dfAll$WHbetween,dfAll$WHwithin,dfAll$WGbetween,
                dfAll$WHwithin,dfAll$Lwithin,dfAll$Lbetween))
png("~/Writing/Figures/NurseLarva/Fig1b.png",width=3000,height=3000,res=300)
p2 <- ggplot(dfAll,aes(x=WHwithin,y=WHbetween))+geom_point()+
  theme_bw()+xlab("within")+ggtitle("WH connection strength")+
  ylab("between")+geom_abline(slope=1)+xlim(0,0.25)+ylim(0,0.25)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))
dev.off()

png("~/Writing/Figures/NurseLarva/Fig1c.png",width=3000,height=3000,res=300)
p3 <- ggplot(dfAll,aes(x=WGwithin,y=WGbetween))+geom_point()+
  theme_bw()+xlab("within")+ggtitle("WG connection strength")+
  ylab("between")+geom_abline(slope=1)+xlim(0,0.25)+ylim(0,0.25)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))
dev.off()

png("~/Writing/Figures/NurseLarva/Fig1d.png",width=3000,height=3000,res=300)
p4 <- ggplot(dfAll,aes(x=Lwithin,y=Lbetween))+geom_point()+
  theme_bw()+xlab("within")+ggtitle("L connection strength")+
  ylab("between")+geom_abline(slope=1)+xlim(0,0.25)+ylim(0,0.25)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))
dev.off()

png("~/Writing/Figures/NurseLarva/Fig1.png",width=3000,height=3000,res=300)
grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)
dev.off()

df$SIL = df$Lbetween - df$Lwithin
df$SIWH = df$WHbetween - df$WHwithin
df$SIWG = df$WGbetween - df$WGwithin
df$SI_Overall = (df$SIL + df$SIWH + df$SIWG)/3
df$SIW = ((df$SIWH+df$SIWG)/2+df$SIL)/2

dfAll = merge(ext,df,by.x="Gene",by.y="gene")
dfAll = merge(dfAll,sn,by.x="Gene",by.y="gene")
dfAll$PS2 = factor(dfAll$PS2, levels = c("cellular","eukaryote","bilaterian","insect","hymenopteran_ant"))

dfAll$midF = (dfAll$f.est+dfAll$BSnIPRE.f)/2
dfAll = dfAll[!is.na(dfAll$PS2),]
d2 <- melt(dfAll[,c(1,30:35,27,48)],id.vars=c("Gene","PS2","midF"))
d2$tissue="Larva"
d2$tissue[grepl("WG",d2$variable)]="WorkerGaster"
d2$tissue[grepl("WH",d2$variable)]="WorkerHead"
d2$connection.type = "within"
d2$connection.type[grepl("between",d2$variable)]="between"
d3 <- dcast(d2, PS2+tissue+Gene+midF ~ connection.type)

d4 <- ddply(dfAll,~PS2,summarise,
            meanSI = mean(SI_Overall),
            c1B = quantile(bootMean(SI_Overall,1000),0.025),
            c2B = quantile(bootMean(SI_Overall,1000),0.975),
            meanF = mean(midF),
            c1W = quantile(bootMean(midF,1000),0.025),
            c2W = quantile(bootMean(midF,1000),0.975))

p3 <- ggplot(d4,aes(x=meanSI,y=meanF,color=PS2))+
  geom_point(size=3)+
  geom_errorbarh(aes(xmin=c1B,xmax=c2B))+
  geom_errorbar(aes(ymin=c1W,ymax=c2W))+
  theme_bw()+xlab("Social Index")+ylab("Neutrality")+
  scale_color_manual(values=cbbPalette)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))
d3$SI = d3$between - d3$within
d3 = d3[!is.na(d3$PS2),]
d3$PS2 = factor(d3$PS2,levels = c("cellular","eukaryote","bilaterian","insect","hymenopteran_ant"))
png("~/Writing/Figures/NurseLarva/Fig2.png",width=3000,height=3000,res=300)
p1 <- ggplot(d3[!is.na(d3$PS2),],aes(x=tissue,y=SI,fill=PS2))+geom_boxplot()+
  theme_bw()+ylab("sociality index")+
  scale_fill_manual(values=cbbPalette)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))
dev.off()

d4 <- ddply(dfAll,~PS2,summarise,
            meanWG = mean(SIWG),
            c1B = quantile(bootMean(SIWG,1000),0.025),
            c2B = quantile(bootMean(SIWG,1000),0.975),
            meanWH = mean(SIWH),
            c1W = quantile(bootMean(SIWH,1000),0.025),
            c2W = quantile(bootMean(SIWH,1000),0.975))

png("~/Writing/Figures/NurseLarva/Fig2b.png",width=3000,height=3000,res=300)
p2 <- ggplot(d4[!is.na(d4$PS2),],aes(x=meanWG,y=meanWH,color=PS2))+
  geom_point(size=3)+
  geom_errorbarh(aes(xmin=c1B,xmax=c2B))+
  geom_errorbar(aes(ymin=c1W,ymax=c2W))+
  theme_bw()+xlab("WG sociality")+ylab("WH sociality")+
  scale_color_manual(values=cbbPalette)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))
dev.off()

png("~/Writing/Figures/NurseLarva/Fig2both.png",width=3000,height=3000,res=300)
grid.arrange(p1,p2)
dev.off()

png("~/Writing/Figures/NurseLarva/Fig3a.png",width=3000,height=3000,res=300)
p1 <- ggplot(dfAll[!is.na(dfAll$PS2),],aes(x=SIWH,y=midF,color=PS2))+geom_point()+theme_bw()+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),
        legend.position = "none")+
  scale_color_manual(values=cbbPalette)+
  ggtitle("cor = 0.11, p = 0.007")
dev.off()

png("~/Writing/Figures/NurseLarva/Fig3b.png",width=3000,height=3000,res=300)
p2 <- ggplot(dfAll[!is.na(dfAll$PS2),],aes(x=SIWG,y=midF,color=PS2))+geom_point()+theme_bw()+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),legend.position = "top")+
  scale_color_manual(values=cbbPalette)+
  ggtitle("cor = 0.17, p < 0.001")
dev.off()

png("~/Writing/Figures/NurseLarva/Fig3c.png",width=3000,height=3000,res=300)
p3 <- ggplot(dfAll[!is.na(dfAll$PS2),],aes(x=SIL,y=midF,color=PS2))+geom_point()+theme_bw()+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20),
        legend.position = "none")+
  scale_color_manual(values=cbbPalette)+
  ggtitle("cor = 0.123, p = 0.002")
dev.off()

png("~/Writing/Figures/NurseLarva/Fig3all.png",width=3000,height=1500,res=300)
grid.arrange(p1,p2,p3,nrow=1)
dev.off()
####
##Fig 4: relationship between f and d
d = droplevels(dfAll[!is.na(dfAll$PS2),])
d$PS2 = factor(d$PS2, levels = c("cellular","eukaryote","bilaterian","insect","hymenopteran_ant"))
corMat = pMat = matrix(nrow=6,ncol=11)

for (j in 1:11){
  for (i in 1:5){
    d2 = d[d$PS2==levels(d$PS2)[i],]
    res= cor.test(d2$midF,d2[,c(29+j)])
    corMat[i,j] = res$estimate
    pMat[i,j] = res$p.value
  }
  res = cor.test(d$midF,d[,c(29+j)])
  corMat[6,j] = res$estimate
  pMat[6,j] = res$p.value
}

colnames(pMat) = colnames(corMat) = colnames(d2)[c(30:40)]
rownames(pMat) = rownames(corMat) = c(levels(d2$PS2),"Overall")

##Quick heatmap
pMat = t(pMat[,-c(11)])
corMat = t(corMat[,-c(11)])
corM = melt(corMat)
corP = melt(pMat)
corM$p.value = corP$value
corM$text = apply(round(corM[,c(3,4)],3),1,paste,collapse="\n")
corM$text[corM$p.value < 0.001]=gsub("\n0","\n<0.001",corM$text[corM$p.value < 0.001])
png("~/Writing/Figures/NurseLarva/Fig4heatmap.png",width=3000,height=3000,res=300)
ggplot(corM,aes(x=Var2,y=Var1))+geom_tile(aes(fill=value))+
  geom_text(aes(label=text))+
  scale_fill_gradient2(name="Pearson cor",low=muted("blue"),high=muted("red"),mid="white")+
  xlab("phylostrata")+ylab("connectivity measurement")+
  ggtitle("Correlation of f with connectivity")+
  theme_bw()+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))
dev.off()

######
##Fig 5: add in queen expression
load("~/Dropbox/monomorium nurses/data.processed/cleandata.RData")
fpkm = log(fpkm + sqrt(fpkm ^ 2 + 1)) #hyperbolic sine transformation to normalize gene expression data

meanQH = rowSums(fpkm[,grepl("QH",colnames(fpkm))])
meanQG = rowSums(fpkm[,grepl("QG",colnames(fpkm))])

expr = data.frame(QH=meanQH,QG=meanQG,Gene=names(meanQH))
d = merge(dfAll,expr,by="Gene")
d = d[!is.na(d$PS2),]
d$PS2 = factor(d$PS2, levels = c("cellular","eukaryote","bilaterian","insect","hymenopteran_ant"))

png("~/Writing/Figures/NurseLarva/Fig5a.png",width=3000,height=3000,res=300)
p1 <- ggplot(d,aes(x=QH,y=SIWH,color=PS2))+geom_point(size=2)+
  scale_color_manual(values=cbbPalette)+
  theme_bw()+
  ylab("nurse head sociality index")+xlab("queen head normalized expression")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))
dev.off()

png("~/Writing/Figures/NurseLarva/Fig5b.png",width=3000,height=3000,res=300)
p2 <- ggplot(d,aes(x=QG,y=SIWG,color=PS2))+geom_point(size = 2)+
  scale_color_manual(values=cbbPalette)+
  theme_bw()+
  ylab("nurse gaster sociality index")+xlab("queen gaster normalized expression")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))
dev.off()

png("~/Writing/Figures/NurseLarva/Fig5all.png",width=3000,height=3000,res=300)
grid.arrange(p1,p2,ncol=1)
dev.off()

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

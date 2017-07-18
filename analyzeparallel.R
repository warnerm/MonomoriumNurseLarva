cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


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

df <- read.csv("~/Data/TopExprWorkerNetNetSocialityDF.csv")

###Fig 4a
d2 <- melt(dfAll[,c(1,30:35)],id.vars="Gene")
d2$tissue="Larva"
d2$tissue[grepl("WG",d2$variable)]="WorkerGaster"
d2$tissue[grepl("WH",d2$variable)]="WorkerHead"
d2$connection.type = "within tissue"
d2$connection.type[grepl("between",d2$variable)]="between tissue"

png("~/Writing/Figures/NurseLarva/Fig1a.png",width=3000,height=3000,res=300)
ggplot(d2,aes(x=tissue,y=value,fill=connection.type))+geom_boxplot()+
  theme_bw()+xlab("tissue")+ylab("connection strength")+
  scale_fill_grey(start = 0.5, end = 0.8,name="connection type")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))
dev.off()

png("~/Writing/Figures/NurseLarva/Fig1b.png",width=3000,height=3000,res=300)
ggplot(dfAll,aes(x=WHwithin,y=WHbetween))+geom_point()+
  theme_bw()+xlab("within")+ggtitle("WH connection strength")+
  ylab("between")+geom_abline(slope=1)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))
dev.off()

png("~/Writing/Figures/NurseLarva/Fig1c.png",width=3000,height=3000,res=300)
ggplot(dfAll,aes(x=WGwithin,y=WGbetween))+geom_point()+
  theme_bw()+xlab("within")+ggtitle("WG connection strength")+
  ylab("between")+geom_abline(slope=1)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))
dev.off()

png("~/Writing/Figures/NurseLarva/Fig1d.png",width=3000,height=3000,res=300)
ggplot(dfAll,aes(x=Lwithin,y=Lbetween))+geom_point()+
  theme_bw()+xlab("within")+ggtitle("L connection strength")+
  ylab("between")+geom_abline(slope=1)+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))
dev.off()

df$SIL = df$Lbetween - df$Lwithin
df$SIWH = df$WHbetween - df$WHwithin
df$SIWG = df$WGbetween - df$WGwithin
df$SI_Overall = (df$SIL + df$SIWH + df$SIWG)/3
df$SIW = ((df$SIWH+df$SIWG)/2+df$SIL)/2

dfAll = merge(ext,df,by.x="Gene",by.y="gene")
dfAll = merge(dfAll,sn,by.x="Gene",by.y="gene")

dfAll$midF = (dfAll$f.est+dfAll$BSnIPRE.f)/2

d2 <- melt(dfAll[,c(1,30:35,27,48)],id.vars=c("Gene","PS2","midF"))
d2$tissue="Larva"
d2$tissue[grepl("WG",d2$variable)]="WorkerGaster"
d2$tissue[grepl("WH",d2$variable)]="WorkerHead"
d2$connection.type = "within"
d2$connection.type[grepl("between",d2$variable)]="between"
d3 <- dcast(d2, PS2+tissue+Gene+midF ~ connection.type)


d3$SI = d3$between - d3$within
d3 = d3[!is.na(d3$PS2),]
d3$PS2 = factor(d3$PS2,levels = c("cellular","bilaterian","eukaryote","insect","hymenopteran_ant"))
png("~/Writing/Figures/NurseLarva/Fig2.png",width=3000,height=3000,res=300)
ggplot(d3[!is.na(d3$PS2),],aes(x=tissue,y=SI,color=PS2))+geom_boxplot()+
  theme_bw()+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))
dev.off()
  
lm <- glm(SIWG ~ PS2, data=dfAll) 
summary(glht(lm, mcp(PS2="Tukey"))) 


ggplot(dfAll,aes(x=SIWH,y=SIWG,color=PS2))+geom_point()

d4 <- ddply(dfAll,~PS2,summarise,
            meanWG = mean(SIWG),
            c1B = quantile(bootMean(SIWG,1000),0.025),
            c2B = quantile(bootMean(SIWG,1000),0.975),
            meanWH = mean(SIWH),
            c1W = quantile(bootMean(SIWH,1000),0.025),
            c2W = quantile(bootMean(SIWH,1000),0.975))

png("~/Writing/Figures/NurseLarva/Fig2b.png",width=3000,height=3000,res=300)
ggplot(d4[!is.na(d4$PS2),],aes(x=meanWG,y=meanWH,color=PS2))+
  geom_point(size=3)+
  geom_errorbarh(aes(xmin=c1B,xmax=c2B))+
  geom_errorbar(aes(ymin=c1W,ymax=c2W))+
  theme_bw()+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))
dev.off()

png("~/Writing/Figures/NurseLarva/Fig3a.png",width=3000,height=3000,res=300)
ggplot(dfAll[!is.na(dfAll$PS2),],aes(x=SIWH,y=midF,color=PS2))+geom_point()+theme_bw()+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))+
  ggtitle("cor = 0.11, p = 0.007")
dev.off()

png("~/Writing/Figures/NurseLarva/Fig3b.png",width=3000,height=3000,res=300)
ggplot(dfAll[!is.na(dfAll$PS2),],aes(x=SIWG,y=midF,color=PS2))+geom_point()+theme_bw()+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))+
  ggtitle("cor = 0.17, p < 0.001")
dev.off()

png("~/Writing/Figures/NurseLarva/Fig3c.png",width=3000,height=3000,res=300)
ggplot(dfAll[!is.na(dfAll$PS2),],aes(x=SIL,y=midF,color=PS2))+geom_point()+theme_bw()+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=20))+
  ggtitle("cor = 0.123, p = 0.002")
dev.off()



png("~/Writing/Figures/NurseLarva/Fig4b.png",width=3000,height=3000,res=300)
ggplot(dfAll[!is.na(dfAll$PS2),],aes(x=SIWG,y=SIWH))+geom_point(size=3)+
  theme_bw()+xlab("Sociality Index (WG)")+ylab("Sociality Index (WH)")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20))
dev.off()

bootMean <- function(data,boots){
  res = c()
  for (i in 1:boots){
    d = sample(data,length(data),replace=TRUE)
    res = c(res,mean(d,na.rm=TRUE))
  }
  return(res)
}

SumSI <- function(dfAll){
  d = ddply(dfAll,~PS2,
            meanSI = mean(bootMean(eval(parse(text=SI)),1000)),summarise,
            c1SI = quantile(bootMean(eval(parse(text=SI)),1000),0.025),
            c2SI = quantile(bootMean(eval(parse(text=SI)),1000),0.975),
            meanF = mean(midF),
            c1F = quantile(bootMean(midF,1000),0.025),
            c2F = quantile(bootMean(midF,1000),0.975))
  return(d)
}

SI = "SI_Overall" ##ddply can't evaluate a passed variable, so have to define it globally
d <- SumSI(dfAll)[1:5,] #remove "NA" PS2

png("~/Writing/Figures/NurseLarva/Fig4c.png",width=3000,height=3000,res=300)
ggplot(d[1:5,],aes(x=meanSI,y=meanF,color=PS2))+geom_point(size=4,shape=1)+
  geom_errorbarh(aes(xmin=c1SI,xmax=c2SI),size=1.5)+
  geom_errorbar(aes(ymin=c1F,ymax=c2F),size=1.5)+
  ylab("% neutral loci")+
  xlab("sociality index, workers")+
  theme_bw()+
  scale_colour_manual(values=cbbPalette,name="PS2")+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=20),
        legend.text=element_text(size=15),legend.title=element_text(size=20))+
  geom_point(data=dfAll,aes(x=SIW,y=midF),size=1)
dev.off()

png("~/Writing/Figures/NurseLarva/Fig4d.png",width=3000,height=3000,res=300)
ggplot(d[1:5,],aes(x=meanSI,y=meanF,color=PS2))+geom_point(size=8,shape=1,stroke=2)+
  geom_errorbarh(aes(xmin=c1SI,xmax=c2SI),size=1)+
  geom_errorbar(aes(ymin=c1F,ymax=c2F),size=1)+
  ylab("mean % neutral loci")+
  xlab("mean sociality index, workers")+
  theme_bw()+
  scale_colour_manual(values=cbbPalette,name="PS2")+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=20),
        legend.text=element_text(size=15),legend.title=element_text(size=20))
dev.off()

####EdgeR
library(edgeR)
load("~/Dropbox/monomorium nurses/data.processed/cleandata.RData")
EdgeR <- function(data,design,coef){
  data <- DGEList(counts=data)
  data <- calcNormFactors(data)
  dat <- estimateGLMTrendedDisp(data, design)
  dat <- estimateGLMTagwiseDisp(dat, design)
  fit <- glmFit(dat,design)
  diff <- glmLRT(fit, coef=coef) 
  out <- topTags(diff,n=Inf,adjust.method="BH")$table 	###Calculate FDR
  return(out)
}

logFCconn <- function(code){
  df <- counts[,grepl(code,colnames(counts))]
  f <- factors[rownames(factors) %in% colnames(df),]
  design <- model.matrix(~Caste,data=droplevels(f))
  a <- EdgeR(df,design,2)
  a$Gene = rownames(a)
  d <- merge(dfAll,a,by="Gene")
  ggplot(d,aes(x=SIWH,y=logFC))+geom_point()
}
fpkm = log(fpkm + sqrt(fpkm ^ 2 + 1)) #hyperbolic sine transformation to normalize gene expression data

meanQH = rowSums(fpkm[,grepl("QH",colnames(fpkm))])
meanQG = rowSums(fpkm[,grepl("QG",colnames(fpkm))])

expr = data.frame(QH=meanQH,QG=meanQG,Gene=names(meanQH))
d = merge(dfAll,expr,by="Gene")
lm <- glm(log(f.est) ~ PS2 + QH + QG + WHwithin + WGwithin+
            WHbetween + WGwithin +
            Lwithin + Lbetween,data=d)

lm <- glm(log(f.est) ~ PS2 + QH + QG + SIWH + SIWG + SIL,data=d)
drop1(lm,.~.,test="Chi") 
summary(glht(lm, mcp(PS2="Tukey"))) 

lm <- glm(log(f.est) ~  SIWH + SIWG,data=d)
lm <- glm(SIWG ~ PS2,data=d)

dfAllH <- merge(dfAll,df,)

df <- counts[,grepl("QG|WG",colnames(counts))]
f <- factors[rownames(factors) %in% colnames(df),]
design <- model.matrix(~Caste,data=droplevels(f))
b <- EdgeR(df,design,2)

d2 <- melt(dfAll[,c(1,30:35,27,48)],id.vars=c("Gene","PS2","midF"))
d2$tissue="Larva"
d2$tissue[grepl("WG",d2$variable)]="WorkerGaster"
d2$tissue[grepl("WH",d2$variable)]="WorkerHead"
d2$connection.type = "within"
d2$connection.type[grepl("between",d2$variable)]="between"
d3 <- dcast(d2, PS2+tissue+Gene+midF ~ connection.type)

ggplot(d3, aes(x=within,y=between,color=PS2,shape=tissue))+
  geom_point()

d4 <- ddply(d3,~PS2+tissue,summarise,
            meanB = mean(between),
            c1B = quantile(bootMean(between,1000),0.025),
            c2B = quantile(bootMean(between,1000),0.975),
            meanW = mean(within),
            c1W = quantile(bootMean(within,1000),0.025),
            c2W = quantile(bootMean(within,1000),0.975))

ggplot(d4[!is.na(d4$PS2),],aes(x=meanW,y=meanB,color=PS2,shape=tissue))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=c1B,ymax=c2B))+
  geom_errorbarh(aes(xmin=c1W,xmax=c2W))+
  theme_bw()

d3$SI = d3$between - d3$within
ggplot(d3[!is.na(d3$PS2),],aes(x=tissue,y=SI,color=PS2))+geom_boxplot()
ggplot(d3[!is.na(d3$PS2),],aes(x=SI,y=midF,color=PS2,shape=tissue))+geom_point()

d4 <- ddply(d3,~PS2+tissue,summarise,
            meanSI = mean(SI),
            c1B = quantile(bootMean(SI,1000),0.025),
            c2B = quantile(bootMean(SI,1000),0.975),
            meanF = mean(midF),
            c1W = quantile(bootMean(midF,1000),0.025),
            c2W = quantile(bootMean(midF,1000),0.975))

ggplot(d4[!is.na(d4$PS2),],aes(x=meanSI,y=meanF,color=PS2,shape=tissue))+
  geom_point(size=3)+
  geom_errorbarh(aes(xmin=c1B,xmax=c2B))+
  geom_errorbar(aes(ymin=c1W,ymax=c2W))+
  theme_bw()

dOld = d[d$PS2=="cellular",]



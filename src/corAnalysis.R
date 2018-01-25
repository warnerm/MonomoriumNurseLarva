#File takes output of 'corApproach.R' and 'getGeneType.R' for further analysis
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

PS_palette = c("white","#deebf7","#9ecae1","#4292c6","#2171b5","#084594")
tissue_palette = c("#9e9ac8","#6a51a3","cornflowerblue","royalblue4")

f.est <- read.csv("~/Data/Nurse_Larva/MKtestConstraintOneAlpha.csv")
colnames(f.est) = c("Gene","f")
ext <- read.csv("~/Downloads/msx123_Supp (1)/MpharAnn.csv") #load in MBE results

#General function for saving plots
saveP <- function(p,name,w,h){
  filename = paste("~/Writing/Figures/NurseLarva/corApproach/",name,".png",sep="")
  png(filename,height = h,width = w,res =300)
  grid.arrange(p)
  dev.off()
  filename = paste("~/Writing/Figures/NurseLarva/corApproach/",name,".pdf",sep="")
  pdf(filename)
  grid.arrange(p)
  dev.off()
}

theme_all = theme(axis.text=element_text(size=20),
                  axis.title=element_text(size=25),
                  legend.text=element_text(size=20),
                  legend.title=element_text(size=25),
                  plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
                  axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)),
                  axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))

theme4 = theme(axis.text=element_text(size=13),
                  axis.title=element_text(size=17),
                  legend.text=element_text(size=13),
                  legend.title=element_text(size=17),
                  plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
                  axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)),
                  axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))

load("socDet.RData")
load("corResults.RData")
#Number of gene-gene social connections
corStat <- function(df){
  N = nrow(df)*ncol(df)
  nS = sum(df=="negSocial")
  pS = sum(df=="posSocial")
  nC = sum(df=="negControl")
  pC = sum(df=="posControl")
  d = c(N = N, nS = nS, pS = pS, nC = nC, pC = pC,tS = nS + pS, tC = nC + pC)
  rates = d[2:7]/d[1]
  names(rates) = paste("rate",names(rates),sep="")
  return(c(d,rates))
}


#Number of significant correlation with larval profiles of nurse DEGs by larval stage
for (i in 1:length(socDet)){
  socDet[[i]]$geneType = apply(socDet[[i]][,c('Pvalue','Excess')],1,function (x) {
    if (!is.na(x['Pvalue'])){
      if (x['Pvalue'] < 0.05/nrow(socDet[[i]])){
        if (x['Excess'] > 0){
          return("social")
        } else {
          return("control")
        }
      } 
    }
    return("NS")
  })
}
numbers <- lapply(socDet,function (x) c(sum(x$geneType!="social"),sum(x$geneType=="social")))

# numbers <- numbers[!grepl('LARV',names(numbers))] #Remove larval numbers
# num = ldply(numbers)
# m <- melt(num,id.vars=".id")
# levels(m$variable) = c('not correlated','correlated')
# 
# p <- ggplot(m,aes(x=.id,y=value,fill=variable))+
#   geom_bar(stat="identity")+
#   scale_fill_manual(name="expression profile correlated with larvae?",values = Soc_palette)+theme_bw()+
#   ylab("number of genes")+
#   xlab("tissue")
# 
# png("~/Writing/Figures/NurseLarva/corApproach/NumGenesBinom.png",width=4000,height=4000,res=300)
# p + theme_all+theme(legend.position = 'bottom',plot.title = element_text(size=20))+
#   ggtitle('correlation of larval and nurse expression profiles,\nfor genes differentially expressed in nurses by larval stage')
# dev.off()

#Overlap of social genes for each of the three

#Numbers of negative and positive connections
PlotPosNeg <- function(x,label){
  p <- ggplot(x[x$geneType=="social",],aes(x=negSocial,y=posSocial))+
    geom_point(size=1.5,alpha = 0.5)+theme_bw()+
    ylab("number positive")+xlab("number negative")+
    theme4+annotate("text",x = max(x$negSocial)/2, y = max(x$posSocial),label=label,size=12)
  return(p)
}

#Do genes mostly exhibit positive or negative correlations?
# sigPosNeg <- function(x){
#   numSig <- sum(apply(x[x$geneType=="social",],1,function(i){
#     if (binom.test(x = as.integer(i['posSocial']), n = as.integer(i['posSocial']) + as.integer(i['negSocial']),p = 0.5)$p.value < 0.05) 1
#     else 0
#   }))
#   return(c(numSig,nrow(x),numSig/nrow(x[x$geneType=="social",])))
# }
# 
# binom.test(x = as.integer(i['posSocial']), n = as.integer(i['posSocial'] + i['negSocial']),p = 0.5)
# 
# numDiff <- lapply(socDet,sigPosNeg)

label = c("focal nurse head","focal nurse abdomen","larvae (head)","larvae (abdomen)")
plots <- mapply(PlotPosNeg,socDet[c(3,4,7,8)],label,SIMPLIFY = FALSE)
a = do.call("grid.arrange", c(plots, ncol=2))
saveP(a,"PosNegNumbers",4000,4000)
# 
# ########
# ##Various ways to compare overlap in sociality definition of genes
# ########
# socGenes <- lapply(socDet,function(x) x$Gene[x$geneType=="social"])
# 
# GeneOverlap <- function(g1,g2){
#   over = sum(g1 %in% g2)
#   d = rbind(c(over,length(g1) - over),
#             c(length(g2) - over,nrow(fpkm) - length(g2) -length(g1) + over))
#   test = fisher.test(d)
#   return(list(d,test))
# }
# 
# a = lapply(c(3,4,7,8),function(i) lapply(c(3,4,7,8),function(j){
#   if (i!=j) GeneOverlap(socGenes[[i]],socGenes[[j]])
#   else NA
# }))
# 
# ###Correlation between number of social connections for a given gene
# corSoc <- function(i,j){
#   d1 = socDet[[i]]
#   d2 = socDet[[j]]
#   d1 = d1[d1$Gene %in% d2$Gene,]
#   d2 = d2[d2$Gene %in% d1$Gene,]
#   return(cor.test(d1$posSocial + d1$negSocial,d2$posSocial+d2$negSocial,method = "spearman"))
# }
# 
# socCor <- lapply(c(3,4,7,8),function(i){
#   lapply(c(3,4,7,8),function(j){
#     if (i!=j) corSoc(i,j)
#     else NA
#   })
# })
# 
# ##Venn diagram of social genes
# names <- c("CH","CG","LARV_CH","LARV_CG")
# soc <- lapply(names, function(x) socDet[[x]]$Gene[socDet[[x]]$geneType=="social"])
# soc2 = list(soc[[1]],soc[[2]],soc[[3]][soc[[3]] %in% soc[[4]]])
# 
# venn.diagram(soc2,category.names = names[1:3],filename="~/Writing/Figures/NurseLarva/corApproach/Venn.png")

#Using connection strengths calculated by Genie
makePlot <- function(nurse){
  df <- read.csv(paste(nurse,"GenieTabConn.csv",sep=""))
  df$tissue = gsub("_.*","",df$Gene)
  df$Gene = gsub(".*_","",df$Gene)
  df = df[,-c(1)]
  p1 <- ggplot(df,aes(x = reg_within, y = reg_between, color = tissue))+
    geom_point(alpha = 0.3)+
    geom_smooth(method = "lm")+ggtitle(paste(nurse, "regulatory"))
  p2 <- ggplot(df,aes(x = targ_within, y = targ_between, color = tissue))+
    geom_point(alpha = 0.3)+
    geom_smooth(method = "lm")+ggtitle(paste(nurse, "target"))
  df$code = nurse
  return(list(p1,p2,df))
}

codes <- c("CH","CG","RH","RG")
plotsG <- lapply(codes,makePlot)

allD <- ldply(lapply(seq(1,length(codes)), function(i) plotsG[[i]][[3]]),data.frame)
m <- melt(allD,id.vars = c("Gene","tissue","code"))

m$variable = factor(m$variable,levels = levels(m$variable)[c(1,3,2,4)])
m = droplevels(m[grepl("reg",m$variable),])
m$tissueC = apply(m[,c('tissue','code')],1,paste,collapse="_")
levels(m$variable) = c("within-tissue","social")
m$tissueC = as.factor(m$tissueC)
levels(m$tissueC) = c("larva (abd)","larva (head)","larva (random abd)","larva (random head)","nurse (abd)","nurse (head)","nurse (random abd)","nurse (random head)")

p1 <- ggplot(m[grepl("C",m$code),],aes(x = variable, y = value, fill = tissueC))+
  geom_boxplot(notch = TRUE)+
  theme_bw()+
  scale_y_log10()+
  ylab("total connectivity")+
  xlab("connection type")+
  scale_fill_manual(values = tissue_palette,name = "tissue")+theme4

png("~/Writing/Figures/NurseLarva/corApproach/connStrength.png",height = 4000,width = 4000,res =300)
p1 + theme_all
dev.off()

socDir <- function(code,df){
  socW = socDet[[code]]
  socL = socDet[[paste("LARV_",code,sep="")]]
  l = df[[3]][df[[3]]$tissue == "larv",]
  n = df[[3]][df[[3]]$tissue == "nurse",]
  nurse = merge(socW,n)
  larv = merge(socL,l)
  cN = cor(nurse[,c(2:7,9:12)])
  cL = cor(larv[,c(2:7,9:12)])
  return(rbind(nurse,larv))
}

dirCH = socDir("CH",plotsG[[1]])
dirCG = socDir("CG",plotsG[[2]])
dir = rbind(dirCH,dirCG)
dir$tissueC = as.factor(apply(dir[,c('tissue','code')],1,paste,collapse="_"))
levels(dir$tissueC) = c("larva (abd)","larva (head)","nurse (abd)","nurse (head)")

regPlots <- lapply(levels(dir$tissueC),function(x){
  list(
    ggplot(dir[dir$tissueC==x,],aes(x = posSocial+1, y = reg_between))+
      geom_point(alpha = 0.5)+
      theme_bw()+theme4+scale_x_log10()+
      ggtitle(x)+ylab("social regulatory strength")+
      xlab("positive social conns")+
      geom_smooth(method = 'lm',se=FALSE),
    ggplot(dir[dir$tissueC==x,],aes(x = negSocial+1, y = reg_between))+
      geom_point(alpha = 0.5)+
      theme_bw()+theme4+scale_x_log10()+
      ggtitle(x)+ylab("social regulatory strength")+
      xlab("negative social conns")+
      geom_smooth(method = 'lm',se=FALSE)
  )
})

targPlots <- lapply(levels(dir$tissueC),function(x){
  list(
    ggplot(dir[dir$tissueC==x,],aes(x = posSocial+1, y = targ_between))+
      geom_point(alpha = 0.5)+
      theme_bw()+theme4+scale_x_log10()+
      ggtitle(x)+ylab("social target strength")+
      xlab("positive social conns")+
      geom_smooth(method = 'lm',se=FALSE),
    ggplot(dir[dir$tissueC==x,],aes(x = negSocial+1, y = targ_between))+
      geom_point(alpha = 0.5)+
      theme_bw()+theme4+scale_x_log10()+
      ggtitle(x)+ylab("social target strength")+
      xlab("negative social conns")+
      geom_smooth(method = 'lm',se=FALSE)
  )
})

png("~/Writing/Figures/NurseLarva/corApproach/regConns.png",height = 4000,width = 2000,res =300)
do.call(grid.arrange,c(unlist(regPlots,recursive=FALSE),ncol=2))
dev.off()

png("~/Writing/Figures/NurseLarva/corApproach/targConns.png",height = 4000,width = 2000,res =300)
do.call(grid.arrange,c(unlist(targPlots,recursive=FALSE),ncol=2))
dev.off()

corHeat <- function(d,tissue){
  d = d[[3]][d[[3]]$tissue == tissue,]
  c = cor(d[,c(2:5)])
  m = melt(c)
  lev = levels(m$Var1)
  m$Var1 = factor(m$Var1,levels = lev[c(1,3,2,4)])
  m$Var2 = factor(m$Var2,levels=rev(levels(m$Var1)))
  levels(m$Var1) = c("within-tissue regulatory","within-tissue target","social regulatory","social target")
  levels(m$Var2) = c("social target","social regulatory","within-tissue target","within-tissue regulatory")
  
  #Make the matrix upper-triangular
  m$value[as.integer(m$Var2) + as.integer(m$Var1) < 5] = NA
  p <- ggplot(m,aes(x = Var2, y = Var1, fill = value))+
    geom_tile()+
    xlab("")+
    ylab("")+theme_bw()+
    ggtitle(tissue)+theme4+
    scale_fill_gradient2(limits = c(-1,1),low = muted("red"), mid = "white", high = muted("blue"))+
    theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))
  return(p)
}
plotsCor <- lapply(plotsG[c(1,2)], function(i) list(corHeat(i,"nurse"),corHeat(i,"larv")))

png("~/Writing/Figures/NurseLarva/corApproach/connHeats.png",height = 3000,width = 4000,res =300)
grid.arrange(plotsCor[[1]][[1]]+ggtitle('nurse (head)'),
             plotsCor[[1]][[2]]+ggtitle('larva (head)'),
             plotsCor[[2]][[1]]+ggtitle('nurse (abd)'),
             plotsCor[[2]][[2]]+ggtitle('larva (abd)'),ncol = 2)
dev.off()

allDf <- merge(dir,f.est,by = "Gene",all.x = TRUE)
allDf$reg_diff = allDf$reg_between - allDf$reg_within
allDf$targ_diff = allDf$targ_between - allDf$targ_within
allDf$tissueC = as.factor(apply(allDf[,c('tissue','code')],1,paste,collapse="_"))
levels(allDf$tissueC) = c("larva (abd)","larva (head)","nurse (abd)","nurse (head)")

p1 <- ggplot(allDf,aes(x = reg_between, y = 1 - log(f),color = tissueC))+
  geom_point(alpha = 0.5)+
  theme_bw()+
  ylab("constraint")+
  xlab("social regulatory strength")+
  scale_color_discrete(name = 'tissue')+
  geom_smooth(method = 'lm',se =FALSE)+theme_all

p2 <- ggplot(allDf,aes(x = targ_between, y = 1 - log(f),color = tissueC))+
  geom_point(alpha = 0.5)+
  theme_bw()+
  ylab("constraint")+
  xlab("social target strength")+
  scale_color_discrete(name = 'tissue')+
  geom_smooth(method = 'lm',se =FALSE)+theme_all

png("~/Writing/Figures/NurseLarva/corApproach/constraint_social.png",height = 4000,width = 3000,res =300)
grid.arrange(p1,p2)
dev.off()

p1 <- ggplot(allDf,aes(x = reg_within, y = 1- log(f),color = tissueC))+
  geom_point(alpha = 0.5)+
  theme_bw()+
  ylab("constraint")+
  xlab("within-tissue regulatory strength")+
  scale_color_discrete(name = 'tissue')+
  geom_smooth(method = 'lm',se =FALSE)+theme_all

p2 <- ggplot(allDf,aes(x = targ_within, y = 1-log(f),color = tissueC))+
  geom_point(alpha = 0.5)+
  theme_bw()+
  ylab("constraint")+
  xlab("within-tissue target strength")+
  scale_color_discrete(name = 'tissue')+
  geom_smooth(method = 'lm',se =FALSE)+theme_all

png("~/Writing/Figures/NurseLarva/corApproach/constraint_within.png",height = 4000,width = 3000,res =300)
grid.arrange(p1,p2)
dev.off()

ggplot(allDf,aes(x = tissueC, y = log(f)))+
  geom_boxplot(notch = TRUE)

lm1 <- glm(log(f) ~ tissueC, data = allDf)
av <- aov(lm1)
summary(av)
summary(lm1)
tukey.test <- TukeyHSD(av,"tissueC")


lm1 <- glm(reg_diff ~ PS2, data = allP[allP$code == "CG" & allP$tissue == "nurse",])
av <- aov(lm1)
summary(av)
tukey.test <- TukeyHSD(av,"PS2")

socF <- function(code){
  fSoc = f.est
  fSoc$social = "non-social"
  fSoc$social[fSoc$Gene %in% socDet[[code]]$Gene[socDet[[code]]$geneType=="social"]]="social"
  fSoc$code = code
  return(fSoc)
}

Fsoc <- ldply(lapply(names(socDet),socF))
ggplot(Fsoc,aes(x = code, y = log(f), fill = social))+
  geom_boxplot(notch = TRUE)

lm <- glm(log(f) ~ social, data = Fsoc[Fsoc$code=="CH",])
av <- aov(lm)
summary(av)

code = "LARV_RH"
wilcox.test(Fsoc$f[Fsoc$code==code & Fsoc$social!="social"],Fsoc$f[Fsoc$code==code &Fsoc$social=="social"])

##########
##Phylostrata
##########
inf = table(ext$Raw.PS)
i = data.frame(name = names(inf),number = inf)
a=taxonomy("Monomorium pharaonis")
a$n = seq(1,nrow(a))
ps = merge(a,i, by = "name",all.x = TRUE,all.y=TRUE,sort=FALSE)
ps = ps[order(ps$n),]
all <- merge(allDf,ext,by="Gene",all.x=TRUE)
allP = all[!is.na(all$Raw.PS),]

allP$PS2 = factor(allP$PS2,levels = c("cellular",'eukaryote','bilaterian','insect','hymenopteran_ant'))
allP$PS2[allP$Raw.PS=="Neoptera"]="insect"

p1 <- ggplot(allP,aes(x = tissueC, y = reg_within,fill=PS2))+
  geom_boxplot(notch = TRUE)+
  theme_bw()+theme_all+ylab("within-tissue regulatory strength")+
  xlab("tissue")

p2 <- ggplot(allP,aes(x = tissueC, y = targ_within,fill=PS2))+
  geom_boxplot(notch = TRUE)+
  theme_bw()+theme_all+ylab("within-tissue target strength")+
  xlab("tissue")

png("~/Writing/Figures/NurseLarva/corApproach/psWithin.png",height = 4000,width = 4000,res =300)
grid.arrange(p1,p2)
dev.off()

p1 <- ggplot(allP,aes(x = tissueC, y = reg_between,fill=PS2))+
  geom_boxplot(notch = TRUE)+
  theme_bw()+theme_all+ylab("social regulatory strength")+
  xlab("tissue")

p2 <- ggplot(allP,aes(x = tissueC, y = targ_between,fill=PS2))+
  geom_boxplot(notch = TRUE)+
  theme_bw()+theme_all+ylab("social target strength")+
  xlab("tissue")

png("~/Writing/Figures/NurseLarva/corApproach/psSocial.png",height = 4000,width = 4000,res =300)
grid.arrange(p1,p2)
dev.off()


sum <- ddply(allP,.(Gene,PS2,f.est),summarize,
             regB = mean(reg_between),
             targB = mean(targ_between),
             regW = mean(reg_within),
             targW = mean(targ_within),
             N = length(reg_between))

lm <- glm(targB ~ PS2,data = sum)
drop1(lm,.~.,test="Chi") 
av <- aov(lm)
TukeyHSD(av,"PS2")

lm <- glm(1-log(f.est) ~ reg_within + reg_between + targ_within + targ_between + PS2, data = allP)
lm <- glm(1-log(f.est) ~ regB + targB + regW + targW, data = sum)


ggplot(allP, aes( x = targ_between - targ_within, y = 1 -f.est))+geom_point()+geom_smooth(method = "lm")


drop1(lm,.~.,test="Chi") 

av <- aov(lm)
summary(av)

sum$regDiff = sum$regB - sum$regW
sum$targDiff = sum$targB - sum$targW

sumM = melt(sum[,c(1,2,4,5,6,7)],id.vars = c("Gene","PS2"))
sumM$variable = factor(sumM$variable,levels = levels(sumM$variable)[c(3,4,1,2)])
levels(sumM$variable) = c("within-tissue regulatory","within-tissue target","social regulatory","social target")
png("~/Writing/Figures/NurseLarva/corApproach/psPooled2.png",height = 4000,width = 4000,res =300)
ggplot(sumM,aes(x = variable, y = value, fill = PS2))+geom_boxplot(notch = TRUE)+
  ylab("overall connection strength")+
  xlab("Type")+theme_bw()+theme_all+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

sumM = melt(sum[,c(1,2,9,10)],id.vars = c("Gene","PS2"))
levels(sumM$variable) = c("regulatory","target")
png("~/Writing/Figures/NurseLarva/corApproach/psPooled.png",height = 4000,width = 4000,res =300)
ggplot(sumM,aes(x = variable, y = value, fill = PS2))+geom_boxplot(notch = TRUE)+
  ylab("social minus within-tissue strength")+
  xlab("connection type")+theme_bw()+theme_all+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

png("~/Writing/Figures/NurseLarva/corApproach/psConstraint.png",height = 4000,width = 4000,res =300)
ggplot(sum,aes(x = PS2, y = 1 - log(f.est),fill = PS2))+
  geom_boxplot(notch = TRUE)+
  theme_bw()+
  theme_all+
  ylab("constraint")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

lm1 <- glm(log(f.est) ~ regW + targW + regB + targB, data = sum)

av <- aov(lm1)
summary(av)
summary(lm1)
TukeyHSD(av,"PS2")



###Module composition of social genes/social connections

#Load larval module definitions
df <- read.table("findK_clusterW_L.txt",header=TRUE)
mods = t(df[11,]) #12 modules results in a high SIL

fpkm <- read.csv("fpkm.csv") #genes are indexed as in fpkm file
rownames(mods)=fpkm$X

#List of mediods
mediods = unique(as.vector(mods))

#List of genes contained in each module
larvMembership <- lapply(mediods,function(x) rownames(mods)[mods == x])

#From assignNurseMods.py, files containing module assignments for nurse genes

#Get module membership of nurse genes
getModMemb <- function(code){
  #First row is true values, subsequent values are derived by scrambling stage codes
  dfH <- read.table(paste("~/Data/Nurse_Larva/sortMods",code,".txt",sep=""),head=TRUE)
  colnames(dfH) = fpkm$X
  
  #True assignments
  Ht <- t(dfH[1,])
  names(Ht) = names(dfH)
  Hlist <- lapply(mediods,function(x) names(Ht)[Ht == x])
  Hd = socMod(code,Ht)
  return(list(Hd,Hlist))
}

socMod <- function(code,modN){
  soc = socDet[[code]]
  dT = data.frame(modN)
  dT$Gene = rownames(dT)
  names(dT)[1] = "Module"
  d = merge(soc,dT,by="Gene",all.y = TRUE)
  d$geneType[is.na(d$geneType)]="NS"
  d$geneType[d$geneType!="social"]="NS"
  return(d)
}


mosaicMod <- function(modN){
  mPallete <- c("red","red","lightcoral","white","white","steelblue1","blue","blue")
  resid = chisq.residuals(table(modN$Module,modN$geneType),std=TRUE)
  r = c(matrix(c(resid[,1], resid[,2]), 2, byrow = T))
  rf = rep(1,length(r))
  rf[r<8]=2
  rf[r<4]=3
  rf[r<2]=4
  rf[r<0]=5
  rf[r< -2]=6
  rf[r< -4]=7
  rf[r< -8]=8
  sT = ddply(modN,.(Module),summarize,
             Freq = length(Module))
  s = ddply(modN,.(Module,geneType),summarize,
            Freq = length(Module))
  s$Mod = apply(s,1,function(x) modO$order[modO$N.mods==x['Module']])
  s$pearson = rf
  s$pearson = factor(s$pearson,levels=c("1","2","3","4","5","6","7","8"))
  levels(s$pearson) = c("8","4","2",""," ","-2","-4","-8")
  
  p <- ggplot(data = s)+
    geom_mosaic(aes(weight = Freq, x = product(geneType,Mod), fill = pearson))+
    xlab("module number")+
    scale_fill_manual(values = mPallete)+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
  return(p)
}

codes <- c("CH","CG","RH","RG")
ModMemb <- lapply(codes,getModMemb)
Lmod <- lapply(codes,function(x) socMod(paste("LARV_",x,sep=""),mods))
LmodCH = socMod("LARV_CH",mods)
LmodCG = socMod("LARV_CG",mods)
chi <- lapply(ModMemb,function(x) chisq.test(table(x[[1]]$Module,x[[1]]$geneType)))

#Order modules by number of larval genes in each
tbl = table(mods)
modO = data.frame(N = tbl)
modO = modO[order(modO$N.Freq,decreasing=TRUE),]
modO$order = seq(1,nrow(modO))
mosaics = lapply(ModMemb[1:3],function(x) mosaicMod(x[[1]]))
mosaicL = lapply(list(LmodCH,LmodCG),mosaicMod)

png("~/Writing/Figures/NurseLarva/corApproach/mosaics.png",height = 3000,width = 4400,res =300)
grid.arrange(mosaics[[1]]+ggtitle("nurse (head)"),
             mosaics[[2]]+ggtitle("nurse (abd)"),
             mosaicL[[1]]+ggtitle("larva (head)"),
             mosaicL[[2]]+ggtitle("larva (abd)"))
dev.off()

modTbl <- lapply(list(ModMemb[[1]][[1]],ModMemb[[2]][[1]],LmodCH,LmodCG),
                      function (x) table(x$geneType,x$Module))

socRate <- ldply(lapply(modTbl,function(x) x[2,]/(x[2,]+x[1,])))

rownames(socRate) = c("CH","CG","LARV_CH","LARV_CG")
c = cor(t(socRate),method = "spearman")
m = melt(c)

png("~/Writing/Figures/NurseLarva/corApproach/mosaicCor.png",height = 3000,width = 4400,res =300)
ggplot(m,aes(x = Var2, y = Var1, fill = value))+
  geom_tile()+
  xlab("")+
  ylab("")+theme_bw()+theme4+
  scale_fill_gradient2(limits = c(-1,1),low = muted("red"), mid = "white", high = muted("blue"))+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))
dev.off()

#Mean connection strength in each module
modGconn <- function(modN,g){
  g = gAll[[g]]
  d = merge(modN,g,by="Gene")
  d$Mod = apply(d,1,function(x) modO$order[modO$N.mods==x['Module']])
  dM = melt(d[,c(1,10:13,16)],id.vars = c("Gene","Mod"))
  p <- ggplot(dM,aes(x = variable, fill = as.factor(Mod), y = value))+
    geom_boxplot(notch = TRUE)+
    theme_bw()
  return(p)
}

allModConn <- mapply(modGconn,modN = list(ModMemb[[1]][[1]],ModMemb[[2]][[1]],LmodCH,LmodCG),g = c(1,3,2,4),SIMPLIFY = FALSE)
png("~/Writing/Figures/NurseLarva/corApproach/modConnH.png",height = 4000,width = 4400,res =300)
grid.arrange(allModConn[[1]]+ggtitle("nurse (head)"),
             allModConn[[3]]+ggtitle("larva (head)"))
dev.off()

png("~/Writing/Figures/NurseLarva/corApproach/modConnG.png",height = 4000,width = 4400,res =300)
grid.arrange(allModConn[[2]]+ggtitle("nurse (abd)"),
             allModConn[[4]]+ggtitle("larva (abd)"))
dev.off()
#Number and rate of pos/neg connections in each module
plotModConn <- function(modM){
  d = modM[modM$geneType=="social",c("Gene","posSocial","negSocial","Module")]
  m = melt(d,id.vars = c("Gene","Module"))
  m$modN = apply(m,1,function(x) modO$order[modO$N.mods==x['Module']])
  s = ddply(d,.(Module),summarize,
            Nneg = sum(negSocial),
            Npos = sum(posSocial))
  sm = melt(s,id.vars = "Module")
  sm$modN = apply(sm,1,function(x) modO$order[modO$N.mods==x['Module']])
  p1 <- ggplot(m,aes(x = as.factor(modN),y = value, fill = variable))+
    geom_boxplot(notch = TRUE)
  p2 <- ggplot(sm,aes(x = as.factor(modN),y = value, fill = variable))+
    geom_bar(stat = "identity")
  return(list(p1,p2))
}

modConn <- lapply(modMall,plotModConn)

########
##GSEA/GO term analysis
########

#Full list of module membership
modMall <- c(lapply(ModMemb,function(x) x[[1]]),Lmod)

#Same order (alternating nurse and larvae as below)
modMall = modMall[c(1,5,2,6,3,7,4,8)]

#Mediod order
mediods = modO$N.mods
  
#Get list of all genie statistics
genieAll <- lapply(plotsG,function(x){
  list(x[[3]][x[[3]]$tissue=="nurse",],
       x[[3]][x[[3]]$tissue=="larv",])
})
gAll = unlist(genieAll,recursive=FALSE)

  
##Set-up GO analysis
#Previously derived GO annotations
go <- read.csv("~/Writing/Data/NurseSpecialization_transcriptomicData/GOannotation.csv")
new <- list()
for (gene in unique(go$gene)){
  d = go[go$gene %in% gene,]
  new[[gene]]=as.character(d$GO)
}

#Not necessary for GSEA, but necessary for identifying significant genes
selectBool <- function(score){
  return(score == 1)
}

selectConn <- function(score){
  return(score > quantile(score,0.9))
}

#GO analysis. Takes in vector of all genes, with a score (use 0 or 1)
GOfunc <- function(d){
  GOdata <- new("topGOdata",
                description="Simple session",ontology="BP",
                allGenes=d,geneSel=selectBool,
                nodeSize = 10,
                annot=annFUN.gene2GO,gene2GO=new)
  resultFS <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  allRes <- GenTable(GOdata,FS=resultFS,pvalCutOff = 0.05)
  return(allRes) 
}

#GSEA analysis. Takes in vector of all genes, with a score for social connection strength
GSEAfunc <- function(d){
  GOdata <- new("topGOdata",
                description="Simple session",ontology="BP",
                allGenes=d,geneSel=selectConn,
                nodeSize = 10,
                annot=annFUN.gene2GO,gene2GO=new)
  
  #Use scoreOrder = "decreasing" because the higher connection strengths are what we are after
  resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks",scoreOrder="decreasing")
  allRes <- GenTable(GOdata,KS=resultKS)
  return(allRes[allRes$KS < 0.05,]) 
}

#Generate input to topGO tests
GOinput <- function(genes){
  d = rep(0,nrow(fpkm))
  names(d) = fpkm$X
  d[names(d) %in% genes]=1
  return(d)
}

#Get top result based on p value
GOtop <- function(res){
  res = res[res$Significant > res$Expected,]
  return(res[1,2])
}

#GO terms by module, highlighting social modules
#Return top GO term for each module 
GOmod <- function(modM){
  GOres <- lapply(mediods,function(x) GOfunc(GOinput(modM$Gene[modM$Module==x])))
  topRes <- ldply(lapply(GOres,GOtop))
  return(list(topRes,GOres))
}
  
#GO terms by module (only social genes)
#Return matrix with top GO term for each module for each tissue, where only social genes are included
GOmodSoc <- function(modM){
  GOres <- lapply(mediods,function(x) GOfunc(GOinput(modM$Gene[modM$Module==x&modM$geneType=="social"])))
  topRes <- ldply(lapply(GOres,GOtop))
  return(list(topRes,GOres))
}  

#Overall GSEA, by reg_between, reg_within, targ_between, and targ_within
#Return matrix with each of the four stats, each tissue
gsea_overall <- function(gS){
  gS[,6] = gS[,3] - gS[,2]
  gS[,7] = gS[,5] - gS[,4]
  gseaInput <- lapply(seq(2,7),function(x) gS[,x])
  for (i in 1:length(gseaInput)){
    names(gseaInput[[i]]) = gS$Gene
  }
  gseaRes <- lapply(gseaInput, GSEAfunc)
  topRes <- ldply(lapply(gseaRes,GOtop))
  return(list(topRes,gseaRes))
}



GO_Mods <- lapply(modMall,GOmod)
GO_socMods <- lapply(modMall,GOmodSoc)
gseaTissue <- lapply(gAll[c(1:4)],gsea_overall)

allGS <- do.call(cbind,lapply(gseaTissue,function(x) x[[1]]))
colnames(allGS) = c("nurse head","larva (head)","nurse abdomen","larva (abdomen)")
allGS$'connection type' = c("within-tissue regulatory","social regulatory","within-tissue target",
                    "social target","difference (regulatory)","difference (target)")
allGS = allGS[,c(5,1:4)]


tt3 <- ttheme_minimal(
  core=list(
    fg_params=list(fontsize=8)),
  colhead=list(fg_params=list(fontface="bold",fontsize=14)))

#Make table grob out of results
resT <- tableGrob(allGS,theme=tt3,rows = NULL)
grid.arrange(resT)
g <- gtable_add_grob(resT,
                     grobs = segmentsGrob( # line across the bottom
                       x0 = unit(0,"npc"),
                       y0 = unit(0,"npc"),
                       x1 = unit(1,"npc"),
                       y1 = unit(0,"npc"),
                       gp = gpar(lwd = 3.0)),
                     t = 1, b = 1, l = 1, r = 5)
separators <- replicate(ncol(g) - 1,
                        segmentsGrob(x1 = unit(0, "npc"), gp=gpar(lty=2)),
                        simplify=FALSE)
## add vertical lines on the left side of columns (after 2nd)
g <- gtable::gtable_add_grob(g, grobs = separators,
                             t = 2, b = nrow(g), l = seq_len(ncol(g)-1)+1)
p <- grid.arrange(g)
ggsave(p,file="~/Writing/Figures/NurseLarva/corApproach/TableGO.png",width=10,height=4,dpi=300)


#######
##Looking at top genie connection pairs
#######
perc.rank <- function(x) trunc(rank(x))/length(x)

topConn <- function(gS){
  g = gS[order(gS$reg_between,decreasing = TRUE),]
  gExt <- merge(g,ext[,c(1:3,19,20,25,27)],by = "Gene",sort = FALSE,all.x = TRUE)
  gExt <- merge(gExt,f.est,by.x = "Gene",by.y = "gene",sort = FALSE,all.x = TRUE)
  gExt <- gExt[order(gExt$reg_between,decreasing = TRUE),]
  gExt <- within(gExt, fR <- perc.rank(f.est))
  return(gExt)
}

topGenie <- lapply(gAll,topConn)

i = 4
a = topGenie[[i]][,c(2,3,8,9,11,12,13,14,15)]
a = a[order(a$reg_between - a$reg_within,decreasing = TRUE),]





getSocMod <- function(dfH){
  #'Bootstrapped' assignments
  dfH = dfH[-c(1),]
  
  #Identify null distribution of numbers of genes per module
  occ <- apply(dfH,1,function(x) unlist(lapply(mediods,function(y) length(colnames(dfH)[x == y]))))
  
  #null hypothesis mean based on resampling 
  null_mean <- as.data.frame(t(apply(occ,1,function(x) c(mean = mean(x),ci95 = quantile(x,0.95)))))
  colnames(null_mean)[2] = "ci95"
  
  #Social modules
  nSoc <- unlist(lapply(mediods,function(y) length(colnames(Ht)[Ht == y])))
  binom <- apply(cbind(nSoc,m = null_mean$mean),1,function(y) binom.test(x = y['nSoc'], p = y['m']/length(Ht),n = length(Ht),alternative = "greater")$p.value)
  socMod <- mediods[binom < 0.05]
  socList <- lapply(socMod,function(x) rownames(Ht)[Ht == x])
  null_mean$trueN = nSoc
  
  names(socList) = socMod
  return(list(socList,null_mean,binom))
}

#Plot null and true results for numbers of genes in each module
plotNurseMod <- function(d){
  null_mean=d[[2]]
  null_mean$Module = as.factor(seq(1,nrow(null_mean),by=1))
  m <- melt(null_mean,id.vars=c("ci95","Module"))
  m$ci95[m$variable=="trueN"]=NA #No confidence interval for true values
  levels(m$variable) = c("mean of bootstraps","true value")
  limits <- aes(ymin = 5,ymax = ci95)
  p <- ggplot(m,aes(x=Module,y=value,fill=variable))+
    geom_errorbar(limits,position = position_dodge(),size=0.5)+ #plotting error bars first makes the bottom ones disappear 
    geom_bar(stat="identity",position=position_dodge())+
    scale_fill_manual(name="expression profile correlated with larvae?",values = Soc_palette)+theme_bw()+
    ylab("number of genes")+
    xlab("tissue")
  return(p)
}

#Get social module, stats on bootstrapping
codes <- c('CH','CG','RH','RG','QCH','QCG')
nurseMods <- lapply(codes,function(x) getSocMod(x))
modPlots <- lapply(nurseMods, plotNurseMod)






go_mods = lapply(mList,GOmod)








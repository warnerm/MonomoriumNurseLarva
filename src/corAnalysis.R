#File takes output of 'corApproach.R' and 'getGeneType.R' for further analysis
setwd("~/Data/Nurse_Larva/")

library(plyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(edgeR)
library(topGO)

PS_palette = c("white","#deebf7","#9ecae1","#4292c6","#2171b5","#084594")
Soc_palette = c("#dadaeb","#9e9ac8","#6a51a3")


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

#Number of significant correlation with larval profiles of nurse DEGs by larval stage
numbers <- lapply(socDet,function (x) c(sum(x$geneType!="social"),sum(x$geneType=="social")))
numbers <- numbers[!grepl('LARV',names(numbers))] #Remove larval numbers
num = ldply(numbers)
m <- melt(num,id.vars=".id")
levels(m$variable) = c('not correlated','correlated')

p <- ggplot(m,aes(x=.id,y=value,fill=variable))+
  geom_bar(stat="identity")+
  scale_fill_manual(name="expression profile correlated with larvae?",values = Soc_palette)+theme_bw()+
  ylab("number of genes")+
  xlab("tissue")

png("~/Writing/Figures/NurseLarva/corApproach/NumGenesBinom.png",width=4000,height=4000,res=300)
p + theme_all+theme(legend.position = 'bottom',plot.title = element_text(size=20))+
  ggtitle('correlation of larval and nurse expression profiles,\nfor genes differentially expressed in nurses by larval stage')
dev.off()

#Numbers of negative and positive connections
PlotPosNeg <- function(x,label){
  p <- ggplot(x,aes(x=posSocial,y=negSocial))+
    geom_point(size=3)+theme_bw()+
    xlab("")+ylab("")+
    theme4+annotate("text",x = max(x$posSocial), y = max(x$negSocial),label=label,size=12)
  return(p)
}

label = names(socDet)
plots <- mapply(PlotPosNeg,socDet,label,SIMPLIFY = FALSE)
png("~/Writing/Figures/NurseLarva/corApproach/PosNegNumbers.png",width=4000,height=8000,res=300)
do.call("grid.arrange", c(plots, ncol=2))
dev.off()

#Previously derived GO annotations
go <- read.csv("~/Writing/Data/NurseSpecialization_transcriptomicData/GOannotation.csv")
new <- list()
for (gene in unique(go$gene)){
  d = go[go$gene %in% gene,]
  new[[gene]]=as.character(d$GO)
}

#Define genes as social or not
#Not necessary for GSEA, but necessary for identifying significant genes
selectFDR <- function(score){
  return(score < 0.05)
}

#GO analysis for categories
GSEA <- function(socialP){
  GOdata <- new("topGOdata",
                description="Simple session",ontology="BP",
                allGenes=socialP,geneSel=selectFDR,
                nodeSize = 10,
                annot=annFUN.gene2GO,gene2GO=new)
  resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks",scoreOrder="increasing")
  allRes <- GenTable(GOdata,KS=resultKS)
  return(allRes[allRes$Significant> allRes$Expected,])
}

getGSEA <- function(x){
  vec = x$social-x$control
  names(vec) = x$Gene
  tabl = GSEA(vec)
  return(tabl$Term[tabl$KS < 0.05])
}

GOres <- lapply(socP,GSEA)

GOres <- lapply(socDet,getGSEA)
res = sapply(GOres,'[',1:max) #keep 10 GO terms
res <- as.data.frame(res)
colnames(res) = names(socDet)

tt3 <- ttheme_minimal(
  core=list(
    fg_params=list(fontsize=8)),
  colhead=list(fg_params=list(fontface="bold",fontsize=14)))

#Make table grob out of results
resT <- tableGrob(res,theme=tt3,rows=NULL)
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
ggsave(p,file="~/Writing/Figures/NurseLarva/corApproach/TableGO.png",width=10,height=4,dpi=300)

#Evolutionary and functional differences between social and non-social genes

#What is the PS distribution of social genes vs NS? vs rest of genome?
#Relationship between number of social links and PS?
#GSEA for social genes


ext <- read.csv("~/Downloads/msx123_Supp (1)/MpharAnn.csv") #load in MBE results

getPS <- function(d){
  d = merge(d,ext,by="Gene",all.y=TRUE)
  d$PS2=as.character(d$PS2)
  d$PS2[d$Raw.PS=="Monomorium pharaonis"]="M. pharaonis"
  d=droplevels(d[!is.na(d$PS2),])
  d$PS2 = factor(d$PS2,levels = c(
    "cellular","eukaryote","bilaterian","insect","hymenopteran_ant","M. pharaonis"
  ))
  levels(d$PS2)[5]="hymenopteran"
  d$soc = 0
  d$soc[d$geneType=="social"]=1
  return(d[,c("Gene","negSocial","geneType","posSocial","negControl","posControl","social","control","PS2","SwissProt","UniProt","soc")])
}

graphPS <- function(d){
  p1 <- ggplot(d[],aes(x=PS2,fill = geneType))+
    geom_bar(stat="count")+theme_bw()

  return(p1)
}

PS_sum <- lapply(socDet,getPS)
d = PS_sum[[1]][order(PS_sum[[1]]$negSocial/(PS_sum[[1]]$posSocial+PS_sum[[1]]$posControl),decreasing=TRUE),]

testPS <- function(d){
  d = d[!is.na(d$negSocial),]
  soc = d$geneType == "social"
  res = matrix(nrow=length(levels(d$PS2)),ncol=2)
  for (i in 1:length(levels(d$PS2))){
    ps = d$PS2 == levels(d$PS2)[i]
    tbl <- rbind(c(sum(soc & ps),sum(soc & !ps)),c(sum(!soc & ps),sum(!soc & !ps)))
    test = chisq.test(tbl)
    res[i,1] = test$statistic
    res[i,2] = test$p.value
  }
  colnames(res) = c("X2","P")
  rownames(res) = levels(d$PS2)
  return(res)
}

test <- lapply(PS_sum,testPS)

lm <- glm(social ~ PS2,data = d,family="binomial")
summary(glht(lm, mcp(PS2="Tukey"))) 


graphs <- lapply(PS_sum,graphPS)


#Graph-based questions:
#What is the connectivity distribution of social genes in non-social networks?
#Using Genie, are social links typically larva to nurse or vice versa
#Clustering of larva genes based on nurse social connections? How well does this clustering match clustering from non-social larva only network?
#Are nurse social genes found in particular modules?

#Note that we probably want to ask all these questions for larval social genes as well

load("~/Dropbox/monomorium nurses/data.processed/cleandata.RData")
tissueNet <- function(code,pwr){
  fpkm <- fpkm[,grepl(code,colnames(fpkm))]
  conn = cor(t(fpkm))
  diag(conn) <- 0
  conns <- rowSums(abs(conn)^pwr)
  return(conns)
}

CH <- tissueNet("CH",6)
CG <- tissueNet("CG",6)
W <- tissueNet("W_L",6)

addConn <- function(d,conn){
  conn = data.frame(Gene = names(conn),kTotal = conn)
  dC = merge(d,conn,by="Gene")
  return(dC)
}

QCH = addConn(PS_sum[[1]],CH)
QCG = addConn(PS_sum[[2]],CG)
QCH_L = addConn(PS_sum[[1]],W)
QCG_L = addConn(PS_sum[[2]],W)

ggplot(QCG_L,aes(x=as.factor(soc),y=kTotal))+geom_boxplot(notch=TRUE)

fpkm$ID = rownames(fpkm)
fpkm = fpkm[,c(ncol(fpkm),1:(ncol(fpkm) -1 ))]
write.table(fpkm[,grepl("ID",colnames(fpkm)) | grepl("W_L",colnames(fpkm))],file="~/Data/Nurse_Larva/fpkm_larva.txt",row.names=FALSE,col.names=TRUE,sep="\t")
write.csv(fpkm[,grepl("CH",colnames(fpkm))],file="~/Data/Nurse_Larva/fpkm_head.csv")
write.csv(fpkm[,grepl("CG",colnames(fpkm))],file="~/Data/Nurse_Larva/fpkm_gaster.csv")

d <- read.table("~/Data/Nurse_Larva/Results_18_Dec_17_6/Clusters_Objects.tsv",header=TRUE,sep="\t")







#Calculate within and between module connectivity for each gene within the larval network. Could also define this as the distance from the centroid

#Relationship between number of social connections (or binomial p value, or social/non-social, or #social - control) and connectivity within the larval network 

#Are 'social' genes non-randomly distributed across the network?

########################
###Do larval social links form clusters? Are they comparable to clusters formed in networks using just larvae samples?
########################

########################
###What about nurse genes? Are social genes non-randomly distributed in nurse-built networks? Relationship between # social connections and within/between connectivity 
########################

########################
####Adding in Genie information- do links tend to be larva -> nurse or vise-versa? How does this relate to pos/neg connection type?
####Ideally, it would be nice to do Genie on a subset of genes
########################

########################
####It could be interesting to perform heirarchical clustering just on the CorResults matrices...in this way, simultaneously cluster
####rows and columns, and are able to ask about how complementary functions form. Within these blocks, could use Genie to assess directionality
########################

########################
####Somewhere down the line, we want to bring this back to development: how do nurse inputs affect larval development?
####How do nurse and larval modules interact over time?

########################
####Note that we also want to ask about function, PS, and constraint/alpha

#Load dataframe identifying social genes
load("~/Data/Nurse_Larva/socDet.RData")

#Load larval module definitions
df <- read.table("~/Data/Nurse_Larva/findK_cluster.txt",header=TRUE)
mods = t(df[11,]) #12 modules results in a high SIL
fpkm <- read.csv("~/Data/Nurse_Larva/fpkm.csv") #genes are indexed as in fpkm file
rownames(mods)=fpkm$X
mediods = unique(as.vector(mods))
mList <- lapply(mediods,function(x) rownames(mods)[mods == x])

#From assignNurseMods.py, files containing module assignments for nurse genes
#socNum is for accessing socDet for specific tissues
getSocMod <- function(code,socNum){
  #First row is true values, subsequent values are derived by scrambling stage codes
  dfH <- read.table(paste("~/Data/Nurse_Larva/sortMods",code,".txt",sep=""),head=TRUE)
  names(dfH) = fpkm$X
  
  #Filter for social genes
  soc = socDet[[socNum]]
  dfH = dfH[,colnames(dfH) %in% soc$Gene[soc$geneType=='social']]
  
  #True assignments
  Ht <- t(dfH[1,])
  
  #'Bootstrapped' assignments
  dfH = dfH[-c(1),]
  
  #Identify null distribution of numbers of genes per module
  occ <- apply(dfH,1,function(x) unlist(lapply(mediods,function(y) length(colnames(dfH)[x == y]))))
  
  #null hypothesis mean based on resampling 
  null_mean <- apply(occ,1,function(x) mean(x))
  
  #Social modules
  nSoc <- unlist(lapply(mediods,function(y) length(colnames(Ht)[Ht == y])))
  binom <- apply(cbind(nSoc,null_mean),1,function(y) binom.test(x = y['nSoc'], p = y['null_mean']/length(Ht),n = length(Ht),alternative = "greater")$p.value)
  socMod <- mediods[binom < 0.05]
  socList <- lapply(socMod,function(x) rownames(Ht)[Ht == x])
  names(socList) = socMod
  return(list(socList,nSoc,null_mean,binom))
}

HsocMod <- getSocMod('CH',1)
GsocMod <- getSocMod('CG',2)

#Not necessary for GSEA, but necessary for identifying significant genes
selectFDR <- function(score){
  return(score == 1)
}

#GO analysis for categories
GSEA <- function(d){
  GOdata <- new("topGOdata",
                description="Simple session",ontology="BP",
                allGenes=d,geneSel=selectFDR,
                nodeSize = 10,
                annot=annFUN.gene2GO,gene2GO=new)
  resultFS <- runTest(GOdata, algorithm = "classic", statistic = "fisher",scoreOrder="increasing")
  allRes <- GenTable(GOdata,FS=resultFS)
  return(allRes) #only return enriched categories
}

GOmod <- function(genes){
  d = rep(0,nrow(fpkm))
  names(d) = fpkm$X
  d[names(d) %in% genes]=1
  gs <- GSEA(d)
  return(gs)
}

go_mods = lapply(mList,GOmod)

##############
##Representation of "social" genes among modules
##############
load("~/Data/Nurse_Larva/socDet.RData")
larvSoc = socDet[5:6]
names(larvSoc)=c("head","abdomen")
getMod <- function(soc,mods){
  for (i in 1:length(mods)){
    soc$mod[soc$Gene %in% mods[[i]]]=i
  }
  return(soc)
}

larvSoc <- lapply(larvSoc,getMod,mods=mList)
lSoc = data.frame(Gene=larvSoc[[1]]$Gene,typeHead = larvSoc[[1]]$geneType,typeGaster = larvSoc[[2]]$geneType,
                  Module = larvSoc[[1]]$mod,headSocConn = (larvSoc[[1]]$social - larvSoc[[1]]$control),
                  gasterSocConn = (larvSoc[[2]]$social - larvSoc[[2]]$control))
lSoc$type = "not_soc"
lSoc$type[lSoc$typeHead=="social"&lSoc$typeGaster!="social"]="head_soc"
lSoc$type[lSoc$typeHead!="social"&lSoc$typeGaster=="social"]="abd_soc"
lSoc$type[lSoc$typeHead=="social"&lSoc$typeGaster=="social"]="both_soc"
lSoc$type=factor(lSoc$type,levels=c("not_soc","head_soc",
                                    "abd_soc","both_soc"))

png("~/Writing/Figures/NurseLarva/GeneCountMods.png")
ggplot(lSoc,aes(x=as.factor(Module),fill=type))+
  geom_bar(stat="count")+
  theme_bw()+
  xlab("Module")+
  theme_all
dev.off()

lSoc$tH = lSoc$tG = 0 
lSoc$tH[lSoc$typeHead=="social"]=1
lSoc$tG[lSoc$typeGaster=="social"]=1

t = table(lSoc$Module,lSoc$tH)
chiResH=lapply(seq(1,12,by=1),function(i) fisher.test(rbind(c(t[i,2],t[i,1]),c(colSums(t[-i,])[2],colSums(t[-i,])[1])),
                    alternative="greater"))

t = table(lSoc$Module,lSoc$tG)
chiResG=lapply(seq(1,12,by=1),function(i) fisher.test(rbind(c(t[i,2],t[i,1]),c(colSums(t[-i,])[2],colSums(t[-i,])[1])),
                                                     alternative="greater"))

lSoc$Module=as.factor(lSoc$Module)
p1 <- ggplot(lSoc,aes(x=Module,y=headSocConn))+
  geom_boxplot()+theme_bw()+
  xlab("Module")+
  ylab("Number social connections")+
  ggtitle("head")+
  scale_fill_discrete(name="connection type")

p2 <- ggplot(lSoc,aes(x=Module,y=gasterSocConn))+
  geom_boxplot()+theme_bw()+
  xlab("Module")+
  ggtitle("abdomen")+
  ylab("Number social connections")+
  scale_fill_discrete(name="connection type")

png("~/Writing/Figures/NurseLarva/ModSocConnections.png")
grid.arrange(p1,p2)
dev.off()

#######
##Calculate connectivity within larval networks
##Separately, calculate connectivity (geometric mean of pearson cor?)
##As well as distance from medoid center
##They should be highly correlated...
#######
rownames(fpkm) = fpkm$X
fpkm = fpkm[,-c(1)]
conn = abs(cor(t(fpkm)))
meds = rownames(fpkm)[mediods+1]

#Calculate within-connectivity as distance to medoid
lSoc$withinConn = apply(lSoc[,c('Gene','Module')],1,
                        function(x) conn[x['Gene'],meds[as.integer(x['Module'])]])

lSocH = lSoc[lSoc$tH==1,]
lSocG = lSoc[lSoc$tG==1,]
#lSocH = merge(lSocH,socDet[[5]],by="Gene",all.y=FALSE)

ggplot(lSoc,aes(x=typeHead,y=withinConn))+
  geom_boxplot(notch=TRUE)

ggplot(lSoc,aes(x=typeGaster,y=withinConn))+
  geom_boxplot(notch=TRUE)

keepH = c(1,2,3,11)
ggplot(lSocH,aes(x=withinConn,y=headSocConn))+
  geom_point()+
  theme_bw()+
  xlab("distance to module center")+
  ylab("number of social connections")+
  ggtitle("nurse head")+
  geom_smooth(method="lm",se=FALSE)

keepG= 11
ggplot(lSocG,aes(x=withinConn,y=gasterSocConn))+
  geom_point()+
  theme_bw()+
  xlab("distance to module center")+
  ylab("number of social connections")+
  ggtitle("nurse abdomen")+
  geom_smooth(method="lm",se=FALSE)

png("~/Writing/Figures/NurseLarva/SocConnConnectivity.png")
grid.arrange(p1,p2)
dev.off()

allF <- merge(f.est,lSoc,by.x='gene',by.y="Gene")

ggplot(allF,aes(x=withinConn,y=withinConn,fill=PS2))+
  geom_boxplot()

ggplot(allPS,aes(x=Module,y=withinConn,fill=PS2))+
  geom_boxplot()

ggplot(allPS[allPS$tH==1,],aes(x=PS2,y=headSocConn))+
  geom_boxplot(notch=TRUE)

ggplot(allPS[allPS$tG==1,],aes(x=PS2,y=gasterSocConn))+
  geom_boxplot(notch=TRUE)

ggplot(allPS[allPS$tG==1,],aes(x=PS2,y=withinConn))+
  geom_boxplot(notch=TRUE)

df <- t(df)

same = matrix(nrow=nrow(df),ncol=nrow(df))
#Tabulate number of times genes are found in the same module
for (i in 1:nrow(df)){
  for (j in 1:nrow(df)){
    same[i,j] = sum(df[i,]==df[j,])
  }
}

diag(same) = 0 #Not relevant here

#One strategy could be to now use heirarchical clustering, 
#Where we calculate correlation of module definition
#Then basically decrease cut height until we have the optimal number of modules. 
#This will leave out some genes, which we can put in a "grey" module or just force to the closest module

#identify clusters that are found every time
clusters = list()
index 
df_temp = df
skip = c()
for (i in 1:nrow(df)){
  for (j in 1:nrow(df)){
    if (j %in% skip | i %in% skip){
      next;
    }
    #Wait until we find a pair of seed clusters
    if (same[i,j]==3){
      keep = seq(1,nrow(same),by=1)[same[i,]==3] #Collect all genes found in the same module
      skip = c(skip,keep) #In the future, skip these
      clusters = c(clusters,list(c(i,keep)))
    }
  }
}


table(df[1,])
d <- read.csv("~/Data/Nurse_Larva/fpkm.csv")
d = d[,grepl("W_L",colnames(d))]
d = d[1:50,]





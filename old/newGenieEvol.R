getTopGenes <- function(CV1,CV2,N){
  top1 = names(CV1)[order(CV1,decreasing=TRUE)][1:floor(N/2)] #First select top N/2 genes from larvae and nurse tissue
  top2 = names(CV2)[order(CV2,decreasing=TRUE)][1:floor(N/2)]
  keep = c(top1,top2)
  keep = keep[!duplicated(keep)] #Some top genes may be duplicated
  samp = 1
  while (length(keep) < N){ #Take turns adding genes until we have N genes
    if (samp %% 2 == 1){
      CV1 = CV1[!names(CV1) %in% keep]
      keep = c(keep,names(CV1)[order(CV1,decreasing=TRUE)][1])
    } else {
      CV2 = CV2[!names(CV2) %in% keep]
      keep = c(keep,names(CV2)[order(CV2,decreasing=TRUE)][1])
    }
    samp = samp + 1
  }
  return(keep)
}

#Select top N genes by CV for analysis; label genes by where they are expressed
selectGene <- function(N,codes,names){
  load("~/Dropbox/monomorium nurses/data.processed/cleandata.RData")
  f1 <- fpkm[,grepl(codes[1],colnames(fpkm))]
  f2 <- fpkm[,grepl(codes[2],colnames(fpkm))]
  CV1 = apply(f1,1, sd)/apply(f1,1,mean)
  CV2 = apply(f2,1, sd)/apply(f2,1,mean)
  keep = getTopGenes(CV1,CV2,N)
  f1 = f1[keep,]
  f2 = f2[keep,]
  rownames(f1) = paste(rownames(f1),names[1],sep="_")
  rownames(f2) = paste(rownames(f2),names[2],sep="_")
  colnames(f1) = gsub("[A-Z]+_.*","",colnames(f1))
  colnames(f2) = gsub("[A-Z]+_.*","",colnames(f2))
  f2 = f2[,colnames(f2) %in% colnames(f1)]
  f1 = f1[,colnames(f1) %in% colnames(f2)]
  f = rbind(f1,f2)
  f = log(f + sqrt(f ^ 2 + 1)) #Variance stablizing transformation
  return(f)
}

#Run Genie algorithm on dataframe
runGenie <- function(run,N){
  setwd("~/Downloads/GENIE3_R_C_wrapper") #Have to switch directories because there are .so files we need
  source("~/Downloads/GENIE3_R_C_wrapper/GENIE3.R")
  res = list()
  i=1
  while (TRUE){
    x <- try(GENIE3(as.matrix(f)))
    if (!inherits(x,"try-error")){ #Sometimes the data get weird, just try again!
      res[[i]]=x
      i = i+1
      if (i > N) break;
    }
  }
  res1 = Reduce("+",res)/length(res)
  return(res1)
}

#tablulate connection strengths based on 
tabConn <- function(net){
  d = get.link.list(net)
  d = d[d$weight > 0,]
  d$tissue1 = gsub('.*_','',d$regulatory.gene)
  d$tissue2 = gsub('.*_','',d$target.gene)
  conn <- ddply(d,~regulatory.gene+tissue1+tissue2,summarize,
                N = length(weight),
                Mean = mean(weight),
                Max = max(weight),
                k = sum(weight))
  return(conn)
}

#Trim edges from focal networks until only social connections that are 'significant', i.e. significantly higher than random nurse networks
#Returns F, with the non-significant connections removed
trimEdges <- function(Foc,Rand,pval){
  #First, just cut out everything but the top 20%
  Foc = Foc[Foc$weight > quantile(Foc$weight,0.8),]
  Rand = Rand[Rand$weight > quantile(Rand$weight,0.8),]
  #Using ranks makes translating between the networks easier
  Foc$rank = seq(1,nrow(Foc),by=1)
  Rand$rank = seq(1,nrow(Rand),by=1)
  
  #Get dfs with just the Larv-Nurse or Nurse-Larv connections
  FLsoc = Foc[grepl("Larv",Foc$regulatory.gene) & grepl("Nurse",Foc$target.gene),]
  RLsoc = Rand[grepl("Larv",Rand$regulatory.gene) & grepl("Nurse",Rand$target.gene),]
  FNsoc = Foc[grepl("Nurse",Foc$regulatory.gene) & grepl("Larv",Foc$target.gene),]
  RNsoc = Rand[grepl("Nurse",Rand$regulatory.gene) & grepl("Larv",Rand$target.gene),]
  
  #Get genes that would be in the top 5% of the random nurse distribution
  cut = quantile(RLsoc$rank,0.05)
  keepL = FLsoc[FLsoc$rank < cut,]  
  cut = quantile(RNsoc$rank,0.05)
  keepN = FNsoc[FNsoc$rank < cut,]
  
  return(list(keepL,keepN))
}

#Make matrix where many connections are 0 (those not in the top ranks of random social connections from trimEdges())
genMat <- function(soc,full){
  fL = full[grepl("Larv",full$regulatory.gene) & grepl("Larv",full$target.gene),]
  fL = fL[fL$weight > min(soc[[1]]$weight),] # keep same number as the number of social connections from larvae
  
  fN = full[grepl("Nurse",full$regulatory.gene) & grepl("Nurse",full$target.gene),]
  fN = fN[fN$weight > min(soc[[2]]$weight),] # keep same number as the number of social connections from nurses
  
  allConn = rbind(soc[[1]][,-c(4)],soc[[2]][,-c(4)],fL,fN)
  mat = acast(allConn,regulatory.gene ~ target.gene,value.var="weight",fill=0) #Make matrix from available connections

  return(mat)
}

#Use sparse matrix to generate modules, according to the method of Clauset, Newman, and Moore 2004 'Detecting community structure in very large networks'
makeModules <- function(mat){
  d = get.link.list(mat)
  d = d[d$weight > 0,]
  d$G3 = apply(d[,c(1,2)],1,paste,collapse='-')
  d$G4 = apply(d[,c(2,1)],1,paste,collapse='-')
  dT = d
  g = c()
  rem = c()
  for (i in 1:nrow(d)){
    g = c(g,as.character(d$G3[i]))
    if (d$G4[i] %in% g){
      rem = c(rem,i) #Remove connections which have duplicate recipricals 
    }
  }
  dT = dT[-rem,]
  e = as.matrix(dT[,c(1,2)])
  graph <- graph_from_edgelist(e,directed=FALSE)
  fc <- cluster_edge_betweenness(graph,weights = dT$weight)
  layout <-layout.fruchterman.reingold(graph)
  
  a = dist(QH)
  hc = hclust(a)
  plot(fc, graph,  layout=layout, vertex.label=NA, vertex.size=5, edge.arrow.size=.2)
  
}

codes = c("QW","QCH")
names = c("WorkLarvQR","NurseH")
f = selectGene(500,codes,names)
genes = unique(gsub("_[A-Za-z]+$","",rownames(f)))

QH = runGenie(f,2)
H = get.link.list(QH)


fc <- cluster_edge_betweenness(graph,weights=H$weight)

codes = c("QW","RH")
names = c("WorkLarvQR","RNurseH")
f = selectGene(500,codes,names)
RH = runGenie(f,2)
R = get.link.list(RH)

Hsoc <- trimEdges(H,R) #Trim edges based on random network ranks
Hmat <- genMat(Hsoc,H) #Make a matrix (with many zeros) of connections between genes, removing those with connection strengths lower than those in the social matrices

Hconn <- tabConn(Hmat)
Hconn$gene = gsub('_.*','',Hconn$regulatory.gene)
Hext <- merge(Hconn,ext,by.x="gene",by.y="Gene",all.x=TRUE)
Hext$PS2 = as.character(Hext$PS2)
Hext$PS2[Hext$Raw.PS=="Monomorium pharaonis"]='M. pharaonis'
Hext$tissueAll = apply(Hext[,c("tissue1","tissue2")],1,paste,collapse='_')

ggplot(Hext,aes(x=tissueAll,y=k,fill=PS2))+geom_boxplot()



m = unique(c(levels(d$regulatory.gene),levels(d$target.gene)))
mat = matrix(0,nrow=length(m),ncol=length(m))
for (i in 1:length(m)){
  for (j in 1:length(m)){
    dt = d$weight[d$regulatory.gene %in% m[i] & d$target.gene %in% m[j]]
    if (length(dt) == 1){
      mat[i,j]=dt
    }
  }
}

library(modMax)
s = spectralOptimization(mat)

QH2 = runGenie(f)
Hconn = tabConn(QH,"WorkLarvQR","NurseH")

codes = c("QW","QCG")
names = c("WorkLarvQR","NurseG")
f = selectGene(500,codes,names)
genes = unique(gsub("_[A-Za-z]+$","",rownames(f)))
QG = runGenie(f)
Gconn = tabConn(QG,"WorkLarvQR","NurseG")

df = merge(Hconn,Gconn,by="gene",all.x=TRUE,all.y=TRUE)
df$SIL_WH=df[,3] - df[,2]
df$SIL_WG=df[,7] - df[,6]
df$SIWH=df[,4]- df[,5]
df$SIWG=df[,8] - df[,9]
df$SIO = (df$SIWH + df$SIWG+(df$SIL_WH+df$SIL_WG)/2)/3
f.est <- read.csv("~/Data/MKtestConstraintOneAlpha.csv")
colnames(f.est) = c("gene","f.est")
dfAll <- merge(df,f.est,by="gene")
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

ggplot(dfAll,aes(x=PS2,y=SIO))+geom_boxplot(notch=TRUE)

d = melt(dfAll[,c(1,27,28:35)],id.vars=c("Gene","PS2"))
ggplot(d,aes(x=variable,y=value,fill=PS2))+geom_boxplot(notch=TRUE)

load("~/Dropbox/monomorium nurses/data.processed/cleandata.RData")

#Get matrices of expression data for larvae and a nurse tissue, to test for correlation
#Filter for top N genes by CV
OrderedMatrices <- function(larv,nurse,N){
  fL = fpkm[,grepl(larv,colnames(fpkm))]
  fN = fpkm[,grepl(nurse,colnames(fpkm))]
  
  colnames(fL) = gsub('[A-Z]+_.*','',colnames(fL))
  colnames(fN) = gsub('[A-Z]+_.*','',colnames(fN))
  
  fL = fL[,colnames(fL) %in% colnames(fN)]
  fN = fN[,colnames(fN) %in% colnames(fL)] 
  
  return(list(fL,fN))
}

#Calculate CV after dropping most variable sample
CVnoOut <- function(vect){
  m = mean(vect)
  d = vect - m
  vect = vect[abs(d) != max(abs(d))]
  return(sd(vect)/mean(vect))
}

#Filter for the 1000 most variable (by average CV)
filterCV <- function(d,N){
  CV1 = as.vector(apply(d[[1]],1, CVnoOut))
  CV2 = as.vector(apply(d[[2]],1, sd)/apply(d[[2]],1,mean))
  CV = (CV1 + CV2)/2
  names(CV) = rownames(d[[1]])
  CV = CV[order(CV,decreasing=TRUE)]
  a = lapply(d,function(d) d[rownames(d) %in% names(CV)[1:N],])
  return(a)
}

#Get p values of pearson correlation
calcPearsonP <- function(d1,d2){
  pVals = matrix(nrow=nrow(d1),ncol=nrow(d2))
  for (i in 1:nrow(d1)){
    for (j in 1:nrow(d2)){
      test = cor.test(t(d1[i,]),t(d2[j,]))
      pVals[i,j] = test$p.value
    }
  }
  rownames(pVals)=rownames(d1)
  colnames(pVals)=rownames(d2)
  return(pVals)
}

#gene-by-gene correlation values comparing nurse to larva trajectories
getCorVec <- function(larv,nurse,N){
  d = OrderedMatrices(larv,nurse)
  d = filterCV(d,N) #Could also filter by geometric mean p value from ~ stage + colony done separately
  rownames(d[[2]])=paste(rownames(d[[2]]),"_nurse",sep="")
  mat = cor(t(d[[1]]),t(d[[2]]))
  #matP = calcPearsonP(d[[1]],d[[2]])
  return(list(mat))
}

allTestDF <- function(N){
  QCH = getCorVec('QW','QCH',5000)
  QCG = getCorVec('QW','QCG',5000)
  RH = getCorVec('QW','RH',5000)
  RG = getCorVec('QW','RG',5000)
  LCH = getCorVec('LW','LCH',N)
  LCG = getCorVec('LW','LCG',N)
  return(list(QCH,QCG,RH,RG,LCH,LCG))
}

numHigh <- function(res,r){
  res = as.vector(res[[1]])
  return(length(res[abs(res) > r&!is.na(res)]))
}

QCH = getCorVec('QW','QCH',10970)
QCG = getCorVec('QW','QCG',10970)
RH = getCorVec('QW','RH',10970)
RG = getCorVec('QW','RG',10970)
LCH = getCorVec('LW','LCH',10970)
LCG = getCorVec('LW','LCG',10970)

Res = list(QCH,QCG,RH,RG,LCH,LCG)
names(Res)=c("QCH","QCG","RH","RG","LCH","LCG")
Res_cor = matrix(ncol=3,nrow=24)
i=1
for (rcut in c(0.7,0.8,0.9,0.95)){
  for (j in 1:length(Res)){
    Res_cor[i,1]=names(Res)[j]
    Res_cor[i,2]=rcut
    Res_cor[i,3] =numHigh(Res[[j]],rcut)
    i=i+1
  }
}

Res_cor = as.data.frame(Res_cor)
colnames(Res_cor) = c("sample","pearson.Cor","N")
Res_cor$tissue='head'
Res_cor$tissue[grepl("G",Res_cor$sample)]='abdomen'
Res_cor$pearson.Cor=as.factor(Res_cor$pearson.Cor)
Res_cor$sample = factor(Res_cor$sample,levels=c("RH","QCH","LCH",
                                         "RG","QCG","LCG"))
levels(Res_cor$sample) = c("rand nurse head","QR nurse head","QL nurse head",
                          "rand nurse abd","QR nurse abd","QL nurse abd")
png("~/Documents/quickplot.png")
ggplot(Res_cor,aes(x=factor(pearson.Cor),fill=sample,y=as.numeric(as.character(N))))+
  geom_bar(stat='identity',position = position_dodge())+
  scale_y_log10()+
  theme_bw()+xlab("minimum pairwise pearson correlation")+
  ylab("Count")+theme_all+theme(legend.position = "right")
dev.off()

ggplot(Res_cor[Res_cor$pearson.Cor!='0.7',],aes(x=factor(pearson.Cor),fill=sample,y=as.numeric(as.character(N))))+
  geom_bar(stat='identity',position = position_dodge())+
  theme_bw()+xlab("minimum pairwise pearson correlation")+
  xlab("Count")

p = as.vector(QCH[[2]])
r = as.vector(QCH[[1]])
keep = p < 0.05
rSig = r[keep]

p = as.vector(RH[[2]])
r = as.vector(RH[[1]])
keep = p < 0.05
rSigR = r[keep]
#Ns = c(100,500,1000,2000,5000,10970)
#p = lapply(Ns,allTestDF)
allTestFull <- allTestDF(10970)




themes = function(N) theme_bw()+ggtitle(paste('top',N,'genes',sep=" "))
lapply(p,function(x) print(x[['df']]))
                          
for (N in Ns){
  p[[N]] <- allTestDF(N)
  p[[N]] <- p[[N]]+)
}

lapply(p[Ns],print)


mean(abs(LCG),na.rm=TRUE)

wilcox.test(abs(QCH),abs(RH),alternative='greater')
wilcox.test(abs(QCG),abs(RG),alternative='greater')

a = data.frame(A=c(1,3,4),B=c(4,3,2),C=c(5,1,2))
b = data.frame(B=c(1,3,4),A=c(4,3,2),C=c(2,6,2))


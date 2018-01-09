library(plyr)
library(rslurm)

setwd("~/Nurse_Larva")
load("cleandata.RData")

#Return list of expression matrices for correlation analysis
formatExpr <- function(codes){
  d = lapply(codes,subExpr) #get expression matrix for each of three sample types
  for (i in 1:2){
    #Needed to align random nurses to QR larvae
    colnames(d[[i]]) = gsub("R","Q",colnames(d[[i]]))
  }
  dF = alignStage(list(d[[1]],d[[2]]))
  return(dF) #Returns data for foc/larv comparison and rand/larv comparison
}

#only keeps relevant genes and lines up larvae and nurses by stage
subExpr <- function(code){
  f = fpkm[grepl(code,colnames(fpkm))]
  #hyperbolic sine transformation
  f = log(f + sqrt(f^2 + 1))
  colnames(f) = substr(colnames(f),start=1,stop=4) #Has colony, stage and queen presence information
  return(f)
}

#Take formatted expression data and remove instances of colonies that don't match (missing data)
alignStage <- function(d){
  d[[1]]=d[[1]][,colnames(d[[1]]) %in% colnames(d[[2]])]
  d[[2]]=d[[2]][,colnames(d[[2]]) %in% colnames(d[[1]])]
  return(d)
}

getModList <- function(codes){
  df <- read.table(paste("sortMods",codes[1],".txt",sep=""),head=TRUE)
  mods <- t(df[1,])
  return(mods)
}

#For each gene in input frame, get connectivities to other genes
getConns <- function(i){
  c = abs(cor(t(d[[1]][i,]),t(d[[2]])))^6
  kTotal = sum(c)
  kMod = sum(c[mods==mods[i]])
  return(data.frame(Gene = rownames(d[[1]])[i],kTotal = kTotal, kMod = kMod))
}

#Function computes WGCNA-like connectivity, signed (power = 6) and unsigned (power = 12)
#Computes connectivity between a pair of expression matrices
WGCNAconn <- function(codes){
  
  #Get module definitions for nurse genes
  mods <<- getModList(codes)
  d <<- formatExpr(codes)
  
  df = data.frame(i = seq(1,nrow(d[[1]]),by=1))
  sjob <- slurm_apply(getConns, df, jobname = 'getConn',
                      nodes = 4, cpus_per_node = 20, add_objects=c("d","mods","unsignedPWR"),submit = TRUE)
  res <- get_slurm_out(sjob,outtype='raw') #get output as lists
  res <- ldply(res)
  
  write.csv(res,file=paste(codes[1],codes[2],"connFrame.csv",sep=""))
}

connCode <- list(
  c('CH','W_L'),
  c('CG','W_L'),
  c('RH','W_L'),
  c('RG','W_L'),
  c('QCH','W_L'),
  c('QCG','W_L'),
  c('CH','CH'),
  c('CG','CG'),
  c('RH','RH'),
  c('RG','RG'),
  c('QCH','QCH'),
  c('QCG','QCG')
)

unsignedPWR = 6
lapply(connCode,WGCNAconn)



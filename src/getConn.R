library(plyr)
setwd("~/Data/Nurse_Larva")
load("~/Dropbox/monomorium nurses/data.processed/cleandata.RData")
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

#Take correlation matrix and calculate connectivity
getConn <- function(cor,mods){
  if (is.null(mods)) mods = rep(1,nrow(cor)) #Define everything as in same module
  res = lapply(seq(1,nrow(cor)), function(i) sum(cor[i,mods==mods[i]])/ncol(cor[,mods==mods[i]]))
  res = unlist(res)
  return(res)
}

#Function computes WGCNA-like connectivity, signed (power = 6) and unsigned (power = 12)
#Computes connectivity between a pair of expression matrices
WGCNAconn <- function(codes,unsignedPWR,signedPWR){
  d <- formatExpr(codes)
  
  nGene = nrow(d[[1]])
  dM <- rbind(d[[1]],d[[2]])
  corMat <- cor(t(dM))
  
  #Only want off diagonal square matrix, which corresponds to between tissue correlation with rows as code1 and columns as code 2
  corMat <- corMat[1:nGene,(nGene+1):ncol(corMat)]
  
  #Necessary because R adds '1' to duplicate rownames
  colnames(corMat) = rownames(d[[2]])
  corUnsigned = abs(corMat^unsignedPWR)
  corSigned = corMat^signedPWR
  
  #Get module definitions for nurse genes
  mods <- getModList(codes)
  
  #Compute within module and total network connectivity
  wMu <- getConn(corUnsigned,mods)
  wMs <- getConn(corSigned,mods)
  tMu <- getConn(corUnsigned,mods=NULL)
  tMs <- getConn(corSigned,mods=NULL)
  kDat <- data.frame(withinMod_unsigned = wMu, withinMod_signed = wMs,
                     total_unsigned = tMu, total_signed = tMs)
  return(kDat)
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

conns <- lapply(connCode,WGCNAconn,unsignedPWR=6,signedPWR=11)
save(conns,file = "connMeas.RData")

setwd("~/GitHub/MonomoriumNurseLarva/")
codes = "CH"
df <- read.table(paste("sortMods",codes[1],".txt",sep=""),head=TRUE)



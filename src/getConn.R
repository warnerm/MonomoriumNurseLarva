library(plyr)
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

#Function computes WGCNA-like connectivity, signed (power = 6) and unsigned (power = 12)
#Computes connectivity between a pair of expression matrices
WGCNAconn <- function(codes){
  d <- formatExpr(codes)
  
  nGene = nrow(d[[1]])
  dM <- rbind(d[[1]],d[[2]])
  corMat <- cor(t(dM))
  
  #Only want off diagonal square matrix, which corresponds to between tissue correlation with rows as code1 and columns as code 2
  corMat <- corMat[1:nGene,(nGene+1):ncol(corMat)]
  
  #Necessary because R adds '1' to duplicate rownames
  colnames(corMat) = rownames(d[[2]])
  corUnsigned = corMat^4
  #corSigned = corMat^11
  return(corUnsigned)
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

conns <- lapply(connCode,WGCNAconn)
save(conns,file = "connMeas.RData")
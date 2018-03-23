library(WGCNA)
library(plyr)
library(DESeq2)
load("~/Dropbox/monomorium nurses/data.processed/cleandata.RData")
enableWGCNAThreads()
formatExpr <- function(codes){
  d = lapply(codes,subExpr) #get expression matrix for each of three sample types
  for (i in 1:length(d)){
    #Needed to align random nurses to QR larvae
    colnames(d[[i]]) = gsub("R","Q",colnames(d[[i]]))
  }
  dF = alignStage(d)
  return(dF) #Returns data for foc/larv comparison and rand/larv comparison
}

#only keeps relevant genes and lines up larvae and nurses by stage
subExpr <- function(code){
  f = fpkm[grepl(code,colnames(fpkm))]
  colnames(f) = substr(colnames(f),start=1,stop=4) #Has colony, stage and queen presence information
  return(f)
}

#Take formatted expression data and remove instances of colonies that don't match (missing data)
alignStage <- function(d){
  all_name = colnames(d[[1]])
  for (i in 2:length(d)){
    all_name = all_name[all_name %in% colnames(d[[i]])]
  }
  
  d = lapply(d,function(x) x[,colnames(x) %in% all_name])
  
  return(d)
}

metaExpr <- function(codes){
  expr = formatExpr(codes)
  for (i in 1:length(expr)){
    expr[[i]]$Gene <- paste(codes[i],rownames(expr[[i]]),sep="_")
  }
  allExpr <- ldply(expr)
  rownames(allExpr) = allExpr$Gene
  allExpr = allExpr[,colnames(allExpr)!=("Gene")]
}

plotSFT <- function(datExpr,name){
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  cex1 = 0.9
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  
  # Mean connectivity as a function of the soft-thresholding power
  png(paste("~/Writing/Figures/NurseLarva/sft_metasample_",name,".png",sep=""))
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  dev.off()
}

d <- metaExpr(codes = c('CH','CG','W_L'))

#Try log2 normalization
d_log2 <-log(d+1,base=2)

#Try rlog transformation
d_rlog <- rlog(floor(as.matrix(d)))

#variance stabilizing transformation
d_vst <- varianceStabilizingTransformation(floor(as.matrix(d)))

#inverse hyperbolic sine transformation
d_hst <- log(d+sqrt(d*d+1))

plotSFT(t(d_log2),"log2_focal")
plotSFT(t(d_rlog),"rlog_focal")
plotSFT(t(d_vst),"vst_focal")
plotSFT(t(d_hst),"hst_focal")

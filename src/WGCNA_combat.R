library(WGCNA)
library(plyr)
library(DESeq2)
library(sva)
load("cleandata.RData")
enableWGCNAThreads(nThreads = 10)
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
  abline(h=0.80,col="blue")
  
  dev.off()
}
quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

runWGCNA <- function(datExpr,softPower,name){
  adjacency = adjacency(datExpr, power = softPower);
  
  TOM = TOMsimilarity(adjacency);
  dissTOM = 1-TOM
  geneTree = hclust(as.dist(dissTOM), method = "average");
  minModuleSize = 30;
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize);
  dynamicColors = labels2colors(dynamicMods)
  MEList = moduleEigengenes(datExpr, colors = dynamicColors)
  MEs = MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs);
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average");
  MEDissThres = 0.25
  # Call an automatic merging function
  merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors;
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs;
  moduleColors = mergedColors
  # Construct numerical labels corresponding to the colors
  colorOrder = c("grey", standardColors(50));
  moduleLabels = match(moduleColors, colorOrder)-1;
  MEs = mergedMEs;
  save(MEs,moduleLabels,moduleColors,geneTree,file = paste("Network",name,".RData",sep=""))
}

datExpr <- fpkm[,grepl("CH|CG|W_L",colnames(fpkm))]
d <- datExpr
tissue = rep('head',ncol(datExpr))
tissue[grepl("G",colnames(datExpr))]='gaster'
tissue[grepl("L",colnames(datExpr))]='larva'
colony = rep(1,ncol(datExpr))
colony[grepl('X2',colnames(datExpr))]=2
colony[grepl('X3',colnames(datExpr))]=3

#Try rlog transformation
d_rlog <- rlog(floor(as.matrix(d)))

#variance stabilizing transformation
d_vst <- varianceStabilizingTransformation(floor(as.matrix(d)))

#inverse hyperbolic sine transformation
d_hst <- log(d+sqrt(d*d+1))

plotSFT(t(quantile_normalisation(ComBat(d_rlog,batch = tissue))),"rlog_combat")
plotSFT(t(quantile_normalisation(ComBat(d_vst,batch = tissue))),"vst_combat")
plotSFT(t(quantile_normalisation(ComBat(d_hst,batch = tissue))),"hst_combat")


datExpr <- fpkm[,grepl("RH|RG|W_L",colnames(fpkm))]
d <- datExpr
tissue = rep('head',ncol(datExpr))
tissue[grepl("G",colnames(datExpr))]='gaster'
tissue[grepl("L",colnames(datExpr))]='larva'
colony = rep(1,ncol(datExpr))
colony[grepl('X2',colnames(datExpr))]=2
colony[grepl('X3',colnames(datExpr))]=3

d_hst <- log(d+sqrt(d*d+1))
d_rlog <- rlog(floor(as.matrix(d)))
d_vst <- varianceStabilizingTransformation(floor(as.matrix(d)))
plotSFT(t(quantile_normalisation(ComBat(d_rlog,batch = tissue))),"rlog_combat_random")
plotSFT(t(quantile_normalisation(ComBat(d_vst,batch = tissue))),"vst_combat_random")
plotSFT(t(quantile_normalisation(ComBat(d_hst,batch = tissue))),"hst_combat_random")
# 
# ###Conclusion: vst, softpower = 10 for both
# datExpr <- fpkm[,grepl("CH|CG|W_L",colnames(fpkm))]
# d <- datExpr
# tissue = rep('head',ncol(datExpr))
# tissue[grepl("G",colnames(datExpr))]='gaster'
# tissue[grepl("L",colnames(datExpr))]='larva'
# d_vst <- varianceStabilizingTransformation(floor(as.matrix(d)))
# datExpr <- t(quantile_normalisation(ComBat(d_vst,batch = tissue)))
# runWGCNA(datExpr,10,"focal_combat")
# 
# datExpr <- fpkm[,grepl("RH|RG|W_L",colnames(fpkm))]
# d <- datExpr
# tissue = rep('head',ncol(datExpr))
# tissue[grepl("G",colnames(datExpr))]='gaster'
# tissue[grepl("L",colnames(datExpr))]='larva'
# d_vst <- varianceStabilizingTransformation(floor(as.matrix(d)))
# datExpr <- t(quantile_normalisation(ComBat(d_vst,batch = tissue)))
# runWGCNA(datExpr,10,"random_combat")
# 







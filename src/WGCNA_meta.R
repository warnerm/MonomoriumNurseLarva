library(WGCNA)
library(plyr)
library(DESeq2)
load("cleandata.RData")
enableWGCNAThreads(nThreads = 20)
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

d <- metaExpr(codes = c('CH','CG','W_L'))

#According to WGCNA_sft.R, inverse hyperbolic sine transformation performs the best by scale-free fit. power = 7
d_hst <- log(d+sqrt(d*d+1))
d_hst = quantile_normalisation(d_hst)
datExpr = t(d_hst)
softPower = 7;

runWGCNA(datExpr,softPower,"focal_metasample")

d <- metaExpr(codes = c('RH','RG','W_L'))

#According to WGCNA_sft.R, inverse hyperbolic sine transformation performs the best by scale-free fit. power = 9 for random nurses (fewer samples)
d_hst <- log(d+sqrt(d*d+1))
d_hst = quantile_normalisation(d_hst)
datExpr = t(d_hst)
softPower = 9;
runWGCNA(datExpr,softPower,"random_metasample")


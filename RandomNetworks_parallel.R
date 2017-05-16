##RandomNetworks_parallel.R 
##this file produces many networks with a randomly chosen set of genes
##And then tabulates connection strength between chosen genes across all networks

library(parallel)


load("~/Dropbox/monomorium nurses/data.processed/cleandata.RData")
fpkm = log(fpkm + sqrt(fpkm ^ 2 + 1)) #hyperbolic sine transformation to normalize gene expression data

#Create many random networks for a given sample set
RandomNetworks <- function(codes,names,boots,nGene){
  input <- GetExpr(codes,names)
  resList <- parallelGenie(codes,names,boots,nGene,input)
  # connStrength <- tabulateResList(resList)
  # return(connStrength)
}

#Parallel wrapper for genie function
parallelGenie <- function(codes,names,boots,nGene,input){
  nReps = floor(boots/bootsPerCore)
  no_cores <- detectCores() - 1
  
  # Initiate cluster
  cl <- makeCluster(no_cores)
  clusterExport(cl = cl, 
                varlist = c(
                  "codes","names","boots","nGene",
                  "input","runGenie","bootsPerCore")) ##Must export these variables for parLapply to see them
  
  # In parallel, go through all permutations
  AllResults = parLapply(cl,1:nReps,function(k) runGenie(codes,names,nGene,input,k))
  stopCluster(cl)
  return(AllResults)
}

#Run bootsPerCore number of Genie runs on each thread
#This is better than full parallelization so we can deal with the try error
runGenie <- function(codes,names,nGene,input,run){
  Results = list()
  setwd("~/Downloads/GENIE3_R_C_wrapper")
  source("~/Downloads/GENIE3_R_C_wrapper/GENIE3.R")
  library(plyr)
  i = 1
  while (i <= bootsPerCore){
    Genes = sample(unique(input$gene),nGene,replace=F) #Pick random genes
    Ginput = input[input$gene %in% Genes,] #Get random gene expression data
    rownames(Ginput) = Ginput$tissue_gene
    x <- try(GENIE3(as.matrix(Ginput[,c(1:5)])))
    if (!inherits(x,"try-error")){ #Sometimes the data get weird, just try again!
      Results[[i]] = get.link.list(x)
      i = i+1
    }
  }
  Results = ldply(Results)
  return(Results)
}

#Take results of many runs and tabulate them
tabulateResList <- function(resList){
  res <- ldply(resList,data.frame)
  d <- ddply(res,c("target.gene","regulatory.gene"),
             Conn = mean(weight))
  d = d[order(d$weight,decreasing=TRUE),]
  return(d)
}

#Get expression over time for given set of samples; label genes by their sample association
GetExpr <- function(codes,names){
  d = list()
  for (i in 1:length(codes)){
    d[[i]] = ExprTime(codes[i],names[i])
  }
  data = ldply(d,data.frame)
  
  data$tissue_gene=with(data,paste(Samp,gene,sep="_"))
  return(data)
}

#Gets expression over time for a given sample
ExprTime <- function(code,name){
  expr = matrix(nrow=nrow(fpkm))
  for (i in 1:5){
    e=rowSums(fpkm[,grepl(code,colnames(fpkm)) & factors$stage %in% i])/ncol(fpkm[,grepl(code,colnames(fpkm)) & factors$stage %in% i])
    expr=cbind(expr,e)
  }
  expr = expr[,-c(1)]
  colnames(expr)=paste("Samp",seq(1,5,by=1),sep="")
  expr = as.data.frame(expr)
  expr$Samp = rep(name,nrow(expr))
  expr$gene=rownames(fpkm)
  return(expr)
}

boots = 1000000
nGene = 10
bootsPerCore = 100
codes = c("W.*_L","C.*WH","C.*WG")
names = c("WorkLarv","WorkNurseH","WorkNurseG")
Wconns = RandomNetworks(codes,names,boots,nGene)
save(Wconns,file="WconnsParallel.RData")


codes = c("LS|1LW","XH|1LCH","XG|1LCG")
names = c("SexLarv","SexNurseH","SexNurseG")
Sconns = RandomNetworks(codes,names,boots,nGene)
save(Sconns,file="SconnsParallel.RData")


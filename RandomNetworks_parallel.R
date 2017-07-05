##RandomNetworks_parallel.R 
##this file produces many networks with a randomly chosen set of genes
##And then tabulates connection strength between chosen genes across all networks

library(parallel)
library(plyr)


load("~/cleandata.RData")

##Sort for 1000 most connected genes to reduce dataset size
corM = cor(t(fpkm))
conn = colSums(abs(corM))
conn = conn[order(conn,decreasing=TRUE)]
keep = names(conn)[1:1000]
fpkm = fpkm[rownames(fpkm) %in% keep,]

fpkm = log(fpkm + sqrt(fpkm ^ 2 + 1)) #hyperbolic sine transformation to normalize gene expression data

#Create many random networks for a given sample set
RandomNetworks <- function(){
  parallelGenie()
  return()
}




#Parallel wrapper for genie function
parallelGenie <- function(){
  nReps = floor(boots/bootsPerCore)
 
  # Initiate cluster; this only works on linux
  cl <- makePSOCKcluster(40,
                         master=system("hostname -i", intern=TRUE))
  
  clusterExport(cl = cl, c(
                  "name","nGene","input","runGenie","bootsPerCore")) ##Must export these variables for parLapply to see them
  
  # In parallel, go through all permutations
  p <- parLapply(cl,1:nReps, function(k) {
    runGenie(k)
  })
  stopCluster(cl)
  return()
}

#Run bootsPerCore number of Genie runs on each thread
#This is better than full parallelization so we can deal with the try error
runGenie <- function(run){
  setwd("~/GENIE3_R_C_wrapper")
  source("~/GENIE3_R_C_wrapper/GENIE3.R")
  Results = list()
  library(plyr)
  i = 1
  while (i <= bootsPerCore){
    Genes = sample(rownames(input),nGene,replace=F) #Pick random genes
    Ginput = input[rownames(input) %in% Genes,] #Get random gene expression data
    x <- try(GENIE3(as.matrix(Ginput)))
    if (!inherits(x,"try-error")){ #Sometimes the data get weird, just try again!
      Results[[i]] = get.link.list(x)
      i = i+1
    }
  }
  Results = ldply(Results)
  save(Results,file=paste("~/Results/",name,run,"GenieParallel.RData",sep=""))
  return()
}

setInput <- function(codes,names){
  AllExpr = list()
  for (i in 1:5){
    nSamp = getNsamp(codes,i)
    expr = list()
    for (j in 1:length(codes)){
      f = factors[grepl(codes[j],rownames(factors))&factors$stage==i,]
      samps = sample(rownames(f),nSamp,replace=FALSE) ##Want same number of samples per sample type
      expr[[j]] = fpkm[,samps]
      expr[[j]]$tissue=names[j]
      expr[[j]]$gene = rownames(expr[[j]])
      expr[[j]]$tissue_gene = with(expr[[j]],paste(tissue,gene,sep="_"))
      colnames(expr[[j]]) = c(paste("Stage",i,"_",seq(1,nSamp),sep=""),"tissue","gene","tissue_gene")
    }
    AllExpr[[i]] = ldply(expr,data.frame)
    rownames(AllExpr[[i]]) = AllExpr[[i]]$tissue_gene
    AllExpr[[i]] = AllExpr[[i]][,c(1:nSamp)]
  }
  input <- do.call(cbind,AllExpr)
  return(input)
}

getNsamp <- function(codes,stage){
  nSamp = c()
  for (code in codes){
    f = factors[grepl(code,rownames(factors))&factors$stage==stage,]
    nSamp=c(nSamp,nrow(f))
  }
  return(min(nSamp))
}

boots = 1000000
nGene = 10
bootsPerCore = 500
d <- fpkm
input = d[,grepl("LS|W.*_L",colnames(d))]
rownames(input) = rownames(fpkm)
name = "Larva"
#RandomNetworks()

codes = c("W.*_L","C.*WH","C.*WG")
names = c("WorkLarv","WorkNurseH","WorkNurseG")
input = setInput(codes,names)
name = "WorkerNet"
#RandomNetworks()

codes = c("1LW|LS","1LCH|XH","1LCG|XG")
names = c("SexLarv","SexNurseH","SexNurseG")
input = setInput(codes,names)
name = "SexualNet"
RandomNetworks()

# input <- GetExpr(codes,names)
# worker = RandomNetworks("Worker")
# save(worker,file="WorkerGenieParallel.RData")
# 
# codes = c("LS|1LW","XH|1LCH","XG|1LCG")
# names = c("SexLarv","SexNurseH","SexNurseG")
# sexual = RandomNetworks("Sexual")
# save(sexual,file="SexualGenieParallel.RData")


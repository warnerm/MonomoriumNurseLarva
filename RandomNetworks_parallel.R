##RandomNetworks_parallel.R 
##this file produces many networks with a randomly chosen set of genes
##And then tabulates connection strength between chosen genes across all networks

library(parallel)
library(plyr)


load("~/cleandata.RData")
fpkm = log(fpkm + sqrt(fpkm ^ 2 + 1)) #hyperbolic sine transformation to normalize gene expression data

#Create many random networks for a given sample set
RandomNetworks <- function(){
  file <- parallelGenie()
  return(file)
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
  parLapply(cl,1:nReps, function(k) {
    runGenie(k)
  })
  stopCluster(cl)
  return()
}

#Run bootsPerCore number of Genie runs on each thread
#This is better than full parallelization so we can deal with the try error
runGenie <- function(run){
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
  save(Results,file=paste(run,name,"~/Results/GenieParallel.RData"))
  return()
}

source("~/GENIE3_R_C_wrapper/GENIE3.R")

boots = 10
nGene = 10
bootsPerCore = 2
d <- fpkm
input = d[,grepl("LS|W.*_L",colnames(d))]
rownames(input) = rownames(fpkm)
name = "Larva"
runGenie()
#RandomNetworks()

# codes = c("W.*_L","C.*WH","C.*WG")
# names = c("WorkLarv","WorkNurseH","WorkNurseG")
# input <- GetExpr(codes,names)
# worker = RandomNetworks("Worker")
# save(worker,file="WorkerGenieParallel.RData")
# 
# codes = c("LS|1LW","XH|1LCH","XG|1LCG")
# names = c("SexLarv","SexNurseH","SexNurseG")
# sexual = RandomNetworks("Sexual")
# save(sexual,file="SexualGenieParallel.RData")


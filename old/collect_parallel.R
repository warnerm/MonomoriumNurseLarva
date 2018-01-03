#Collects RandomNetworks_parallel results, summarizes data frame to get mean connection weight for pairings
library(rslurm)
library(plyr)
args <- commandArgs(TRUE)

name <- args[1]
#Load in results of networks with top 1000 connected genes
load(paste(name,"InitialData.RData",sep="_"))
genes <- rownames(fpkm)
files <- list.files(pattern=paste(name,".*results.csv",sep="")) #List all output genie files
files = data.frame(file=as.character(files),name=name)
fun <- function(file,name){
  load(paste("../",as.character(name),"_InitialData.RData",sep=""))
  genes <- rownames(fpkm)
  library(plyr)
  d <- read.csv(paste("../",as.character(file),sep=""))
  res <- vector("list",length(genes))
  names(res) = genes
  for (gene in genes){
    res[[gene]]=vector("list",9)
    names(res[[gene]])=apply(expand.grid(names,names),1,function(x) paste(x[1],x[2],sep=".")) #List all 9 combinations (L-L, L-WH, etc)
    d2 = d[gsub(".*_","",d$regulatory.gene) %in% gene,] #subset dataframe for rows in which the regulatory gene is present
    if (nrow(d2)==0){
      next;
    }
    for (j in 1:nrow(d2)){
      samp = paste(gsub("_.*","",d2$regulatory.gene[j]),gsub("_.*","",d2$target.gene[j]),sep=".") #"Sample" is the regulatory connection between tissues...i.e. NurseHeadH.WorkLarv is WH->L
      res[[gene]][[samp]]=c(res[[gene]][[samp]],d2$weight[j]) #Connection strength from genie
    }
  }
  res
}

sjob <- slurm_apply(fun, files, jobname = 'collect_parGenie',
                    nodes = 4, cpus_per_node = 15, submit = TRUE)

res <- get_slurm_out(sjob,outtype='raw') #get output as lists
save(res,file="test.RData")

resL = vector("list",length=nrow(fpkm))
#Combine all lists so that for each gene we have a vector of regulatory observations
for (gene in genes){
  resL[[gene]]=list()
  for (col in apply(expand.grid(names,names),1,function(x) paste(x[1],x[2],sep="."))){
    for (i in 1:length(res)){
      resL[[gene]][[col]]=c(resL[[gene]][[col]],res[[i]][[gene]][[col]])
    }
  }
}

resL2 <- vector("list",length=nrow(fpkm))
names(resL2) = rownames(fpkm)
b=FALSE
for (gene in rownames(fpkm)){
  resL2[[gene]]['gene']=gene
  for (col in apply(expand.grid(names,names),1,function(x) paste(x[1],x[2],sep="."))){
    resL2[[gene]][col] = mean(resL[[gene]][[col]],na.rm=TRUE)
  }
}

resL3 <- ldply(resL2)
write.csv(resL3,file=paste(name,"connStrengths.csv",sep="_"))

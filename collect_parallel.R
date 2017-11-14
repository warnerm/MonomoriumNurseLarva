#Collects RandomNetworks_parallel results, summarizes data frame to get mean connection weight for pairings
library(rslurm)
library(plyr)
#args <- commandArgs(TRUE)

#name <- args[1]
#Load in results of networks with top 1000 connected genes
#load(paste(name,"InitialData.RData",sep="_")) #load initial codes, names, and fpkm
load(paste("worker_500000_1000_InitialData.RData",sep="_"))
genes <- rownames(fpkm)
name='worker_500000_1000'
files <- list.files(pattern=paste(name,".*csv",sep="")) #List all output genie files
files = data.frame(file=as.character(files[1:6]))
WorkerRes = list()
fun <- function(file){
  library(plyr)
  load(paste("worker_500000_1000_InitialData.RData",sep="_")) #load initial codes, names, and fpkm
  genes <- rownames(fpkm)
  print(as.character(file))
  d <- read.csv(as.character(file))
  d <- d[1:500,]
  res <- vector("list",length(genes))
  names(res) = genes
  for (i in 1:length(res)){
    res[[i]]=vector("list",9)
    names(res[[i]])=apply(expand.grid(names,names),1,function(x) paste(x[1],x[2],sep=".")) #List all 9 combinations (L-L, L-WH, etc)
    d2 = d[gsub(".*_","",d$regulatory.gene) %in% genes[i],] #subset dataframe for rows in which the regulatory gene is present
    for (j in 1:nrow(d2)){
      samp = paste(gsub("_.*","",d2$regulatory.gene[j]),gsub("_.*","",d2$target.gene[j]),sep=".") #"Sample" is the regulatory connection between tissues...i.e. NurseHeadH.WorkLarv is WH->L
      res[[i]][[samp]]=c(res[[i]][[samp]],d2$weight[j]) #Connection strength from genie
    }
  }
  return(res)
}


sjob <- slurm_apply(fun, files, jobname = 'collect_parGenie',
                    add_objects = c("names","genes"),
                    nodes = 2, cpus_per_node = 3, submit = TRUE)

res <- get_slurm_out(sjob,outtype='raw') #get output as lists
save(res,file="test.RData")

resL <- do.call(Map,c(c,res)) #Combine all lists so that for each gene we have a vector of regulatory observations

resL2 <- vector("list",length=nrow(fpkm))
names(resL2) = rownames(fpkm)
for (gene in rownames(fpkm)){
  resL2[[gene]]=matrix(nrow=1,ncol=10)
  colnames(resL2[[gene]])=c("gene",apply(expand.grid(names,names),1,function(x) paste(x[1],x[2],sep=".")))
  resL2[[gene]]['gene']=gene
  for (col in colnames(resL2[[gene]])[2:10]){
    resL2[[gene]][col] = mean(resL[[gene]][[col]])
  }
}

resL3 <- ldply(resL2,matrix)
write.csv(resL3,file=paste(name,"connStrengths.csv",sep="_"))

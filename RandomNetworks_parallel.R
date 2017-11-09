##RandomNetworks_parallel.R 
##this file produces many networks with a randomly chosen set of genes
##And then tabulates connection strength between chosen genes across all networks

library(parallel)
library(plyr)

load(paste(name,"InitialData.RData",sep="_")) #load initial codes, names, and fpkm

nGene = 10
bootsPerCore = 500

#This is better than full parallelization so we can deal with the try error
runGenie <- function(run){
  setwd("~/GENIE3_R_C_wrapper")
  source("~/GENIE3_R_C_wrapper/GENIE3.R")
  Results = list()
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
  df <- read.csv(paste(name,"results.csv",sep="_"))
  df = df[,-c(1)]
  df <- rbind(df,Results)
  write.csv(df,file=paste(name,"results.csv",sep="_"))
  return()
}

#subset fpkm based on which samples requested, how many genes
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

input <- setInput(codes,names)
RandomNetworks()

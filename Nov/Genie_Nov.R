library(plyr)
library(preprocessCore)
library(rslurm)
load("cleandata.RData")
load("genes_for_genie.RData")

#Return list of expression matrices for correlation analysis
formatExpr <- function(dN,dL){
  d = lapply(c(dN,dL),subExpr) #get expression matrix for each of three sample types
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
  colnames(f) = substr(colnames(f),start=1,stop=4) #Has colony, stage and queen presence information
  return(f)
}

#Take formatted expression data and remove instances of colonies that don't match (missing data)
alignStage <- function(d){
  d[[1]]=d[[1]][,colnames(d[[1]]) %in% colnames(d[[2]])]
  d[[2]]=d[[2]][,colnames(d[[2]]) %in% colnames(d[[1]])]
  return(d)
}

#Tabulate connection strengths
tabGenie <- function(geneList,samp){
  Ngene = unique(geneList$regulatory.gene[geneList$reg==samp])
  n = data.frame(Gene = Ngene,reg_within = 0,reg_between = 0,targ_within=0,targ_between=0)
  for (i in 1:length(Ngene)){
    n$reg_within[i] = sum(geneList$weight[geneList$regulatory.gene==Ngene[i] & geneList$targ==samp])
    n$reg_between[i] = sum(geneList$weight[geneList$regulatory.gene==Ngene[i] & geneList$targ!=samp])
    n$targ_within[i] = sum(geneList$weight[geneList$target.gene==Ngene[i] & geneList$reg==samp])
    n$targ_between[i] = sum(geneList$weight[geneList$target.gene==Ngene[i] & geneList$reg!=samp])
  }
  return(n)
}

selectSocial <- function(codeN,gN,gL,boots=100){
  setwd("~/GENIE3_R_C_wrapper") #Have to switch directories because there are .so files we need
  source("~/GENIE3_R_C_wrapper/GENIE3.R")
  expr = formatExpr(codeN,'W_L')
  keep = list(gN,gL)
  expr <- lapply(c(1,2),function(i) expr[[i]][rownames(expr[[i]]) %in% keep[[i]],])
  rownames(expr[[1]]) = paste("nurse",rownames(expr[[1]]),sep="_")
  rownames(expr[[2]]) = paste('larv',rownames(expr[[2]]),sep="_")
  
  #Global variable for slurm
  allExpr <<- rbind(expr[[1]],expr[[2]])
  df <- data.frame(run=seq(1,boots,by=1))
  sjob <- slurm_apply(runGenie, df, jobname = 'parGenie',
                      nodes = 1, cpus_per_node = 20, add_objects=c("allExpr"),submit = TRUE)
  res <- get_slurm_out(sjob,outtype='raw') #get output as lists
  Averaged <- Reduce("+", res) / length(res)
  geneList <- get.link.list(Averaged)
  geneList$reg = gsub("_LOC.*","",geneList$regulatory.gene)
  geneList$targ = gsub("_LOC.*","",geneList$target.gene)
  nurse <- tabGenie(geneList,"nurse")
  larv <- tabGenie(geneList,"larv")
  results <- rbind(nurse,larv)
  write.csv(Averaged,paste("~/Nurse_Larva/DEC18",codeN,"GenieMat.csv",sep=""))
  write.csv(results,file=paste("~/Nurse_Larva/DEC18",codeN,"GenieTabConn.csv",sep=""))
}

runGenie <- function(run){
  setwd("~/GENIE3_R_C_wrapper") #Have to switch directories because there are .so files we need
  source("~/GENIE3_R_C_wrapper/GENIE3.R")
  while (TRUE){
    x <- try(GENIE3(as.matrix(allExpr)))
    if (!inherits(x,"try-error")){ #Sometimes the data get weird, just try again!
      return(x)
    }
  }
}

#fpkm_norm = as.data.frame(normalize.quantiles(as.matrix(fpkm)))
#colnames(fpkm_norm) = colnames(fpkm)
#rownames(fpkm_norm) = rownames(fpkm)
fpkm = log(fpkm+sqrt(fpkm^2+1))
selectSocial("CH",keepH,keepL)
selectSocial("CG",keepG,keepL)
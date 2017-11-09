args <- commandArgs(TRUE)
#First arg is fpkm, second is factors file, third is sample subset, fourth is number of bootstraps, fifth is number of genes
fpkm <- read.csv(args[1])
factors <- read.csv(args[2])
samp <- args[3]
boots <- args[4]
totGene <- args[5]

##Sort for N highest expressed genes to reduce dataset size
sortData <- function(N,fpkm){
  rownames(fpkm) = fpkm[,1]
  fpkm <- fpkm[,-c(1)]
  rowS = rowSums(fpkm)
  keep = rowS[order(rowS,decreasing=TRUE)]
  keep = names(keep)[1:N]
  fpkm = fpkm[keep,]
  fpkm = log(fpkm + sqrt(fpkm ^ 2 + 1)) #hyperbolic sine transformation to normalize gene expression data
  
  return(fpkm)
}

#Derive sample codes and names based on the samp command-line argument
deriveCodes <- function(samp){
  if (samp=="worker"){
    codes = c("W.*_L","C.*WH","C.*WG")
    names = c("WorkLarv","WorkNurseH","WorkNurseG")
  } else if (samp=="random"){
    codes = c("QW.*_L","R.*WH","R.*WG")
    names = c("WorkLarvQR","RandNurseH","RandNurseG")
  } else if (samp=="forager"){
    codes = c("W.*_L","F.*WH","F.*WG")
    names = c("WorkeLarv","ForagerH","ForagerG")
  } else if (samp=="workerQR"){
    codes = c("QW.*_L","QCH","QCG")
    names = c("WorkLarvQR","WorkNurseHQR","WorkNurseGQR")
  }
  return(list(codes,names))
}

fpkm <- sortData(totGene,fpkm)
c = deriveCodes(samp)
codes = c[[1]]
names = c[[2]]
name = paste(samp,boots,totGene,sep="_")
save(codes,names,fpkm,factors,file=paste(name,"InitialData.RData",sep="_"))

#intialize dataframes for which to store results
results <- matrix(nrow=1,ncol=3)
colnames(results) = c("regulatory.gene","target.gene","weight")
write.csv(results,file=paste(name,"results.csv",sep="_"))








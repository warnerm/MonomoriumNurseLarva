library(rslurm)

#Get top genes by CV of larvae and a nurse tissue, keeping the top N/2 of each by CV
getTopGenes <- function(CV1,CV2,N){
  top1 = names(CV1)[order(CV1,decreasing=TRUE)][1:floor(N/2)] #First select top N/2 genes from larvae and nurse tissue
  top2 = names(CV2)[order(CV2,decreasing=TRUE)][1:floor(N/2)]
  keep = c(top1,top2)
  keep = keep[!duplicated(keep)] #Some top genes may be duplicated
  samp = 1
  while (length(keep) < N){ #Take turns adding genes until we have N genes
    if (samp %% 2 == 1){
      CV1 = CV1[!names(CV1) %in% keep]
      keep = c(keep,names(CV1)[order(CV1,decreasing=TRUE)][1])
    } else {
      CV2 = CV2[!names(CV2) %in% keep]
      keep = c(keep,names(CV2)[order(CV2,decreasing=TRUE)][1])
    }
    samp = samp + 1
  }
  return(keep)
}

#Select top N genes by CV for analysis; label genes by where they are expressed
selectGene <- function(N,codes,names){
  fpkm <- read.csv("~/Nurse_Larva/fpkm.csv")
  rownames(fpkm) = fpkm[,1]
  fpkm <- fpkm[,-c(1)]
  fpkm = log(fpkm + sqrt(fpkm ^ 2 + 1))
  f1 <- fpkm[,grepl(codes[1],colnames(fpkm))]
  f2 <- fpkm[,grepl(codes[2],colnames(fpkm))]
  CV1 = apply(f1,1, sd)/apply(f1,1,mean)
  CV2 = apply(f2,1, sd)/apply(f2,1,mean)
  keep = getTopGenes(CV1,CV2,N)
  f1 = f1[keep,]
  f2 = f2[keep,]
  rownames(f1) = paste(rownames(f1),names[1],sep="")
  rownames(f2) = paste(rownames(f2),names[2],sep="")
  colnames(f1) = gsub("[A-Z]+_.*","",colnames(f1))
  colnames(f2) = gsub("[A-Z]+_.*","",colnames(f2))
  f2 = f2[,colnames(f2) %in% colnames(f1)]
  f1 = f1[,colnames(f1) %in% colnames(f2)]
  return(rbind(f1,f2))
}

#Run Genie algorithm on dataframe
runGenie <- function(run){
  setwd("~/GENIE3_R_C_wrapper") #Have to switch directories because there are .so files we need
  source("~/GENIE3_R_C_wrapper/GENIE3.R")
  done = 0
  while (TRUE){
    x <- try(GENIE3(as.matrix(f)))
    if (!inherits(x,"try-error")){ #Sometimes the data get weird, just try again!
      return(x)
    }
  }
}

#Using slurm, run Genie repeatedly on a range of N values on cluster
allBoots <- function(boots){
  res = list()
  runs = data.frame(run = seq(1,N,by=1))
  sjob <- slurm_apply(runGenie, runs, jobname = 'parGenie',
                      nodes = 4, cpus_per_node = 20, add_objects="f",submit = TRUE)
  res <- get_slurm_out(sjob,outtype='raw') #get output as lists
  return(res)
}

codes = c("QW","QCH")
names = c("WorkLarvQR","NurseH")
Ns <- c(100,500,1000,2500,5000)
boots = 10
workQRH = list()
for (N in Ns){
  f = selectGene(N,codes,names)
  workQRH[[N]] = allBoots(boots)
}
save(workQRH,file="GenieWorkQRH.RData")

codes = c("QW","QCG")
names = c("WorkLarvQR","NurseG")
workQRG = list()
for (N in Ns){
  f = selectGene(N,codes,names)
  workQRG[[N]] = allBoots(boots)
}
save(workQRG,file="GenieWorkQRG.RData")

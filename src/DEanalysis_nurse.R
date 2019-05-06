library(plyr)
library(edgeR)
library(reshape2)

#Load counts, factors
ext <- read.csv("~/Writing/Data/NurseSpecialization_transcriptomicData/MpharAnn.csv") #load in MBE results
load("~/Dropbox/monomorium nurses/data.processed/ps_genelevelJuly29.RData")
load("~/Dropbox/monomorium nurses/data.processed/cleandata.RData")
a = TAIgene$Mphar_E5 
ext <- merge(ext,a,by.x="Gene",by.y="gene")

#Filter for genes with phylostrata calls (they are long enough)
keep = ext$Gene[!is.na(ext$ps)]
counts = counts[rownames(counts) %in% keep,] 

####Inital steps of standard edgeR analysis, as per manual
EdgeR <- function(data,design,coef){
  data <- DGEList(counts=data)
  data <- calcNormFactors(data)
  dat <- estimateGLMTrendedDisp(data, design)
  dat <- estimateGLMTagwiseDisp(dat, design)
  fit <- glmFit(dat,design)
  diff <- glmLRT(fit, coef=coef) 
  out <- topTags(diff,n=Inf,adjust.method="BH")$table 	###Calculate FDR
  return(out)
}

#Test for genes with P value < 0.05 as candidates
stageGenes <- function(code){
  f = droplevels(factors[grepl(code,rownames(factors)),])
  design <- model.matrix(~stage+colony,data=f)
  nStage = length(levels(f$stage))
  d <- EdgeR(counts[,colnames(counts) %in% rownames(f)],design,2:nStage)
  return(list(rownames(d)[d$PValue < 0.05],d))
}

#Tests to identify genes with P < 0.05 by stage as candidates for further analysis
tests <- c("CH","CG","RH","RG","W_L")
DEgene <- lapply(tests,stageGenes)
names(DEgene) = tests

save(DEgene,file = "~/Data/Nurse_Larva/DEstage_results.RData")
lapply(DEgene,function(x) length(x[[1]]))

#Only keep genes which aren't DE for random nurses
keepH = DEgene[[1]][[1]][!DEgene[[1]][[1]] %in% DEgene[[3]][[1]]]
keepG = DEgene[[2]][[1]][!DEgene[[2]][[1]] %in% DEgene[[4]][[1]]]

Htop <- DEgene[[1]][[2]][rownames(DEgene[[1]][[2]]) %in% keepH,]
Gtop <- DEgene[[2]][[2]][rownames(DEgene[[2]][[2]]) %in% keepG,]

#Keep top 1000 genes for regulatory network reconstruction
keepH <- rownames(Htop)[1:1000]
keepG <- rownames(Gtop)[1:1000]
keepL <- rownames(DEgene[[5]][[2]])[1:1000]

save(keepH,keepG,keepL,file = "genes_for_genie.RData")






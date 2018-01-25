library(plyr)
library(edgeR)
library(reshape2)
library(rslurm)

#Load counts, factors
load("cleandata.RData")

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
  return(rownames(d)[d$PValue < 0.05])
}

#Tests to identify genes with P < 0.05 by stage as candidates for further analysis
tests <- c("QCH","QCG","CH","CG","RH","RG","W_L")
DEgene <- lapply(tests,stageGenes)
names(DEgene) = tests

#Identify genes with pearson correlations different from random genes
#'Social' genes have pearson cor CI non-overlapping with random nurse cor CI AND more abs(cor_focal) > abs(cor_random)
#'Control' genes have pearson cor CI non-overlapping AND abs(cor_rand) > abs(cor_focal)
testCor <- function(){
  Lgenes = rownames(fpkm)
  df = data.frame(geneL = Lgenes)
  sjob <- slurm_apply(getSocConns, df, jobname = 'parCor',
                     nodes = 4, cpus_per_node = 20, add_objects=c("data","genes","getCorCI","conType","compareCor"),submit = TRUE)
  res <- get_slurm_out(sjob,outtype='raw') #get output as lists
  res <- ldply(res,rbind) #Each job is a row (i.e. connections with nurse genes for each larval gene)
  colnames(res) = genes
  rownames(res) = Lgenes
  return(res)
}

#Function to parallelize to get social connections
getSocConns <- function(geneL){
  res = vector()
  for (j in 1:length(genes)){
    res[j] = conType(genes[j],as.character(geneL))
  }
  return(res)
}

#Return list of expression matrices for correlation analysis
formatExpr <- function(larv,foc,rand){
  d = lapply(c(larv,foc,rand),subExpr) #get expression matrix for each of three sample types
  for (i in 2:3){
    #Needed to align random nurses to QR larvae
    colnames(d[[i]]) = gsub("R","Q",colnames(d[[i]]))
  }
  dF = alignStage(list(d[[1]],d[[2]]))
  dR = alignStage(list(d[[1]],d[[3]]))
  return(list(dF,dR)) #Returns data for foc/larv comparison and rand/larv comparison
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

getCorCI <- function(d,GeneN,GeneL){
  dL = d[[1]][rownames(d[[1]])==GeneL,]
  dN = d[[2]][rownames(d[[2]])==GeneN,]
  if (sd(dN) == 0) return(c(0,0)) #some genes not expressed at all
  test = cor.test(t(dL),t(dN))
  return(test$conf.int[1:2])
}

#Test if pearson correlation CI differ for focal and random nurses
conType <- function(GeneN,GeneL){
  corFoc <- getCorCI(data[[1]],GeneN,GeneL)
  corRand <- getCorCI(data[[2]],GeneN,GeneL)
  return(compareCor(corFoc,corRand))
}

#Compare correlations between focal vs random nurses and larvae
compareCor <- function(cF,cR){
  if (abs(min(cF)) > abs(max(cR))){ #We know we have significant social gene; provided both values of cF have same sign
    if (cF[1] > 0 & cF[2] > 0) return("posSocial")
    else if (cF[1] < 0 & cF[2] < 0) return("negSocial")
  } else if (abs(max(cF)) < abs(min(cR))){ #We know we have significant control gene; provided both values of cR have same sign
    if (cR[1] > 0 & cR[2] > 0) return("posControl")
    else if (cR[1] < 0 & cR[2] < 0) return("negControl")
  }
  #If none of these gates are met, it is a non-significant difference
  return("NS")
}

data <- formatExpr("QW","QCH","RH")
genes <- DEgene[[1]]
corQCH <- testCor()
data <- formatExpr("QW","QCG","RG")
genes <- DEgene[[2]]
corQCG <- testCor()
data <- formatExpr("W_L","CH","RH")
genes <- DEgene[[3]]
corCH <- testCor()
data <- formatExpr("W_L","CG","RG")
genes <- DEgene[[4]]
corCG <- testCor()

#By treating focal nurses as 'random' and random as focal (by switching their order),
#we can make another test that focal nurses indeed exhibit more correlations
data <- formatExpr("QW","RH","QCH")
genes <- DEgene[[5]]
corRH <- testCor()
data <- formatExpr("QW","RG","QCG")
genes <- DEgene[[6]]
corRG <- testCor()

save(corQCH,corQCG,corCH,corCG,corRH,corRG,DEgene,file='CorResults.RData')





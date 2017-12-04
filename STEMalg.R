#The following algorithm is based on Ernst et al. 2005 Bioinformatics "Clustering short time series gene expression data"

#This file sorts genes into dummy profiles separately for nurses and larvae.

library(rslurm)
library(plyr)
library(combinat)

###Step 1: intialize dummy profiles (according to Ernst et al)
SetProfiles <- function(timepoints){
  profiles <- matrix(ncol=timepoints)
  profiles = as.data.frame(profiles)
  colnames(profiles) = paste("Stage",seq(1,timepoints,1))
  posMoves = c(-1,0,1)
  for (i2 in -1:1){
    for (i3 in -2:2){
      for (i4 in -3:3){
        if (timepoints == 5){
          for (i5 in -4:4){
            if ((i2 - i3) %in% posMoves & (i3 - i4) %in% posMoves & (i4 - i5) %in% posMoves){
              profiles = rbind(profiles,c(0,i2,i3,i4,i5))
            }
          }
        } else if (timepoints == 4){
          if ((i2 - i3) %in% posMoves & (i3 - i4) %in% posMoves){
            profiles = rbind(profiles,c(0,i2,i3,i4))
          }
        }
        
      }
    }
  }
  profiles = profiles[-c(1),]
  profiles = profiles[rowSums(abs(profiles))!=0,]
  return(profiles)
}

###As input to the algorithm, we need the mean expression across samples of a given type at each developmental stage
CountsbyStage <- function(code){
  load("cleandata.RData")
  counts <- counts[,grep(code,colnames(counts))]
  factors <- factors[grep(code,rownames(factors)),]
  StageExpr <- matrix(nrow=nrow(counts),ncol=timepoints)
  stages = seq(6-timepoints,5,1) ##for sexuals, start at 2
  for (j in stages){
    st = j-5+timepoints ##Converts actual stage to relative (so sexuals go 1:4 instead of 2:5)
    for (i in 1:nrow(counts)){
      StageExpr[i,st] = mean(as.numeric(counts[i,factors$stage==j])) ##Calculate initial expression
    }
  }
  for (j in 2:timepoints){
    for (row in 1:nrow(StageExpr)){
      StageExpr[row,j]=log2((StageExpr[row,j]+1)/(StageExpr[row,1]+1)) ##Expression at stage X is the log2 ratio of expression at stage X to expression at intitial stage
    }
  }
  StageExpr[,1]=0
  rownames(StageExpr) = rownames(counts)
  return(StageExpr)
}

####Sort genes into the profile they mostly closely resemble. 
SortGenes <- function(){
  GeneMembership <- c()
  for (i in 1:nrow(StageExpr)){
    member = GeneMember(StageExpr[i,]) ##For each gene, return number of profile it belongs to.
    GeneMembership = c(GeneMembership,member)
  }
  names(GeneMembership) = rownames(StageExpr)
  Membership = data.frame(Module=seq(1,nrow(profiles),1),Number=1)
  for (i in 1:nrow(Membership)){
    Membership$Number[i] = sum(GeneMembership==i,na.rm=T) ##Returns number of genes belonging to a given module
  }
  return(list(Membership,GeneMembership))
}

##Function identifies, on a gene-by-gene basis, the closest profile match
GeneMember <- function(expr){
  diss = 1
  member = NA
  for (j in 1:nrow(profiles)){
    dist = 1 - cor(as.numeric(profiles[j,]),as.numeric(expr)) ## 1 - pearson correlation
    if (is.na(dist)){
      next;
    }
    if (dist < diss){
      diss = dist ##update distance, membership if correlation is stronger than previously tested profiles
      member = j
    }
  }
  return(member)
}

##Permute stage labels for each sample. Calculates null distribution for membership of each profile given resampled data
NullMembership <- function(){ 
  p = data.frame(k=seq(1,nrow(stagePerm),by=1))
  
  # In parallel, go through all permutations
  sjob <- slurm_apply(NullSort, p, jobname = 'parSTEM',
                      nodes = 4, cpus_per_node = 10, 
                      add_objects=c("NullSort","profiles","GeneMember","timepoints","StageExpr"),
                      submit = TRUE)
  res <- get_slurm_out(sjob,outtype='raw') #get output as lists

  GeneMembership = unlist(res) ##Collapse all results to one list
  Expected = data.frame(Module=seq(1,nrow(profiles),1),Number=1)
  for (i in 1:nrow(Expected)){
    Expected$Number[i] = sum(GeneMembership==i,na.rm=T)/nrow(stagePerm) ##GeneMembership stores results for all permutations, so find the mean number of genes for each profile
  }
  return(Expected)
}

##Assign genes to profiles in null distributions
NullSort <- function(k){
  GeneMembership = c()
  for (i in 1:nrow(StageExpr)){ ##One gene at a time
    expr = StageExpr[i,]
    member = GeneMember(expr) ##Sort genes using permuted data
    GeneMembership = c(GeneMembership,member)
  }
  return(GeneMembership)
}


##Based on real results and null distributions, identify profiles to retain as significant/dominant
SigProfiles <- function(Expected,Membership){
  nGene = nrow(StageExpr)
  sigProf = allprob = c()
  for (i in 1:nrow(Membership)){
    pi = Expected$Number[i]/nGene ##Null hypothesis mean based on resampled data
    prob = binom.test(Membership$Number[i],nGene,pi,alternative="greater")$p.value
    allprob = c(allprob,prob)
    if (prob < 0.1/(nrow(Membership)-1)){  ##Retain profiles with FDR < 0.1
      sigProf = c(sigProf,i)
    }
  }
  return(list(sigProf,allprob))
}

##Full algorithm
Alg <- function(code,name){
  Sort <- SortGenes() ##Sort genes into profiles
  Expected <- NullMembership() ##Generate null distribution
  sigProf <- SigProfiles(Expected,Sort[[1]])[[1]] ##find significant profiles
  GeneMembership <- Sort[[2]]
  GeneMembership[!GeneMembership %in% sigProf] = NA ##Set membership to NA if the module is non-significant
  results = list(GeneMembership=GeneMembership,Expected=Expected)
  save(results,file=paste(name,"STEMmembership.RData"))
}

profiles <- SetProfiles(5) ##Initialize the same profiles for all samples with 5 stages
timepoints = 5
stages = seq(6-timepoints,5,1)
stagePerm = ldply(permn(stages)) ##All possible permutations of stages
codes = c("W.*_L","QW","CH","CG","QCH","QCG","R.*_WH","R.*_WG")
names = c("WLarv","WlarvQR","NurseH","NurseG",
          "NurseHQR","NurseGQR","RNurseH","RNurseG")
names(names) = codes
for (code in codes){
  StageExpr <- CountsbyStage(code)
  Alg(names[code])
}
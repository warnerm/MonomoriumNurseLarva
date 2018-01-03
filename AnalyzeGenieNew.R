library(gtools)
library(rslurm)

#Get correlation vector for a single path. Returns vector of length boots - 1
getCor <- function(i){
  p = sample(1:boots,boots,replace=FALSE)
  corL = c()
  meanVal = dN[[p[1]]]
  diag(meanVal)=NA
  meanValN = meanVal
  nGene = nrow(meanVal)
  for (j in p[-c(1)]){
    for (row in 1:nGene){ #Update mean conn strengths for each gene-gene pair
      for (col in 1:nGene){
        if (row==col){ #Skip these; they will end up being NA
          next;
        }
        meanValN[row,col] = (meanVal[row,col]*(j-1) + dN[[j]][row,col])/2 #Average of previous mean and new results, weighted by number of values that went into mean calculation   
      }
    }
    m1 = c(meanVal)[!is.na(c(meanVal))]
    m2 = c(meanValN)[!is.na(c(meanValN))]
    c = cor(c(m1),c(m2))
    corL = c(corL,c) #Calculate average correlation to previous mean value
  }
  return(corL)
}

#Get vector of average correlations of connections
getCorTrend <- function(){
  b = data.frame(i = seq(1,nTry,by=1))
  sjob <- slurm_apply(getCor, b, jobname = 'dGenie',
                      nodes = 4, cpus_per_node = 20, add_objects=c("dN","boots"),submit = TRUE)
  res <- get_slurm_out(sjob,outtype='raw') #get output as lists
  res2 = do.call(cbind,res)
  #res = list()
  #for (i in 1:nTry){
  #  res[[i]]=getCor(i)
  #}
  m = c1 = c2 = c()
  for (i in 1:(boots-1)){
    m[i]=rowSums(res2)[i]/ncol(res2)
    c1[i] = quantile(res2[i,],0.025)
    c2[i] = quantile(res2[i,],0.975)
  }
  corL = list(MEAN = m, C1 = c1, C2 = c2)
  return(corL)
}

runSlurm <- function(N){
  res = list()
  runs = data.frame(run = seq(1,N,by=1))
  sjob <- slurm_apply(runGenie, runs, jobname = 'parGenie',
                      nodes = 4, cpus_per_node = 20, add_objects="f",submit = TRUE)
  res <- get_slurm_out(sjob,outtype='raw') #get output as lists
  return(res)
}

nTry = 100
v = seq(1,boots,by=1)
perm = permutations(n=boots,r=boots,v=v,repeats.allowed = FALSE)
Ns <- c(100,500,1000,2500,5000)
boots = 10

corVecH = list()
load("GenieWorkQRH.RData")
d = workQRH
for (N in Ns){
  dN = d[[N]]
  corVecH[[N]]=getCorTrend()
}
save(corVecH,file="CorVecH.RData")

corVecG = list()
load("GenieWorkQRG.RData")
d = workQRG
for (N in Ns){
  dN = d[[N]]
  corVecG[[N]]=getCorTrend()
}
save(corVecG,file="CorVecG.RData")




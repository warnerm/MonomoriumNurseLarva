library(gtools)

#Get correlation vector for a single path. Returns vector of length boots - 1
getCor <- function(i){
  p = perm[i,]
  corL = c()
  meanVal = dN[[i]]
  nGene = nrow(meanVal)
  for (j in 2:boots){
    corL = c(corL,sum(colSums(cor(meanVal,dN[[j]])))/(nGene*nGene)) #Calculate average correlation to previous mean value
    for (row in 1:nGene){ #Update mean conn strengths for each gene-gene pair
      for (col in 1:nGene){
        meanVal[row,col] = (meanVal[row,col]*(j-1) + dN[[j]][row,col])/2 #Average of previous mean and new results, weighted by number of values that went into mean calculation   
      }
    }  
  }
  return(corL)
}

#Get vector of average correlations of connections
getCorTrend <- function(){
  b = data.frame(i = seq(1,boots,by=1))
  sjob <- slurm_apply(getCor, b, jobname = 'dGenie',
                      nodes = 4, cpus_per_node = 20, add_objects="dN",submit = TRUE)
  res <- get_slurm_out(sjob,outtype='raw') #get output as lists
  res2 = do.call(cbind,res)
  m = c1 = c2 = c()
  for (i in 1:(boots-1)){
    m[i]=colSums(res2)[i]/nrow(res2)
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

boots = 100
v = seq(1,boots,by=1)
perm = permutations(n=boots,r=boots,v=v,repeats.allowed = FALSE)
Ns <- c(10,30,50,100,250,500,750,1000)
Ns <- c(10,30,50)


corVecH = list()
load("GenieWorkQRHs.RData")
d = workQRH
for (N in Ns){
  dN = d[[N]]
  corVecH[[N]]=getCorTrend()
}
save(corVecH,file="CorVecH.RData")

load("GenieWorkQRGs.RData")
d = workQRG
for (N in Ns){
  dN = d[[N]]
  corVecG[[N]]=getCorTrend()
}
save(corVecG,file="CorVecG.RData")




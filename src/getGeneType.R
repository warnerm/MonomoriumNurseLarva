library(reshape)
library(reshape2)
library(plyr)

load("CorResults.RData") #After running 'corApproach.R' on cluster
#Tabulate type (posSocial, negSocial, etc of social interactions)
tabType <- function(d){
  d = as.data.frame(t(d))
  d$Gene = rownames(d)
  m <- melt(d,id.vars = "Gene")
  m2 <- cast(m,Gene~value,length)
  return(m2)
}

#Binomial test to see if gene is social
socBinom <- function(i){
  test = binom.test(m2$social[i],(m2$social[i]+m2$control[i]),p = 0.5)
  m2$geneTypePval[i] = test$p.value
  if (test$p.value < 0.05){ #Two sided test
    if (test$estimate > 0.5) m2$geneType[i] = return(c(test$p.value,"social"))
    else m2$geneType[i] = return(c(test$p.value,"control"))
  }
  return(c(test$p.value,'NS'))
}

#Determine whether or not genes are 'social' by lumping positive and negative connections
#Given that some correlation disparity happened, use binomial test to say whether we see more social connections than control
socGene <- function(m2){
  m2$social = m2$posSocial + m2$negSocial
  m2$control = m2$posControl + m2$negControl
  df = data.frame(i = seq(1,nrow(m2),by=1))
  sjob <- slurm_apply(socBinom, df, jobname = 'geneType',
                      nodes = 4, cpus_per_node = 20, add_objects=c("m2"),submit = TRUE)
  res <- get_slurm_out(sjob,outtype='raw') #get output as lists
  res = ldply(res,rbind) #Each job is a row (i.e. connections with nurse genes for each larval gene)
  colnames(res) = c('p.value','geneType')
  return(cbind(m2,res))
}

allDat <- list(corQCH,corQCG,corCH,corCG,t(corQCH),t(corQCG),t(corCH),t(corCG),corRH,corRG) #Last two get Larva-Nurse social genes
table <- lapply(allDat,tabType)
socDet <- lapply(table,socGene)
names(socDet) = c("QCH","QCG",'CH','CG','LARV_QCH','LARV_QCG','LARV_CH','LARV_CG','RH','RG')
save(socDet,file="socDet.RData")




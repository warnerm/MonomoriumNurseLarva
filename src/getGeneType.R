library(reshape)
library(reshape2)
library(plyr)
library(rslurm)

setwd("~/Nurse_Larva")
load("CorResults.RData") #After running 'corApproach.R' on cluster

getSocConns <- function(i){
  social = sum(grepl("Social",t(d[i,])))
  control = sum(grepl("Control",t(d[i,])))
  if (social + control > 0) pvalue = binom.test(social,social+control,p = 0.5)$p.value
  else pvalue = NA
  return(data.frame(Gene = rownames(d)[i],Pvalue = pvalue, Excess = social - control,
             posSocial = sum(grepl("posSocial",d[i,])),negSocial = sum(grepl("posSocial",d[i,])),
             posControl = sum(grepl("posControl",d[i,])),negControl = sum(grepl("posControl",d[i,]))))
}

#Tabulate type (posSocial, negSocial, etc of social interactions)
tabType <- function(df){
  d <<- as.data.frame(t(df)) #Global variable for slurm
  var <- data.frame(i = seq(1,nrow(d),by=1))
  sjob <- slurm_apply(getSocConns, var, jobname = 'parGetType',
                      nodes = 4, cpus_per_node = 20, add_objects=c("d"),submit = TRUE)
  res <- get_slurm_out(sjob,outtype='raw') #get output as lists
  res <- ldply(res,rbind) #Each job is a row (i.e. connections with nurse genes for each larval gene)
  return(res)
}

allDat <- list(corQCH,corQCG,corCH,corCG,t(corQCH),t(corQCG),t(corCH),t(corCG),corRH,corRG) #Last two get Larva-Nurse social genes
socDet <- lapply(allDat,tabType)
names(socDet) = c("QCH","QCG",'CH','CG','LARV_QCH','LARV_QCG','LARV_CH','LARV_CG','RH','RG')
save(socDet,file="socDet.RData")





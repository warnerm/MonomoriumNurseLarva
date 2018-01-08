library(reshape)
library(reshape2)
library(plyr)
#library(rslurm)

setwd("~/Data/Nurse_Larva")
load("CorResults.RData") #After running 'corApproach.R' on cluster
#Tabulate type (posSocial, negSocial, etc of social interactions)
tabType <- function(df){
  d <<- as.data.frame(t(df)) #Global variable for slurm
  results <- list()
  for (i in 1:nrow(d)){
    social = sum(grepl("Social",t(d[i,])))
    control = sum(grepl("Control",t(d[i,])))
    if (social + control > 0) pvalue = binom.test(social,social+control,p = 0.5)$p.value
    else pvalue = NA
    results[[i]] = data.frame(Gene = rownames(d)[i],Pvalue = pvalue, Excess = social - control,
                     posSocial = sum(grepl("posSocial",d[i,])),negSocial = sum(grepl("posSocial",d[i,])),
                     posControl = sum(grepl("posControl",d[i,])),negControl = sum(grepl("posControl",d[i,])))
  }
  results = ldply(results)
  return(results)
}

allDat <- list(corQCH,corQCG,corCH,corCG,t(corQCH),t(corQCG),t(corCH),t(corCG),corRH,corRG) #Last two get Larva-Nurse social genes
socDet <- lapply(allDat,tabType)
names(socDet) = c("QCH","QCG",'CH','CG','LARV_QCH','LARV_QCG','LARV_CH','LARV_CG','RH','RG')
save(socDet,file="socDet.RData")

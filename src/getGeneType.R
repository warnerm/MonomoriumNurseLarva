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

#Determine whether or not genes are 'social' by lumping positive and negative connections
#Given that some correlation disparity happened, use binomial test to say whether we see more social connections than control
socGene <- function(m2){
  m2$social = m2$posSocial + m2$negSocial
  m2$control = m2$posControl + m2$negControl
  m2$geneType = "NS"
  m2$geneTypePval = 1
  for (i in 1:nrow(m2)){
    if (m2$social[i] > 0){
      test = binom.test(m2$social[i],(m2$social[i]+m2$control[i]),p = 0.5)
      m2$geneTypePval[i] = test$p.value
      if (test$p.value < 0.05){ #Two sided test
        if (test$estimate > 0.5) m2$geneType[i] = "social"
        else m2$geneType[i] = "control"
      }
    }
  }
  return(m2)
}

allDat <- list(corQCH,corQCG,corCH,corCG,t(corQCH),t(corQCG)) #Last two get Larva-Nurse social genes
table <- lapply(allDat,tabType)
socDet <- lapply(table,socGene)
save(socDet,file="socDet.RData")




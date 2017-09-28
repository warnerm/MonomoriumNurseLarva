#Collects RandomNetworks_parallel results, summarizes data frame to get mean connection weight for pairings
library(plyr)

#Load in results of networks with top 1000 connected genes
setwd("~/Results/")

collate <- function(name){
  files <- list.files(pattern=paste(name,".*RData",sep=""))
  for (i in 1:length(files)){
    load(files[i])
    WorkerRes[[i]] = ddply(Results,~regulatory.gene + target.gene,summarize,parallel=TRUE,
                           N = length(weight),
                           meanW = mean(weight))
  }
  WorkerRes = ldply(WorkerRes,data.frame)
  #Calculate pairwise mean connection values
  WRsum = ddply(WorkerRes,~regulatory.gene + target.gene,summarize,parallel=TRUE,
                N = sum(N),
                meanW = mean(meanW))
  WRsum = WRsum[order(WRsum$meanW,decreasing=TRUE),]
  stat <- c(mean(WRsum$N),ci1=quantile(WRsum$N,0.025),ci2=quantile(WRsum$N,0.975))
  write.csv(WRsum,"CompiledConnections.csv")
  write.table(stat,"MeanNumberConns.txt")

  WRsoc = data.frame(gene = unique(gsub(".*_","",WRsum$regulatory.gene)))
  WRsoc$Lwithin=WRsoc$Lbetween=WRsoc$WHwithin=WRsoc$WHbetween=WRsoc$WGwithin=WRsoc$WGbetween=WRsoc$nNet=0
  for (i in 1:nrow(WRsoc)){
    d = WRsum[grepl(WRsoc$gene[i],WRsum$regulatory.gene),]
    WRsoc$Lwithin[i] = mean(d$meanW[grepl("Larv",d$target.gene)&grepl("Larv",d$regulatory.gene)])
    WRsoc$Lbetween[i] = mean(d$meanW[(!grepl("Larv",d$target.gene))&grepl("Larv",d$regulatory.gene)])
    WRsoc$WHwithin[i] = mean(d$meanW[grepl("NurseH",d$target.gene)&grepl("NurseH",d$regulatory.gene)])
    WRsoc$WHbetween[i] = mean(d$meanW[(grepl("Larv",d$target.gene))&grepl("NurseH",d$regulatory.gene)])
    WRsoc$WGwithin[i] = mean(d$meanW[grepl("NurseG",d$target.gene)&grepl("NurseG",d$regulatory.gene)])
    WRsoc$WGbetween[i] = mean(d$meanW[(grepl("Larv",d$target.gene))&grepl("NurseG",d$regulatory.gene)])
    WRsoc$WG.WH[i] = mean(d$meanW[(grepl("NurseH",d$target.gene))&grepl("NurseG",d$regulatory.gene)])
    WRsoc$WH.WG[i] = mean(d$meanW[(grepl("NurseG",d$target.gene))&grepl("NurseH",d$regulatory.gene)])
    WRsoc$nNet[i] = nrow(d)  
  }
  write.csv(WRsoc,paste(name,"SocialityDF.csv",sep=""))
}

collate("TopExprWorkerNet")


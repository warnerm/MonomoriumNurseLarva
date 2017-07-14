#Collects RandomNetworks_parallel results, summarizes data frame to get mean connection weight for pairings
library(plyr)

#Load in results of networks with top 1000 connected genes
setwd("~/Results")

collate <- function(name){
  files <- list.files(pattern=paste(name,".*RData",sep=""))
  load(files[1])
  WorkerRes = Results
  for (i in 2:length(files)){
    load(files[i])
    WorkerRes = rbind(WorkerRes,Results)
  }
  print(head(WorkerRes))
  #Calculate pairwise mean connection values
  WRsum = ddply(WorkerRes,~regulatory.gene + target.gene,summarize,
                N = length(weight),
                meanW = mean(weight))
  WRsum = WRsum[order(WRsum$meanW,decreasing=TRUE),]
  WRsum$regT = WRsum$targT = "Larva"
  WRsum$regT[!grepl("Larv",WRsum$regulatory.gene)]="Worker"
  WRsum$targT[!grepl("Larv",WRsum$target.gene)]="Worker"
  betWR = WRsum[WRsum$regT!=WRsum$targT,]
  betWR = betWR[c(1:10000),]
  write.csv(betWR,file=paste(name,"TopConns.csv",sep=""))
  print(head(WRsum))
  WRsum$targReg =with(WRsum,paste0(regulatory.gene,target.gene))
  #Calculate "socialility index" as mean conn outside individual - mean conn inside tissue
  
  write.csv(WRsum,file=paste(name,"WRsum.csv",sep=""))
  WRsoc = data.frame(gene = unique(gsub(".*_","",WRsum$regulatory.gene)))
  WRsoc$Lwithin=WRsoc$Lbetween=WRsoc$WHwithin=WRsoc$WHbetween=WRsoc$WGwithin=WRsoc$WGbetween=0
  WRsocMax = WRsoc
  for (i in 1:nrow(WRsoc)){
    d = WRsum[grepl(WRsoc$gene[i],WRsum$regulatory.gene),]
    WRsoc$Lwithin[i] = mean(d$meanW[grepl("Larv",d$target.gene)&grepl("Larv",d$regulatory.gene)])
    WRsoc$Lbetween[i] = mean(d$meanW[(!grepl("Larv",d$target.gene))&grepl("Larv",d$regulatory.gene)])
    WRsoc$WHwithin[i] = mean(d$meanW[grepl("NurseH",d$target.gene)&grepl("NurseH",d$regulatory.gene)])
    WRsoc$WHbetween[i] = mean(d$meanW[(grepl("Larv",d$target.gene))&grepl("NurseH",d$regulatory.gene)])
    WRsoc$WGwithin[i] = mean(d$meanW[grepl("NurseG",d$target.gene)&grepl("NurseG",d$regulatory.gene)])
    WRsoc$WGbetween[i] = mean(d$meanW[(grepl("Larv",d$target.gene))&grepl("NurseG",d$regulatory.gene)])
    WRsocMax$Lwithin[i] = max(d$meanW[grepl("Larv",d$target.gene)&grepl("Larv",d$regulatory.gene)])
    WRsocMax$Lbetween[i] = max(d$meanW[(!grepl("Larv",d$target.gene))&grepl("Larv",d$regulatory.gene)])
    WRsocMax$WHwithin[i] = max(d$meanW[grepl("NurseH",d$target.gene)&grepl("NurseH",d$regulatory.gene)])
    WRsocMax$WHbetween[i] = max(d$meanW[(grepl("Larv",d$target.gene))&grepl("NurseH",d$regulatory.gene)])
    WRsocMax$WGwithin[i] = max(d$meanW[grepl("NurseG",d$target.gene)&grepl("NurseG",d$regulatory.gene)])
    WRsocMax$WGbetween[i] = max(d$meanW[(grepl("Larv",d$target.gene))&grepl("NurseG",d$regulatory.gene)])
  }
  print(head(WRsocMax))
  
  write.csv(WRsoc,paste(name,"NetSocialityDF.csv",sep=""))
  write.csv(WRsocMax,paste(name,"NetSocialityDFmax.csv",sep=""))
}

collate("TopExprWorkerNet")
collate("TopExprSexualNet")
collate("TopConnWorkerNet")
collate("TopConnSexualNet")



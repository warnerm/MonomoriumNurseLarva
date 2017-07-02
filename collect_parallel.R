#Collects RandomNetworks_parallel results, summarizes data frame to get mean connection weight for pairings
library(plyr)

#Load in results of networks with top 1000 connected genes
setwd("~/Results/Worker_Results")
files <- list.files()
load(files[1])
WorkerRes = Results
for (i in 2:length(files)){
  load(files[i])
  WorkerRes = rbind(WorkerRes,Results)
}

#Calculate pairwise mean connection values
WRsum = ddply(WorkerRes,~regulatory.gene + target.gene,summarize,
              N = length(weight),
              meanW = mean(weight))

WRsum$targReg = apply(WRsum[,("regulatory.gene","target.gene")],1,paste,sep="-")

#Calculate "socialility index" as mean conn outside individual - mean conn inside tissue

WRsoc = data.frame(gene = unique(gsub(".*_","",)))
WRsoc$Lwithin=WRsoc$Lbetween=WRsoc$WHwithin=WRsoc$WHbetween=WRsoc$WGwithin=WRsoc$WGbetween=0
for (i in 1:nrow(WRsoc)){
  d = WRsum[grepl(WRsoc$gene[i],WRsum$regulatory.gene),]
  WRsoc$Lwithin[i] = mean(d$meanW[grepl("Larv",d$target.gene)&grepl("Larv",d$regulatory.gene)])
  WRsoc$Lbetween[i] = mean(d$meanW[(!grepl("Larv",d$target.gene))&grepl("Larv",d$regulatory.gene)])
  WRsoc$WHwithin[i] = mean(d$meanW[grepl("WorkNurseH",d$target.gene)&grepl("WorkNurseH",d$regulatory.gene)])
  WRsoc$WHbetween[i] = mean(d$meanW[(grepl("Larv",d$target.gene))&grepl("WorkNurseH",d$regulatory.gene)])
  WRsoc$WGwithin[i] = mean(d$meanW[grepl("WorkNurseG",d$target.gene)&grepl("WorkNurseG",d$regulatory.gene)])
  WRsoc$WGbetween[i] = mean(d$meanW[(grepl("Larv",d$target.gene))&grepl("WorkNurseG",d$regulatory.gene)])
}

write.csv(WRsoc,"WorkerNetSocialityDF.csv")

#Keep top 5000 connections for further processing
keep = WRsum$targReg[order(WRsum$meanW,decreasing=TRUE)][1:5000]
WRsumTop = WRsum[WRsum$targReg %in% keep,]

write.csv(WRsumTop,file="TopConns.csv")\




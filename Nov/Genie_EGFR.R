
ann <- read.csv("~/Downloads/msx123_Supp (1)/MpharAnn.csv") #load in MBE results


jelly = as.character(ann$Gene[grepl("jelly",ann$SwissProt)])

egfr_genes2 <- c("LOC105828088","LOC105830494","LOC105834310",
                 "LOC105837907","LOC105838568","LOC105832621",
                 "LOC105838372","LOC105829781")

gl <- "LOC105830675"
spitz <- "LOC105833934"

go <- read.table("~/GitHub/devnetwork/data/dmel_ann.txt",sep="\t",header = FALSE,stringsAsFactors = FALSE)
load("~/GitHub/devnetwork/results/collectedPhylo.RData")
g = go[go$V2=="GO:0007173",] #EGFR signaling
egfr_genes <- as.character(ENDogg$gene_Mphar[ENDogg$gene_Dmel %in% g$V1])

test_genes <- unique(c(egfr_genes,egfr_genes2,gl,spitz,jelly))

makePlot <- function(nurse){
  df <- read.csv(paste("~/GitHub/MonomoriumNurseLarva/results/DEC16",nurse,"GenieTabConn.csv",sep=""))
  df$tissue = gsub("_.*","",df$Gene)
  df$Gene = gsub(".*_","",df$Gene)
  df$code = nurse
  df = df[,-c(1)]
  return(df)
}

codes <- c("CH","CG")
plotsG <- lapply(codes,makePlot)

allD <- ldply(plotsG,data.frame)
allDf <- merge(allD,f.est,by = "Gene",all.x = TRUE)

allDf$reg_diff = allDf$reg_between - allDf$reg_within
allDf$targ_diff = allDf$targ_between - allDf$targ_within
allDf$tissueC = "larva \u2192 nurse head"
allDf$tissueC[allDf$code=="CG" & allDf$tissue == "larv"] = "larva \u2192 nurse abdomen"
allDf$tissueC[allDf$code=="CH" & allDf$tissue == "nurse"] = "nurse head \u2192 larva"
allDf$tissueC[allDf$code=="CG" & allDf$tissue == "nurse"] = "nurse abdomen \u2192 larva"
allDf$tissueC = factor(allDf$tissueC,levels = c("larva \u2192 nurse head","larva \u2192 nurse abdomen","nurse head \u2192 larva","nurse abdomen \u2192 larva"))

# 
# ch <- allD[allD$tissue=="nurse" & allD$code=="CH",]
# ch$regR = rank(ch$reg_between)/nrow(ch)
# ch = merge(ch,ann[,c("Gene","SwissProt")],by="Gene")
# ch = ch[order(ch$regR,decreasing = T),]
# chT = ch[ch$Gene %in% test_genes,]
# 
# lh <- allD[allD$tissue=="larv" & allD$code=="CH",]
# lh$targR = rank(lh$targ_between)/nrow(lh)
# lh$regR = rank(lh$reg_between)/nrow(lh)
# 
# lh = merge(lh,ann[,c("Gene","SwissProt")],by="Gene")
# lh = lh[order(lh$targR,decreasing = T),]
# lhT = lh[lh$Gene %in% test_genes,] 
# 
# net_genes_nurse <- gl
# 
# netF = fpkm[rownames(fpkm) %in% test_genes,grepl("W_L|CH",colnames(fpkm))]
# netL = netF[,grepl("W_L",colnames(netF))]
# netN = netF[,grepl("CH",colnames(netF))]
# netN = netN[rownames(netN) %in% net_genes_nurse,]
# rownames(netL) = paste("Larv_",rownames(netL),sep="")
# rownames(netN) = paste("Nurse_",rownames(netN),sep="")
# 
# 
# #Take formatted expression data and remove instances of colonies that don't match (missing data)
# alignStage <- function(d){
#   colnames(d[[1]]) = substr(colnames(d[[1]]),start=1,stop=4) #Has colony, stage and queen presence information
#   colnames(d[[2]]) = substr(colnames(d[[2]]),start=1,stop=4) #Has colony, stage and queen presence information
#   d[[1]]=d[[1]][,colnames(d[[1]]) %in% colnames(d[[2]])]
#   d[[2]]=d[[2]][,colnames(d[[2]]) %in% colnames(d[[1]])]
#   return(d)
# }
# 
# expr <- alignStage(list(netN,netL))
# allE = rbind(expr[[1]],expr[[2]])
# 
# runGenie <- function(expr){
#   setwd("~/GENIE3_R_C_wrapper") #Have to switch directories because there are .so files we need
#   source("~/GENIE3_R_C_wrapper/GENIE3.R")
#   while (TRUE){
#     x <- try(GENIE3(as.matrix(log(expr+sqrt(expr+1)))))
#     if (!inherits(x,"try-error")){ #Sometimes the data get weird, just try again!
#       return(x)
#     }
#   }
# }
# 
# results <- runGenie(allE)
# links <- get.link.list(results)
# soc <- links[grepl("Nurse",links$regulatory.gene) & grepl("Larv",links$target.gene),]
# soc$target.gene = gsub("Larv_","",soc$target.gene)
# soc = merge(soc,ann[c("Gene","SwissProt")],by.x = "target.gene",by.y="Gene")
# soc = merge(soc,ENDogg,by.x="target.gene",by.y="gene_Mphar",all.x=T)
# soc = soc[order(soc$weight,decreasing=T),]


plot2 <- function(e){
  p <- ggplot(e,aes(x=stage,y=expression,color=sample))+
    geom_point()+
    #ylim(1,3.5)+
    geom_smooth(aes(group=sample,fill=sample),se=FALSE)+
    xlab("larval developmental stage")+
    #scale_color_manual(values = tissue_palette[c(1,3,4)])+
    #scale_fill_manual(values = tissue_palette[c(1,3,4)])+
    theme(legend.position = "none",
          panel.background = element_blank(),
          axis.line.x = element_line(color='black'),
          axis.line.y = element_line(color='black'),
          axis.title = element_text(size = 20, face = "bold"),
          axis.text = element_text(size = 15),legend.title = element_blank(),legend.text = element_text(size=15))
  
  return(p)
}


nPlot <- function(genes,codes,samps){
  nS = length(codes)
  f = lapply(codes,function(c) factors[grepl(c,rownames(factors)),])
  e = lapply(seq(1:nS),function(i) as.numeric(log(fpkm[genes[i],colnames(fpkm) %in% rownames(f[[i]])])))
  eA = lapply(seq(1:nS),function(i){
    x = e[[i]] - mean(e[[i]][f[[i]]$stage==1])
    return(x)
  })
  eD = lapply(seq(1:nS),function(i) data.frame(expression = eA[[i]],stage=f[[i]]$stage,sample=samps[i]))
  
  e = ldply(eD,data.frame)

  levels(e$stage) = c("L1","L2","L3","L4","L5")
  p <- plot2(e)
  return(p)
}

plB <- nPlot(c("LOC105830675","LOC105837907"),c("CH","W_L"),c("giant-lens (nurse head)","eps8 (larva)"))
#plB <- nPlot(c("LOC105830675","LOC105830675"),c("1LCH|XH","1LW|LS"),c("giant-lens (nurse head)","giant-lens (larva)"))

plB <- plB+theme(legend.position = c(0.2,0.2))


plB <- nPlot(c("LOC105837907","LOC105837907"),c("1LW|LS","W_L"),c("eps8 (sex larva)","eps8 (worker larva)"))

plC <- nPlot(c("LOC105833314","LOC105838382"),c("W_L","W_L"),c("importin-7 (worker larva)","asteroid (worker larva)"))
plC <- plC+theme(legend.position = c(0.6,0.8))
plD <- nPlot(c("LOC105836930","LOC105829992","LOC105831738"),rep("CH",3),c("mrjp-1 (nurse head)","mrjp-2 (nurse head)","mrjp-3 (nurse head)"))
plD <- plD+theme(legend.position = c(0.2,0.2))

tog <- nPlot(c("LOC105833314","LOC105838382","LOC105837907"),
             c("W_L","W_L","W_L"),c("importin-7 (larva)","asteroid (larva)","eps8 (larva)"))


tog <- nPlot(c("LOC105833314","LOC105838382","LOC105837907"),
             rep("1LW|LS",3),c("importin-7 (larva)","asteroid (larva)","eps8 (larva)"))

tog <- tog+theme(legend.position = c(0.2,0.2))
p <- arrangeGrob(plB,tog,nrow=1)
ggsave(p,file="~/GitHub/MonomoriumNurseLarva/Figures/glEGFR_time.png",height=6,width=12,dpi=300)


f2d <- ggplot(allDf, aes(x = reg_within/max(allDf$reg_within), y = reg_between/max(allDf$reg_within), color = tissueC))+
  geom_point(alpha = 0.2)+
  geom_smooth(method = "lm",se = FALSE,size=1)+apatheme+
  ylab("social connectivity")+
  xlab("within-tissue connectivity")+
  scale_color_manual(values = tissue_palette,name = bquote(underline("connection type")))+
  theme(legend.position = c(0.8,0.85),
        axis.line.x = element_line(color='black'),
        axis.line.y = element_line(color='black'),
        legend.title = element_text(size = 17,face="bold",margin = margin(t=0,b=-5,r=0,l=0)),
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 15),
        legend.background = element_blank())

ogg <- read.csv("~/Writing/Data/NurseSpecialization_transcriptomicData/ThreeWayOGGMap.csv") #Import 3-way OGG map
df <- read.csv("~/Writing/Data/NurseSpecialization_transcriptomicData/Dmel_secreted.csv") #Derived from http://www.flyrnai.org/tools/glad/web/
df = df[,c(1,2)]
colnames(df)[2] = "Flybase"
key <- read.table("~/Writing/Data/NurseSpecialization_transcriptomicData/DmelKey.txt") #Generated from Drosophila melanogaster gff file
oggK = merge(ogg,key,by.x="gene_Dmel_FLYBASE",by.y="V1")
oggK = oggK[!duplicated(oggK$OGG),]
oggK$secreted = "no"
oggK$secreted[oggK$V2 %in% df$Flybase] = "yes"
oggK$secreted = factor(oggK$secreted,levels = c("yes","no"))
dS = merge(all,oggK,by.x="Gene",by.y="gene_Mphar")

dS$regR = rank(dS$reg_between)/nrow(dS)

dS1 = dS[dS$tissue=="nurse" & dS$code=="CH",]
dS2 = dS[dS$tissue=="nurse" & dS$code=="CG",]
wilcox.test(dS1$reg_between[dS1$secreted=="yes"],dS1$reg_between[dS1$secreted=="no"])
wilcox.test(dS2$reg_between[dS2$secreted=="yes"],dS2$reg_between[dS2$secreted=="no"])
lm <- glm(reg_between ~ secreted,data=dS2)

wilcox.test(dS1$targ_between[dS1$secreted=="yes"],dS1$targ_between[dS1$secreted=="no"])
wilcox.test(dS2$targ_between[dS2$secreted=="yes"],dS2$targ_between[dS2$secreted=="no"])
dS1 = dS1[order(dS1$regR,decreasing = T),]
dS2 = dS2[order(dS2$regR,decreasing = T),]
dSe1 = dS1[dS1$secreted=="yes" & dS1$SwissProt!="-",]
dSe2 = dS2[dS2$secreted=="yes"& dS2$SwissProt!="-",]

topG = data.frame(NH=dSe1$SwissProt[1:10],NG=dSe2$SwissProt[1:10])
colnames(topG) = c("nurse head","nurse abdomen")
makeTbl2(topG,"secGenes",8)

levels(dS$secreted) = c(" secreted   "," not secreted   ")

p <- ggplot(dS[dS$tissue=="nurse",],aes(x=tissueC,y=reg_between,fill=secreted))+
  geom_boxplot(outlier.shape = NA)+
  #geom_violin()+
  apatheme+
  ylab("social connectivity")+
  xlab("tissue")+
  coord_cartesian(ylim=c(0,0.3))+
  theme(legend.position="top",
        legend.title = element_blank(),
        axis.line = element_line(color="black"))+
  scale_fill_manual(values=c("grey30","grey60"))+
  annotate("text",x=1,y=0.28,label="*",size=12)+
  annotate("text",x=2,y=0.3,label="ns",size=12)


p <- ggplot(dS[dS$tissue=="nurse",],aes(x=tissueC,y=reg_within,fill=secreted))+
  geom_boxplot(outlier.shape = NA)+
  #geom_violin()+
  apatheme+
  ylab("social connectivity")+
  xlab("tissue")+
  #coord_cartesian(ylim=c(0,0.3))+
  theme(legend.position="top",
        legend.title = element_blank(),
        axis.line = element_line(color="black"))+
  scale_fill_manual(values=c("grey30","grey60"))+
  annotate("text",x=1,y=0.28,label="*",size=12)+
  annotate("text",x=2,y=0.3,label="ns",size=12)

pA <- arrangeGrob(p,plB,nrow=1)
ggsave(pA,file="~/GitHub/MonomoriumNurseLarva/Figures/Fig3_all.png",height=6,width=12,dpi=300)


# 
# plE <- nPlot(c("LOC105828088","LOC105833314","LOC105838382","LOC105840195","LOC105829992","LOC105831738"),
#              rep("1LW|LS",6),
#              c("giant-lens","importin-7","asteroid","mrjp-1","mrjp-2","mrjp-3"))+theme(legend.position=c(0.2,0.3))
# 
# plF <- nPlot(c("LOC105828088","LOC105833314","LOC105838382","LOC105840195","LOC105829992","LOC105831738"),
#              rep("W_L",6),
#              c("giant-lens","importin-7","asteroid","mrjp-1","mrjp-2","mrjp-3"))+theme(legend.position=c(0.2,0.3))
# 


getGene <- function(g){
  return(ann$SwissProt[ann$Gene==g])
}








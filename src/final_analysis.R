setwd("~/GitHub/MonomoriumNurseLarva/")
library(plyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(edgeR)
library(topGO)
library(multcomp)
library(gtable)
library(magrittr)

#Generate means and bootstrapped confidence intervals
bootCI <- function(v,boots){
  v = v[!is.na(v)]
  bV = unlist(lapply(seq(1,boots),function(x){
    mean(v[sample(x = seq(1,length(v)),size = length(v),replace=TRUE)])
  }))
  c1 = quantile(bV, 0.025)
  c2 = quantile(bV, 0.975)
  m = mean(v)
  return(data.frame(mean = m, c1 = c1, c2 = c2))
}

#Function to make, save table for supplemental files
makeTbl2 <- function(tbl,name,font){
  
  tt3 <- ttheme_minimal(
    core=list(
      fg_params=list(fontsize=font)),
    colhead=list(fg_params=list(fontface="bold",fontsize=10)),
    rowhead=list(fg_params=list(fontsize=10)))
  
  #Make table grob out of results
  resT <- tableGrob(tbl,theme=tt3,rows=NULL)
  grid.arrange(resT)
  g <- gtable_add_grob(resT,
                       grobs = segmentsGrob( # line across the bottom
                         x0 = unit(0,"npc"),
                         y0 = unit(0,"npc"),
                         x1 = unit(1,"npc"),
                         y1 = unit(0,"npc"),
                         gp = gpar(lwd = 3.0)),
                       t = 1, b = 1, l = 1, r = ncol(tbl))
  separators <- replicate(ncol(g) - 1,
                          segmentsGrob(x1 = unit(0, "npc"), gp=gpar(lty=2)),
                          simplify=FALSE)
  ## add vertical lines on the left side of columns (after 2nd)
  g <- gtable::gtable_add_grob(g, grobs = separators,
                               t = 2, b = nrow(g), l = seq_len(ncol(g)-1)+1)
  p <- grid.arrange(g)
  ggsave(p,file=paste("Figures/",name,".png",sep=""),width=10,height=4,dpi=300)
  return(p)
  
}

#Function for collecting genie results
collectRes <- function(nurse){
  df <- read.csv(paste("~/GitHub/MonomoriumNurseLarva/results/DEC16",nurse,"GenieTabConn.csv",sep=""))
  df$tissue = gsub("_.*","",df$Gene)
  df$Gene = gsub(".*_","",df$Gene)
  df$code = nurse
  df = df[,-c(1)]
  return(df)
}

#Generate one legend for a set of figures
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right","top")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none",
                                            axis.line.x = element_line(color='black'),
                                            axis.line.y = element_line(color='black'),
                                            plot.margin = unit(c(0.25,0.58,0.5,0.5),"cm")))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)),
                     "top" = arrangeGrob(legend,
                                         do.call(arrangeGrob, gl),
                                         ncol = 1,
                                         heights = unit.c(lheight,unit(1, "npc") - lheight))
                     
  )
  
  return(combined)
  
}

#Load counts, fpkm data from Warner et al 2017 (MBE)
load("data/cleandata.RData")

apatheme=theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        text=element_text(family='Arial'),
        legend.position=c("top"),
        axis.line.x = element_line(color='black'),
        axis.line.y = element_line(color='black'),
        axis.text=element_text(size=13),
        axis.title=element_text(size=17,face="bold"),
        legend.text=element_text(size=13),
        axis.title.y = element_text(margin = unit(c(0,0.5,0,0),"cm")),
        axis.title.x = element_text(margin = unit(c(0,0.5,0,0),"cm")))


PS_palette = c("#deebf7","#9ecae1","#4292c6","#2171b5","#084594","darkblue","black")
tissue_palette = c("cornflowerblue","royalblue4","firebrick1","firebrick4")


#Constraint, calculated in Warner et al 2019 (in prep) using MKtest 2.0
f.est <- read.csv("data/MKtestConstraintOneAlpha.csv")
colnames(f.est) = c("Gene","f")

#Previously generated M. pharaonis annotation (Warner et al. 2017 MBE)
ext <- read.csv("data/MpharAnn.csv") #load in MBE results

#Collect Genie results
codes <- c("CH","CG")
plotsG <- lapply(codes,collectRes)

allD <- ldply(plotsG,data.frame)

allDf <- merge(allD,f.est,by = "Gene",all.x = TRUE)

#Label tissue combinations, with arrows for plotting
allDf$tissueC = "larva \u2192 nurse head"
allDf$tissueC[allDf$code=="CG" & allDf$tissue == "larv"] = "larva \u2192 nurse abdomen"
allDf$tissueC[allDf$code=="CH" & allDf$tissue == "nurse"] = "nurse head \u2192 larva"
allDf$tissueC[allDf$code=="CG" & allDf$tissue == "nurse"] = "nurse abdomen \u2192 larva"
allDf$tissueC = factor(allDf$tissueC,levels = c("larva \u2192 nurse head","larva \u2192 nurse abdomen","nurse head \u2192 larva","nurse abdomen \u2192 larva"))

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

ggsave(f2d,file="~/GitHub/MonomoriumNurseLarva/Figures/FigS4.tiff",height=6,width=6.5,dpi=300)

df <- allDf[,c(1:3,9)]

#Normalize connectivity to top value of within-tissue connectivity (as plotted above)
df$reg_between = df$reg_between/max(df$reg_within)
df$reg_within = df$reg_within/max(df$reg_within)
colnames(df)[c(2:4)] = c("within-tissue connectivity","social connectivity","connection type")
write.table(df,file = "data_files/FigS4.txt",row.names=F)

allDf$SI = allDf$reg_between-allDf$reg_within #Definition of sociality index

#Correlation of selective constraint (1-f) and connectivity
cor.test(1-allDf$f[allDf$tissue=="nurse" & allDf$code == "CH"],allDf$reg_within[allDf$tissue=="nurse" & allDf$code == "CH"],method = "spearman")
cor.test(1-allDf$f[allDf$tissue=="nurse" & allDf$code == "CH"],allDf$reg_between[allDf$tissue=="nurse" & allDf$code == "CH"],method = "spearman")
cor.test(1-allDf$f[allDf$tissue=="nurse" & allDf$code == "CH"],allDf$SI[allDf$tissue=="nurse" & allDf$code == "CH"],method = "spearman")
cor.test(1-allDf$f[allDf$tissue=="nurse" & allDf$code == "CG"],allDf$reg_within[allDf$tissue=="nurse" & allDf$code == "CG"],method = "spearman")
cor.test(1-allDf$f[allDf$tissue=="nurse" & allDf$code == "CG"],allDf$reg_between[allDf$tissue=="nurse" & allDf$code == "CG"],method = "spearman")
cor.test(1-allDf$f[allDf$tissue=="nurse" & allDf$code == "CG"],allDf$SI[allDf$tissue=="nurse" & allDf$code == "CG"],method = "spearman")

#Get just nurse->larva results
allDN = droplevels(allDf[allDf$tissue!="larv",])
allDN$SI = allDN$reg_between- allDN$reg_within

#Rename for plotting aesthetics
levels(allDN$tissueC) = c("nurse head ","nurse abdomen ")

nurse_pallete=tissue_palette[c(3,4)]

f3a <- ggplot(allDN,aes( x = reg_within/max(allDN$reg_within), y = 1 - f, color = tissueC))+
  geom_point(alpha = 0.5)+
  geom_smooth(method = "lm",se=FALSE)+
  xlim(0,1)+ylim(0,1)+
  apatheme+
  ylab("constraint")+  
  scale_color_manual(values = nurse_pallete)+
  annotate("text",x = .98,y=0.25,label='bold("A")',parse=TRUE,size=8,color="gray29")+
  guides(shape = guide_legend(override.aes = list(size=3)))+
  xlab("within-tissue connectivity")+
  theme(legend.title = element_blank(),
        axis.text = element_text(size=10),
        axis.title = element_text(size = 14),
        plot.margin = unit(c(-0.5,0,0,0),"cm"),
        legend.justification = "center",
        axis.line = element_line(color = "black"))
  

f3b <- ggplot(allDN,aes( x = reg_between/max(allDN$reg_between), y = 1 - f, color = tissueC))+
  geom_point(alpha = 0.5)+
  geom_smooth(method = "lm",se=FALSE)+
  xlim(0,1)+ylim(0,1)+
  xlab("social connectivity")+
  apatheme+
  ylab("constraint")+  
  theme(legend.title = element_blank(),
        axis.text = element_text(size=10),
        axis.title = element_text(size = 14),
        axis.line = element_line(color = "black"))+
  scale_color_manual(values = nurse_pallete)+
  annotate("text",x = .98,y=0.25,label='bold("B")',parse=TRUE,size=8,color="gray29")

f3c <- ggplot(allDN,aes( x = SI, y = 1-f, color = tissueC))+
  geom_point(alpha = 0.5)+
  apatheme+
  ylim(0,1)+
  ylab("constraint")+
  geom_smooth(method = "lm",se=FALSE)+
  xlab("sociality index")+
  scale_x_continuous(limits = c(-1.25,0.25))+
  theme(legend.title = element_blank(),
        axis.text = element_text(size=10),
        axis.title = element_text(size = 14),
        axis.line = element_line(color = "black"))+
  scale_color_manual(values = nurse_pallete)+
  annotate("text",x = 0.23,y=0.25,label='bold("C")',parse=TRUE,size=8,color="gray29")


#Previously estimated phylostratigraphic categories (Warner et al. 2019 Nat Comm)
load("data/collectedPhylo.RData")
newP = merge(allDN,Aps,by="Gene",all.x=T)
df = newP[,c(1:3,10,6,8,14)]
levels(df$psName)[1] = "ancient"

#Convert to constraint (1-f)
df$f = 1-df$f
colnames(df)[c(2:6)] = c("within-tissue connectivity","social connectivity","sociality index","tissue","constraint")
df$tissue = as.factor(df$tissue)
levels(df$tissue) = c("nurse abdomen","nurse head")
write.table(df,file = "data_files/Fig4.txt",row.names=F)

#Change categorization for plotting 
newP$psName[newP$psName=="novel"]="ant"
newP$psName[newP$psName=="aculeata"] = "hymenoptera"
newP = droplevels(newP)
levels(newP$psName)[1] = "ancient"

#Statistics for sociality index and phylostrata
lm <- glm(SI ~ psName,data=newP)
summary(glht(lm,mcp(psName="Tukey")))
drop1(lm,test="Chi")

#Both f and ps affect sociality index
lm <- glm(SI ~ psName + f + tissueC, data = newP)
drop1(lm,.~.,test="Chi") 

#Summarize phylostrata across genes (i.e. average between nurse head and abdomen)
pm = ddply(newP, ~ Gene + psName,summarize,
           f = mean(f,na.rm = TRUE),
           SI = mean(SI),
           regW = mean(reg_within),
           regB = mean(reg_between)
           )
pm = pm[!is.na(pm$f),]

#Bootstrap confidence intervals for sociality index & constraint
fM = lapply(levels(pm$psName), function(x){
  cbind(bootCI(1-pm$f[pm$psName == x],1000),bootCI(pm$SI[pm$psName == x],1000))
})

fM = ldply(fM,data.frame)
colnames(fM)[c(1,4)] = c("f","SI")
fM$phylostrata = factor(levels(newP$psName),levels = levels(newP$psName))

f4b <- ggplot(fM,aes(x=SI,y=f))+
  geom_point(aes(fill=phylostrata),color="black",size = 6,pch=21)+
  geom_errorbar(aes(ymin = c1, ymax = c2,color=phylostrata),size=1,width=0.005)+
  geom_errorbarh(aes(xmin = c1.1, xmax = c2.1,color=phylostrata),size=1,height=0.005)+
  scale_fill_manual(values = PS_palette[c(2,3,4,6)])+
  scale_color_manual(values = PS_palette[c(2,3,4,6)])+
  apatheme+
  ylab("constraint")+
  theme(legend.position = "none",
        axis.line.x = element_line(color='black'),
        axis.line.y = element_line(color='black'),
        axis.title.y = element_text(margin = unit(c(0,0.3,0,0),"cm")),
        axis.title.x = element_text(margin = unit(c(0.3,0,0,0),"cm")))+
  xlab("sociality index")+
  annotate("text",x = -0.285, y = 0.818, label = 'italic("ant (31)")',parse=TRUE,size = 4,angle = 45)+
  annotate("text",x = -0.395, y = 0.82, label = 'italic("hymenopteran (67)")',parse=TRUE,size = 4,angle = 45)+
  annotate("text",x = -0.33, y = 0.873, label = 'italic("insect (275)")',parse=TRUE,size = 4,angle = 45)+
  annotate("text",x = -0.42, y = 0.89, label = 'italic("ancient (1120)")',parse=TRUE,size = 4,angle=45)

tiff("~/GitHub/MonomoriumNurseLarva/Figures/Fig4.tiff",height = 1558,width = 2250,res =260)
grid.arrange(grid_arrange_shared_legend(f3a,f3b,f3c,ncol=1,nrow=3,position = "top"),
             f4b+theme(plot.margin = unit(c(0.5,0.5,1,0.9),"cm"),
                       axis.line = element_line(color = "black"))+
               annotate("text",x = -0.425, y = 0.78, label = 'bold("D")',parse = TRUE,size = 8,color="gray29"),
             nrow=1,widths = c(0.35,0.65))
grid.text("more social",x = 0.62, y = 0.04, gp = gpar(fontface="italic",fontsize=18,col="red"))
grid.lines(arrow = arrow(angle = 30, length = unit(0.2, "inches"),
                 ends = "last", type = "open"),
           x=unit(c(0.7,0.9),"npc"),y=unit(c(0.04,0.04),"npc"),gp=gpar(lwd=2.5,lty="solid",col="red"))
dev.off()


#########
##Secreted proteins and social connectivty
#########

#Load in orthology map (Walsh et al. 2018 Animal Behavior)
ogg <- read.csv("data/ThreeWayOGGMap.csv") #Import 3-way OGG map
df <- read.csv("data/Dmel_secreted.csv") #Derived from http://www.flyrnai.org/tools/glad/web/. This is the list of secreted proteins
df = df[,c(1,2)]
colnames(df)[2] = "Flybase"
key <- read.table("data/DmelKey.txt") #Generated from Drosophila melanogaster gff file
oggK = merge(ogg,key,by.x="gene_Dmel_FLYBASE",by.y="V1")
oggK = oggK[!duplicated(oggK$OGG),]
oggK$secreted = "no"
oggK$secreted[oggK$V2 %in% df$Flybase] = "yes"
oggK$secreted = factor(oggK$secreted,levels = c("yes","no"))
dS = merge(allDN,oggK,by.x="Gene",by.y="gene_Mphar")
dS = merge(dS,ext,by="Gene")

dS$regR = rank(dS$reg_between)/nrow(dS)

dS1 = dS[dS$tissue=="nurse" & dS$code=="CH",]
dS2 = dS[dS$tissue=="nurse" & dS$code=="CG",]

#social connectivity based on whether or not protein is secreted
wilcox.test(dS1$reg_between[dS1$secreted=="yes"],dS1$reg_between[dS1$secreted=="no"])
wilcox.test(dS2$reg_between[dS2$secreted=="yes"],dS2$reg_between[dS2$secreted=="no"])

dS1 = dS1[order(dS1$regR,decreasing = T),]
dS2 = dS2[order(dS2$regR,decreasing = T),]

#Identify top secreted proteins for table
dSe1 = dS1[dS1$secreted=="yes" & dS1$SwissProt!="-",]
dSe2 = dS2[dS2$secreted=="yes"& dS2$SwissProt!="-",]

topG = data.frame(tissue = rep(c("nurse head","nurse abdomen"),each=10),social= c(dSe1$reg_between[1:10],dSe2$reg_between[1:10]),
                  Gene = c(as.character(dSe1$SwissProt[1:10]),as.character(dSe2$SwissProt[1:10])))
colnames(topG) = c("tissue","social connectivity","gene")
p <- makeTbl2(topG,"TableS6",8)
ggsave(p,file="Figures/TableS6.png",height=6,width=6.5,dpi=300)

#Rename for plotting aesthetics
levels(dS$secreted) = c(" secreted   "," not secreted   ")

df = dS[c("code","secreted","reg_between")]
colnames(df) = c("tissue","secreted in D. melanogaster?","social connectivity")
df$tissue = as.factor(df$tissue)
levels(df$tissue) = c("nurse abdomen","nurse head")
write.table(df,file="data_files/Fig3a.txt",row.names = F)

p <- ggplot(dS[dS$tissue=="nurse",],aes(x=tissueC,y=reg_between,fill=secreted))+
  geom_boxplot(outlier.shape = NA)+
  apatheme+
  ylab("social connectivity")+
  xlab("tissue")+
  coord_cartesian(ylim=c(0,0.3))+
  theme(legend.position = "top",
        panel.background = element_blank(),
        axis.line.x = element_line(color='black'),
        axis.line.y = element_line(color='black'),
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 15),legend.title = element_blank(),
        legend.text = element_text(size=15))+
  scale_fill_manual(values=c("grey30","grey60"))+
  annotate("text",x=1,y=0.28,label="P = 0.025",size=4)+
  annotate("text",x=2,y=0.28,label="P = 0.067",size=4)+
  annotate("text",x=0.6,y=0.29,label="A",size=12,fontface="bold",color="gray29")


#Plot genes across time for giant-lens, etc
plot2 <- function(e){
  p <- ggplot(e,aes(x=stage,y=expression,color=sample))+
    geom_point(size=2,alpha=0.5)+
    #ylim(1,3.5)+
    geom_smooth(aes(group=sample,fill=sample),se=FALSE,size=1.5)+
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

#Summarize expression data across time for plotting
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
  return(list(p,e))
}

#Get expression profiles for statistics
getExpr <- function(genes,codes){
  nS = length(codes)
  f = lapply(codes,function(c) factors[grepl(c,rownames(factors)),])
  e = lapply(seq(1:nS),function(i) as.numeric(log(fpkm[genes[i],colnames(fpkm) %in% rownames(f[[i]])])))
  
  eD = lapply(seq(1:nS),function(i) data.frame(expression = e[[i]],stage=f[[i]]$stage,colony=f[[i]]$colony,qp = f[[i]]$queen_presence,sample=codes[i]))
  
  eD = lapply(eD,function(x){
    x$col = apply(x[,c("stage","colony","qp")],1,paste,collapse='_')
    return(x)
  }) 
  return(eD)
}

plB <- nPlot(c("LOC105830675","LOC105837907"),c("CH","W_L"),c(" giant-lens (nurse head)   "," eps8 (larva)"))

d = plB[[2]]
write.table(d,"data_files/Fig3b.txt",row.names=F)
dM = dcast(d, expression ~ sample,value.var="expression")

#Expression is negatively correlated
cor.test(dM$expression[c(2:5)],dM$expression[c(6:9)],method="spearman")

plB <- plB[[1]]+theme(legend.position = "top")

plB <- plB + theme(legend.margin = margin(t=5,l=30,b=5,r=10))+
  scale_color_manual(values=tissue_palette[c(3,1)])+
  annotate("text",x=1,y=0.933,label="B",size=12,fontface="bold",color="gray29")

pA <- arrangeGrob(p,plB,nrow=1)
ggsave(pA,file="~/GitHub/MonomoriumNurseLarva/Figures/Fig3.tiff",height=4,width=12,dpi=300)


#Statistical significance for Fig 3b
eD <- getExpr(c("LOC105830675","LOC105837907"),c("CH","W_L"))
F3b = merge(eD[[1]],eD[[2]],by = "col")
cor.test(F3b$expression.x,F3b$expression.y,method="spearman") #Correlation of each worker-larva pair

eS = lapply(eD,function(x){
  e1 = ddply(x, ~stage, summarize,
        val = mean(expression)
        )
  e1$val - e1$val[1] #normalized to first stage
})


#Statistics for Fig S5
eD <- getExpr(c("LOC105830675","LOC105830675"),c("CH","W_L"))
wil = lapply(seq(1:5),function(i){
  wilcox.test(eD[[1]]$expression[eD])
})

plg <- nPlot(c("LOC105830675","LOC105830675"),c("CH","W_L"),c(" giant-lens (nurse head)   "," giant-lens (larva)   "))
d = plg[[2]]

write.table(d,file="data_files/FigS5.txt",row.names = F)
levels(d$sample) = c("nurse","larva")
wil = lapply(seq(1:5),function(i){
  wilcox.test(d$expression[d$sample=="nurse" & d$stage==paste("L",i,sep="")],
              d$expression[d$sample=="larva" & d$stage==paste("L",i,sep="")])
})
plg <- plg[[1]]+ theme(legend.position="top",legend.margin = margin(t=5,l=30,b=5,r=10))+
  scale_color_manual(values=tissue_palette[c(3,1)])+
  annotate("text",label = c("ns","ns","ns","**","**"),x=seq(1,5),y=-4,size=8)
ggsave(plg,file="~/GitHub/MonomoriumNurseLarva/Figures/FigS5.tiff",height=6,width=6,dpi=300)



#Fig S6
pls <- nPlot(c("LOC105837907","LOC105837907"),c("1LW|LS","W_L"),c(" eps8 (sex larva)   "," eps8 (worker larva)   "))

d = pls[[2]]
write.table(d,file="data_files/FigS6.txt",row.names = F)

levels(d$sample) = c("nurse","larva")
levels(d$stage) = c(1,2,3,4,5)

#Make stage an ordinal variable
d$stage = as.ordered(d$stage)

lm <- glm(expression ~ stage*sample,data=d)
drop1(lm,test="Chi")

pls <- pls[[1]] + theme(legend.position="top",legend.margin = margin(t=5,l=30,b=5,r=10))+
  scale_color_manual(values=tissue_palette[c(1,2)])
ggsave(pls,file="~/GitHub/MonomoriumNurseLarva/Figures/FigS6.tiff",height=6,width=6,dpi=300)

############
#####Gene Set Enrichment Analysis (GSEA)
###########
#Load d. melanogaster GO annotation
go <- read.csv("data/GOannotation.csv")
new <- list()
for (gene in unique(go$gene)){
  d = go[go$gene %in% gene,]
  new[[gene]]=as.character(d$GO)
}

#placeholder function
selectConn <- function(score){
  return(score > quantile(score,0.9))
}

#GSEA analysis. Takes in vector of all genes, with a score for social connection strength
GSEAfunc <- function(d){
  GOdata <- new("topGOdata",
                description="Simple session",ontology="BP",
                allGenes=d,geneSel=selectConn,
                nodeSize = 10,
                annot=annFUN.gene2GO,gene2GO=new)
  
  #Use scoreOrder = "decreasing" because the higher connection strengths are what we are after
  resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks",scoreOrder="decreasing")
  allRes <- GenTable(GOdata,KS=resultKS,numChar=100)
  return(allRes[allRes$KS < 0.05,]) 
}

NH = allDN$reg_between[allDN$tissueC=="nurse head "]
names(NH) = allDN$Gene[allDN$tissueC=="nurse head "]
GsNH = GSEAfunc(NH)

NG = allDN$reg_between[allDN$tissueC=="nurse abdomen "]
names(NG) = allDN$Gene[allDN$tissueC=="nurse abdomen "]
GsNG = GSEAfunc(NG)


colnames(GsNH)[6] = colnames(GsNG)[6] = "P-value"
makeTbl2(GsNH,"TableS3",10)
makeTbl2(GsNG,"TableS4",10)

#List of species used
df <- read.table("data/phylostratigraphy_codes.txt",sep="\t")
names(df) = c("Species","code","NCBI Taxonomy ID")
df = df[,-c(2)]
write.table(df,"data/phylostratigraphy_species.txt",row.names = FALSE)

#making file for P-values, sociality index
load("~/GitHub/MonomoriumNurseLarva/results/DEstage_results.RData")
Pvals = lapply(names(DEgene),function(x){
  data.frame(Gene=rownames(DEgene[[x]][[2]]),Pvalue=DEgene[[x]][[2]]$PValue,TissueM = x)
})

a = allDf[,c(1,3,6,7,11)]
aL = a[a$tissue=="larv",]
aCH = a[a$code=="CH" & a$tissue=="nurse",]
aCG = a[a$code=="CG" & a$tissue=="nurse",]
aL = merge(aL,Pvals[[5]],by="Gene",all.y=T)
aCH = merge(aCH,Pvals[[1]],by="Gene",all.y=T)
aCG = merge(aCG,Pvals[[2]],by="Gene",all.y=T)
a = do.call(rbind,list(aL,aCH,aCG))
a = a[a$Pvalue < 0.05,]
levels(a$TissueM) = c("worker larva","nurse head","nurse abdomen")
a$used = "yes"
a$used[is.na(a$reg_between)] = "no"
newP = merge(a,Aps,by="Gene",all.x=T)
all = merge(newP,ext,by="Gene",all.x=T)
all = merge(all,f.est,by="Gene",all.x=T)

all = all[,c(1,2,5:8,12,13,42)]
all = all[c(1,5,4,6,3,2,7:9)]
all$psName=as.character(all$psName)
all$psName[all$psName=="old"]="ancient"
all$f = 1-all$f #make it constraint (1-f)
colnames(all) = c("Gene","differential expression tissue","differential expression P-value","used in regulatory network reconstruction?",
                  "estimated regulatory direction","social connectivity","phylostrata category","SwissProt ID","constraint")
write.csv(all,file = "~/GitHub/MonomoriumNurseLarva/results/GeneResults.csv",row.names=F)

nurse = all[grepl("nurse",all$`differential expression tissue`),]
nurse = nurse[nurse$`used in regulatory network reconstruction?`=="yes",]
nurse = nurse[nurse$`social connectivity` > quantile(nurse$`social connectivity`,0.99),c(1,5,6,8)]
nurse = nurse[order(nurse$`social connectivity`,decreasing = TRUE),]
nurse$`social connectivity` = round(nurse$`social connectivity`,3)

p <- makeTbl2(nurse,"TableS5",10)
ggsave(p,file="~/GitHub/MonomoriumNurseLarva/Figures/TableS5.tiff",height=8,width=10,dpi=300)


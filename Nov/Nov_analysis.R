setwd("~/Data/Nurse_Larva/")

library(plyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(edgeR)
library(topGO)
library(multcomp)
library(gtable)

load("~/Dropbox/monomorium nurses/data.processed/cleandata.RData")

theme4 = theme(axis.text=element_text(size=13),
               axis.title=element_text(size=17),
               legend.text=element_text(size=13),
               legend.title=element_text(size=17),
               plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
               axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)),
               axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))

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
        axis.title.y = element_text(margin = margin(t=0,r=10,b=0,l=0)),
        axis.title.x = element_text(margin = margin(t=10,r=0,b=0,l=0)))

apatheme_f3=theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        text=element_text(family='Arial'),
        legend.title=element_blank(),
        legend.position=c("top"),
        axis.line.x = element_line(color='black'),
        axis.line.y = element_line(color='black'),
        axis.text=element_text(size=13),
        axis.title=element_text(size=17),
        legend.text=element_text(size=13))



PS_palette = c("#deebf7","#9ecae1","#4292c6","#2171b5","#084594","darkblue","black")
tissue_palette = c("cornflowerblue","royalblue4","firebrick1","firebrick4")


#Constraint
mp <- read.csv("~/GitHub/popgenAM/results/Mphar.substitutions.csv")
mp1 = mp[mp$FN > 0 & mp$FS > 0 & mp$PS > 0,]
mC <- read.table("~/GitHub/popgenAM/results/MKtest_globalAlpha_locusF_Mphar")
mp1 = cbind(mp1,f=as.numeric(as.character((t(mC[2,4:(ncol(mC) - 1)])))))
f.est = mp1[,c(1,9)]

ext <- read.csv("~/Downloads/msx123_Supp (1)/MpharAnn.csv") #load in MBE results

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

mD <- melt(allDf,id.vars=c("Gene","tissue","code","f"))

mD$tissue_code = as.factor(apply(mD[,c(2,3)],1,paste,collapse="_"))

cT <- lapply(levels(mD$tissue_code),function(x){
  lapply(levels(mD$variable),function(y){
    
    #print(x);print(y)
    cor.test(1-mD$f[mD$tissue_code==x&mD$variable==y],mD$value[mD$tissue_code==x&mD$variable==y],method="spearman")
  })
})


allDf$reg_diff = allDf$reg_between - allDf$reg_within
allDf$targ_diff = allDf$targ_between - allDf$targ_within
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

ggsave(f2d,file="~/GitHub/MonomoriumNurseLarva/Figures/FigS4.png",height=6,width=6.5,dpi=300)


##Statistics for regulatory network topology
wilcox.test(allDf$reg_within[allDf$tissue=="larv"],allDf$reg_within[allDf$tissue=="nurse"],alternative = "greater")
wilcox.test(allDf$reg_between[allDf$tissue=="larv"],allDf$reg_between[allDf$tissue=="nurse"],alternative = "greater")

wilcox.test(allDf$reg_within[allDf$tissue=="larv" & allDf$code == "CH"],allDf$reg_between[allDf$tissue=="larv" & allDf$code == "CH"],alternative = "greater")
wilcox.test(allDf$reg_within[allDf$tissue=="larv" & allDf$code == "CG"],allDf$reg_between[allDf$tissue=="larv" & allDf$code == "CG"],alternative = "greater")
wilcox.test(allDf$reg_within[allDf$tissue=="nurse" & allDf$code == "CH"],allDf$reg_between[allDf$tissue=="nurse" & allDf$code == "CH"],alternative = "greater")
wilcox.test(allDf$reg_within[allDf$tissue=="nurse" & allDf$code == "CG"],allDf$reg_between[allDf$tissue=="nurse" & allDf$code == "CG"],alternative = "greater")

allDf$SI = allDf$reg_between-allDf$reg_within
cor.test(allDf$reg_within[allDf$tissue=="larv" & allDf$code == "CH"],allDf$reg_between[allDf$tissue=="larv" & allDf$code == "CH"],method = "spearman")
cor.test(allDf$reg_within[allDf$tissue=="larv" & allDf$code == "CG"],allDf$reg_between[allDf$tissue=="larv" & allDf$code == "CG"],method = "spearman")
cor.test(allDf$reg_within[allDf$tissue=="nurse" & allDf$code == "CH"],allDf$reg_between[allDf$tissue=="nurse" & allDf$code == "CH"],method = "spearman")
cor.test(allDf$reg_within[allDf$tissue=="nurse" & allDf$code == "CG"],allDf$reg_between[allDf$tissue=="nurse" & allDf$code == "CG"],method = "spearman")
cor.test(1-allDf$f[allDf$tissue=="nurse" & allDf$code == "CH"],allDf$reg_within[allDf$tissue=="nurse" & allDf$code == "CH"],method = "spearman")
cor.test(1-allDf$f[allDf$tissue=="nurse" & allDf$code == "CH"],allDf$reg_between[allDf$tissue=="nurse" & allDf$code == "CH"],method = "spearman")
cor.test(1-allDf$f[allDf$tissue=="nurse" & allDf$code == "CH"],allDf$SI[allDf$tissue=="nurse" & allDf$code == "CH"],method = "spearman")
cor.test(1-allDf$f[allDf$tissue=="nurse" & allDf$code == "CG"],allDf$reg_within[allDf$tissue=="nurse" & allDf$code == "CG"],method = "spearman")
cor.test(1-allDf$f[allDf$tissue=="nurse" & allDf$code == "CG"],allDf$reg_between[allDf$tissue=="nurse" & allDf$code == "CG"],method = "spearman")
cor.test(1-allDf$f[allDf$tissue=="nurse" & allDf$code == "CG"],allDf$SI[allDf$tissue=="nurse" & allDf$code == "CG"],method = "spearman")




allDN = droplevels(allDf[allDf$tissue!="larv",])
levels(allDN$tissueC) = c("nurse head ","nurse abdomen ")
allDN$SI = allDN$reg_between- allDN$reg_within
allDN$SI_rank = rank(allDN$SI)
allDN$fRank = rank(allDN$f)

nurse_pallete=tissue_palette[c(3,4)]

f3a <- ggplot(allDN,aes( x = reg_within/max(allDN$reg_within), y = 1 - f, color = tissueC))+
  geom_point(alpha = 0.5)+
  geom_smooth(method = "lm",se=FALSE)+
  xlim(0,1)+ylim(0,1)+
  apatheme+
  ylab("constraint")+  
  xlab("within-tissue connectivity")+
  theme(legend.title = element_blank(),
        axis.text = element_text(size=10),
        axis.title = element_text(size = 14),
        legend.margin = margin(c(0,0,0,0)),
        axis.line = element_line(color = "black"))+
  scale_color_manual(values = nurse_pallete)+
  annotate("text",x = .98,y=0.25,label='bold("A")',parse=TRUE,size=8,color="gray29")+
  guides(shape = guide_legend(override.aes = list(size=3)))

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
                                            plot.margin = unit(c(0.5,0.55,0.5,0.5),"cm")))
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



##########
###Phylostrata
##########
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

load("~/GitHub/devnetwork/results/collectedPhylo.RData")
newP = merge(allDN,Aps,by="Gene")
newP$psName[newP$psName=="novel"]="ant"
newP$psName[newP$psName=="aculeata"] = "hymenoptera"
newP = droplevels(newP)
levels(newP$psName)[1] = "ancient"

lm <- glm(SI ~ psName,data=newP)

drop1(lm,test="Chi")

pm = ddply(newP, ~ Gene + psName,summarize,
           f = mean(f,na.rm = TRUE),
           SI = mean(SI),
           regW = mean(reg_within),
           regB = mean(reg_between)
           )
pm = pm[!is.na(pm$f),]

fM = lapply(levels(pm$psName), function(x){
  cbind(bootCI(1-pm$f[pm$psName == x],1000),bootCI(pm$SI[pm$psName == x],1000))
})

fM = ldply(fM,data.frame)
colnames(fM)[c(1,4)] = c("f","SI")
fM$phylostrata = factor(levels(newP$psName),levels = levels(newP$psName))
levels(fM$phylostrata) = sapply(levels(fM$phylostrata),function(x) paste(x," (",sum(pm$psName==x),")",sep=""))

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
        axis.title.y = element_text(margin = margin(t = 0,r=12,b=0,l=0)),
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 20, l = 0)))+
  xlab("sociality index")+
  annotate("text",x = -0.3, y = 0.58, label = 'italic("ant (26)")',parse=TRUE,size = 4,angle = 45)+
  annotate("text",x = -0.37, y = 0.61, label = 'italic("hymenopteran (59)")',parse=TRUE,size = 4,angle = 45)+
  annotate("text",x = -0.315, y = 0.8, label = 'italic("insect (246)")',parse=TRUE,size = 4,angle = 45)+
  annotate("text",x = -0.4, y = 0.82, label = 'italic("ancient (1029)")',parse=TRUE,size = 4,angle=45)

png("~/GitHub/MonomoriumNurseLarva/Figures/Fig4.png",height = 1800,width = 2600,res =300)
grid.arrange(grid_arrange_shared_legend(f3a,f3b,f3c,ncol=1,nrow=3,position = "top"),
             f4b+theme(plot.margin = unit(c(0.5,0.5,0.5,0.9),"cm"),
                       axis.line = element_line(color = "black"))+
               annotate("text",x = -0.425, y = 0.5, label = 'bold("D")',parse = TRUE,size = 8,color="gray29"),
             nrow=1,widths = c(0.35,0.65))
grid.text("more social",x = 0.62, y = 0.04, gp = gpar(fontface="italic",fontsize=18,col="red"))
grid.lines(arrow = arrow(angle = 30, length = unit(0.2, "inches"),
                 ends = "last", type = "open"),
           x=unit(c(0.7,0.9),"npc"),y=unit(c(0.04,0.04),"npc"),gp=gpar(lwd=2.5,lty="solid",col="red"))
dev.off()

###########
##Statistics for evolutionary rates
###########
lm <- glm(SI ~ psName + f + tissueC, data = newP)
drop1(lm,.~.,test="Chi") 



cor.test(allDf$reg_within[allDf$tissue=="nurse" & allDf$code == "CH"],1-allDf$f[allDf$tissue=="nurse" & allDf$code == "CH"],method = "pearson")
cor.test(allDf$reg_within[allDf$tissue=="nurse" & allDf$code == "CG"],1-allDf$f[allDf$tissue=="nurse" & allDf$code == "CG"],method = "pearson")
cor.test(allDf$reg_between[allDf$tissue=="nurse" & allDf$code == "CH"],1-allDf$f[allDf$tissue=="nurse" & allDf$code == "CH"],method = "pearson")
cor.test(allDf$reg_between[allDf$tissue=="nurse" & allDf$code == "CG"],1-allDf$f[allDf$tissue=="nurse" & allDf$code == "CG"],method = "pearson")
cor.test(allDf$reg_diff[allDf$tissue=="nurse" & allDf$code == "CH"],1-allDf$f[allDf$tissue=="nurse" & allDf$code == "CH"],method = "pearson")
cor.test(allDf$reg_diff[allDf$tissue=="nurse" & allDf$code == "CG"],1-allDf$f[allDf$tissue=="nurse" & allDf$code == "CG"],method = "pearson")


lm <- glm(1-f ~ reg_within + reg_between, data = allDf[allDf$tissue=="nurse" & allDf$code == "CH",])
drop1(lm,.~.,test="Chi") 

lm <- glm(1-f ~ reg_within + reg_between, data = allDf[allDf$tissue=="nurse" & allDf$code == "CG",])
drop1(lm,.~.,test="Chi") 

lm <- glm(1-f ~ reg_between, data = allP[allP$code=="CH",])
drop1(lm,.~.,test="Chi") 

lm <- glm(SI ~ PS4, data = pm)
drop1(lm,.~.,test="Chi") 

pm$psCat = "old"
pm$psCat[pm$PS4=="M. pharaonis" | pm$PS4 == "hymenopteran" | pm$PS4 == "formicidae"] = "young"
wilcox.test(pm$SI[pm$psCat=="young"],pm$SI[pm$psCat=="old"],alternative="greater")

lm <- glm(1-f ~ PS4, data = pm)
drop1(lm,.~.,test="Chi") 

wilcox.test(pm$f[pm$psCat=="young"],pm$f[pm$psCat=="old"],alternative="greater")

lm <- glm(1-f ~ PS4 + SI, data = pm)
drop1(lm,.~.,test="Chi") 


########
##Top genes
########
all = merge(allDf,ext,by="Gene")
all = all[all$SwissProt!="-",]
lT = all[all$tissue=="larv",]
lTb = lT[order(lT$reg_between,decreasing=TRUE),c("Gene","tissueC","SwissProt")]
lTs = lT[order(lT$targ_between,decreasing=TRUE),c("Gene","tissueC","SwissProt")]

nT = all[all$tissue=="nurse",]
nTb = nT[order(nT$reg_between,decreasing=TRUE),c("Gene","tissueC","SwissProt")]
nTs = nT[order(nT$reg_diff,decreasing=TRUE),c("Gene","tissueC","SwissProt")]

colnames(lTb) = colnames(lTs) = colnames(nTb) = colnames(nTs) = c("Gene","Connection Type","SwissProt Annotation")

makeTbl <- function(tbl,name,font){
  
  tt3 <- ttheme_minimal(
    core=list(
      fg_params=list(fontsize=font)),
    colhead=list(fg_params=list(fontface="bold",fontsize=10)),
    rowhead=list(fg_params=list(fontsize=10)))
  
  #Make table grob out of results
  resT <- tableGrob(tbl,theme=tt3,rows = NULL)
  grid.arrange(resT)
  g <- gtable_add_grob(resT,
                       grobs = segmentsGrob( # line across the bottom
                         x0 = unit(0,"npc"),
                         y0 = unit(0,"npc"),
                         x1 = unit(1,"npc"),
                         y1 = unit(0,"npc"),
                         gp = gpar(lwd = 3.0)),
                       t = 1, b = 1, l = 1, r = 3)
  separators <- replicate(ncol(g) - 1,
                          segmentsGrob(x1 = unit(0, "npc"), gp=gpar(lty=2)),
                          simplify=FALSE)
  ## add vertical lines on the left side of columns (after 2nd)
  g <- gtable::gtable_add_grob(g, grobs = separators,
                               t = 2, b = nrow(g), l = seq_len(ncol(g)-1)+1)
  separators2 <- replicate(nrow(g) - 2,
                          segmentsGrob(y1 = unit(0, "npc"), gp=gpar(lty=1)),
                          simplify=FALSE)
  g <- gtable::gtable_add_grob(g, grobs = separators2,
                               t = seq_len(nrow(g) - 2) + 1, l = 1, r = ncol(g))
  p <- grid.arrange(g)
  ggsave(p,file=paste("~/GitHub/MonomoriumNurseLarva/Figures/",name,".png",sep=""),width=8,height=8,dpi=300)
  
}

makeTbl(nTb[1:20,],"nurseTopB",9)
makeTbl(nTs[1:20,],"nurseTopS",9)
makeTbl(lTb[1:20,],"larvTopB",9)
makeTbl(lTs[1:30,],"larvTopS",9)


########
##Secreted stuff
#######

#########
##Test if DE genes tend to secreted in Drosophila melanogaster

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

plB <- nPlot(c("LOC105830675","LOC105837907"),c("CH","W_L"),c(" giant-lens (nurse head)   "," eps8 (larva)"))

d = plB[[2]]
dM = dcast(d, expression ~ sample,value.var="expression")
cor.test(dM$expression[c(2:5)],dM$expression[c(6:9)],method="spearman")

plB <- plB[[1]]+theme(legend.position = "top")

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
  annotate("text",x=1,y=0.28,label="P = 0.011",size=4)+
  annotate("text",x=2,y=0.28,label="P = 0.094",size=4)+
  annotate("text",x=0.6,y=0.29,label="A",size=12,fontface="bold",color="gray29")

plB <- plB + theme(legend.margin = margin(t=5,l=30,b=5,r=10))+
  scale_color_manual(values=tissue_palette[c(3,1)])+
  annotate("text",x=1,y=0.933,label="B",size=12,fontface="bold",color="gray29")

pA <- arrangeGrob(p,plB,nrow=1)
ggsave(pA,file="~/GitHub/MonomoriumNurseLarva/Figures/Fig3_all.png",height=4,width=12,dpi=300)

#####GSEA
go <- read.csv("~/Writing/Data/NurseSpecialization_transcriptomicData/GOannotation.csv")
new <- list()
for (gene in unique(go$gene)){
  d = go[go$gene %in% gene,]
  new[[gene]]=as.character(d$GO)
}

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
  ggsave(p,file=paste("~/GitHub/MonomoriumNurseLarva/Figures/",name,".png",sep=""),width=10,height=4,dpi=300)
  
}

colnames(GsNH)[6] = colnames(GsNG)[6] = "P-value"
makeTbl2(GsNH,"NHGO",10)
makeTbl2(GsNG,"NGGO",10)


#spitz is LOC105833934
#EGFR substrate 8 is LOC105837907
fCH = factors[grepl("CH",rownames(factors)),]
eCH = log(fpkm["LOC105830675",colnames(fpkm) %in% rownames(fCH)])
fLW = factors[grepl("W_L",rownames(factors)),]
eLW = log(fpkm["LOC105837907",colnames(fpkm) %in% rownames(fLW)])

expr = data.frame(gL = t(eCH),stage = fCH$stage)
expr2 = data.frame(gL = t(eLW),stage = fLW$stage)
colnames(expr)[1] = colnames(expr2)[1] = "expression"
expr$sample = "giant-lens  (nurse head)"
expr2$sample = "eps8  (larva)"
e = rbind(expr,expr2)
#e$stage = as.numeric(as.character(e$stage))
levels(e$stage) = c("L1","L2","L3","L4","L5")

p <- ggplot(e,aes(x=stage,y=expression,color=sample))+
  geom_point()+
  ylim(1,3.5)+
  geom_smooth(aes(group=sample,fill=sample),color="black",se=FALSE)+
  xlab("larval developmental stage")+
  scale_color_manual(values = tissue_palette[c(1,3)])+
  scale_fill_manual(values = tissue_palette[c(1,3)])+
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line.x = element_line(color='black'),
        axis.line.y = element_line(color='black'),
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 15))+
  annotate("text",x=4.1,y=3.2,label=parse(text=paste0("italic('eps8')~(larva)")),size=6)+
  annotate("text",x=4.3,y=1.2,label=parse(text=paste0("italic('giant-lens')~(nurse)")),size=6)

ggsave(p,file="~/GitHub/MonomoriumNurseLarva/Figures/gl_egfr.png",height=6,width=6.5,dpi=600)


p <- ggplot(e[grepl("nurse",e$sample),],aes(x=stage,y=expression,color=sample))+
  geom_point()+
  ylim(1,3.5)+
  geom_smooth(aes(group=sample,fill=sample),color="black",se=FALSE)+
  xlab("larval developmental stage")+
  scale_color_manual(values = tissue_palette[c(3)])+
  scale_fill_manual(values = tissue_palette[c(3)])+
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line.x = element_line(color='black'),
        axis.line.y = element_line(color='black'),
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 15))+
  annotate("text",x=4.3,y=1.2,label=parse(text=paste0("italic('giant-lens')~(nurse)")),size=6)

ggsave(p,file="~/GitHub/MonomoriumNurseLarva/Figures/gl_egfr_justnurse.png",height=6,width=6.5,dpi=600)

##With sexual larvae
fCH = factors[grepl("CH",rownames(factors)),]
eCH = log(fpkm["LOC105830675",colnames(fpkm) %in% rownames(fCH)])
fLW = factors[grepl("1LW|LS",rownames(factors)),]
eLW = log(fpkm["LOC105837907",colnames(fpkm) %in% rownames(fLW)])
fLW2 = factors[grepl("LS",rownames(factors)),]
eLW2 = log(fpkm["LOC105837907",colnames(fpkm) %in% rownames(fLW2)])

expr = data.frame(gL = t(eCH),stage = fCH$stage)
expr2 = data.frame(gL = t(eLW),stage = fLW$stage)
expr3 = data.frame(gL = t(eLW2),stage = fLW2$stage)
colnames(expr)[1] = colnames(expr2)[1] = colnames(expr3)[1] = "expression"
expr$sample = "giant-lens  (nurse head)"
expr2$sample = "eps8  (larva)"
expr3$sample = "eps8 (repr. larva)"
e = rbind(expr,expr2)
e2 = rbind(e,expr3)
#e$stage = as.numeric(as.character(e$stage))
levels(e$stage) = c("L1","L2","L3","L4","L5")

p <- ggplot(e,aes(x=stage,y=expression,color=sample))+
  geom_point()+
  geom_smooth(aes(group=sample,fill=sample),color="black")+
  xlab("larval developmental stage")+
  scale_color_manual(values = tissue_palette[c(1,3,4)])+
  scale_fill_manual(values = tissue_palette[c(1,3,4)])+
  theme(legend.position = "none",
        panel.background = element_blank(),
        axis.line.x = element_line(color='black'),
        axis.line.y = element_line(color='black'),
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 15))+
  annotate("text",x=4.1,y=3.3,label=parse(text=paste0("italic('eps8')~('reproductive larva')")),size=6)+
  annotate("text",x=4.3,y=1.2,label=parse(text=paste0("italic('giant-lens')~(nurse)")),size=6)

ggsave(p,file="~/GitHub/MonomoriumNurseLarva/Figures/gl_egfr_sex.png",height=6,width=6.5,dpi=600)

theme_black = function(base_size = 12, base_family = "") {
  
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    
    theme(
      # Specify axis options
      axis.line = element_blank(),  
      axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.ticks = element_line(color = "white", size  =  0.2),  
      axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),  
      legend.key = element_rect(color = "white",  fill = "black"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*0.8, color = "white"),  
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA),  
      panel.border = element_rect(fill = NA, color = "white"),  
      panel.grid.major = element_line(color = "grey35"),  
      panel.grid.minor = element_line(color = "grey20"),  
      panel.margin = unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = base_size*0.8, color = "white"),  
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),  
      plot.title = element_text(size = base_size*1.2, color = "white"),  
      plot.margin = unit(rep(1, 4), "lines")
      
    )
  
}

pal = c('royalblue1','#FF3030')

p <- ggplot(e,aes(x=stage,y=expression,color=sample,fill=sample))+
  geom_point(size=2.5)+
  apatheme+
  ylim(1,3.5)+
  geom_smooth(aes(group=sample,fill=sample),color="white")+
  xlab("larval developmental stage")+
  scale_color_manual(values = pal)+
  scale_fill_manual(values = pal)+
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, color = "white"),  
        axis.ticks = element_line(color='white'),
        panel.background = element_rect(fill = "black", color  =  NA),  
        axis.title = element_text(size = 25, face = "bold",color='white'),
        axis.text = element_text(size = 20,color='white'),
        plot.background = element_rect(color = "black", fill = "black"))+
  annotate("text",x=4.5,y=3.2,label=parse(text=paste0("italic('eps8')~(larva)")),size=8,color='white')+
  annotate("text",x=4.5,y=1.2,label=parse(text=paste0("italic('giant-lens')~(nurse)")),size=8,color='white')

ggsave(p,file="~/GitHub/MonomoriumNurseLarva/Figures/F4b.png",height=6,width=7,dpi=600)


pal = c('royalblue1','#FF3030')

p <- ggplot(e[grepl("nurse",e$sample),],aes(x=stage,y=expression,color=sample,fill=sample))+
  geom_point(size=2.5)+
  apatheme+
  ylim(1,3.5)+
  geom_smooth(aes(group=sample,fill=sample),color="white")+
  xlab("larval developmental stage")+
  scale_color_manual(values = pal[2])+
  scale_fill_manual(values = pal[2])+
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA, color = "white"),  
        axis.ticks = element_line(color='white'),
        panel.background = element_rect(fill = "black", color  =  NA),  
        axis.title = element_text(size = 25, face = "bold",color='white'),
        axis.text = element_text(size = 20,color='white'),
        plot.background = element_rect(color = "black", fill = "black"))+
  annotate("text",x=4.5,y=1.2,label=parse(text=paste0("italic('giant-lens')~(nurse)")),size=8,color='white')

ggsave(p,file="~/GitHub/MonomoriumNurseLarva/Figures/F4b2.png",height=6,width=7,dpi=600)

a = merge(antSB,allDf,by="Gene")
ggplot(a,aes(x=reg_between,y=cb,color=tissueC))+
  geom_point()+
  geom_smooth()

aN = a[a$tissue=="nurse" & a$code=="CH",]
cor.test(aN$reg_between,aN$cb,method="spearman")

jelly = as.character(ann$gene[grepl("jelly",ann$DescriptionSP)])
aN = allDf[allDf$tissue=="nurse",]
aN$rRank = rank(aN$reg_between)
e = aN[aN$Gene %in% a$gene,]
e2 = allDf[allDf$Gene %in% a$gene,]

a = ann[grepl("growth factor",ann$DescriptionSP),c("gene","DescriptionSP")]
egfr_genes2 <- c("LOC105828088","LOC105830494","LOC105834310",
                "LOC105837907","LOC105838568","LOC105832621",
                "LOC105838372","LOC105829781")

gl <- "LOC105830675"
spitz <- "LOC105833934"

go <- read.table("~/GitHub/devnetwork/data/dmel_ann.txt",sep="\t",header = FALSE,stringsAsFactors = FALSE)
load("~/GitHub/devnetwork/results/collectedPhylo.RData")
g = go[go$V2=="GO:0007173",]
egfr_genes <- as.character(ENDogg$gene_Mphar[ENDogg$gene_Dmel %in% g$V1])
a = ann[ann$gene %in% egfr_genes,c("gene","DescriptionSP")]

netF = fpkm[rownames(fpkm) %in% unique(c(egfr_genes,egfr_genes2,gl,spitz,jelly)),grepl("W_L|CH",colnames(fpkm))]
netL = netF[,grepl("W_L",colnames(netF))]
netN = netF[,grepl("CH",colnames(netF))]
rownames(netL) = paste("Larv_",rownames(netL),sep="")
rownames(netN) = paste("Nurse_",rownames(netN),sep="")
genIn = rbind(netL,netN)


a = ann[ann$gene %in% rownames(netF),c("gene","DescriptionSP")]

procoll ="LOC105833984"

df <- read.csv("~/GitHub/MonomoriumNurseLarva/results/DEC16CGGenieMat.csv")
rownames(df) = df[,1]
d = as.vector(df[rownames(df)==paste("nurse",procoll,sep="_"),])
d=d[-c(1)]
d=t(d)
d = d[grepl("larv",rownames(d)),]
top = names(d)[d > quantile(d,0.9)]
topL = gsub("larv_","",top)
d = data.frame(Gene = topL,num=seq(1,100))
dA = merge(d,ann[,c("gene","DescriptionSP")],by.x="Gene",by.y='gene')
a = ENDogg[ENDogg$gene_Mphar %in% topL,]
a[a$gene_Mphar %in% unique(c(egfr_genes,egfr_genes2,gl,spitz,jelly)),]


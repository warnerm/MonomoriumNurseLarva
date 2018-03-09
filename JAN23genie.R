setwd("~/Data/Nurse_Larva/")

library(plyr)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(edgeR)
library(topGO)
library(scales)
library(ggmosaic)
library(questionr)
library(VennDiagram)
library(multcomp)
library(boot)
library(gtable)


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
        axis.title=element_text(size=17),
        legend.text=element_text(size=13))

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

f.est <- read.csv("~/Data/Nurse_Larva/MKtestConstraintOneAlpha.csv")
colnames(f.est) = c("Gene","f")
ext <- read.csv("~/Downloads/msx123_Supp (1)/MpharAnn.csv") #load in MBE results

makePlot <- function(nurse){
  df <- read.csv(paste("FEB26/JAN23",nurse,"GenieTabConn.csv",sep=""))
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

f2b <- ggplot(allDf,aes (x = tissueC, y = reg_within/max(allDf$reg_within)))+
  geom_boxplot(notch = TRUE)+
  ylab("within-tissue regulatory strength")+
  xlab("connection type")+
  apatheme+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1))+
  theme(legend.position = "none",
        plot.margin = unit(c(0.5,2.5,0.5,0),"cm"),
        axis.text.x = element_text(angle = 315, hjust = 0, vjust = 1),
        axis.line = element_line(color="black"))

f2c <- ggplot(allDf,aes (x = tissueC, y = reg_between/max(allDf$reg_within)))+
  geom_boxplot(notch = TRUE)+
  ylab("social regulatory strength")+
  xlab("connection type")+
  apatheme+
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 315, hjust = 0, vjust = 1),
        plot.margin = unit(c(0.5,2.5,0.5,0),"cm"),
        axis.line = element_line(color="black"))

f2d <- ggplot(allDf, aes(x = reg_within/max(allDf$reg_within), y = reg_between/max(allDf$reg_within), color = tissueC))+
  geom_point(alpha = 0.2)+
  geom_smooth(method = "lm",se = FALSE,size=1)+apatheme+
  ylab("social regulatory strength")+
  xlab("within-tissue regulatory strength")+
  scale_color_manual(values = tissue_palette,name = bquote(underline("connection type")))

png("~/Writing/Figures/NurseLarva/corApproach/fig3.png",height = 2000,width = 2000,res =300)
f2d+theme(legend.position = c(0.8,0.9),
          axis.line.x = element_line(color='black'),
          axis.line.y = element_line(color='black'),
          legend.title = element_text(size = 17),
          legend.background = element_blank())
dev.off()

png("~/Writing/Figures/NurseLarva/corApproach/figS2.png",height = 2000,width = 3000,res =300)
grid.arrange(f2b+ylim(0,1)+annotate("text",x=3.5,y=0.95,label='bold("a")',parse=TRUE,size=9,color="gray29"),
             f2c+ylim(0,0.4)+annotate("text",x=3.5,y=0.95*0.4,label='bold("b")',parse=TRUE,size=9,color="gray29"),ncol = 2)
dev.off()

##Statistics for regulatory network topology
wilcox.test(allDf$reg_within[allDf$tissue=="larv"],allDf$reg_within[allDf$tissue=="nurse"],alternative = "greater")
wilcox.test(allDf$reg_between[allDf$tissue=="larv"],allDf$reg_between[allDf$tissue=="nurse"],alternative = "greater")

wilcox.test(allDf$reg_within[allDf$tissue=="larv" & allDf$code == "CH"],allDf$reg_between[allDf$tissue=="larv" & allDf$code == "CH"],alternative = "greater")
wilcox.test(allDf$reg_within[allDf$tissue=="larv" & allDf$code == "CG"],allDf$reg_between[allDf$tissue=="larv" & allDf$code == "CG"],alternative = "greater")
wilcox.test(allDf$reg_within[allDf$tissue=="nurse" & allDf$code == "CH"],allDf$reg_between[allDf$tissue=="nurse" & allDf$code == "CH"],alternative = "greater")
wilcox.test(allDf$reg_within[allDf$tissue=="nurse" & allDf$code == "CG"],allDf$reg_between[allDf$tissue=="nurse" & allDf$code == "CG"],alternative = "greater")

cor.test(allDf$reg_within[allDf$tissue=="larv" & allDf$code == "CH"],allDf$reg_between[allDf$tissue=="larv" & allDf$code == "CH"],method = "pearson")
cor.test(allDf$reg_within[allDf$tissue=="larv" & allDf$code == "CG"],allDf$reg_between[allDf$tissue=="larv" & allDf$code == "CG"],method = "pearson")
cor.test(allDf$reg_within[allDf$tissue=="nurse" & allDf$code == "CH"],allDf$reg_between[allDf$tissue=="nurse" & allDf$code == "CH"],method = "pearson")
cor.test(allDf$reg_within[allDf$tissue=="nurse" & allDf$code == "CG"],allDf$reg_between[allDf$tissue=="nurse" & allDf$code == "CG"],method = "pearson")


allDN = droplevels(allDf[allDf$tissue!="larv",])
levels(allDN$tissueC) = c("nurse head","nurse abdomen")
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
  xlab("within-tissue regulatory strength")+
  theme(legend.title = element_blank(),
        axis.text = element_text(size=10),
        axis.title = element_text(size = 14),
        legend.margin = margin(c(0,0,0,0)))+
  scale_color_manual(values = nurse_pallete)+
  annotate("text",x = 1,y=0.25,label='bold("a")',parse=TRUE,size=8,color="gray29")+
  guides(shape = guide_legend(override.aes = list(size=3)))

f3b <- ggplot(allDN,aes( x = reg_between/max(allDN$reg_between), y = 1 - f, color = tissueC))+
  geom_point(alpha = 0.5)+
  geom_smooth(method = "lm",se=FALSE)+
  xlim(0,1)+ylim(0,1)+
  xlab("social regulatory strength")+
  apatheme+
  ylab("constraint")+  
  theme(legend.title = element_blank(),
        axis.text = element_text(size=10),
        axis.title = element_text(size = 14))+
  scale_color_manual(values = nurse_pallete)+
  annotate("text",x = 1,y=0.25,label='bold("b")',parse=TRUE,size=8,color="gray29")

f3c <- ggplot(allDN,aes( x = SI, y = 1-f, color = tissueC))+
  geom_point(alpha = 0.5)+
  apatheme+
  ylab("constraint")+
  geom_smooth(method = "lm",se=FALSE)+
  xlab("sociality index")+
  scale_x_continuous(limits = c(-1.25,0.25))+
  theme(legend.title = element_blank(),
        axis.text = element_text(size=10),
        axis.title = element_text(size = 14))+  
  scale_color_manual(values = nurse_pallete)+
  annotate("text",x = 0.25,y=0.25,label='bold("c")',parse=TRUE,size=8,color="gray29")

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
                                            plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")))
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

all <- merge(allDN,ext,by="Gene",all.x=TRUE)
allP = all[!is.na(all$Raw.PS),]
allP$PS2 = factor(allP$PS2,levels = c("cellular",'eukaryote','bilaterian','insect','hymenopteran_ant','M.pharaonis'))
allP$PS2[allP$Raw.PS=="Neoptera"]="insect"
allP$PS3 = as.character(allP$PS2)
allP$PS3[allP$Raw.PS=="Monomorium pharaonis"]="M. pharaonis"
allP$PS4 = allP$PS3
allP$PS4[allP$Raw.PS=="Myrmicinae"|allP$Raw.PS=="Formicidae"]="formicidae"
allP$PS3 = factor(allP$PS3,levels = c("cellular",'eukaryote','bilaterian','insect','hymenopteran_ant','M. pharaonis'))
allP$PS4 = factor(allP$PS4,levels = c("cellular",'eukaryote','bilaterian','insect','hymenopteran_ant','formicidae','M. pharaonis'))
levels(allP$PS4)[5]="hymenopteran"


grid.arrange(grid_arrange_shared_legend(f3a,f3b,f3c,ncol=3,nrow=1,position = "top"),
             left=textGrob("constraint",gp=gpar(fontsize=20,fontface="bold"),rot =90),
             bottom = textGrob("regulatory strength\n",gp=gpar(fontsize=20,fontface="bold")))

pm = ddply(allP, ~ Gene + PS4,summarize,
           f = mean(f,na.rm = TRUE),
           SI = mean(SI)
           )

fM = lapply(levels(pm$PS4), function(x){
  cbind(bootCI(1-pm$f[pm$PS4 == x],1000),bootCI(pm$SI[pm$PS4 == x],1000))
})

fM = ldply(fM,data.frame)
colnames(fM)[c(1,4)] = c("f","SI")
fM$phylostrata = factor(levels(allP$PS4),levels = levels(allP$PS4))
levels(fM$phylostrata) = sapply(levels(fM$phylostrata),function(x) paste(x," (",sum(pm$PS4==x),")",sep=""))

f4b <- ggplot(fM,aes(x=SI,y=f))+
  geom_point(aes(fill=phylostrata),color="black",size = 6,pch=21)+
  geom_errorbar(aes(ymin = c1, ymax = c2,color=phylostrata),size=1)+
  geom_errorbarh(aes(xmin = c1.1, xmax = c2.1,color=phylostrata),size=1)+
  scale_fill_manual(values = PS_palette)+
  scale_color_manual(values = PS_palette)+
  apatheme+
  ylab("constraint")+
  theme(legend.position = "none",
        axis.line.x = element_line(color='black'),
        axis.line.y = element_line(color='black'),
        axis.title.y = element_text(margin = margin(t = 0,r=12,b=0,l=0)),
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 20, l = 0)))+
  xlab("sociality index")+
  ylim(0.675,0.93)+
  annotate("text",x = -0.25, y = 0.835, label = 'italic("M. pharaonis (18)")',parse=TRUE,size = 3.5,angle = 45)+
  annotate("text",x = -0.328, y = 0.78, label = 'italic("ant (27)")',parse=TRUE,size = 3.5,angle = 45)+
  annotate("text",x = -0.375, y = 0.825, label = 'italic("hymenopteran (96)")',parse=TRUE,size = 3.5,angle = 45)+
  annotate("text",x = -0.391, y = 0.854, label = 'italic("insect (132)")',parse=TRUE,size = 3.5,angle = 45)+
  annotate("text",x = -0.36, y = 0.915, label = 'italic("cellular (680)")',parse=TRUE,size = 3.5,angle=45)+
  annotate("text",x = -0.395, y = 0.92, label = 'italic("bilaterian (301)")',parse=TRUE,size = 3.5,angle = 45)+
  annotate("text",x = -0.44, y = 0.845, label = 'italic("eukaryote (395)")',parse=TRUE,size = 3.5,angle = 45)
  

png("~/Writing/Figures/NurseLarva/corApproach/fig4.png",height = 2000,width = 3000,res =300)
grid.arrange(grid_arrange_shared_legend(f3a,f3b,f3c,ncol=1,nrow=3,position = "top"),
             f4b+theme(plot.margin = unit(c(0.5,0.5,0.5,1.25),"cm"),
                       axis.line = element_line(color = "black"))+
               annotate("text",x = -0.2, y = 0.73, label = 'bold("d")',parse = TRUE,size = 8,color="gray29"),
             nrow=1,widths = c(0.35,0.65))
grid.text("more social",x = 0.62, y = 0.04, gp = gpar(fontface="italic",fontsize=18,col="red"))
grid.lines(arrow = arrow(angle = 30, length = unit(0.2, "inches"),
                 ends = "last", type = "open"),
           x=unit(c(0.7,0.9),"npc"),y=unit(c(0.04,0.04),"npc"),gp=gpar(lwd=2.5,lty="solid",col="red"))
dev.off()

###########
##Statistics for evolutionary rates
###########
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

lm <- glm(BSnIPRE.est ~ reg_within + reg_between, data = allP[allP$code=="CH",])
drop1(lm,.~.,test="Chi") 

lm <- glm(BSnIPRE.est ~ reg_within + reg_between, data = allP[allP$code=="CG",])
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

pirate <- function(ps,column,tissueSub = NULL){
  if (!is.null(tissueSub)){
    allP = allP[allP$code==tissueSub,]
  }
  d = droplevels(allP[,c(ps,column,"code")])
  colnames(d)=c("phylostrata","stat","tissue")
  d$tissue = as.factor(d$tissue)
  levels(d$tissue) = c("abdomen","head")
  #Bootstrap for confidence intervals
  sum = ldply(lapply(levels(d$phylostrata),function(x) bootCI(d$stat[d$phylostrata==x],1000)))
  sum$ps = levels(d$phylostrata)
  p1 <- ggplot(data = sum, aes(x = ps, y = mean))+
    geom_violin(data = d,aes(x = phylostrata, y = stat))+
    geom_jitter(data = d, aes(x = phylostrata, y = stat,color=tissue), shape = 1, width = .1)+
    geom_point(size = 3)+
    geom_errorbar(aes(ymax = c2, ymin = c1))+
    apatheme+
    scale_color_manual(values = c("firebrick","dodgerblue1"))+
    ggtitle(column)+
    xlab("")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.title = element_blank(),
          legend.text = element_text(size=13))+
    guides(color = guide_legend(override.aes = list(size=3)))
  return(p1)
}

p1 <- pirate("PS4","reg_within")+ggtitle("")+ylab("within-tissue regulatory strength")+
  annotate("text",x = 6.5, y = 1.25, label = 'bold("a")',parse = TRUE,color="gray29",size=9)
p2 <- pirate("PS4","reg_between")+ggtitle("")+ylab("social regulatory strength")+
  annotate("text",x = 6.5, y = 0.45, label = 'bold("b")',parse = TRUE,color="gray29",size=9)
p3 <- pirate("PS4","reg_diff")+ggtitle("")+ylab("sociality index")+
  annotate("text",x = 6.5, y = 0.2, label = 'bold("c")',parse = TRUE,color="gray29",size=9)

png("~/Writing/Figures/NurseLarva/corApproach/pirates.png",height = 4000,width = 1500,res=300)
grid.arrange(grid_arrange_shared_legend(p1,p2,p3,ncol=1,nrow=3,position = "top"),
             bottom = textGrob("phylostrata\n",gp=gpar(fontsize=18,fontface="bold")))
dev.off()

########
##Top genes
########
all = merge(allDf,ext,by="Gene")
all = all[all$SwissProt!="-",]
lT = all[all$tissue=="larv",]
lTb = lT[order(lT$reg_between,decreasing=TRUE),c("Gene","tissueC","SwissProt")]
lTs = lT[order(lT$reg_diff,decreasing=TRUE),c("Gene","tissueC","SwissProt")]

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
  ggsave(p,file=paste("~/Writing/Figures/NurseLarva/corApproach/",name,".png",sep=""),width=8,height=8,dpi=300)
  
}

makeTbl(nTb[1:20,],"nurseTopB",9)
makeTbl(nTs[1:20,],"nurseTopS",9)
makeTbl(lTb[1:20,],"larvTopB",9)
makeTbl(lTs[1:20,],"larvTopS",9)


#####GSEA
genieAll <- lapply(plotsG,function(x){
  list(x[x$tissue=="nurse",],
       x[x$tissue=="larv",])
})
gAll = unlist(genieAll,recursive=FALSE)
gAll = gAll[c(1,3,2,4)]

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
  allRes <- GenTable(GOdata,KS=resultKS)
  return(allRes[allRes$KS < 0.05,]) 
}

#Get top result based on p value
GOtop <- function(res){
  res = res[res$Significant > res$Expected,]
  return(res[1,2])
}

#Overall GSEA, by reg_between, reg_within, targ_between, and targ_within
#Return matrix with each of the four stats, each tissue
gsea_overall <- function(gS){
  gS[,8] = gS[,3] - gS[,2]
  gseaInput <- lapply(c(2,3,8),function(x) gS[,x])
  for (i in 1:length(gseaInput)){
    names(gseaInput[[i]]) = gS$Gene
  }
  gseaRes <- lapply(gseaInput, GSEAfunc)
  topRes <- ldply(lapply(gseaRes,GOtop))
  return(list(topRes,gseaRes))
}

gseaTissue <- lapply(gAll,gsea_overall)

allGS <- do.call(cbind,lapply(gseaTissue,function(x) x[[2]][[3]]$Term[1:5]))
colnames(allGS) = c("nurse head \u2192 larva","nurse abdomen \u2192 larva","larva \u2192 nurse head","larva \u2192 nurse abdomen")
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
                       t = 1, b = 1, l = 1, r = 4)
  separators <- replicate(ncol(g) - 1,
                          segmentsGrob(x1 = unit(0, "npc"), gp=gpar(lty=2)),
                          simplify=FALSE)
  ## add vertical lines on the left side of columns (after 2nd)
  g <- gtable::gtable_add_grob(g, grobs = separators,
                               t = 2, b = nrow(g), l = seq_len(ncol(g)-1)+1)
  p <- grid.arrange(g)
  ggsave(p,file=paste("~/Writing/Figures/NurseLarva/corApproach/",name,".png",sep=""),width=10,height=4,dpi=300)
  
}

makeTbl2(allGS,"GOtable",8)


load("~/Dropbox/monomorium nurses/data.processed/cleandata.RData")
tissueExpr <- function(code){
  f = fpkm[,grepl(code,colnames(fpkm))]
  meanF = rowSums(f)/ncol(f)
  geomMeanF = apply(f,1,function(x) exp(mean(log(x))))
  d = data.frame(Gene = rownames(f),mean = meanF,geom_mean = geomMeanF)
  d$code = code
  return(d)
}

headExpr <- tissueExpr("CH")
abdExpr <- tissueExpr("CG")
QH <- tissueExpr("QH")
QG <- tissueExpr("QG")
QH$code = "CH"
QG$code = "CG"

expr <- rbind(QH,QG)


allPe = merge(allP,expr,by = c("Gene","code"))

allPeH = allPe[allPe$code=="CH",]
allPeG = allPe[allPe$code=="CG",]
cor.test(allPeH$reg_within,allPeH$geom_mean,method = "pearson")
cor.test(allPeH$reg_between,allPeH$geom_mean,method = "pearson")
cor.test(allPeH$reg_diff,allPeH$geom_mean,method = "pearson")
cor.test(allPeG$reg_within,allPeG$geom_mean,method = "pearson")
cor.test(allPeG$reg_between,allPeG$geom_mean,method = "pearson")
cor.test(allPeG$reg_diff,allPeG$geom_mean,method = "pearson")

lm <- glm(1-log(f) ~ reg_between + reg_within + PS2 + geom_mean,data = allPeH)
av <- aov(lm)
summary(av)

lm <- glm(BSnIPRE.est ~ reg_within + reg_between + geom_mean + PS2, data = allPeG)
drop1(lm,.~.,test="Chi") 

ggplot(allPe,aes (x = code, y = mean+1, fill = PS2))+
  geom_boxplot(notch = TRUE)+scale_y_log10()

png("~/Writing/Figures/NurseLarva/corApproach/exprRegulation.png",height = 4000,width = 4000,res =300)
ggplot(allPe[allPe$code=="CH",],aes (x = geom_mean+1, y = reg_within, color = PS2))+
  geom_point()+scale_x_log10()+geom_smooth(method = "lm",se=FALSE,size=3)+
  ylab("within-tissue reg strength")+
  xlab("geometric mean expression")+theme_bw()+theme_all+
  scale_color_manual(values = PS_palette)
dev.off()

ggplot(allPe[allPe$code=="CG",],aes (x = geom_mean+1, y = reg_within, color = PS2))+
  geom_point()+scale_x_log10()+geom_smooth(method = "lm",se=FALSE)


ggplot(allPe[allPe$code=="CH",],aes (x = geom_mean+1, y = log(f), color = PS2))+
  geom_point()+scale_x_log10()+geom_smooth(method = "lm",se=FALSE)

lm <- glm(reg_between ~ geom_mean*PS2,data = allPe[allPe$code=="CG",])
av <- aov(lm)
summary(av)
TukeyHSD(av,"PS2")




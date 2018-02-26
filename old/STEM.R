#First, use STEMalg.R to sort genes into profiles.

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


makeTbl <- function(tbl,name,font){
  
  tt3 <- ttheme_minimal(
    core=list(
      fg_params=list(fontsize=font)),
    colhead=list(fg_params=list(fontface="bold",fontsize=10)),
    rowhead=list(fg_params=list(fontsize=10)))
  
  #Make table grob out of results
  resT <- tableGrob(tbl,theme=tt3)
  grid.arrange(resT)
  g <- gtable_add_grob(resT,
                       grobs = segmentsGrob( # line across the bottom
                         x0 = unit(0,"npc"),
                         y0 = unit(0,"npc"),
                         x1 = unit(1,"npc"),
                         y1 = unit(0,"npc"),
                         gp = gpar(lwd = 3.0)),
                       t = 1, b = 1, l = 1, r = 5)
  separators <- replicate(ncol(g) - 1,
                          segmentsGrob(x1 = unit(0, "npc"), gp=gpar(lty=2)),
                          simplify=FALSE)
  ## add vertical lines on the left side of columns (after 2nd)
  g <- gtable::gtable_add_grob(g, grobs = separators,
                               t = 2, b = nrow(g), l = seq_len(ncol(g)-1)+1)
  p <- grid.arrange(g)
  ggsave(p,file=paste("~/Writing/Figures/NurseLarva/corApproach/",name,".png",sep=""),width=10,height=4,dpi=300)
  
}

#After transferring files:
files <- dir("~/Data/Nurse_Larva/FEB26/",pattern=".*membership.RData")

#trim file names to get the sample information
names <- lapply(1:length(files),function(x) gsub(" STEMmembership.RData","",files[x]))
STEMdata <- list()

ext <- read.csv("~/Writing/Data/NurseSpecialization_transcriptomicData/MpharAnn.csv") #load in MBE results
load("~/Dropbox/monomorium nurses/data.processed/ps_genelevelJuly29.RData")
a = TAIgene$Mphar_E5 
ext <- merge(ext,a,by.x="Gene",by.y="gene")


#Filter for genes with phylostrata calls (they are long enough)
keep = ext$Gene[!is.na(ext$ps)]

#load all in, store, in STEM data
for (i in 1:length(files)){
  load(paste("~/Data/Nurse_Larva/FEB26/",files[i],sep=""))
  STEMdata[[i]] = results
}
names(STEMdata) = names


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

cbbPalette <- c("#000000", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", 
                "#CC79A7", "#F0E442","#999999","purple2","springgreen3")
pal1 = c("lightgoldenrod1","cornflowerblue","brown1","tan1","black")
pal2 = cbbPalette[c(8,4,2,5,11)]
pal3 = c(pal1,pal2)

SetProfiles <- function(timepoints){
  profiles <- matrix(ncol=timepoints)
  profiles = as.data.frame(profiles)
  colnames(profiles) = paste("Stage",seq(1,timepoints,1))
  posMoves = c(-1,0,1)
  for (i2 in -1:1){
    for (i3 in -2:2){
      for (i4 in -3:3){
        if (timepoints == 5){
          for (i5 in -4:4){
            if ((i2 - i3) %in% posMoves & (i3 - i4) %in% posMoves & (i4 - i5) %in% posMoves){
              profiles = rbind(profiles,c(0,i2,i3,i4,i5))
            }
          }
        } else if (timepoints == 4){
          if ((i2 - i3) %in% posMoves & (i3 - i4) %in% posMoves){
            profiles = rbind(profiles,c(0,i2,i3,i4))
          }
        }
        
      }
    }
  }
  profiles = profiles[-c(1),]
  profiles = profiles[rowSums(abs(profiles))!=0,]
  rownames(profiles) = seq(1,nrow(profiles),by=1)
  return(profiles)
}

profiles5 <- SetProfiles(5)

CountsbyStage <- function(code){
  load("~/Dropbox/monomorium nurses/data.processed/cleandata.RData")
  counts <- counts[rownames(counts) %in% keep,grep(code,colnames(counts))]
  factors <- factors[grep(code,rownames(factors)),]
  StageExpr <- matrix(nrow=nrow(counts),ncol=timepoints)
  stages = seq(6-timepoints,5,1) ##for sexuals, start at 2
  for (j in stages){
    st = j-5+timepoints ##Converts actual stage to relative (so sexuals go 1:4 instead of 2:5)
    for (i in 1:nrow(counts)){
      StageExpr[i,st] = mean(as.numeric(counts[i,factors$stage==j])) ##Calculate initial expression
    }
  }
  for (j in 2:timepoints){
    for (row in 1:nrow(StageExpr)){
      StageExpr[row,j]=log2((StageExpr[row,j]+1)/(StageExpr[row,1]+1)) ##Expression at stage X is the log2 ratio of expression at stage X to expression at intitial stage
    }
  }
  StageExpr[,1]=0
  rownames(StageExpr) = rownames(counts)
  return(StageExpr)
}

###As input to the algorithm, we need the mean expression across samples of a given type at each developmental stage
timepoints = 5


#Plot the actual log2 expression values of genes over time of the largest 5 larval and nurse modules
sharedDf <- function(nurse,larv,lSTEM,nSTEM){
  nE = as.data.frame(CountsbyStage(nurse))
  lE = as.data.frame(CountsbyStage(larv))
  nE$Gene = rownames(nE)
  lE$Gene = rownames(lE)
  nE$Module = STEMdata[[nSTEM]]$GeneMembership
  lE$Module = STEMdata[[lSTEM]]$GeneMembership
  nEm = melt(nE,id.vars = c("Gene","Module"))
  lEm = melt(lE,id.vars = c("Gene","Module"))
  nEm$tissue = "nurse"
  lEm$tissue = "larva"
  all = rbind(nEm,lEm)
  all$variable=as.numeric(all$variable)
  all$modN = 81 - all$Module
  Mods = as.data.frame(table(all$Module,all$tissue))
  colnames(Mods) = c("Module","tissue","nGene")
  nMod = as.data.frame(table(all$modN,all$tissue))
  colnames(nMod) = c("Module","tissue","nGene")
  both = merge(Mods,nMod,by="Module")
  posShared = unique(nE$Module)[unique(nE$Module) %in% unique(lE$Module)]
  negShared = both$Module[both$tissue.x=="larva" & both$tissue.y=="nurse" & both$nGene.x!=0 & both$nGene.y!=0]
  return(list(all,posShared,negShared))
}

Hshare <- sharedDf("CH","W_L","WLarv","NurseH")
RHshare <- sharedDf("RH","QW","WlarvQR","RNurseH")
Gshare <- sharedDf("CG","W_L","WLarv","NurseG")
RGshare <- sharedDf("RG","QW","WlarvQR","RNurseG")
res = list(Hshare,RHshare,Gshare,RGshare)

#Table for numbers of genes and number of modules
nMod <- lapply(c("NurseH","RNurseH","NurseG","RNurseG"),function(x) length(unique(STEMdata[[x]][[1]])))
nPos <- lapply(res,function(x) length(x[[2]][!is.na(x[[2]])]))
nNeg <- lapply(res,function(x) length(x[[3]][!is.na(x[[3]])]))
nGene <- lapply(res,function(x){
  sum(x[[1]]$tissue== "nurse" & x[[1]]$Module %in% c(x[[2]],81-as.integer(as.character(x[[3]]))) & !is.na(x[[1]]$Module))/5
})

tbl = do.call(cbind,lapply(list(nMod,nPos,nNeg,nGene),unlist))
colnames(tbl) = c("total modules","modules positively\nshared with larvae","modules negatively\nshared with larvae","number of genes\nin shared modules")
rownames(tbl) = c("focal nurse head","random nurse head","focal nurse abdomen","random nurse abdomen")

makeTbl(tbl,"STEMtable",8)

#Fisher's exact test
lNum = length(unique(STEMdata[["WLarv"]][[1]])) - 1
geneH = rbind(c(tbl[1,4],tbl[2,4]),c(9344-tbl[1,4],9344-tbl[2,4]))
geneG = rbind(c(tbl[3,4],tbl[4,4]),c(9344-tbl[3,4],9344-tbl[4,4]))

#Number of larva genes
nGene <- lapply(res,function(x){
  sum(x[[1]]$tissue== "larva" & x[[1]]$Module %in% c(x[[2]],as.integer(as.character(x[[3]]))) & !is.na(x[[1]]$Module))/5
})


#All negative
Hmods = as.integer(as.character(Hshare[[3]]))
Hmod_inv = 81 - Hmods
nGene = unlist(lapply(Hmod_inv,function(x){
  sum(Hshare[[1]]$Module==x & Hshare[[1]]$tissue =="nurse",na.rm = TRUE)
}))

lMod = Hmods[order(nGene,decreasing = TRUE)][1:5]
Hmod = Hmod_inv[order(nGene,decreasing = TRUE)][1:5]


#Many positive for abdomens
nGene = unlist(lapply(Gshare[[2]],function(x){
  sum(Gshare[[1]]$Module==x & Gshare[[1]]$tissue =="nurse",na.rm = TRUE)
}))

Gmod = Gshare[[2]][order(nGene,decreasing = TRUE)][1:5]

sharedPlot <- function(mod,tissue,df,pal,flip = FALSE){
  dN = df[(df$Module %in% mod & df$tissue == tissue),]
  dN$Module = factor(dN$Module,levels = mod)
  nTop = profiles5[mod,]
  nTop$profile = levels(dN$Module)
  nTopM = melt(nTop,id.vars = "profile")
  nTopM$variable = as.factor(nTopM$variable)
  nTopM$profile = factor(nTopM$profile,levels = levels(dN$Module))
  nMed = ddply(dN,~ Module + variable,summarize,
               med = median(value))
  levels(dN$Module) = levels(nMed$Module) = sapply(levels(dN$Module),function(x) sum(dN$Module==x)/5)
  if (flip){
    dN$value = - dN$value
    nMed$med = -nMed$med
    nTopM$value = - nTopM$value
  }
  dN$variable=as.factor(dN$variable)
  levels(dN$variable) = c("L1","L2","L3","L4","L5")
  p1 <- ggplot(dN,aes(x = variable, y = value))+
    geom_line(alpha = 0.1,aes(color = Module,group = Gene))+
    apatheme+
    xlab("larval developmental stage")+
    ylab("log2 expression change from initial stage")+
    theme(legend.position = c(0.85,0.83),
          axis.line.x = element_line(color="black"),
          axis.line.y = element_line("black"),
          legend.background = element_blank(),
          legend.title = element_text(size=15),
          legend.title.align = 0.5,
          legend.justification = "center",
          legend.text = element_text(size = 13))+
    scale_color_manual(name = bquote(underline("genes in top 5\nshared profiles")),values=pal)
  p1 <- p1 + geom_line(data = nMed,aes(x = variable,y=med,group=Module),alpha = 1,size = 1.5,color="black")+
    geom_line(data = nMed,aes(x = variable,y=med,group=Module,color = Module),alpha = 1,size = 1.2)
  levels(nTopM$variable) = c("L1","L2","L3","L4","L5")
  p1o <- ggplot(nTopM)+
    geom_line(size = 1.5,data = nTopM, alpha = 0.8,
              aes(x = variable, y = value,color = profile,group=profile))+
    apatheme+
    xlab("")+
    ylab("")+
    scale_color_manual(values=pal)+
    theme(legend.position = "none",
          axis.line.x = element_line(color="black"),
          axis.line.y = element_line("black"),
          plot.background = element_blank())
  return(list(p1,p1o))
}

pieCharts <- function(data,tissue,flip=FALSE){
  mods = c(data[[2]],as.integer(as.character(data[[3]])))
  mods = mods[!is.na(mods)]
  if (flip) mods = 81 - mods
  totShared = sum(data[[1]]$Module %in% mods & data[[1]]$tissue == tissue)/5
  d = data.frame(Prop = c(totShared/9344,(9344-totShared)/9344),Shared = c('a','b'),lev = c('same','same'))
  p <- ggplot(d,aes(x = lev,y = Prop, fill = Shared))+
    geom_bar(stat="identity",width =1)+
    coord_polar(theta = "y")+
    theme_bw()+
    scale_fill_manual(values = c("black","lightgrey"))+
    theme(legend.position = "none",
          axis.line = element_blank(),
          plot.background = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank())
  return(p)
}

Gplot <- sharedPlot(Gmod,"nurse",Gshare[[1]],pal1)
GLplot <- sharedPlot(Gmod,"larva",Gshare[[1]],pal1)
Hplot <- sharedPlot(Hmod,"nurse",Hshare[[1]],pal2,flip = TRUE)
HLplot <- sharedPlot(lMod,"larva",Hshare[[1]],pal2)

pieG <- pieCharts(Gshare,"nurse")
pieGL <- pieCharts(Gshare,"larva")
pieH <- pieCharts(Hshare,"nurse",flip=TRUE)
pieHL <- pieCharts(Hshare,"larva")

vp <- viewport(width = 0.3, height = 0.3, x = 0.25, y = 0.25)

library(cowplot)
library(rsvg)
image <- rsvg("~/Downloads/smaller larva with shadow.svg")
mult=dim(image)[1]/dim(image)[2]


png("~/Writing/Figures/NurseLarva/corApproach/allTogether.png",height = 4000, width = 6000, res = 300)
ggdraw()+
  draw_plot(Hplot[[1]]+scale_y_continuous(limits = c(-1.5,1.5))+
              theme(legend.position = c(0.7,0.2),
                    axis.line = element_line(color = "black"))+
              annotate("text",x = 1.2,y=1.2,label="b",size=14,color="gray29")+
              ylab("-log2 expression change from initial stage"), x = 0.44, y = 0.55, width = 0.4, height = 0.45)+
  draw_plot(Hplot[[2]]+scale_y_continuous(limits = c(-2,2)), x = 0.48, y = 0.59, width = 0.15, height = 0.175)+
  draw_plot(pieH,x = 0.705, y = 0.61,width=0.1,height=0.1)+
  draw_plot(HLplot[[1]]+annotate("text",x = 1.2,y=1.6,label="a",size=14,color="gray29")+
              scale_y_continuous(limits = c(-2,2))+
              theme(legend.position = c(0.7,0.2),
                    axis.line = element_line(color = "black")),
            x = 0, y = 0.5, width = 0.4, height = 0.45)+
  draw_plot(pieHL,x = 0.26, y = 0.56,width=0.1,height=0.1)+
  draw_plot(HLplot[[2]]+scale_y_continuous(limits = c(-2,2)), x = 0.04, y = 0.54, width = 0.15, height = 0.175)+
  draw_plot(Gplot[[1]]+annotate("text",x = 1.2,y=.4,label="d",size=14,color="gray29")+
              scale_y_continuous(limits = c(-1.5,0.5))+
              theme(legend.position = c(0.85,0.85),
                    axis.line = element_line(color = "black")),
            x = 0.55, y = 0.04, width = 0.4, height = 0.45)+
  draw_plot(pieG,x = 0.87, y = 0.36,width=0.1,height=0.1)+
  draw_plot(Gplot[[2]],  x = 0.59, y = 0.08, width = 0.15, height = 0.175)+
  draw_plot(GLplot[[1]]+
              annotate("text",x = 1.2,y=.8,label="c",size=14,color="gray29")+
              scale_y_continuous(limits = c(-2.5,1))+
              theme(legend.position = c(0.8,0.78),
                    axis.line = element_line(color = "black")),
            x = 0.08, y = 0, width = 0.4, height = 0.45)+
  draw_plot(pieGL,x = 0.375, y = 0.293,width=0.1,height=0.1)+
  draw_plot(Gplot[[2]], x = 0.12, y = 0.04, width = 0.15, height = 0.175)+
  draw_plot(rasterGrob(image,x = 0, y = 0),x=.45,y=0.56,width=0.3/mult,height=0.3)
grid.lines(x=unit(c(0.2,0.36),"npc"),y=unit(c(0.5,0.49),"npc"),gp=gpar(lwd=2,lty="dashed"))
grid.lines(x=unit(c(0.40,0.42),"npc"),y=unit(c(0.85,0.5),"npc"),gp=gpar(lwd=2,lty="dashed"))
grid.lines(x=unit(c(0.36,0.25),"npc"),y=unit(c(0.44,0.41),"npc"),gp=gpar(lwd=2,lty="dashed"))
grid.lines(x=unit(c(0.45,0.48),"npc"),y=unit(c(0.42,0.22),"npc"),gp=gpar(lwd=2,lty="dashed"))
grid.lines(x=unit(c(0.53,0.55),"npc"),y=unit(c(0.49,0.3),"npc"),gp=gpar(lwd=2,lty="dashed"))
grid.lines(x=unit(c(0.56,0.65),"npc"),y=unit(c(0.53,0.49),"npc"),gp=gpar(lwd=2,lty="dashed"))
grid.lines(x=unit(c(0.455,0.46),"npc"),y=unit(c(0.51,0.62),"npc"),gp=gpar(lwd=2,lty="dashed"))
grid.lines(x=unit(c(0.47,0.5),"npc"),y=unit(c(0.515,0.57),"npc"),gp=gpar(lwd=2,lty="dashed"))
dev.off()







cbbPalette <- c("#000000", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", 
                "#CC79A7", "#F0E442","#999999","purple2","springgreen3")
pal1 = c("lightgoldenrod1","cornflowerblue","brown1","tan1","black")
pal2 = cbbPalette[c(8,4,2,5,11)]
pal3 = c(pal1,pal2)

png("~/Downloads/test.png")
plot(x = 1:length(pal1), y = rep(1, length(pal1)), pch = 19, cex = 5, 
     col = pal1, yaxt = "n", bty = "n", xaxt = "n", xlab = "", ylab = "")
dev.off()

plot(x = 1:length(cbbPalette), y = rep(1, length(cbbPalette)), pch = 19, cex = 5, 
     col = cbbPalette, yaxt = "n", bty = "n", xaxt = "n", xlab = "", ylab = "")

png("~/Downloads/test.png")
plot(x = 1:length(pal3), y = rep(1, length(pal3)), pch = 19, cex = 5, 
     col = pal3, yaxt = "n", bty = "n", xaxt = "n", xlab = "", ylab = "")
dev.off()


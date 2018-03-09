
############
##Making supplemental figure one
############
fontsize = 15
MakeBox <- function(center){
  grid.rect(x=unit(center[1],"npc"),y=unit(center[2],"npc"),
            width = unit(0.125,"npc"),height = unit(0.28,"npc"),
            gp=gpar(fill="transparent",lwd=5,col="red"))
}

noReprText <- function(center){
  grid.text("larva\ns-s nurse\nrandom nurse",x=center[1],y=center[2],gp=gpar(fontsize=fontsize,fontface="bold"),just="center")
}

ReprText <- function(center){
  grid.text("larva\ns-s nurse",x=center[1],y=center[2],gp=gpar(fontsize=fontsize,fontface="bold"),just="center")
}

StageText <- function(center,label){
  grid.text(label,x=center[1],y=center[2],gp=gpar(fontsize=17,fontface="bold"),just="center")
}
labels = c("L1","L2","L3","L4","L5")

svg("~/Writing/Figures/NurseLarva/corApproach/FigS1_boxes.svg",height=3.5,width=12)
par(mar=c(0,0,0,0))
grid.lines(x=unit(c(0.09,0.13),"npc"),y=unit(c(0.5,0.7),"npc"),gp=gpar(lwd=2,lty="solid"))
grid.lines(x=unit(c(0.09,0.13),"npc"),y=unit(c(0.5,0.3),"npc"),gp=gpar(lwd=2,lty="solid"))

for (i in 1:5){
  StageText(c(0.05+(i*0.15),0.9),labels[i])
}

for (i in 0:4){
  MakeBox(c(0.2+(i*0.15),0.3))
  MakeBox(c(0.2+(i*0.15),0.7))
}

for (i in 0:4){
  noReprText(c(0.2+(i*0.15),0.3))
  ReprText(c(0.2+(i*0.15),0.7))
}

grid.text("queen-present",x = 0.937, y = 0.3,gp=gpar(fontsize=17,fontface="italic"),just="center")
grid.text("queen-absent",x = 0.937, y = 0.7,gp=gpar(fontsize=17,fontface="italic"),just="center")
grid.text("mixed\ncolony\nsource\n",x = 0.05, y = 0.47,gp = gpar(fontsize=17),just="center")
grid.text("Larval Stage:",x = 0.1, y = 0.9,gp=gpar(fontsize=18,fontface="bold"),just="center")
grid.rect(x=unit(0.05,"npc"),y=unit(0.50,"npc"),
          width = unit(0.08,"npc"),height = unit(0.34,"npc"),
          gp=gpar(fill="transparent",lwd=4,col="black"))
dev.off()



##########
##Differential expression based on caste fed
##########
f = droplevels(factors[grepl("XH|[^1]LCH",factors$sample.id),]) #skip the first stage (before caste can be identified)
f$stage = as.factor(f$stage)
design <- model.matrix(~type+stage+colony,data=f) #additive model, controlling for replicate and developmental stage
c = counts[,colnames(counts) %in% f$sample.id]
outH <- EdgeR(c,design,2)

f = droplevels(factors[grepl("XG|[^1]LCG",factors$sample.id),])
design <- model.matrix(~type+stage+colony,data=f)
c = counts[,colnames(counts) %in% f$sample.id]
outG <- EdgeR(c,design,2)

#There are 0 DE genes (FDR < 0.05) for gasters
Hgenes <- outH[outH$FDR < 0.05,] #Note all are upregulated, i.e. upregulated in reproductive nurses
ext <- read.csv("MpharAnn.csv")
ann <- ext[ext$Gene %in% rownames(Hgenes),] #Identify a few interesting DEGs




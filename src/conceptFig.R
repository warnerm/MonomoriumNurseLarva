########
##conceptual figure
library(rsvg)
image <- rsvg("~/Downloads/smaller larva with shadow.svg")
mult=dim(image)[1]/dim(image)[2]

#yT is the center of the cluster, layD is the number of points on each side of the grid,
#nGene is number of genes/cluster, rad is proportional to the size of the circle, in terms of proportion of grid
randCluster <- function(yT,layD,nGene,rad){
  rad = layD*rad
  ret = c()
  angle.inc = 2*pi/nGene
  for (i in 1:nGene){
    yadj = floor(sin(i*angle.inc)*rad+0.5)
    xadj = floor(cos(i*angle.inc)*rad+0.5)
    adj = yT+xadj+yadj*layD
    ret = c(ret,adj)
  }
  return(ret)
}


library(igraph)
nGene = 5
nGT = nGene*3
Soc_palette = c("#dadaeb","#9e9ac8","#6a51a3")


links <- data.frame(tissue1=c(rep(rep(c("larv","wh","wg"),each=nGene),nGT),"not","not"),
                    tissue2=c(rep(c("larv","wh","wg"),each=nGene*nGT),"not","not"),
                    ID1=paste("gene",c(rep(1:nGT,each=nGT),nGT+1,nGT+2)),
                    ID2=paste("gene",c(rep(seq(1:nGT),nGT),nGT+1,nGT+2)),weight=c(rep(rnorm(nGT*nGT,0,1)),5,5),
                    edgeCol = "gold")
links$edgeCol=as.character(links$edgeCol)
links$edgeCol[links$tissue1==links$tissue2]="#6a51a3"
links$edgeCol[links$tissue1!=links$tissue2]="#dadaeb"
links$edgeCol[(links$tissue1=="wh"&links$tissue2=="wg")|(links$tissue1=="wg"&links$tissue2=="wh")]="#9e9ac8"
links$weight[links$tissue1==links$tissue2]=rnorm(sum(links$tissue1==links$tissue2),0.9,1)

#Remove wh-wg connections
links = links[!(links$tissue1=="wh" & links$tissue2 == "wg") & !(links$tissue1=="wg" & links$tissue2 == "wh"),]

##Trim number of connections
keep1 = links[sample(rownames(links[links$edgeCol=="#6a51a3",]),24),]
keep2 = links[sample(rownames(links[links$edgeCol=="#dadaeb",]),5),]
#keep3 = links[sample(rownames(links[links$edgeCol=="#9e9ac8",]),5),]
keep = rbind(keep1,keep2)
keep = keep[keep$ID1!=keep$ID2,]
missing = unique(links$ID1)[!unique(links$ID1) %in% c(as.character(keep$ID1),as.character(keep$ID2))]
for (gene in missing){
  keep = rbind(keep,links[sample(rownames(links[links$ID1==gene,]),1),])
}

keep$lty = "solid"
keep$lty[keep$edgeCol=="#dadaeb"]="dotted"
#keep$lty[keep$edgeCol=="#9e9ac8"]="dashed"
keep$curved= FALSE
keep$curved[keep$edgeCol!="#6a51a3"]=0.2
keep = keep[,c(3:8)]

nodes <- paste("gene",seq(1:(nGT+2)))
net <- graph_from_data_frame(d=keep, vertices=nodes, directed=T) 

E(net)$width = rep(3,nrow(keep))

V(net)$color = rep("black",nGT+2)
E(net)$edge.color = rep("black",nrow(keep))

net <- simplify(net, remove.multiple = F, remove.loops = T) 


######
##Layout on lattice
######

layD = 100
lay <- layout_on_grid(make_lattice( c(layD,layD)))
lcenter = floor(layD/2)*layD+floor(layD/6)+4
whcenter = floor(3*layD/4)*layD+floor(3*layD/6)
wgcenter = floor(3*layD/5)*layD+layD+floor(5*layD/6)+1

larv = randCluster(lcenter,layD,nGene,0.1)
wh = randCluster(whcenter,layD,nGene,0.1)
wg = randCluster(wgcenter,layD,nGene,0.1)

layout =lay[c(larv,wh,wg,floor(layD/2),layD*layD-floor(layD/2)),]
V(net)$size=rep(4,nGT+2)
V(net)$color=c(rep(c("red","darkblue","lightblue"),each=nGene),"black","black")
V(net)$shape="square"

plot.new()
a=heightDetails(grid.text("within-tissue\n",x=0.15,y=0.88,gp=gpar(fontsize=20),just="left"))
lineh = as.numeric(substr(a,1,6))/8.3 ##8.3 is number of inches the figure is tall

svg("~/Github/MonomoriumNurseLarva/Figures/conceptFig.svg")
par(mar=c(0,0,0,0))
plot(net,edge.color=E(net)$edge.color,edge.arrow.size=0.8,vertex.label=NA,
     layout=layout,edge.curved=E(net)$curved)
grid.raster(image,x=.45,y=0.45,width=0.8/mult,height=0.8)
grid.curve(0.13,0.28,0.04,0.43,curvature=0.2,gp=gpar(lwd=1.5,col="darkgrey"))
grid.curve(0.2,0.32,0.26,0.48,curvature=-0.1,gp=gpar(lwd=1,col="darkgrey"))
grid.curve(0.45,0.31,0.37,0.67,curvature=0.2,gp=gpar(lwd=1.5,col="darkgrey"))
grid.curve(0.53,0.33,0.58,0.7,curvature=-0.2,gp=gpar(lwd=1,col="darkgrey"))
grid.curve(0.75,0.39,0.7,0.55,curvature=0.2,gp=gpar(lwd=1.5,col="darkgrey"))
grid.curve(0.85,0.37,0.9,0.5,curvature=-0.2,gp=gpar(lwd=1,col="darkgrey"))
grid.rect(x=0.45,y=0.05,width=0.1,height=0.05,gp=gpar(col="white"))
grid.rect(x=0.45,y=0.95,width=0.1,height=0.05,gp=gpar(col="white"))
grid.text("within-tissue\nsocial",x=0.15,y=0.86+lineh/2,gp=gpar(fontsize=20),just="left")
grid.lines(x=unit(c(0.04,0.13),"npc"),y=unit(c(0.86+lineh,0.86+lineh),"npc"),gp=gpar(lwd=1.5,lty="solid"))
grid.lines(x=unit(c(0.04,0.13),"npc"),y=unit(c(0.86-lineh/4,0.86-lineh/4),"npc"),gp=gpar(lwd=1.5,lty="dotted"))
grid.text(expression(underline("connection type")),x=0.15,y=0.96,just="left",gp=gpar(fontsize=22,fontface="bold"))
dev.off()
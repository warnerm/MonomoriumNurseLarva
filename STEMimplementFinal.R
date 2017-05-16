library(combinat)
library(plyr)
library(edgeR)
library(parallel)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(GSEABase)
library(GOstats)
library(GSEABase)
load("~/Dropbox/monomorium nurses/data.processed/cleandata.RData")
load("~/Dropbox/monomorium nurses/data.processed/GODB.RData")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#load annotation data, and also dnds, because the gene_id order is the same in dnds and annotation
dnds=read.table("~/Dropbox/monomorium nurses/data/dnds.txt",head=F)   #dim 22476
colnames(dnds)=c("gene_id","dnds")
ann=read.table("~/Dropbox/monomorium nurses/data/monomorium.annotation.txt",sep="\t",head=T,quote="\"",stringsAsFactors = FALSE)   # note, problem caused by "'" in file

#add gene_id column to annotation based on the same order in dnds file
ann$gene_id=dnds$gene_id
ann$HSPEvalueTR=as.numeric(ann$HSPEvalueTR)

#order list based on TR Evalue
ann=ann[order(ann$HSPEvalueTR),]

#remove duplicates so only isoform with highest HSPEvalueTR remaining; 12648 gene_id now
ann=ann[!duplicated(ann$gene_id),]   
annNames = ann[,c(10,23,43)]


#######################
###Part 1: Differential Expression Analysis
#######################



####Inital steps of standard edgeR analysis, as per manual
EdgeR <- function(data,design,coef){
	data <- DGEList(counts=data)
	data <- calcNormFactors(data)
	dat <- estimateGLMTrendedDisp(data, design)
	dat <- estimateGLMTagwiseDisp(dat, design)
	fit <- glmFit(dat,design)
	diff <- glmLRT(fit, coef=coef) 
	out <- topTags(diff,n=Inf,adjust.method="BH")$table 	###Calculate FDR
	return(out)
}

##Calculate number of DEs by queen presence at each stage
QPStageDE <- function(code){
	load("~/Dropbox/monomorium nurses/data.processed/cleandata.RData")
	counts <- counts[,grep(code,colnames(counts))]
	factors <- droplevels(factors[grep(code,rownames(factors)),])
	numDE = c()
	DEres = list()
	for (i in 1:5){
		f = factors[factors$stage==i,]
		c = counts[,colnames(counts) %in% rownames(f)]
		design <- model.matrix(~queen_presence,data=f)
		x <- EdgeR(c,design,2)
		numDE[i] = sum(x$FDR < 0.05)
		DEres[[i]] = x
	}

	design <- model.matrix(~queen_presence+colony,data=factors)
	overall <- EdgeR(counts,design,2)
	
	return(list(numDE,DEres,overall))
}


CH = QPStageDE("C.*WH")
CG = QPStageDE("CG")
FH = QPStageDE("F.*WH")
FG = QPStageDE("F.*WG")
W = QPStageDE("W.*_L")
LP = QPStageDE("P.*_L")

#########################
###Correlation over time
########################


##Main function to calculate stage-by-stage correlations and make heatmap
CorByStage <- function(code,name){
	load("~/Dropbox/monomorium nurses/data.processed/cleandata.RData")
	counts <- counts[,grep(code,colnames(counts))]
	fpkm <- fpkm[,grep(code,colnames(fpkm))]
	factors <- droplevels(factors[grep(code,rownames(factors)),])
	matCor <- CorMat(counts,fpkm,factors,Full=TRUE)
	p = CorHeat(matCor,code,name)
	return(list(matCor,p))
}

##Make heatmap of stage-by-stage correlations
CorHeat <- function(mat,code,name){
	m = as.data.frame(mat)
	colnames(m) = rownames(m) = c("L1","L2","L3","L4","L5")
	m$StageX = factor(rownames(m),levels=colnames(m))
	m1 = melt(m,id.vars="StageX")
	p <- ggplot(m1,aes(variable,StageX))+
		geom_tile(aes(fill=value))+
		xlab("Stage")+ylab("Stage")+
		scale_y_discrete(limits=rev(levels(m1$StageX)))+
		ggtitle(name)
	return(p)
}

##Make matrix of correlations
CorMat <- function(counts,fpkm,factors,Full,code){
	nStage = length(levels(factors$stage))
	mat = matrix(nrow = nStage, ncol = nStage)
	for (i in 1:nStage){
		for (j in i:nStage){
			i1 = 5 - nStage + i
			j1 = 5 - nStage + j
			if (Full){
				datA = fpkm[,factors$stage==i]
				datB = fpkm[,factors$stage==j]
				mat[i,j] = sum(rowSums(as.data.frame(cor(datA,datB))))/(ncol(datA) * ncol(datB))
			} else {
				if (i == j){
					mat[i,j] = 0
				}  else {
					mat[i,j] = numDE(counts,factors,i,j)
				}
			}
		}
	}
	return(mat)
}


LW = CorByStage("LW")
QW = CorByStage("QW")
LS = CorByStage("LS|1LW")
LFH = CorByStage("LFH|FLH")
QFH = CorByStage("QFH|FQH") ##Not enough sampling..
LFG = CorByStage("LFG|FLG")
QFG = CorByStage("QFG|FQG")
LCH = CorByStage("LCH")
QCH = CorByStage("QCH")
XH = CorByStage("XH|1LCH")
LCG = CorByStage("LCG")
QCG = CorByStage("QCG")
XG = CorByStage("XG|1LCG")
RH = CorByStage("RH")
RG = CorByStage("RG")

FH = CorByStage("F.*WH")
FG = CorByStage("F.*WG")
LCH2 = CorByStage("LCH")
QCH2 = CorByStage("QCH")
XH2 = CorByStage("XH|1LCH")

LW = CorByStage("LW","WorkerLarvae")
LS = CorByStage("LS|1LW","SexualLarvae")
LCH = CorByStage("LCH","QL WorkerNurseHead")
XH = CorByStage("XH|1LCH","SexualNurseHead")
LCG = CorByStage("LCG","QL WorkerNurseGaster")
XG = CorByStage("XG|1LCG","SexualNurseGaster")
QCH = CorByStage("QCH","QR WorkerNurseHead")
QCG = CorByStage("QCG","QR WorkerNurseGaster")
RH = CorByStage("RH","RandomNurseHead")
RG = CorByStage("RG","RandomNurseGaster")


adjN <- function(plot){
	p <- plot + scale_fill_continuous(limits=c(0.825, 1), breaks=seq(0.825,1,by=0.025))
	return(p)
}

adjL <- function(plot){
	p <- plot + scale_fill_continuous(limits=c(0.2, 1), breaks=seq(0.2,1,by=0.1))
	return(p)
}

png("AllHeats.png",width=3000,height=3000,res=300)
grid.arrange(adjL(LW[[3]]),adjL(LS[[3]]),adjN(LCH[[3]]),adjN(LCG[[3]]),adjN(XH[[3]]),adjN(XG[[3]]),nrow=3,ncol=2)
dev.off()

pdf("AllHeats.png")
grid.arrange(LW[[3]],LS[[3]],LCH[[3]],LCG[[3]],XH[[3]],XG[[3]],nrow=3,ncol=2)
dev.off()

png("AllHeatsQR.png",width=3000,height=3000,res=300)
grid.arrange(adjN(QCH[[3]]),adjN(QCG[[3]]),adjN(RH[[3]]),adjN(RG[[3]]),nrow=2,ncol=2)
dev.off()

CVbyStage <- function(code){
	load("~/Dropbox/monomorium nurses/data.processed/cleandata.RData")
	fpkm <- fpkm[,grep(code,colnames(fpkm))]
	factors <- droplevels(factors[grep(code,rownames(factors)),])
	CV = c()
	for (i in levels(factors$stage)){
		s_d = sd(rowSums(fpkm[,factors$stage==i]))
		m = mean(rowSums(fpkm[,factors$stage==i]))
		CV = c(CV,s_d/m)
	}
	return(CV)
}

CLW = CVbyStage("LW")
CLCG = CVbyStage("LCG")

CLS = CVbyStage("LS")
CLCH = CVbyStage("LCH")
CXH = CVbyStage("XH|1LCH")
CQCH = CVbyStage("QCH")


#########################
###Part 2: Dominant Patterns
#########################
#First, use STEMalg.R to sort genes into profiles.

#After transferring files:
files <- dir("STEMgroup/",pattern=".*membership.RData")

#trim file names to get the sample information
names <- lapply(1:length(files),function(x) gsub(" STEMmembership.RData","",files[x]))
STEMdata <- list()

#load all in, store, in STEM data
for (i in 1:length(files)){
	load(paste("STEMgroup/",files[i],sep=""))
	STEMdata[[i]] = results
}
names(STEMdata) = names

profiles4 <- read.csv("STEMprofiles4stage.csv")
profiles5 <- read.csv("STEMprofiles5stage.csv")

##Get shared significant profiles for two sets of STEM results
##Returns two lists. First is positive sharing, second is negative (in which we look for profiles in d2 that are the opposite of those in d1)
##Keep in mind this must be done between samples with the same stages
SharedProfiles <- function(d1,d2,nProfile){
	mod1 = GetMods(d1)
	mod2 = GetMods(d2)
	mod1N = nProfile + 1 - mod1 #Because profiles are generated sequentially, the opposite profile is in the opposite position of the list
	Pos = getShared(mod1,mod2,d2)
	Neg = getShared(mod1N,mod2,d2)
	totShared = unique(c(unlist(Pos[[4]]),unlist(Neg[[4]])))
	return(list(Pos,Neg,totShared))
}

getShared <- function(mod1,mod2,d){
	Shared = mod2[mod2 %in% mod1] ##Keeping track of modules from the second sample..
	if (length(Shared) == 0){ 
		return(list(NA,NA,NA,NA))
	} else {
		NGenes = getNgenes(Shared,d) ##We will keep track of number of genes for the second sample
		Genes = list()
		for (i in 1:length(Shared)){
			Genes[[i]] = names(d[[1]])[d[[1]] %in% Shared[i]] ##Gets genes for each profile
		}
		return(list(Shared = Shared,NGenes = NGenes,totGene = sum(NGenes),GeneNames = Genes))
	}
	
}

#Returns significant profiles from gene mSembership data
GetMods <- function(d){
	mod = unique(d[[1]])
	mod = mod[!is.na(mod)]
	return(mod)
}

#Returns number of genes belonging to each significant profile
getNgenes <- function(mods,d){
	num = c()
	for (i in 1:length(mods)){
		num = c(num,sum(d[[1]] %in% mods[i]))
	}
	return(num)	
}

WorkH = SharedProfiles(STEMdata[["WLarv"]],STEMdata[["FocNurseH"]],80)
WorkG = SharedProfiles(STEMdata[["WLarv"]],STEMdata[["FocNurseG"]],80)
RWorkH = SharedProfiles(STEMdata[["WLarv"]],STEMdata[["RNurseH"]],80)
RWorkG = SharedProfiles(STEMdata[["WLarv"]],STEMdata[["RNurseG"]],80)
QLWorkH = SharedProfiles(STEMdata[["WLarv"]],STEMdata[["FocNurseHQL"]],80)
QLWorkG = SharedProfiles(STEMdata[["WLarv"]],STEMdata[["FocNurseGQL"]],80)
QRWorkH = SharedProfiles(STEMdata[["WLarv"]],STEMdata[["FocNurseHQR"]],80)
QRWorkG = SharedProfiles(STEMdata[["WLarv"]],STEMdata[["FocNurseGQR"]],80)

SexH = SharedProfiles(STEMdata[["SLarv"]],STEMdata[["SNurseH"]],26)
SexG = SharedProfiles(STEMdata[["SLarv"]],STEMdata[["SNurseG"]],26)
QLWorkH4stage = SharedProfiles(STEMdata[["WLarvQL4"]],STEMdata[["FocNurseHQL4"]],26)
QLWorkG4stage = SharedProfiles(STEMdata[["WLarvQL4"]],STEMdata[["FocNurseGQL4"]],26)
QRWorkH4stage = SharedProfiles(STEMdata[["WLarvQR4"]],STEMdata[["FocNurseHQR4"]],26)
QRWorkG4stage = SharedProfiles(STEMdata[["WLarvQR4"]],STEMdata[["FocNurseGQR4"]],26)
SexH5 = SharedProfiles(STEMdata[["SexLarv5"]],STEMdata[["SexNurseH5"]],80)
SexG5 = SharedProfiles(STEMdata[["SexLarv5"]],STEMdata[["SexNurseG5"]],80)

CountsbyStage <- function(code,timepoints){
	load("~/Dropbox/monomorium nurses/data.processed/cleandata.RData")
	counts <- counts[,grep(code,colnames(counts))]
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

GenePlots <- function(name,code,mods,timepoints,title=NULL){
	StageExpr <- as.data.frame(CountsbyStage(code,timepoints))

	StageExpr$profile=STEMdata[[name]][[1]]
	StageExpr = StageExpr[StageExpr$profile %in% mods,]
	StageExpr$gene=rownames(StageExpr)
	SE = melt(StageExpr,id.vars=c("profile","gene"))
	colnames(SE) = c("profile","gene","Stage","Expression")
	levels(SE$Stage) = (6-timepoints):5
	SE$profile=factor(SE$profile,levels=mods)

	p <- ggplot(SE,aes(x=Stage,y=Expression,color=profile,group=gene))+
		geom_line(alpha=0.1)+
		theme_bw()+
	  ylim(-4,4)+
		theme(text=element_text(size=16))+
		scale_colour_manual(values=cbbPalette,name="Module")+
		xlab("Larval Developmental Stage")+
		ylab("Relative Expression")+
		theme(legend.position="none")+
		ggtitle(title)
	#ggsave(plot=p,file=paste(name,title,"OverTime.png",sep=""),height=11,width=9,dpi=300)
	return(p)
}

GetTop5mods <- function(name){
	d = as.numeric(names(sort(table(STEMdata[[name]][[1]]),decreasing=TRUE))[1:5])
	return(d)
}

NegMod <- function(mods,nProfile){
	mods = nProfile + 1 - mods
	return(mods)
}


GenePlots("SLarv","LS",GetTop5mods("SLarv"),4)
GenePlots("SNurseG","XG",GetTop5mods("SLarv"),4)


a1 = GenePlots("WLarv","W.*_L",NegMod(WorkH[[2]][[1]],80),5,"A")
a2 = GenePlots("FocNurseH","CH",WorkH[[2]][[1]],5,"C")

WGmods = WorkG[[1]][[1]][order(WorkG[[1]][[2]],decreasing=TRUE)]


a3 = GenePlots("WLarv","W.*_L",WGmods[1:5],5,"B")
a4 = GenePlots("FocNurseG","CG",WGmods[1:5],5,"D")
png("AllSTEMplot.png",width=3000,height=3000,res=300)
grid.arrange(a1,a3,a2,a4,nrow=2,ncol=2)
dev.off()


a1 = GenePlots("WLarv","W.*_L",NegMod(WorkH[[2]][[1]],80),5,"A")
a2 = GenePlots("FocNurseH","CH",WorkH[[2]][[1]],5,"C")

WGmods = WorkG[[1]][[1]][order(WorkG[[1]][[2]],decreasing=TRUE)]


a3 = GenePlots("WLarv","W.*_L",WGmods[1:5],5,"B")
a4 = GenePlots("FocNurseG","CG",WGmods[1:5],5,"D")
png("AllSTEMplot.png",width=3000,height=3000,res=300)
grid.arrange(a1,a3,a2,a4,nrow=2,ncol=2)
dev.off()


##################
##GRN analysis
##################

load("~/Dropbox/monomorium nurses/data.processed/cleandata.RData")
fpkm = log(fpkm + sqrt(fpkm ^ 2 + 1)) #hyperbolic sine transformation to normalize gene expression data
candidateList <- read.csv("~/Dropbox/MonomoriumCandidateGenes.csv")

#Make null and real networks of gene regulation occurring between and within groups of samples using genie 
MakeNets <- function(codes,names,boots){
  setwd("~/Downloads/GENIE3_R_C_wrapper")
  source("~/Downloads/GENIE3_R_C_wrapper/GENIE3.R")
	if (length(codes) != length(names)){
		return("Codes and Names not same length")
	}
	input <- GetExpr(codes,names)
	NullNetStats <- NullNetwork(input,boots) #Generate networks with candidate and random genes to make sure it "works"
	RealNetConns <- RealNet(input,boots) #Generate real networks by combining 1000 runs to be confident in connections
	RandomNetStats <- RandomNet(input,boots,names) #Generate a bunch of random networks and calculate connection strengths between each name
	return(list(NullNetStats,RealNetConns,RandomNetStats,names))
}


#Get data frame with expression, gene name, and annotation for a few candidates
GetExpr <- function(codes,names){
	d = list()
	for (i in 1:length(codes)){
		d[[i]] = ExprTime(codes[i],fpkm,names[i])
	}
	data = ldply(d,data.frame)
	
	annDat <- merge(data,candidateList,by="gene",all.x=TRUE)
	annDat$tissue_gene=with(annDat,paste(Samp,Common.Name,sep="_"))
	return(annDat)
}

#Calculate expression at each time point for a given set of samples
ExprTime <- function(code,fpkm,name){
	expr = matrix(nrow=nrow(fpkm))
	for (i in 1:5){
		e=rowSums(fpkm[,grepl(code,colnames(fpkm)) & factors$stage %in% i])/ncol(fpkm[,grepl(code,colnames(fpkm)) & factors$stage %in% i])
		expr=cbind(expr,e)
	}
	expr = expr[,-c(1)]
	colnames(expr)=paste("Samp",seq(1,5,by=1),sep="")
	expr = as.data.frame(expr)
	expr$Samp = rep(name,nrow(expr))
	expr$gene=rownames(fpkm)
	return(expr)
}

#Calculate connection strengths between candidates and random genes for a number of mock networks
NullNetwork <- function(input,boots){
	Rdata = matrix(nrow=boots,ncol=6)
	i=1
	while (i <= boots){
		Ginput = GenNullExpr(input)
		x <- try(GENIE3(as.matrix(Ginput)))
		if (!inherits(x,"try-error")){
		  Rdata[i,] = CalcConnsNull(get.link.list(x))
		  i = i+1
		}
	}
	colnames(Rdata)=c("ConnCand","ConnRand","ConnBetween","RankCand","RankRand","RankBetween")
	return(as.data.frame(Rdata))
}

#Make data frame of expression to input from genie that combines random genes and candidates
GenNullExpr <- function(input){
	candDat = input[input$Common.Name %in% CandGenes,] #Get expression of candidates
	nullDat = input[!input$gene %in% CandGenes,] #Get expression of all other genes
	NullGenes = sample(unique(nullDat$gene),length(CandGenes),replace=F) #Pick random genes
	NullData = input[nullDat$gene %in% NullGenes,] #Get random gene expression data
	NullData$tissue_gene=with(NullData,paste(tissue_gene,rep(1:length(CandGenes),each=3))) #Rename genes so they aren't identical..we just want to use NA for random genes
	GRNdata = rbind(candDat,NullData) #Combine real and random data
	Ginput = GRNdata[,c(2:6)]
	rownames(Ginput)=GRNdata$tissue_gene
	return(Ginput)
}

#Calculates mean connection strength and median rank of connection strength of connections between candidates, between candidate-random, and between random genes
CalcConnsNull <- function(links){
	Rdata = c()
	links$rank = as.numeric(rownames(links))
	links = links[gsub("_.*","",links$target.gene)==gsub("_.*","",links$regulatory.gene),] ##Sorts for connections within a sample type
	Rdata[1] = mean(links$weight[!grepl("NA",links$regulatory.gene)&!grepl("NA",links$target.gene)])
	Rdata[2] = mean(links$weight[grepl("NA",links$regulatory.gene)&grepl("NA",links$target.gene)])
	Rdata[3] = mean(links$weight[!grepl("NA",links$regulatory.gene)+!grepl("NA",links$target.gene)==1])
	Rdata[4] = median(links$rank[!grepl("NA",links$regulatory.gene)&!grepl("NA",links$target.gene)])
	Rdata[5] = median(links$rank[grepl("NA",links$regulatory.gene)&grepl("NA",links$target.gene)])
	Rdata[6] = median(links$rank[!grepl("NA",links$regulatory.gene)+!grepl("NA",links$target.gene)==1])
	return(Rdata)
}

#Return connection strengths calculated using only candidate genes, averaged out "boots" # of runs
RealNet <- function(input,boots){
	candDat = input[input$Common.Name %in% CandGenes,] #Get expression of candidates
	Ginput = candDat[,c(2:6)]
	rownames(Ginput)=candDat$tissue_gene
	MeanConns = IterateGenie(Ginput,boots)
	return(MeanConns)
}

#Run genie boots # of times, calculate average connection strength of each connection
IterateGenie <- function(Ginput,boots){
	GenieOut <- list()
	i = 0
	while (i <= boots){
	  x <- try(GENIE3(as.matrix(Ginput))) #generates matrix of connection strength
	  if (!inherits(x),"try-error"){
	    GenieOut[[i]] = x
	    i = i+1
	  }
	}

	MeanConns = GenieOut[[1]]

	for (i in 1:nrow(Ginput)){
		for (j in 1:nrow(Ginput)){
			sum = 0
			for (k in 1:boots){
				sum = sum + GenieOut[[k]][i,j]
			}
			MeanConns[i,j]=sum/boots #Calculate mean connection strength of a given interaction
		}
	}
	return(MeanConns)
}


#Calculate connection strengths between random genes for a number of mock networks. Return data from each mock network
RandomNet <- function(input,boots,names){
	SampleConns = list()
	for (i in 1:boots){
		Ginput <- genRandInput(input)
		Gres <- GENIE3(as.matrix(Ginput))
		SampleConns[[i]] = as.data.frame(getSampConns(Gres,names))
	}
	SampleConns = ldply(SampleConns,data.frame)
	return(SampleConns)
}

#Generate input for Genie of random genes
genRandInput <- function(input){
	NullGenes = sample(unique(input$gene),length(CandGenes),replace=F) #Pick random genes
	NullData = input[input$gene %in% NullGenes,] #Get random gene expression data
	NullData$tissue_gene=with(NullData,paste(tissue_gene,rep(1:length(CandGenes),each=3))) #Rename genes so they aren't identical..we just want to use NA for random genes
	Ginput = NullData[,c(2:6)]
	rownames(Ginput)=NullData$tissue_gene
	return(Ginput)
}

#Calculate mean connection strength between samples of a given type for a genie network
getSampConns <- function(Gres,names){
	Conns <- matrix(ncol = 3,nrow = length(names)*length(names))
	for (i in 1:length(names)){
		for (j in 1:length(names)){
			Conns[(i-1)*3+j,1] = names[i] #Regulatory sample
			Conns[(i-1)*3+j,2] = names[j] #Target Sample
			Conns[(i-1)*3+j,3] = mean(Gres[grep(names[i],rownames(Gres)),grep(names[j],colnames(Gres))]) #Regulatory sample is the row, target the column
		}
	}
	colnames(Conns) = c("RegulatorySamp","TargetSamp","weight")
	return(Conns)
}

#plot figures for genie results
genGraphs <- function(nets,name){
	setwd("~/Dropbox/workspace/")
	p1 <- plotRandCand(nets[[1]],name)
	p2 <- plotSampStrengthsReal(nets[[2]],name,nets[[4]])
	nets[[3]]$weight = as.numeric(as.character(nets[[3]]$weight))
	p3 <- plotSampStrength(nets[[3]],paste("RandNet",name))
	return(list(p1,p2,p3))
}

#Plot connection strengths of candidate-candidate, candidate-random, and random-random connections
plotRandCand <- function(net,name){
	dat = data.frame(GenePair=c(rep("Candidate-Candidate",boots),rep("Candidate-Random",boots),rep("Random-Random",boots)),
		MeanConnectionStrength=c(net$ConnCand,net$ConnBetween,net$ConnRand))
	dat$GenePair=as.factor(dat$GenePair)
	dat$MeanConnectionStrength=dat$MeanConnectionStrength/max(dat$MeanConnectionStrength)
	dat$id=rep(1:boots,times=3)
	dt = dcast(data=dat,formula = id~GenePair,value.var="MeanConnectionStrength")

	p <- ggplot(dat,aes(x=GenePair,y=MeanConnectionStrength))+
		geom_boxplot(notch=TRUE)+
			ylab("")+
			xlab("")+
			theme_bw()+
			theme(text=element_text(size=10),
				axis.title=element_text(size=15),
				axis.title.x=element_text(margin = unit(c(0.75,0,0,0), "cm")),
				axis.title.y=element_text(margin=unit(c(0,0.75,0,0), "cm")))
	ggsave(plot=p,file=paste(name,"NullGenieConn.png",sep=""),height=11,width=9,dpi=300)
  return(p)
}

#prepare candidate gene network results for plotting of connection strength based on target/regulatory samples
plotSampStrengthsReal <- function(net,name,SampNames){
	links <- get.link.list(net)
	links$RegulatorySamp = gsub("_.*","",links$regulatory.gene)
	links$TargetSamp = gsub("_.*","",links$target.gene)
	links$TargetSamp = factor(links$TargetSamp,levels=SampNames)
	links$RegulatorySamp = factor(links$RegulatorySamp,levels=SampNames)
	links$weight = links$weight/max(links$weight)

	p <- plotSampStrength(links,paste(name,"CandGenes"))
  return(p)
}

#plot connection strength based on sample combinations
plotSampStrength <- function(links,name){
	levels(links$TargetSamp) = paste(levels(links$TargetSamp)," ")
	links$weight = links$weight/max(links$weight)
	p <- ggplot(links,aes(x=RegulatorySamp,y=weight,fill=TargetSamp))+
		geom_boxplot(notch=TRUE)+
		theme_bw()+
		theme(text=element_text(size=18),axis.text.x=element_text(size=16))+
		scale_fill_manual(values=cbbPalette,name="Target Gene:  ")+
		theme(legend.position="bottom",legend.text=element_text(size=18))+
		theme(axis.text.y=element_text(size=16))+
		ylab("")+
		xlab("")
	ggsave(plot=p,file=paste(name,"SampStrength.png",sep=""),height=11,width=9,dpi=300)
  return(p)
}

CandGenes = c("transformer","Vg2","JHE1","ILP1","vasa","nanos","dsx","IRS","InR1","JHEH1")
boots=1000

codes = c("W.*_L","C.*WH","C.*WG")
names = c("WorkLarv","WorkNurseH","WorkNurseG")
Wnets = MakeNets(codes,names,boots)
Wgraphs <- genGraphs(Wnets,"Worker")

codes = c("LS|1LW","XH|1LCH","XG|1LCG")
names = c("SexLarv","SexNurseH","SexNurseG")
Snets = MakeNets(codes,names,boots)
Sgraphs <- genGraphs(Snets,"Sexual")

codes = c("QW","R.*WH","R.*WG")
names = c("QRLarv","RandNurseH","RandNurseG")
Rnets = MakeNets(codes,names,boots)
Rgraphs <- genGraphs(Rnets,"Random")

save(Wnets,Rnets,Snets,file="GenieNetsStats.RData")






























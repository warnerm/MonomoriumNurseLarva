setwd("~/Dropbox/workspace")
load("~/Dropbox/monomorium nurses/data.processed/cleandata.RData")
library(GOstats)
library(GSEABase)
library(data.table)

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

annGeneInfo = ann[,c(2,37:41,43)]

write.csv(annNames,file="annNames.csv")
#separate domBPId GO terms so each GO term has an individual line

#separate domBPId GO terms so each GO term has an individual line
go=ann[,c("gene_id","domBPId")]
go$domBPId=gsub("]---[",";",go$domBPId,fixed=TRUE)
go=subset(go,domBPId!="-")
go=data.table(go)
go=go[,list(domBPId=unlist(strsplit(domBPId,";"))),by=gene_id]
go$evidence<-"ISS"
colnames(go)=c("gene","GO","evidence")
setcolorder(go,c("GO","evidence","gene"))

#GOstats, define gene universe, etc.
universe <- as.character(unique(go$gene[go$gene %in% rownames(fpkm)]))
goFrame=GOFrame(go[go$gene %in% universe,],organism="Monomorium pharaonis")
goAllFrame=GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

GOstat <- function(genes){##Employ GOstat to get a list of GO terms for a given list of genes
  universe = 
	x = try(hyperGTest(
	 			GSEAGOHyperGParams(name = "",
                                    geneSetCollection=gsc,geneIds = intersect(genes,universe),
                                    universeGeneIds=universe,ontology = "BP",
                                    pvalueCutoff = 0.05,conditional = FALSE,
                                    testDirection = "over")))
	if (inherits(x,"try-error")){
		frame = data.frame(GOBPID=NA,Pvalue=NA,OddsRatio=NA,ExpCount=NA,
			Count=NA,Size=NA,Term=NA)
	} else {
		frame = summary(x)
	}
	return(frame)
}

GOstatUniverse <- function(genes,universe){##Employ GOstat to get a list of GO terms for a given list of genes
  x = try(hyperGTest(
    GSEAGOHyperGParams(name = "",
                       geneSetCollection=gsc,geneIds = intersect(genes,universe),
                       universeGeneIds=universe,ontology = "BP",
                       pvalueCutoff = 0.05,conditional = FALSE,
                       testDirection = "over")))
  if (inherits(x,"try-error")){
    frame = data.frame(GOBPID=NA,Pvalue=NA,OddsRatio=NA,ExpCount=NA,
                       Count=NA,Size=NA,Term=NA)
  } else {
    frame = summary(x)
  }
  return(frame)
}


GOstatmf <- function(genes){##Employ GOstat to get a list of GO terms for a given list of genes
	x = try(hyperGTest(
	 			GSEAGOHyperGParams(name = "",
                                    geneSetCollection=gsc,geneIds = intersect(genes,universe),
                                    universeGeneIds=universe,ontology = "MF",
                                    pvalueCutoff = 0.05,conditional = FALSE,
                                    testDirection = "over")))
	if (inherits(x,"try-error")){
		frame = data.frame(GOBPID=NA,Pvalue=NA,OddsRatio=NA,ExpCount=NA,
			Count=NA,Size=NA,Term=NA)
	} else {
		frame = summary(x)
	}
	return(frame)
}

GOstatcc <- function(genes){##Employ GOstat to get a list of GO terms for a given list of genes
	x = try(hyperGTest(
	 			GSEAGOHyperGParams(name = "",
                                    geneSetCollection=gsc,geneIds = intersect(genes,universe),
                                    universeGeneIds=universe,ontology = "CC",
                                    pvalueCutoff = 0.05,conditional = FALSE,
                                    testDirection = "over")))
	if (inherits(x,"try-error")){
		frame = data.frame(GOBPID=NA,Pvalue=NA,OddsRatio=NA,ExpCount=NA,
			Count=NA,Size=NA,Term=NA)
	} else {
		frame = summary(x)
	}
	return(frame)
}

save(gsc,universe,GOstat,GOstatUniverse,GOstatmf,GOstatcc,file="GODB.RData")





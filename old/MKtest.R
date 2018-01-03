snipre=read.csv("~/Dropbox/monomorium nurses/data/bayesian_results.csv",head=T)  # 10913 rows
snipre1=subset(snipre,PS>0&FR>0&FS>0)
rownames(snipre1)=snipre1$gene
snipre1$dSLs=snipre1$FS/snipre1$Tsil
snipre1$pSLs=snipre1$PS/snipre1$Tsil
m.dist.order <- order(mahalanobis(snipre1[c("dSLs","pSLs")], 
                                  colMeans(snipre1[c("dSLs","pSLs")]), cov(snipre1[c("dSLs","pSLs")])), decreasing=TRUE)
snipre1=snipre1[m.dist.order[16:length(m.dist.order)],]  

MKtest=snipre1[c("FR","Trepl","FS","Tsil","PR","Trepl","PS","Tsil")]  # new order for 3.2
MKtest$n=42
MKtest1=MKtest
MKtest1$Chr=1  # low diversity 3004
MKtest1[rownames(subset(snipre1,pSLs>=quantile(snipre1$pSLs,0.333))),"Chr"]=2  # mid diversity 3003
MKtest1[rownames(subset(snipre1,pSLs>=quantile(snipre1$pSLs,0.666))),"Chr"]=3  # high diversity 3013
write.table(MKtest1,file="~/Downloads/MKtest-3.2bugfree/allGenes.txt",row.names=FALSE,col.names=FALSE)

##After running MKtest
df <- read.table("~/Downloads/MKtest-3.2bugfree/AllGenes")
a = as.matrix(df)
colnames(df) = a[1,]
df = df[-c(1),]
f.est = as.numeric(as.matrix(df[1,grepl("^f",names(df))]))
names(f.est) = rownames(MKtest1)
write.csv(f.est,file="~/Data/MKtestConstraintOneAlpha.csv")

##After running MKtest with gene-specific alpha
df <- read.table("~/Downloads/MKtest-3.2bugfree/AllGenesMultiAlphaMultiConstraint")
a = as.matrix(df)
colnames(df) = a[1,]
df = df[-c(1),]
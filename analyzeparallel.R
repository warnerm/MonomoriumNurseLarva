

ext <- read.csv("~/Downloads/msx123_Supp (1)/External Database S1.csv")

snipre <- read.csv("~/Dropbox/monomorium nurses/data/bayesian_results.csv")
sn = snipre[,c(1,10,11,15,18,21,22,25)]

dnds <- read.table("~/Dropbox/monomorium nurses/data/dnds.txt")
colnames(dnds) = c("gene","dnds")
ann=read.table("~/Dropbox/monomorium nurses/data/monomorium.annotation.txt",sep="\t",head=T,quote="\"",stringsAsFactors = FALSE)   # note, problem caused by "'" in file
ann$gene=dnds$gene
a1 = ann[,c("TransLength","gene")]
ann$HSPEvalueTR=as.numeric(ann$HSPEvalueTR)

#order list based on TR Evalue
ann=ann[order(ann$HSPEvalueTR),]

#remove duplicates so only isoform with highest HSPEvalueTR remaining; 12648 gene_id now
ann=ann[!duplicated(ann$gene),]   

dnds <- merge(dnds,a1,by="gene")
dnds <- dnds[order(dnds$TransLength,decreasing=TRUE),]
dnds <- dnds[!duplicated(dnds$gene),]

f.est <- read.csv("~/Data/MKtestConstraintOneAlpha.csv")
colnames(f.est) = c("gene","f.est")
ext = merge(ext,f.est,by.x="Gene",by.y="gene")

glmConn <- function(name){
  df <- read.csv(paste("~/Data/",name,".csv",sep=""))
  
  df$SIL = df$Lbetween - df$Lwithin
  df$SIWH = df$WHbetween - df$WHwithin
  df$SIWG = df$WGbetween - df$WGwithin
  df$SI_Overall = (df$SIL + df$SIWH + df$SIWG)/3
  
  dfAll = merge(ext,df,by.x="Gene",by.y="gene")
  dfAll = merge(dfAll,sn,by.x="Gene",by.y="gene")
  
  AllLm <- SixLm(dfAll)
  
  AllDf <- parseLm(AllLm)
  
  temp = list()
  for (i in 1:3){
    temp[[i]] =  rbind(AllDf[[i*2-1]],AllDf[[i*2]])
    colnames(temp[[i]]) = paste(colnames(temp[[i]]),gsub("\\..*","",names(AllDf)[i*2]))
  }
  
  df <- do.call(cbind,temp)
  return(list(df,AllLm))
}

parseLm <- function(AllLm){
  AllDf <- list()
  for (i in 1:length(AllLm)){
    df <- summary(AllLm[[i]])$coefficients
    AllDf[[i]] = df[2:nrow(df),c(3,4)]
  }
  names(AllDf) = names(AllLm)
  return(AllDf)
}

SixLm <- function(dfAll){
  lm1 <- glm(BSnIPRE.est.x ~ Lbetween+Lwithin+WHbetween+WHwithin+WGbetween+WGwithin,data=dfAll)
  lm2 <- glm(BSnIPRE.est.x ~ SIL+SIWH+SIWG,data=dfAll)
  lm3 <- glm(log(f.est) ~ Lbetween+Lwithin+WHbetween+WHwithin+WGbetween+WGwithin,data=dfAll)
  lm4 <- glm(log(f.est) ~ SIL+SIWH+SIWG,data=dfAll)
  lm5 <- glm(log(BSnIPRE.f) ~ Lbetween+Lwithin+WHbetween+WHwithin+WGbetween+WGwithin,data=dfAll)
  lm6 <- glm(log(BSnIPRE.f) ~ SIL+SIWH+SIWG,data=dfAll)
  
  AllLm <- list(lm1,lm2,lm3,lm4,lm5,lm6)
  names(AllLm) = c("est.conn","est.soc","MKf.conn","MKf.soc","BSnIPREf.conn","BSnIPREf.soc")
  return(AllLm)
}

resWexpr <- glmConn("TopExprWorkerNetNetSocialityDF")
resSexpr <- glmConn("TopExprSexualNetNetSocialityDF")
resRexpr <- glmConn("TopExprRandomNetNetSocialityDF")

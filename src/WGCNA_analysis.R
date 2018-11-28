library(WGCNA)

load("~/Data/Nurse_Larva/WGCNA_metanet.RData")
load("~/Data/Nurse_Larva/DEstage_results.RData")
for (lev in names(DEgene)){
  DEgene[[lev]][[2]]$Gene = rownames(DEgene[[lev]][[2]])
}

#For metanet, want to find modules with nurse and larval genes that are DE by stage, 
#with the expectation that there will be a positive association (Compared to random nurses)

load("~/Data/Nurse_Larva/Networkfocal_combat.RData")
MEs_foc = MEs
ML_foc = moduleLabels
MC_foc = moduleColors

stageGLM <- function(ME,stage){
  d = data.frame(ME=ME,stage=stage)
  lm <- glm(ME ~ as.numeric(stage),data=d)
  a = drop1(lm,.~.,test="Chi") 
  return(a$`Pr(>Chi)`[2])
}

datExpr <- fpkm[,grepl("CH|CG|W_L",colnames(fpkm))]
stage <- gsub("X[0-9]","",colnames(datExpr))
stage <- gsub("[Q|L].*","",stage)
tissue <- gsub(".*_","",colnames(datExpr))

glmRes <- lapply(c("L","WH","WG"),function(t) apply(MEs_foc[tissue==t,],2,function(x) stageGLM(x,stage[tissue==t])))
sigStage_foc <- lapply(glmRes,function(x) names(x)[x < 0.05])
cor.test(glmRes[[1]],glmRes[[2]],method="spearman")
cor.test(glmRes[[1]],glmRes[[3]],method="spearman")

load("~/Data/Nurse_Larva/Networkrandom_combat.RData")
MEs_rand = MEs
ML_rand = moduleLabels
MC_rand = moduleColors

datExpr <- fpkm[,grepl("RH|RG|W_L",colnames(fpkm))]
stage <- gsub("X[0-9]","",colnames(datExpr))
stage <- gsub("R.*","",stage)
tissue <- gsub(".*_","",colnames(datExpr))

glmRes <- lapply(c("L","WH","WG"),function(t) apply(MEs_rand[tissue==t,],2,function(x) stageGLM(x,stage[tissue==t])))
sigStage_rand <- lapply(glmRes,function(x) names(x)[x < 0.05])
cor.test(glmRes[[1]],glmRes[[2]],method="spearman")
cor.test(glmRes[[1]],glmRes[[3]],method="spearman")
cor.test(glmRes[[2]],glmRes[[3]],method="spearman")


df <- data.frame(colorFoc = MC_foc, colorRand = MC_rand, Gene = rownames(fpkm))
dfDE <- join_all(lapply(DEgene,function(x) data.frame(Gene = rownames(x[[2]]),Pval = x[[2]]$PValue)),by="Gene")
colnames(dfDE)[2:6]= names(DEgene)
df_all <- merge(df,dfDE,by = "Gene")
df_all_cut <- df_all
df_all_cut[,c(4:8)] = apply(df_all_cut[,c(4:8)],2,function(x) x < 0.05)

df_mods <- melt(df_all,id.vars = colnames(df_all_cut)[-c(2,3)])
df_sumP <- ddply(df_mods, ~variable + value,summarize,
                 CH = -sum(log(CH))/length(CH),
                 CG = -sum(log(CG))/length(CH),
                 RH = -sum(log(RH))/length(CH),
                 RG = -sum(log(RG))/length(CH),
                 W_L = -sum(log(W_L))/length(CH))
df_sumPr = df_sumP[df_sumP$variable == "colorRand",]
df_sumPf = df_sumP[df_sumP$variable == "colorFoc",]

cor.test(df_sumPf$CH,df_sumPf$W_L)
cor.test(df_sumPf$CG,df_sumPf$W_L)
cor.test(df_sumPr$RH,df_sumPr$W_L)
cor.test(df_sumPr$RG,df_sumPr$W_L)

cor.test(df_mods$CH,df_mods$W_L)
cor.test(df_mods$CG,df_mods$W_L)
cor.test(df_mods$RH,df_mods$W_L)
cor.test(df_mods$RG,df_mods$W_L)

stageGLM <- function(ME,stage){
  d = data.frame(ME=ME,stage=stage)
  lm <- glm(ME ~ stage,data=d)
}

load("~/Data/Nurse_Larva/Networkfocal_metasample.RData")
d <- metaExpr(codes = c('CH','CG','W_L'))
stage <- gsub("X[0-9]","",colnames(d))
stage <- gsub("[Q|L]","",stage)

df <- data.frame(Gene = gsub(".*_","",rownames(d)), tissue = gsub("_.*","",rownames(d)),module = moduleColors)
df$tissue = as.character(df$tissue)
df$tissue[df$tissue=="W"]="W_L"
df$tissue = as.factor(df$tissue)

df_tissue <- lapply(levels(df$tissue),function(x) df[df$tissue==x,])
names(df_tissue)=levels(df$tissue)
df_tissue <- lapply(levels(df$tissue),function(x) merge(df_tissue[[x]],DEgene[[x]][[2]],by="Gene"))

for (i in 1:3){
  df_tissue[[i]]$DE = "no"
  df_tissue[[i]]$DE[df_tissue[[i]]$PValue < 0.05]="yes"
}

df_all = ldply(df_tissue)

summed <- ddply(df_all,~tissue + module,summarize,
                N = length(DE),
                NDE = sum(DE=="yes"),
                propDE = sum(DE=="yes")/length(DE),
                propStrictDE = sum(FDR < 0.05)/length(DE))

summed$N[is.na(summed$N)]=0
m = melt(summed,id.vars = c("tissue","module"))
s = dcast(m,module+variable~tissue,value.var = "value")
for (i in 3:5){
  s[,i][is.na(s[,i]) & s$variable=="N"]=0
}

stageGLM <- function(ME,stage){
  d = data.frame(ME=ME,stage=stage)
  lm <- glm(ME ~ as.numeric(stage),data=d)
  a = drop1(lm,.~.,test="Chi") 
  return(a$`Pr(>Chi)`[2])
}

tissues = c("CH","CG","W_L")
MEstage<- apply(MEs,2,function(x) stageGLM(x,stage))


cor.test(s$CH[s$variable=="NDE"],s$W_L[s$variable=="NDE"],method = "spearman")


load("~/Data/Nurse_Larva/Networkrandom_metasample.RData")
d <- metaExpr(codes = c('RH','RG','W_L'))
d = d[apply(d,1,sd)!=0,]
stage <- gsub("X[0-9]","",colnames(d))
stage <- gsub("[Q|L]","",stage)

df <- data.frame(Gene = gsub(".*_","",rownames(d)), tissue = gsub("_.*","",rownames(d)),module = moduleColors)
df$tissue = as.character(df$tissue)
df$tissue[df$tissue=="W"]="W_L"
df$tissue = as.factor(df$tissue)

df_tissue <- lapply(levels(df$tissue),function(x) df[df$tissue==x,])
names(df_tissue)=levels(df$tissue)
df_tissue <- lapply(levels(df$tissue),function(x) merge(df_tissue[[x]],DEgene[[x]][[2]],by="Gene"))

for (i in 1:3){
  df_tissue[[i]]$DE = "no"
  df_tissue[[i]]$DE[df_tissue[[i]]$PValue < 0.05]="yes"
}

df_all = ldply(df_tissue)

summed <- ddply(df_all,~tissue + module,summarize,
                N = length(DE),
                NDE = sum(DE=="yes"),
                propDE = sum(DE=="yes")/length(DE),
                propStrictDE = sum(FDR < 0.05)/length(DE))

summed$N[is.na(summed$N)]=0
m = melt(summed,id.vars = c("tissue","module"))
s = dcast(m,module+variable~tissue,value.var = "value")
for (i in 3:5){
  s[,i][is.na(s[,i]) & s$variable=="N"]=0
}

cor.test(s$RH[s$variable=="NDE"],s$W_L[s$variable=="NDE"])


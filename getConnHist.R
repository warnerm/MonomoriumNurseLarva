load("cleandata.RData")
library(rslurm)

#Return list of expression matrices for correlation analysis
formatExpr <- function(dN,dL){
  d = lapply(c(dN,dL),subExpr) #get expression matrix for each of three sample types
  for (i in 1:2){
    #Needed to align random nurses to QR larvae
    colnames(d[[i]]) = gsub("R","Q",colnames(d[[i]]))
  }
  dF = alignStage(list(d[[1]],d[[2]]))
  return(dF) #Returns data for foc/larv comparison and rand/larv comparison
}

#only keeps relevant genes and lines up larvae and nurses by stage
subExpr <- function(code){
  f = fpkm[grepl(code,colnames(fpkm))]
  colnames(f) = substr(colnames(f),start=1,stop=4) #Has colony, stage and queen presence information
  return(f)
}

#Take formatted expression data and remove instances of colonies that don't match (missing data)
alignStage <- function(d){
  d[[1]]=d[[1]][,colnames(d[[1]]) %in% colnames(d[[2]])]
  d[[2]]=d[[2]][,colnames(d[[2]]) %in% colnames(d[[1]])]
  return(d)
}

getConn <- function(i){
  return(sum(abs(unlist(lapply(seq(1,nrow(expr[[1]])),function(j){
    cor(t(expr[[1]][i,]),t(expr[[2]][j,]))^6
  })))))
}

#Get vector of nurse-larva correlations
vecCor <- function(dN,dL){
  expr <<- formatExpr(dN,dL)
  df <- data.frame(i=seq(1,nrow(expr[[1]]),by=1))
  sjob <- slurm_apply(getConn, df, jobname = 'getConn',
                      nodes = 4, cpus_per_node = 20, add_objects=c("expr"),submit = TRUE)
  res <- get_slurm_out(sjob,outtype='raw') #get output as lists
  return(res)
}

vecH <- vecCor("CH","W_L")
vecG <- vecCor("CG","W_L")

vecRH <- vecCor("RH","QW")
vecRG <- vecCor("RG","QW")

vecQCH <- vecCor("QCH","QW")
vecQCG <- vecCor("QCG","QW")

connRes <- list(vecH,vecG,vecRH,vecRG,vecQCH,vecQCG)
save(connRes,file = "connRes.RData")

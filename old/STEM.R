#First, use STEMalg.R to sort genes into profiles.

#After transferring files:
files <- dir("~/Data/Nurse_Larva/",pattern=".*membership.RData")

#trim file names to get the sample information
names <- lapply(1:length(files),function(x) gsub(" STEMmembership.RData","",files[x]))
STEMdata <- list()

#load all in, store, in STEM data
for (i in 1:length(files)){
  load(paste("~/Data/Nurse_Larva/",files[i],sep=""))
  STEMdata[[i]] = results
}
names(STEMdata) = names

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
  return(profiles)
}

profiles5 <- SetProfiles(5)
# 
# ##Get shared significant profiles for two sets of STEM results
# ##Returns two lists. First is positive sharing, second is negative (in which we look for profiles in d2 that are the opposite of those in d1)
# ##Keep in mind this must be done between samples with the same stages
# SharedProfiles <- function(d1,d2,nProfile){
#   mod1 = GetMods(d1)
#   mod2 = GetMods(d2)
#   mod1N = nProfile + 1 - mod1 #Because profiles are generated sequentially, the opposite profile is in the opposite position of the list
#   Pos = getShared(mod1,mod2,d2)
#   Neg = getShared(mod1N,mod2,d2)
#   totShared = unique(c(unlist(Pos[[4]]),unlist(Neg[[4]])))
#   return(list(Pos,Neg,totShared))
# }
# 
# getShared <- function(mod1,mod2,d){
#   Shared = mod2[mod2 %in% mod1] ##Keeping track of modules from the second sample..
#   if (length(Shared) == 0){ 
#     return(list(NA,NA,NA,NA))
#   } else {
#     NGenes = getNgenes(Shared,d) ##We will keep track of number of genes for the second sample
#     Genes = list()
#     for (i in 1:length(Shared)){
#       Genes[[i]] = names(d[[1]])[d[[1]] %in% Shared[i]] ##Gets genes for each profile
#     }
#     return(list(Shared = Shared,NGenes = NGenes,totGene = sum(NGenes),GeneNames = Genes))
#   }
#   
# }
# 
# #Returns significant profiles from gene mSembership data
# GetMods <- function(d){
#   mod = unique(d[[1]])
#   mod = mod[!is.na(mod)]
#   return(mod)
# }
# 
# #Returns number of genes belonging to each significant profile
# getNgenes <- function(mods,d){
#   num = c()
#   for (i in 1:length(mods)){
#     num = c(num,sum(d[[1]] %in% mods[i]))
#   }
#   return(num)	
# }
# 
# WorkH = SharedProfiles(STEMdata[["WLarv"]],STEMdata[["NurseH"]],80)
# WorkG = SharedProfiles(STEMdata[["WLarv"]],STEMdata[["NurseG"]],80)
# RWorkH = SharedProfiles(STEMdata[["WlarvQR"]],STEMdata[["RNurseH"]],80)
# RWorkG = SharedProfiles(STEMdata[["WlarvQR"]],STEMdata[["RNurseG"]],80)
# QRWorkH = SharedProfiles(STEMdata[["WlarvQR"]],STEMdata[["NurseHQR"]],80)
# QRWorkG = SharedProfiles(STEMdata[["WlarvQR"]],STEMdata[["NurseGQR"]],80)
# QLWorkH = SharedProfiles(STEMdata[["WLarvQL"]],STEMdata[["NurseHQL"]],80)
# QLWorkG = SharedProfiles(STEMdata[["WLarvQL"]],STEMdata[["NurseGQL"]],80)


###As input to the algorithm, we need the mean expression across samples of a given type at each developmental stage
timepoints = 5
CountsbyStage <- function(code){
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

sharedPlot <- function(nMod,lMod,df){
  dN = df[(df$Module %in% nMod & df$tissue == "nurse"),]
  dL = df[(df$Module %in% lMod & df$tissue == "larva"),]
  dN$Module = factor(dN$Module,levels = nMod)
  dL$Module = factor(dL$Module, levels = lMod)
  p1 <- ggplot(dN,aes(x = as.integer(as.numeric(variable)), y = value, color = Module,group = Gene))+
    geom_line(alpha = 0.3)+
    theme_bw()+
    xlab("larval developmental stage")+
    ylab("log2 expression change")
  p2 <- ggplot(dL,aes(x = as.integer(as.numeric(variable)), y = value, color = Module,group = Gene))+
    geom_line(alpha = 0.2)+
    theme_bw()+
    xlab("larval developmental stage")+
    ylab("log2 expression change")
  return(list(p1,p2))
}

Hshare <- sharedDf("CH","W_L","WLarv","NurseH")
Gshare <- sharedDf("CG","W_L","WLarv","NurseG")

#All negative
lMod = Hshare[[3]]
Hmod = 81 - as.integer(as.character(Hshare[[3]]))
Hplot <- sharedPlot(Hmod,lMod,Hshare[[1]])

#Many positive for abdomens
nGene = unlist(lapply(Gshare[[2]],function(x){
  sum(Gshare[[1]]$Module==x & Gshare[[1]]$tissue =="nurse",na.rm = TRUE)
}))

Gmod = Gshare[[2]][order(nGene,decreasing = TRUE)][1:5]
Gplot <- sharedPlot(Gmod,Gmod,Gshare[[1]])

png("~/Writing/Figures/NurseLarva/corApproach/fig1.png",height = 3000,width = 4000,res =300)
grid.arrange(Hplot[[1]] + ggtitle("nurse head"),
             Hplot[[2]] + ggtitle("larva shared with nurse head"),
             Gplot[[1]] + ggtitle("nurse abdomen"),
             Gplot[[2]] + ggtitle("larva shared with nures abdomen"))
dev.off()



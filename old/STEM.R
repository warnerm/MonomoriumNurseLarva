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

WorkH = SharedProfiles(STEMdata[["WLarv"]],STEMdata[["NurseH"]],80)
WorkG = SharedProfiles(STEMdata[["WLarv"]],STEMdata[["NurseG"]],80)
RWorkH = SharedProfiles(STEMdata[["WlarvQR"]],STEMdata[["RNurseH"]],80)
RWorkG = SharedProfiles(STEMdata[["WlarvQR"]],STEMdata[["RNurseG"]],80)
QRWorkH = SharedProfiles(STEMdata[["WlarvQR"]],STEMdata[["NurseHQR"]],80)
QRWorkG = SharedProfiles(STEMdata[["WlarvQR"]],STEMdata[["NurseGQR"]],80)
QLWorkH = SharedProfiles(STEMdata[["WLarvQL"]],STEMdata[["NurseHQL"]],80)
QLWorkG = SharedProfiles(STEMdata[["WLarvQL"]],STEMdata[["NurseGQL"]],80)

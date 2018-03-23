library(WGCNA)

load("~/Data/Nurse_Larva/WGCNA_metanet.RData")
load("~/Data/Nurse_Larva/DEstage_results.RData")

#For metanet, want to find modules with nurse and larval genes that are DE by stage, 
#with the expectation that there will be a positive association (Compared to random nurses)
dFE = read.delim("Input-EFSA/Fieldedge.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
ds = dFE[dFE$Species == "LECAcl4plb", ]
write.table(ds, file="Input-EFSA/Fieldedge_LECAcl4plb.txt", sep="\t", row.names=FALSE, quote=FALSE)

#################
dBE2TH_full <- read.csv("Input-EFSA/Plant traits full dataset_BE2TH.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
dBE2TH_full <- dBE2TH_full[dBE2TH_full$Clonal != "unknown", ]
dBE2TH_full <- dBE2TH_full[!is.na(dBE2TH_full$LMR), ]
dBE2TH_full <- dBE2TH_full[!is.na(dBE2TH_full$MaxMass..mg.) , ]
dBE2TH_full <- dBE2TH_full[!is.na(dBE2TH_full$mSeed..mg.), ]
dBE2TH_full <- dBE2TH_full[dBE2TH_full$Note != "Tree", ] # remove Trees
write.csv(dBE2TH_full, file="Input-EFSA/BE2TH_cleaned.csv", row.names=FALSE)

# select biogeography
# "STE_I1", "ANA_I1", "BLS_I1", "ALP_I1", "ATL_I1", "BOR_I1", "CON_I1", "PAN_I1", "MED_I1", "MAC_I1"
grBE2TH = c("STE_I1", "BLS_I1", "ALP_I1", "ATL_I1", "BOR_I1", "CON_I1", "PAN_I1", "MED_I1", "MAC_I1")
ls_dBE2TH <- setNames(
  lapply(grBE2TH, function(g) dBE2TH_full[dBE2TH_full[[g]] == 1, ]),
  grBE2TH
)
lapply(ls_dBE2TH, nrow) # check number of species per group

# dfefsa = read.delim("Model-files/Fieldedge_efsa.txt")
# dfefsa$Species = gsub(" ", "_", dfefsa$Species)
# dfefsa$Species = gsub("\\.", "", dfefsa$Species)
# dfefsa$Species = iconv(dfefsa$Species, "latin1", "ASCII//TRANSLIT", sub="")
# colnames(dfefsa)[colnames(dfefsa) == "SeedMass"] <- "mSeed"
# dfefsa$palat = 0.5 # cannot be 0.0,  exclude after
# dfefsa$Overwintering = NULL
# write.table(dfefsa, "Model-files/Fieldedge_efsa1.txt", col.names=FALSE, row.names=FALSE, sep="\t")
# 
# 
# dfefsa2 = dfefsa[1:52, ]
# write.table(dfefsa2, "Model-files/Fieldedge_efsa2.txt", col.names=FALSE, row.names=FALSE, sep="\t")
# dforigin = read.delim("Model-files/Fieldedge1.txt") ; dforigin1=dforigin
# colnames(dfefsa2)
# colnames(dforigin)
# COL = c("Species", "MaxAge", "AllocSeed", "LMR", "m0", "MaxMass", "mSeed",
#         "Dist", "pEstab", "Gmax", "SLA", "palat", "memo")
# for(c in COL){
#     dforigin1[[c]] = dfefsa2[[c]]
# }
# write.table(dforigin1, "Model-files/Fieldedge2.txt", col.names=FALSE, row.names=FALSE, sep="\t")

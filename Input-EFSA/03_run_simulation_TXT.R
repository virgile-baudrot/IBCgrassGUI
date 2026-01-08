source("Input-EFSA/02_simulation_setting.R")

dfBE2TH <- read.csv("Input-EFSA/BE2TH_cleaned.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

PATH = "Model-files/"
grBE2TH = c("STE_I1", "BLS_I1", "ALP_I1", "ATL_I1", "BOR_I1", "CON_I1", "PAN_I1", "MED_I1", "MAC_I1")
STRESSOR_TYPE = c("low", "high", "medium")
STRESSORS = c(no="NO", pb="biomass", sb="SEbiomass", sv="survival", es="establishment", st="sterility", ns="seednumber")

############################## LOAD ARGUMENTS
# group = "STE_I1"
# stressor_type = "high"
# stressor="NO"

args <- commandArgs(trailingOnly = TRUE)

print(paste("group :", args[1]))
print(paste("stressor_type :", args[2]))
print(paste("stressor :", args[3]))
group =  args[1]
stressor_type =  args[2]
stressor= args[3]
#############################################

idScenario=build_id_scenario(group, stressor_type, stressor)

if(stressor=="NO"){
    CONTROL = TRUE
} else{
    CONTROL = FALSE
}

ModelVersion <- 3
CellNum <- 50 # 173
Tmax <- 9 # if(CONTROL) 50 else 100
InitDuration <- 3 # 35
NamePftFile <- "Fieldedge.txt"
SeedInput <- 10
belowres <- 90
abres <- 100
abampl <- 0
tramp <- 0.1
graz <- 0
cut <- 0 # or 1,2,3 ...
week_start <- 1
HerbDuration <- 10
HerbEffectType <- 1
EffectModel <- 0
scenario <- idScenario
MC <- 0
nMC <- 5

# SAVE CONFIGURATION AND COPY INPUT FILES TO OUTPUT FOLDER
# ID_OUT = formatC(sample(1:1e8,1), width = 9, format = "d", flag = "0")
# PATH_OUTPUT = paste0("sim_", group,"_", stressor, "_", stressor_type, "_", ID_OUT, "/")
# dir.create(PATH_OUTPUT)

# Create Fieldedge.txt for EFSA simulations
df_fieldedge <- build_Fieldedge_BE2TH(
    dfBE2TH, group = group,
    stressor=stressor,
    stressor_type=stressor_type)
# ################ REDUCE
# df_fieldedge <- df_fieldedge[1:3, ]
# ######################
write.table(df_fieldedge, file=paste0(PATH, "Fieldedge.txt"), sep="\t", row.names=FALSE, quote=FALSE)
# write.table(df_fieldedge, file=paste0(PATH_OUTPUT, "Fieldedge.txt"), sep="\t", row.names=FALSE, quote=FALSE)

# # Create AppRate.txt
# df_AppRate <- build_AppRate(application_rate=HerbEffectType, herbicide_duration=HerbDuration)
# write.table(df_AppRate, paste0(PATH,"AppRate.txt"), col.names=FALSE, row.names=FALSE, sep="\t")
# # write.table(df_AppRate, paste0(PATH_OUTPUT,"AppRate.txt"), col.names=FALSE, row.names=FALSE, sep="\t")

# Create HerbFactor.txt for control simulation
df_HerbFact<- build_HerbFact(HerbDuration=HerbDuration, stressor=stressor, stressor_type=stressor_type)
write.table(df_HerbFact, paste0(PATH,"HerbFact.txt"), col.names=TRUE, row.names=FALSE, sep="\t")

# save_configuration(filepath = paste0(PATH_OUTPUT, "simulation_config"),
save_configuration(filepath = paste0(PATH, "simulation_config"),
                   ModelVersion, CellNum, Tmax, InitDuration, NamePftFile,
                   SeedInput, belowres, abres, abampl, tramp, graz, cut,
                   week_start, HerbDuration, HerbEffectType,
                   EffectModel, scenario, MC, nMC
)

## Setup for parallel processing
# no_cores <- max(detectCores()-2, 1)
# cl <- makeCluster(no_cores)
# registerDoParallel(cl)
# setwd('Model-files')
# CellNum <- 4
# Tmax <- 2
# InitDuration <- 1
# EffectModel <- 0
# MC <- 2
# mycall <- paste('./IBCgrassGUI', ModelVersion, CellNum, Tmax, InitDuration,
#                 NamePftFile, SeedInput, belowres, abres, abampl, tramp, graz, cut,
#                 week_start, HerbDuration, HerbEffectType, EffectModel, scenario, MC, sep=" ")
# system(mycall, intern=TRUE)
# setwd('..')

print("Starting herbicide simulations (dose-response)...")
# Run repetitions for treatment
setwd('Model-files')
# Construct and run the command
for(MC in 1:nMC){
    print(paste("Control simulation, MC run:", MC, "over", nMC))
    mycall <- paste('./IBCgrassGUI', ModelVersion, CellNum, Tmax, InitDuration,
                    NamePftFile, SeedInput, belowres, abres, abampl, tramp, graz, cut,
                    week_start, HerbDuration, HerbEffectType, EffectModel, scenario, MC, sep=" ")
    system(mycall, intern=TRUE)
}
setwd('..')

print("Herbicide simulations finished.")
print(paste("Output files (Pt_*.txt, Grd_*.txt) are located in the Model-files/ directory."))
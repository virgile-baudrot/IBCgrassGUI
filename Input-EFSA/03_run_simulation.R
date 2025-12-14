source("Input-EFSA/02_simulation_setting.R")

dfBE2TH <- read.csv("Input-EFSA/BE2TH_cleaned.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

PATH = "Model-files/"
grBE2TH = c("STE_I1", "BLS_I1", "ALP_I1", "ATL_I1", "BOR_I1", "CON_I1", "PAN_I1", "MED_I1", "MAC_I1")
STRESSOR_TYPE = c("low", "high", "medium")
STRESSORS = c(no="NO", pb="biomass", sb="SEbiomass", sv="survival", es="establishment", st="sterility", ns="seednumber")

############################## LOAD ARGUMENTS
# group = "STE_I1"
# stressor_type = "low"
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
CellNum <- 173
Tmax <- if(CONTROL) 50 else 100
InitDuration <- 35
NamePftFile <- "Fieldedge.txt"
SeedInput <- 10
belowres <- 90
abres <- 100
abampl <- 0
tramp <- 0.1
graz <- 0
cut <- 1
week_start <- 1
HerbDuration <- 30
HerbEffectType <- if(CONTROL) 0 else 1
EffectModel <- 2
scenario <- idScenario
MC <- 0
nMC <- 10


# SAVE CONFIGURATION AND COPY INPUT FILES TO OUTPUT FOLDER
ID_OUT = formatC(sample(1:1e8,1), width = 9, format = "d", flag = "0")
PATH_OUTPUT = paste0("sim_", group,"_", stressor, "_", stressor_type, "_", ID_OUT, "/")
dir.create(PATH_OUTPUT)

# Create Fieldedge.txt for EFSA simulations
df_fieldedge <- build_Fieldedge_BE2TH(
    dfBE2TH, group = group,
    stressor=stressor,
    stressor_type=stressor_type)
write.table(df_fieldedge, file=paste0(PATH, "Fieldedge.txt"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(df_fieldedge, file=paste0(PATH_OUTPUT, "Fieldedge.txt"), sep="\t", row.names=FALSE, quote=FALSE)

# Create AppRate.txt for control simulation
df_AppRate <- build_AppRate(application_rate=HerbEffectType, herbicide_duration=HerbDuration)
write.table(df_AppRate, paste0(PATH,"AppRate.txt"), col.names=FALSE, row.names=FALSE, sep="\t")
write.table(df_AppRate, paste0(PATH_OUTPUT,"AppRate.txt"), col.names=FALSE, row.names=FALSE, sep="\t")

save_configuration(filepath = paste0(PATH_OUTPUT, "simulation_config"),
  ModelVersion, CellNum, Tmax, InitDuration, NamePftFile,
  SeedInput, belowres, abres, abampl, tramp, graz, cut,
  week_start, HerbDuration, HerbEffectType,
  EffectModel, scenario, MC, nMC
)

## Setup for parallel processing
# no_cores <- max(detectCores()-2, 1)
# cl <- makeCluster(no_cores)
# registerDoParallel(cl)

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
# foreach(MC = 1:nMC)  %dopar% {
#    print(paste("Control simulation, MC run:", MC, "over", nMC))
#    mycall <- paste('./IBCgrassGUI', ModelVersion, CellNum, Tmax, InitDuration,
#                  NamePftFile, SeedInput, belowres, abres, abampl, tramp, graz, cut,
#                  week_start, HerbDuration, HerbEffectType, EffectModel, scenario, MC, sep=" ")
#    system(mycall, intern=TRUE)
# }
# stopCluster(cl)
setwd('..')

# Clean up the temporary PFT file for the run
Grd_path = list.files(path = PATH, pattern = "^Grd_")
Pt_path = list.files(path = PATH, pattern = "^Pt_")
for(f in Grd_path){
  if(file.exists(f)){
    file.copy(paste0(PATH, Grd_path), paste0(PATH_OUTPUT, Grd_path))
    file.remove(paste0(PATH, Grd_path))
  }
}
for(f in Pt_path){
  if(file.exists(f)){
    file.copy(paste0(PATH, f), paste0(PATH_OUTPUT, f))
    file.remove(paste0(PATH, f))
  }
}
file.remove(paste0(PATH, "Fieldedge.txt"))
file.remove(paste0(PATH, "AppRate.txt"))

print("Herbicide simulations finished.")
print(paste("Output files (Pt_*.txt, Grd_*.txt) are located in the Model-files/ directory."))

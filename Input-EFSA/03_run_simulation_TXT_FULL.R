source("Input-EFSA/02_simulation_setting.R")

dfBE2TH <- read.csv("Input-EFSA/BE2TH_cleaned.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

PATH = "Model-files/"

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
Tmax <- 50 # if(CONTROL) 9 or 50 else 100
InitDuration <- 35 # 35 or 3
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
df_fieldedge_INIT <- build_Fieldedge_BE2TH(
    dfBE2TH, group = group,
    stressor=stressor,
    stressor_type=stressor_type)
#  -------------------------------------------------
# REDUCE AT 20 SPECIES
# --------------------------------------------------
MC = 1
IT = 1
while(MC<6 & IT<21){
  IT=IT+1
#for(MC in 1:nMC){
  df_fieldedge <- df_fieldedge_INIT
  # --------------------------------------------------
  write.table(df_fieldedge, file=paste0(PATH, "Fieldedge.txt"), sep="\t", row.names=FALSE, quote=FALSE)
  
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
  
  print("Starting herbicide simulations (dose-response)...")
  # SET IN GOOD FILE
  setwd('Model-files')
  # Construct and run the command
    print(paste("group", group, "stressor", stressor, "stressor_type", stressor_type))
    print(paste("MC run:", MC, "over", nMC))
    # 600s (10 minutes)
    mycall <- paste('timeout 600 ./IBCgrassGUI', ModelVersion, CellNum, Tmax, InitDuration,
                    NamePftFile, SeedInput, belowres, abres, abampl, tramp, graz, cut,
                    week_start, HerbDuration, HerbEffectType, EffectModel, scenario, MC, sep=" ")
    # RUN
    res <- system(mycall, intern=TRUE, ignore.stderr = TRUE)
    # CHECK
    exit_status <- attr(res, "status")
    if (!is.null(exit_status) && exit_status == 124) {
      print(paste("Run", MC, "stop after 60 seconds."))
    } else{
      MC=MC+1
    }
  setwd('..')
}


print("Herbicide simulations finished.")
print(paste("Output files (Pt_*.txt, Grd_*.txt) are located in the Model-files/ directory."))
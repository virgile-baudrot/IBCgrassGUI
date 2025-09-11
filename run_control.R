#####
# libraries
#####
library(foreach)
library(doParallel)
library(labeling)

#####
# Notes
#####
# -------------- Make sure to compile IBC before running this code ----------- #
# This script runs the CONTROL simulation (no herbicide effects).

#####
# Load previously saved simulation settings 
#####
# You can use the GUI to generate the SimulationSettings.Rdata,
# or create it manually.
# This script assumes the required input files are in the correct directories.
load("ExampleAnalyses/DoseResponse/HerbicideSettings/SimulationSettings.Rdata") 

#####
# Define Parameters
#####
ModelVersion <- 3
PFTfileName <- get("IBCcommunity", envir=SaveEnvironment)
PFTHerbEffectFile <- "./HerbFact.txt"
AppRateFile <- "./AppRate.txt"
MCruns <- get("IBCrepetition", envir=SaveEnvironment)
GridSize <- get("IBCgridsize", envir=SaveEnvironment)
SeedInput <- get("IBCSeedInput", envir=SaveEnvironment)
belowres <- get("IBCbelres", envir=SaveEnvironment)
abres <- get("IBCabres", envir=SaveEnvironment)
abampl <- get("IBCabampl", envir=SaveEnvironment)
graz <- get("IBCgraz", envir=SaveEnvironment)
tramp <- get("IBCtramp", envir=SaveEnvironment)
cut <- get("IBCcut", envir=SaveEnvironment)
week_start <- get("IBCweekstart", envir=SaveEnvironment)-10
HerbDuration <- get("IBCDuration", envir=SaveEnvironment)
RecovDuration <- get("IBCRecovery", envir=SaveEnvironment)
InitDuration <- get("IBCInit", envir=SaveEnvironment)
Tmax <- InitDuration + HerbDuration + RecovDuration
HerbEff <- get("IBCherbeffect", envir=SaveEnvironment)
if(HerbEff=="txt-file") EffectModel <- 0
if(HerbEff=="dose-response") EffectModel <- 2
Scenarios <- as.numeric(get("IBCScenarios", envir=SaveEnvironment))
nb_data <- as.numeric(get("nb_data", envir=SaveEnvironment))

#####
# Running control simulation
#####
scenario <- 0 # for control runs
# Prepare files
path <- "Model-files/"
write.table(get("IBCcommunityFile", envir=SaveEnvironment), paste(path,PFTfileName, sep=""), sep="\t", quote=F,row.names=F)
file.copy("Input-files/HerbFact.txt", path, overwrite = TRUE)
file.copy("Input-files/AppRate.txt", path, overwrite = TRUE)

# Change directory to run the model
setwd('Model-files')

# Setup for parallel processing
no_cores <- max(detectCores()-2,1)
cl <- makeCluster(no_cores)
registerDoParallel(cl)

# Start control simulations
print("Starting control simulations...")
foreach(MC = 1:MCruns)  %dopar% {
  system(paste('./IBCgrassGUI', ModelVersion, GridSize, Tmax, InitDuration, PFTfileName, SeedInput, belowres, abres, abampl, tramp, graz, cut,
               week_start, HerbDuration, 0, EffectModel, scenario, MC, sep=" "), intern=T)
}
stopCluster(cl)
print("Control simulations finished.")

# Go back to the parent directory
setwd('..')

# Note: The original script moved output files. In this version, output files
# will remain in the 'Model-files' directory. You can manage them as needed.
print(paste("Output files (Pt_*.txt, Grd_*.txt) are located in the", path, "directory."))

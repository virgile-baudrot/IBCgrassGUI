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
# This script runs the HERBICIDE simulation.
# It uses the "dose-response" method as an example.

#####
# Load previously saved simulation settings 
#####
load("ExampleAnalyses/DoseResponse/HerbicideSettings/SimulationSettings.Rdata")

#####
# Define Parameters
#####
ModelVersion <- 3
PFTfileName <- get("IBCcommunity", envir=SaveEnvironment)
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
HerbEff <- "dose-response" # Hardcoded for this script
EffectModel <- 2
Scenarios <- as.numeric(get("IBCScenarios", envir=SaveEnvironment))
nb_data <- as.numeric(get("nb_data", envir=SaveEnvironment))

#####
# Running treatment based on dose responses
#####
PFTfile <- get("IBCcommunityFile", envir=SaveEnvironment)
AppRateScenarios <- data.frame(get("IBCAppRateScenarios", envir=SaveEnvironment))
PFTsensitivity <- get("PFTSensitivityFile", envir=SaveEnvironment)

# Setup for parallel processing
no_cores <- max(detectCores()-2,1)
cl <- makeCluster(no_cores)
registerDoParallel(cl)

path_DR <- "ExampleAnalyses/DoseResponse/HerbicideSettings/"

print("Starting herbicide simulations (dose-response)...")
# Run repetitions for treatment
foreach(MC = 1:MCruns, .export=c("PFTfile", "PFTsensitivity", "PFTfileName", "EffectModel",
                                 "ModelVersion", "belowres", "abres", "abampl", "Tmax", "InitDuration", "GridSize", "SeedInput",
                                 "week_start", "HerbDuration", "tramp", "graz", "cut", "Scenarios", "nb_data",
                                 "AppRateScenarios"))  %dopar% {
  
  # This section generates a PFT file with specific sensitivities for each repetition.
  # It reads EC50 and slope values from external files.
  # Ensure these files are present in the root directory when running the script.
  # (e.g., EC50andslope_Biomass.txt, EC50andslope_Survival.txt, etc.)
  
  PFTfile_run <- merge(PFTfile, PFTsensitivity, by="Species")
  PFTfile_run[,28:39] <- 0 # Clear previous values
  
  # Example for Biomass - this logic is repeated for other effects in the original script
  if("EC50andslope_Biomass.txt" %in% list.files(path_DR)){
    DR <- read.table(paste0(path_DR, "EC50andslope_Biomass.txt"), sep="\t", header=T)
    # EC50
    PFTfile_run[PFTfile_run$Sensitivity=="random",28] <- runif(nrow(PFTfile_run[PFTfile_run$Sensitivity=="random",]), min=max(0,(DR[nb_data+1,2]-DR[nb_data+2,2])), max=max(0,(DR[nb_data+1,2]+DR[nb_data+2,2])))
    # slope
    PFTfile_run[PFTfile_run$Sensitivity=="random",29] <- runif(nrow(PFTfile_run[PFTfile_run$Sensitivity=="random",]), min=max(0,(DR[nb_data+1,3]-DR[nb_data+2,3])), max=max(0,(DR[nb_data+1,3]+DR[nb_data+2,3])))
    # ... (add more logic for other sensitivities if needed) ...
  }

  # Save the generated PFT file for this run
  run_pft_filename <- paste0(unlist(strsplit(PFTfileName,".txt")), MC, ".txt")
  write.table(PFTfile_run, file=run_pft_filename, row.names=FALSE, quote=FALSE, sep="\t")

  # Prepare files for C++ executable
  path <- "Model-files/"
  file.copy(run_pft_filename, path, overwrite = TRUE)
  file.copy("Input-files/HerbFact.txt", path, overwrite = TRUE) # Dummy file needed by C++ code

  setwd('Model-files')

  for(scenario in 1:Scenarios){
    # Create AppRate.txt for the current scenario
    write.table(AppRateScenarios[,scenario], "AppRate.txt", col.names=FALSE, row.names=FALSE, sep="\t")

    # Construct and run the command
    mycall <- paste('./IBCgrassGUI', ModelVersion, GridSize, Tmax, InitDuration, 
                    run_pft_filename, SeedInput, belowres, abres, abampl, tramp, graz, cut, 
                    week_start, HerbDuration, 1, EffectModel, scenario, MC, sep=" ")
    system(mycall, intern=TRUE)
  }
  
  setwd('..')
  # Clean up the temporary PFT file for the run
  file.remove(run_pft_filename)
}

stopCluster(cl)
print("Herbicide simulations finished.")
print(paste("Output files (Pt_*.txt, Grd_*.txt) are located in the Model-files/ directory."))

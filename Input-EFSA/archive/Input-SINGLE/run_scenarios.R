# run_scenarios.R

args <- commandArgs(trailingOnly = TRUE) 
if (length(args) == 0) {
  stop("Erreur : Aucun numéro n'a été fourni en argument !", call. = FALSE)
}
n_scen <- as.numeric(args[1])
n_size <- as.numeric(args[2])

library(dplyr)
source("simulation_setting.R")
source("run_ibcgrass.R")

# 1. Load the species data and isolate a single species
species_data <- read.csv("BE2TH_cleaned.csv", stringsAsFactors = FALSE)
single_species <- species_data[1, , drop=FALSE] # Preserve dataframe format

# Define the scenarios and their corresponding 50% effect modifiers
p_effect=0.9
if(n_scen == 100){
  scenarios <- list(
    "control"       = list(effect = "none", value = 0.0),
    "biomass"       = list(effect = "biomass", value = p_effect),
    "SEbiomass"     = list(effect = "SEbiomass", value = p_effect),
    "survival"      = list(effect = "survival", value = p_effect),
    "establishment" = list(effect = "establishment", value = p_effect),
    "sterility"     = list(effect = "sterility", value = p_effect),
    "seednumber"    = list(effect = "seednumber", value = p_effect)
  )
}
n_rep=30
if(n_scen == 200){
  scenarios <- setNames(
      replicate(n_rep, list(effect = "none", value = 0.0), simplify = FALSE),
      paste0("control", "_", 1:n_rep)
    )
}
if(n_scen == 300){
  name_var = "survival"
  scenarios <- setNames(
    c(replicate(n_rep, list(effect = "none", value = 0.0), simplify = FALSE),
      replicate(n_rep, list(effect = name_var, value = p_effect), simplify = FALSE)
    ),
    c(paste0("control", "_", 1:n_rep), paste0(name_var, "_", 1:n_rep))
  )
}

small_landscape <- list(
  ModelVersion = 3,
  CellNum = n_size,      # 900 = 30x30 (dilue l'effet des survivants isolés)
  Tmax = 40,             # 40 ans (1200 semaines)
  InitDuration = 20,     # 20 ans (600 semaines) pour stabiliser la population
  SeedInput = 10,
  belowres = 90,
  abres = 100,
  abampl = 0,
  tramp = 0.1,
  graz = 0,
  cut = 1,
  week_start = 1,
  HerbDuration = 10      # Stress appliqué pendant 10 ans (épuise la banque de graines)
)

# 2. Loop through scenarios, apply effects, and run the model
scenario_id <- 1 + n_scen
MC=1
for (scenario_name in names(scenarios)) {
  cat("\nRunning scenario:", scenario_name, "\n")
  
  current_setup <- scenarios[[scenario_name]]
  is_control <- (current_setup$effect == "none")
  
  # Run the IBC-grass model
  run_ibcgrass(
    species_data = single_species,
    landscape_config = small_landscape,
    scenario_setup = current_setup,
    scenario_id = scenario_id,
    MC = MC,
    is_control = is_control
  )
  if(n_scen != 100){
    MC <- MC+1
  }
  
  # Save the configuration files used for this specific scenario
  # Copy HerbFact.txt and Fieldedge.txt from Model-files and rename them
  if (file.exists("../Model-files/HerbFact.txt")) {
    file.copy(from = "../Model-files/HerbFact.txt", 
              to = paste0("HerbFact_", scenario_name, ".txt"), 
              overwrite = TRUE)
  }
  if (file.exists("../Model-files/Fieldedge.txt")) {
    file.copy(from = "../Model-files/Fieldedge.txt", 
              to = paste0("Fieldedge_", scenario_name, ".txt"), 
              overwrite = TRUE)
  }
  
  # 5. Extract output and save to the Input-SINGLE folder
  # Model generates Pt_scenarioID_MC.txt
  pt_file <- paste0("../Model-files/Pt_", scenario_id, "_1.txt")
  if (file.exists(pt_file)) {
    pt_data <- read.table(pt_file, header = TRUE, sep = "\t")
    pt_file <- paste0("pt_output_", scenario_name, ".csv")
    write.csv(pt_data, pt_file, row.names = FALSE)
  } else {
    warning(paste("Model output not found for scenario:", scenario_name))
  }
  grd_file <- paste0("../Model-files/Grd_", scenario_id, "_1.txt")
  if (file.exists(pt_file)) {
    grd_data <- read.table(grd_file, header = TRUE, sep = "\t")
    grd_file <- paste0("grd_output_", scenario_name, ".csv")
    write.csv(grd_data, grd_file, row.names = FALSE)
  } else {
    warning(paste("Model output not found for scenario:", scenario_name))
  }
  scenario_id <- scenario_id + 1
}

cat("\nAll scenarios completed successfully.\n")
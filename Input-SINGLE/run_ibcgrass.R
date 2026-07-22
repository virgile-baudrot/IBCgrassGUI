# run_ibcgrass.R

run_ibcgrass <- function(species_data, landscape_config, scenario_setup, scenario_id = 1, MC = 1, is_control = TRUE) {
  # 1. Define paths assuming Input-SINGLE and Model-files are siblings
  path <- "../Model-files/"
  pft_filename <- "Fieldedge.txt"
  
  # 2. Write the species data (Fieldedge.txt)
  df_fieldedge <- build_Fieldedge_Single(species_data)
  write.table(df_fieldedge, file = paste0(path, pft_filename), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # 3. Write the HerbFact.txt file for the specific effect
  df_herbfact <- build_HerbFact_fixed(landscape_config$HerbDuration, scenario_setup$effect, scenario_setup$value)
  write.table(df_herbfact, file = paste0(path, "HerbFact.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)
  
  # 4. Determine Herbicide Effect and Model types
  # 0 forces the TXT way. HerbEffectType 0 = control, 1 = effect
  HerbEffectType <- ifelse(is_control, 0, 1)
  EffectModel    <- 0 
  
  # 5. Construct the system call string
  mycall <- paste('timeout 1800 ./IBCgrassGUI', 
                  landscape_config$ModelVersion, landscape_config$CellNum, 
                  landscape_config$Tmax, landscape_config$InitDuration,
                  pft_filename, landscape_config$SeedInput, 
                  landscape_config$belowres, landscape_config$abres, 
                  landscape_config$abampl, landscape_config$tramp, 
                  landscape_config$graz, landscape_config$cut,
                  landscape_config$week_start, landscape_config$HerbDuration, 
                  HerbEffectType, EffectModel, scenario_id, MC, sep = " ")
  
  # 6. Execute the model in the correct directory
  setwd(path)
  cat("Executing system call:\n", mycall, "\n")
  res <- system(mycall, intern = TRUE, ignore.stderr = TRUE)
  
  exit_status <- attr(res, "status")
  if (!is.null(exit_status) && exit_status == 124) {
    warning(paste("Run", MC, "stopped after 180 seconds due to timeout."))
  }
  
  setwd('../Input-SINGLE')
  return(res)
}
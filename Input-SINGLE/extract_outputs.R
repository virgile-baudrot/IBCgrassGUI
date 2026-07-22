# extract_outputs.R

# 1. Récupération des arguments pour identifier le set de scénarios
args <- commandArgs(trailingOnly = TRUE) 
if (length(args) == 0) {
  stop("Erreur : Aucun numéro n'a été fourni en argument !", call. = FALSE)
}
n_scen <- as.numeric(args[1])

# Installation et chargement des librairies
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(tidyr)) install.packages("tidyr")
if (!require(dplyr)) install.packages("dplyr")

library(ggplot2)
library(tidyr)
library(dplyr)

# 2. Définition des plages d'IDs attendus selon la série (100, 200, 300)
# Plus besoin de définir n_rep : on ratisse large en se basant sur le fait que n_rep < 100.
if(n_scen == 100){
  expected_ids <- as.character(101:199)
} else if(n_scen == 200){
  expected_ids <- as.character(201:299)
} else if(n_scen == 300){
  # Pour 300, on a 2 * n_rep (contrôles + traitements). On va donc jusqu'à 499.
  expected_ids <- as.character(301:499)
} else {
  stop("Erreur : n_scen non reconnu.")
}

# ==========================================
# Fonction générique de traitement des outputs bruts (Pt et Grd)
# ==========================================

process_simulation_data <- function(file_prefix, expected_ids, output_rds_path) {
  all_data <- data.frame()
  
  # A. Scanner le dossier Model-files dynamiquement
  pattern <- paste0("^", file_prefix, "_type_.*_Scenario_\\d+_MCrun_\\d+\\.txt$")
  files <- list.files("../Model-files", pattern = pattern, full.names = TRUE)
  
  if(length(files) == 0) {
    cat("Aucun fichier brut trouvé dans Model-files/ pour le préfixe :", file_prefix, "\n")
    return(NULL)
  }
  
  # B. Lire et compiler les fichiers
  fichiers_lus <- 0
  for (file_name in files) {
    bname <- basename(file_name)
    scen_id <- sub(".*_Scenario_(\\d+)_.*", "\\1", bname)
    
    # On ne lit le fichier QUE s'il appartient à la plage demandée
    if(scen_id %in% expected_ids) {
      temp_data <- read.table(file_name, header = TRUE, sep = "\t")
      
      scen_type <- sub(".*_type_([a-zA-Z]+)_.*", "\\1", bname)
      temp_data$Scenario <- scen_id
      temp_data$Scenario_Type <- scen_type
      
      if("Time" %in% colnames(temp_data)) names(temp_data)[names(temp_data) == "Time"] <- "time"
      all_data <- rbind(all_data, temp_data)
      fichiers_lus <- fichiers_lus + 1
    }
  }
  
  if(nrow(all_data) == 0) {
    cat("Aucune donnée correspondant à la série demandée trouvée pour :", file_prefix, "\n")
    return(NULL)
  }
  
  # C. Passage au format Long
  long_data <- all_data %>%
    pivot_longer(
      cols = -any_of(c("time", "Scenario", "Scenario_Type", "PFT", "X", "Y", "x", "y")), 
      names_to = "OutputVariable",
      values_to = "Value"
    )
  
  # D. Calcul des métriques relatives au contrôle
  control_ids <- unique(long_data$Scenario[long_data$Scenario_Type == "control"])
  
  if (length(control_ids) > 0) {
    grouping_cols <- c("time", "OutputVariable")
    if("PFT" %in% colnames(long_data)) grouping_cols <- c(grouping_cols, "PFT")
    if("X" %in% colnames(long_data)) grouping_cols <- c(grouping_cols, "X")
    if("Y" %in% colnames(long_data)) grouping_cols <- c(grouping_cols, "Y")
    
    control_data <- long_data %>%
      filter(Scenario %in% control_ids) %>%
      group_by(across(all_of(grouping_cols))) %>%
      summarise(Control_Value = mean(Value, na.rm = TRUE), .groups = 'drop')
    
    long_data <- long_data %>%
      left_join(control_data, by = grouping_cols) %>%
      mutate(
        Ratio_to_Control = Value / Control_Value,
        Diff_from_Control = Control_Value - Value
      )
  } else {
    cat("\nNote : Aucun contrôle trouvé pour", file_prefix, ". Métriques relatives ignorées.\n")
  }
  
  # E. Sauvegarde
  saveRDS(long_data, output_rds_path)
  cat("=> Fichier combiné sauvegardé :", output_rds_path, "(", fichiers_lus, "fichiers sourcés )\n")
  
  return(long_data)
}

# ==========================================
# Exécution pour les deux types de fichiers
# ==========================================

cat("\n--- Traitement des données de Population (Pt_) ---\n")
pt_long_data <- process_simulation_data("Pt_", expected_ids, "../output-SINGLE/pt_long_data.rds")

cat("\n--- Traitement des données Spatiales (Grd_) ---\n")
grd_long_data <- process_simulation_data("Grd_", expected_ids, "../output-SINGLE/grd_long_data.rds")

cat("\nScript d'extraction terminé avec succès.\n")
source("Input-EFSA/02_simulation_setting.R")

dfBE2TH <- read.csv("Input-EFSA/BE2TH_cleaned.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

PATH = "Model-files/"

############################## LOAD ARGUMENTS
args <- commandArgs(trailingOnly = TRUE)

print(paste("group :", args[1]))
print(paste("stressor_type :", args[2]))
group =  args[1]
stressor_type =  args[2]
#############################################

# Définition de la liste complète des stresseurs
all_stressors <- c("biomass", "SEbiomass", "survival", "establishment", "sterility", "seednumber")
run_stressors <- c("NO", all_stressors)

ModelVersion <- 3
CellNum <- 173 
Tmax <- 50 
InitDuration <- 20 
NamePftFile <- "Fieldedge.txt"
SeedInput <- 10
belowres <- 90
abres <- 100
abampl <- 0
tramp <- 0.1
graz <- 0
cut <- 0 
week_start <- 1
HerbDuration <- 10
HerbEffectType <- 1
EffectModel <- 0
nMC <- 10

# INITIALISATION : On crée la base avec TOUS les stresseurs pour que build_dose_response() 
# assigne la même sensibilité (rxp) dans toutes les colonnes EC50_ de chaque espèce.
df_fieldedge_INIT <- build_Fieldedge_BE2TH(
  dfBE2TH,
  group = group,
  stressor = all_stressors,
  stressor_type = stressor_type)

#  -------------------------------------------------
# BOUCLE DE MONTE CARLO
# --------------------------------------------------
MC = 1
IT = 1

while(MC <= nMC & IT <= 21){
  IT = IT + 1
  
  # Sélection de la communauté de 20 espèces (fixe pour tous les stresseurs de ce MC)
  nselect = sample(1:nrow(df_fieldedge_INIT), 20)
  df_fieldedge_base <- df_fieldedge_INIT[nselect, ]
  
  community_success <- TRUE
  
  # On boucle sur chaque stresseur (Control + 6 effets) pour cette même communauté
  for(current_stressor in run_stressors){
    
    df_fieldedge <- df_fieldedge_base
    
    # On met à 0 les EC50 des autres stresseurs pour isoler l'effet actuel
    for(s in all_stressors){
      if(s != current_stressor){
        df_fieldedge[[paste0("EC50_", s)]] <- 0
      }
    }
    
    write.table(df_fieldedge, file=paste0(PATH, "Fieldedge.txt"), sep="\t", row.names=FALSE, quote=FALSE)
    
    # Création du fichier HerbFact spécifique
    df_HerbFact <- build_HerbFact(HerbDuration=HerbDuration, stressor=current_stressor, stressor_type=stressor_type)
    write.table(df_HerbFact, paste0(PATH,"HerbFact.txt"), col.names=TRUE, row.names=FALSE, sep="\t")
    
    idScenario = build_id_scenario(group, stressor_type, current_stressor)
    
    # On sauvegarde la config en incluant le nom du stresseur pour ne pas écraser les fichiers
    save_configuration(filepath = paste0(PATH, "simulation_config_", current_stressor, "_", MC),
                       ModelVersion, CellNum, Tmax, InitDuration, NamePftFile,
                       SeedInput, belowres, abres, abampl, tramp, graz, cut,
                       week_start, HerbDuration, HerbEffectType,
                       EffectModel, idScenario, MC, nMC)
    
    # EXECUTION DU MODELE
    setwd('Model-files')
    print(paste("Running -> group:", group, "| type:", stressor_type, "| stressor:", current_stressor, "| MC:", MC))
    
    mycall <- paste('timeout 600 ./IBCgrassGUI', ModelVersion, CellNum, Tmax, InitDuration,
                    NamePftFile, SeedInput, belowres, abres, abampl, tramp, graz, cut,
                    week_start, HerbDuration, HerbEffectType, EffectModel, idScenario, MC, sep=" ")
    
    time_taken <- system.time({
      res <- system(mycall, intern=TRUE, ignore.stderr = TRUE)
    })
    exit_status <- attr(res, "status")
    print(paste("Temps d'exécution :", round(time_taken["elapsed"], 2), "secondes"))
    setwd('..')
    
    # GESTION DU BLOCAGE
    # Si UN SEUL stresseur plante pour cette communauté, on l'abandonne entièrement pour garder des données équilibrées
    if (!is.null(exit_status) && exit_status == 124) {
      print(paste("Run", MC, "for stressor", current_stressor, "stopped after 600 seconds. Discarding community."))
      community_success <- FALSE
      break 
    }
  }
  
  # Si la communauté a passé les 7 simulations avec succès, on valide le Run MC
  if(community_success){
    MC = MC + 1
  }
}

print("Herbicide simulations finished for all stressors on generated communities.")
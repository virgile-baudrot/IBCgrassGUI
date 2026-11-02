# 00_test_performance_parallel.R

# 1. CHARGEMENT DES LIBRAIRIES
# (Décommentez la ligne suivante si les packages ne sont pas installés)
# install.packages(c("foreach", "doParallel"))

library(foreach)
library(doParallel)

source("Input-EFSA/02_simulation_setting.R")
dfBE2TH <- read.csv("Input-EFSA/BE2TH_cleaned.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

# --------------------------------------------------
# 2. CONFIGURATION DU PARALLELISME
# --------------------------------------------------
n_trials <- 5 # Nombre d'essais par configuration
num_cores <- min(detectCores() - 1, n_trials) # Utilise au max 5 coeurs
registerDoParallel(cores = num_cores)

print(paste("Démarrage du cluster parallèle avec", num_cores, "coeurs..."))

# --------------------------------------------------
# 3. DEFINITION DES PARAMETRES (Mode Contrôle)
# --------------------------------------------------
group <- "STE_I1"
stressor_type <- "low"
stressor <- "NO"

ModelVersion <- 3
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
scenario <- build_id_scenario(group, stressor_type, stressor)
MC <- 1
nMC <- 1

# Création des données communes pour tous les essais
df_fieldedge_INIT <- build_Fieldedge_BE2TH(dfBE2TH, group = group, stressor = stressor, stressor_type = stressor_type)
nselect <- sample(1:nrow(df_fieldedge_INIT), 20)
df_fieldedge <- df_fieldedge_INIT[nselect, ]

df_HerbFact <- build_HerbFact(HerbDuration=HerbDuration, stressor=stressor, stressor_type=stressor_type)

# --------------------------------------------------
# 4. FONCTION D'EVALUATION PARALLELE
# --------------------------------------------------
evaluate_performance <- function(cell_num, trials) {
  
  print(paste("=> Début de l'évaluation pour CellNum =", cell_num, "(", trials, "essais en parallèle )"))
  
  # La boucle foreach distribue les essais sur les coeurs disponibles
  times <- foreach(i = 1:trials, .combine = c) %dopar% {
    
    # Création d'un environnement de travail isolé pour éviter les conflits d'écriture
    trial_dir <- paste0("temp_eval_cell_", cell_num, "_trial_", i)
    model_dir <- paste0(trial_dir, "/Model-files")
    
    dir.create(trial_dir, showWarnings = FALSE)
    dir.create(model_dir, showWarnings = FALSE)
    
    # Copie de l'exécutable
    file.copy("Model-files/IBCgrassGUI", model_dir, overwrite = TRUE)
    # Rendre le fichier exécutable au cas où les permissions sautent pendant la copie
    system(paste("chmod +x", paste0(model_dir, "/IBCgrassGUI")))
    
    # Ecriture des fichiers d'entrée spécifiques à cet essai
    write.table(df_fieldedge, file=paste0(model_dir, "/Fieldedge.txt"), sep="\t", row.names=FALSE, quote=FALSE)
    write.table(df_HerbFact, paste0(model_dir, "/HerbFact.txt"), col.names=TRUE, row.names=FALSE, sep="\t")
    
    # Commande bash : on se déplace dans le sous-dossier isolé PUIS on lance le modèle
    mycall <- paste('cd', model_dir, '&& ./IBCgrassGUI', ModelVersion, cell_num, Tmax, InitDuration,
                    NamePftFile, SeedInput, belowres, abres, abampl, tramp, graz, cut,
                    week_start, HerbDuration, HerbEffectType, EffectModel, scenario, MC)
    
    # Mesure du temps
    t <- system.time({
      system(mycall, intern=TRUE, ignore.stderr = TRUE)
    })["elapsed"]
    
    # Nettoyage : suppression du dossier temporaire pour libérer l'espace disque
    unlink(trial_dir, recursive = TRUE)
    
    return(as.numeric(t))
  }
  
  # Calcul des statistiques
  res_summary <- data.frame(
    CellNum = cell_num,
    Min_sec = round(min(times), 2),
    Max_sec = round(max(times), 2),
    Mean_sec = round(mean(times), 2)
  )
  
  return(res_summary)
}

# --------------------------------------------------
# 5. LANCES DES TESTS
# --------------------------------------------------
print("========================================")

# Lancement pour les 3 configurations
res_173 <- evaluate_performance(173, n_trials)
res_50  <- evaluate_performance(50, n_trials)
res_10  <- evaluate_performance(10, n_trials)

print("========================================")
print("BILAN DES PERFORMANCES (en secondes)")
print("========================================")

# Compilation des résultats dans un tableau final
final_results <- rbind(res_173, res_50, res_10)
print(final_results)

# Arrêt propre du cluster parallèle
stopImplicitCluster()
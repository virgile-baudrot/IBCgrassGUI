# 1. Installation des packages si vous ne les avez pas déjà
# install.packages(c("tidyverse", "data.table"))
# install.packages("tidyverse")
# install.packages("data.table")

library(tidyverse) # Contient stringr, dplyr, purrr, etc.
library(data.table) # Ultra-rapide pour lire et manipuler de gros fichiers texte

# 2. Définir le chemin racine from IBCgrass
dossier_racine <- "/home/IBCgrass/runs_effect_r20all"

Pt_fichiers <- list.files(
  path = dossier_racine, 
  pattern = "^Pt__.*\\.txt$", 
  recursive = TRUE, 
  full.names = TRUE
)
Grd_fichiers <- list.files(
  path = dossier_racine, 
  pattern = "^Grd__.*\\.txt$", 
  recursive = TRUE, 
  full.names = TRUE
)

# 4. Créer la fonction qui lit un fichier et ajoute les 5 colonnes
lire_et_annoter <- function(chemin_fichier) {
  # 1. SECURITE : Si le fichier est vide (0 octet), on l'ignore silencieusement
  if (file.info(chemin_fichier)$size == 0) {
    return(NULL) # map_dfr va juste passer au suivant sans erreur
  }
  nom_fichier <- basename(chemin_fichier)
  # --- EXTRACTION (Mise à jour pour fonctionner avec Pt__ ET Grd__) ---
  biome <- str_extract(chemin_fichier, "[A-Z]{3}(?=_I1_)")
  stress_level <- str_extract(chemin_fichier, "(?<=_I1_)z\\d{2}")
  # On cherche juste après "__type_" pour que ça marche pour Pt__ et Grd__
  treatment_type <- str_extract(nom_fichier, "(?<=__type_)[a-zA-Z0-9]+(?=_Scenario)")
  stress_type <- str_extract(nom_fichier, "(?<=Scenario_)\\d{3}")
  repetition <- str_extract(nom_fichier, "(?<=MCrun_)\\d+(?=\\.txt)")
  # 2. LECTURE : suppressWarnings() rend fread silencieux face aux "footers" mal formés
  df <- suppressWarnings(fread(chemin_fichier))
  # 3. SECURITE bis : Si fread n'a rien pu lire du tout malgré tout
  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }
  # 4. AJOUT DES COLONNES
  df <- df %>%
    mutate(
      Biome = biome,
      stress_level = stress_level,
      treatment_type = treatment_type,
      stress_type = stress_type,
      repetition = as.integer(repetition)
    )
  return(df)
}
# 5. Exécution : on boucle sur tous les fichiers et on fusionne (merge) le tout
# map_dfr applique la fonction à chaque chemin et "bind_rows" (fusionne) les résultats en un seul tableau
df_Pt_final <- map_dfr(Pt_fichiers, lire_et_annoter)
df_Grd_final <- map_dfr(Grd_fichiers, lire_et_annoter)

# --- VERIFICATION ET SAUVEGARDE ---
head(df_Pt_final)
head(df_Grd_final)
# Optionnel : Sauvegarder le super-dataframe final dans un fichier unique (sans toucher aux originaux)
# fwrite(df_Pt_final, "Pt__all_16.csv")
# fwrite(df_Grd_final, "Grd__all_16.csv")


# fwrite(df_Pt_final, "/stockage/Pt__all_32.csv")
# fwrite(df_Grd_final, "/stockage/Grd__all_32.csv")
# Version recommandée avec compression
fwrite(df_Pt_final, "/stockage/Pt__all_32.csv.gz")
fwrite(df_Grd_final, "/stockage/Grd__all_32.csv.gz")

fwrite(df_Pt_final, "/stockage/Pt__all_16.csv.gz")
fwrite(df_Grd_final, "/stockage/Grd__all_16.csv.gz")


###############################################################################
#
# CLEAN FILES
library(tidyverse)
library(data.table)
df_Pt1 <- read_csv("/stockage/out/Pt__all_32.csv.gz")
df_Pt2 <- read_csv("/stockage/out/Pt__all_16.csv.gz")
df_Pt_final <- bind_rows(df_Pt1, df_Pt2)
# df_Pt_final <- rbindlist(list(df_Pt1, df_Pt2))
fwrite(df_Pt_final, "/stockage/Pt__all_full.csv.gz")


df_Pt_final <- read_csv("/stockage/Pt__all_full.csv.gz")

df_Pt_final_noPFT <- df_Pt_final
df_Pt_final_noPFT$PFT = NULL
# 1. Convertir en data.table sans copier les données (instantané, 0 RAM supplémentaire)
setDT(df_Pt_final_noPFT)
# 2. Définir les variables
cols_to_summarise <- c("Inds", "seedlings", "seeds", "repromass", "cover", "shootmass")
group_vars <- c("PFT", "Biome", "stress_level", "Time", "treatment_type", "stress_type")
group_vars <- c("Biome", "stress_level", "Time", "treatment_type", "stress_type")
# 3. Créer une liste vide pour stocker les résultats temporaires
liste_resultats <- list()
# 4. Calculer le résumé colonne par colonne
for (var in cols_to_summarise) {
  # data.table calcule les stats uniquement sur la colonne ciblée (get(var))
  res_temp <- df_Pt_final_noPFT[, .(
    mean   = as.numeric(mean(get(var), na.rm = TRUE)),
    median = as.numeric(median(get(var), na.rm = TRUE)),
    qinf95 = as.numeric(quantile(get(var), 0.025, na.rm = TRUE)),
    qsup95 = as.numeric(quantile(get(var), 0.975, na.rm = TRUE)),
    min    = as.numeric(min(get(var), na.rm = TRUE)),
    max    = as.numeric(max(get(var), na.rm = TRUE))
  ), by = group_vars]
  # On ajoute la colonne "VarOutput" au petit tableau de résultat
  res_temp[, VarOutput := var]
  # On stocke ce petit résumé dans la liste
  liste_resultats[[var]] <- res_temp
}
# 5. Empiler les petits tableaux pour créer le résumé final
df_Pt_summary <- rbindlist(liste_resultats)
fwrite(df_Pt_summary, "/stockage/df_Pt_summary_noPFT.csv.gz")


df_Grd1 <- read_csv("/stockage/out/Grd__all_32.csv.gz")
df_Grd2 <- read_csv("/stockage/out/Grd__all_16.csv.gz")
df_Grd_final <- bind_rows(df_Grd1, df_Grd2)
fwrite(df_Grd_final, "/stockage/Grd__all_full.csv.gz")

df_Grd_final <- read_csv("/stockage/Grd__all_full.csv.gz")
setDT(df_Grd_final)
# 2. Définir les variables de groupe (pas de PFT ici)
group_vars_grd <- c("Biome", "stress_level", "Time", "treatment_type", "stress_type")
# 3. Identifier automatiquement les colonnes à résumer
# On prend toutes les colonnes, SAUF les groupes et la répétition
cols_to_summarise_grd <- setdiff(names(df_Grd_final), c(group_vars_grd, "repetition"))
# 4. Liste pour stocker les résultats
liste_resultats_grd <- list()
# 5. Boucle de calcul très économe en mémoire
for (var in cols_to_summarise_grd) {
  res_temp <- df_Grd_final[, .(
    mean   = as.numeric(mean(get(var), na.rm = TRUE)),
    median = as.numeric(median(get(var), na.rm = TRUE)),
    qinf95 = as.numeric(quantile(get(var), 0.025, na.rm = TRUE)),
    qsup95 = as.numeric(quantile(get(var), 0.975, na.rm = TRUE)),
    min    = as.numeric(min(get(var), na.rm = TRUE)),
    max    = as.numeric(max(get(var), na.rm = TRUE))
  ), by = group_vars_grd]
  # Ajout du nom de la variable analysée
  res_temp[, VarOutput := var]
  # Sauvegarde temporaire
  liste_resultats_grd[[var]] <- res_temp
}
# 6. Empiler les résultats
df_Grd_summary <- rbindlist(liste_resultats_grd)
# Optionnel : réorganiser les colonnes pour que VarOutput soit au début avec les groupes
setcolorder(df_Grd_summary, c(group_vars_grd, "VarOutput"))
fwrite(df_Grd_summary, "/stockage/df_Grd_summary.csv.gz")
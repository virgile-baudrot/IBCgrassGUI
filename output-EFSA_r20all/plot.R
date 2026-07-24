library(readr)
library(dplyr)
library(ggplot2)
library(zoo)

# Tableaux de correspondance inverses (décodage)
be2th_map <- c(
  "1" = "STE", "2" = "BLS", "3" = "ALP", "4" = "ATL", 
  "5" = "BOR", "6" = "CON", "7" = "PAN", "8" = "MED", "9" = "MAC"
)
stress_level_map <- c(
  "1" = "low", "2" = "high", "3" = "medium",
  "4" = "z01", "5" = "z02", "6" = "z03", "7" = "z04", "8" = "z05", "9" = "z00"
)
stressor_map <- c(
  "0" = "NO", "1" = "biomass", "2" = "SEbiomass",
  "3" = "survival", "4" = "establishment", "5" = "sterility", "6" = "seednumber"
)
plt_var <- function(df, VAR="shannon", YLAB, STRESS, TITLE=""){
  df_sel <- df |>
    filter(VarOutput==VAR, stressor==STRESS)
  ggplot(data=df_sel,
         aes(x=Time/30, y=median_relative,
             ymin=qinf95_relative, ymax=qsup95_relative,
             fill=period)) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    geom_ribbon(alpha=0.7) +
    geom_line() +
    scale_fill_manual(
      values = c("1_stabilize"="lightgreen", "2_treatment"="#DD8888", "3_recovery"="#4488DD"),
      labels = c(
        "1_stabilize" = "0-20yr (stabilize)", 
        "2_treatment" = "20-30yr (treatment)", 
        "3_recovery"  = "30-50yr (recovery)"
      ),
      name = "Phase"
    ) +
    labs(title=paste("Herbicide effect on", TITLE),
         x="Time (years)", y = paste("Ratio Effect/Control -", YLAB)) +
    facet_grid(Biome~stress_level, scale="free")
}
plt_var_roll <- function(df, VAR="shannon", YLAB, STRESS, TITLE=""){
  df_sel <- df |>
    filter(VarOutput==VAR, stressor==STRESS)
  ggplot(data=df_sel,
         aes(x=Time/30, y=median_relative_roll,
             ymin=qinf95_relative_roll, ymax=qsup95_relative_roll,
             fill=period)) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    geom_ribbon(alpha=0.7) +
    geom_line() +
    scale_fill_manual(
      values = c("1_stabilize"="lightgreen", "2_treatment"="#DD8888", "3_recovery"="#4488DD"),
      labels = c(
        "1_stabilize" = "0-20yr (stabilize)", 
        "2_treatment" = "20-30yr (treatment)", 
        "3_recovery"  = "30-50yr (recovery)"
      ),
      name = "Phase"
    ) +
    labs(title=paste("Herbicide effect on", TITLE),
         x="Time (years)", y = paste("Ratio Effect/Control -", YLAB)) +
    facet_grid(Biome~stress_level, scale="free")
}

################################################################################
#
# GRD
#
#
Grd_summary <- read_csv("output-EFSA_r20all/Grd_summary.csv.gz")
Grd_summary <- Grd_summary |> filter(stress_type != 110) |>
  mutate(
    Time_years = Time / 30,
    period = case_when(
      Time_years >= 0  & Time_years < 20 ~ "1_stabilize",
      Time_years >= 20 & Time_years < 30 ~ "2_treatment",
      Time_years >= 30 & Time_years <= 50 ~ "3_recovery",
      TRUE ~ "other"
    )
  )

# Application à Grd_summary
Grd_summary <- Grd_summary %>%
  rename(stress_type_num = stress_type) %>%
  mutate(
    code_str = sprintf("%03d", as.integer(stress_type_num)),
    # Extraire chaque chiffre
    c_hundreds = substr(code_str, 1, 1),
    c_tens     = substr(code_str, 2, 2),
    c_units    = substr(code_str, 3, 3),
    # Remplacer par le nom d'origine
    be2th_decoded        = be2th_map[c_hundreds],
    stress_level_decoded = stress_level_map[c_tens],
    stressor             = stressor_map[c_units]
  ) %>%
  select(-code_str, -c_hundreds, -c_tens, -c_units)

head(Grd_summary)

df_control <- Grd_summary %>%
  filter(stressor == "NO") %>%
  select(Biome, Time, stress_level, VarOutput, period,
         mean_control = mean, median_control = median)

# 2. Joindre avec le tableau complet et calculer le ratio
df_relative <- Grd_summary %>%
  left_join(df_control, by = c("Biome", "Time", "stress_level", "VarOutput", "period")) %>%
  mutate(
    mean_relative = mean / mean_control,
    median_relative = median / median_control,
    qinf95_relative = qinf95 / median_control,
    qsup95_relative = qsup95 / median_control,
  ) %>%
  filter(!is.na(median_relative), !is.infinite(median_relative))


plt_var(df_relative, "shannon", "Shannon Diversity", "survival", "direct mortality")
plt_var(df_relative, "NInd", "Number of Individuals", "survival", "direct mortality")
plt_var(df_relative, "NPFT", "Number of PFT", "survival", "direct mortality")
plt_var(df_relative, "abovemass", "Above Biomass", "survival", "direct mortality")
plt_var(df_relative, "belowmass", "belowmass Biomass", "survival", "direct mortality")

df_sel <- df_relative |>
  filter(VarOutput=="NPFT", stressor=="survival")
ggplot(data=df_sel,
       aes(x=Time/30, y=median_relative, ymin=qinf95_relative, ymax=qsup95_relative,
           fill=period)) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  geom_ribbon(alpha=0.7) +
  geom_line() +
  scale_fill_manual(
    values = c("1_stabilize"="lightgreen", "2_treatment"="#DD8888", "3_recovery"="#4488DD"),
    labels = c(
      "1_stabilize" = "0-20yr (stabilize)", 
      "2_treatment" = "20-30yr (treatment)", 
      "3_recovery"  = "30-50yr (recovery)"
    ),
    name = "Phase"
  ) +
  labs(title="Herbicide effect on direct mortality",
       x="Time (years)", y = "Ratio Effect/Control - Number of Plant Functional Type") +
  facet_grid(Biome~stress_level, scale="free")

df_relative_roll <- df_relative %>%
  arrange(Biome, stress_level, VarOutput, period, stressor, Time) %>%
  group_by(Biome, stress_level, VarOutput, period, stressor) %>%
  mutate(
    across(
      c(mean_relative, median_relative, qinf95_relative, qsup95_relative),
      # Remplacement de rollmean par rollapply
      ~ rollapply(.x, width = 10, FUN = mean, na.rm = TRUE, fill = NA, align = "center"),
      .names = "{.col}_roll"
    )
  ) %>%
  ungroup()

plt_var_roll(df_relative_roll, "shannon", "Shannon Diversity", "survival", "direct mortality")
plt_var_roll(df_relative_roll, "NInd", "Number of Individuals", "survival", "direct mortality")
plt_var_roll(df_relative_roll, "NPFT", "Number of PFT", "survival", "direct mortality")
plt_var_roll(df_relative_roll, "abovemass", "Above Biomass", "survival", "direct mortality")
plt_var_roll(df_relative_roll, "belowmass", "belowmass Biomass", "survival", "direct mortality")

df_sel <- df_relative_roll |>
  filter(VarOutput=="shannon", stressor=="survival")
ggplot(data=df_sel,
       aes(x=Time/30, y=median_relative_roll,
           ymin=qinf95_relative_roll, ymax=qsup95_relative_roll,
           fill=period)) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  geom_ribbon(alpha=0.7) +
  geom_line() +
  scale_fill_manual(
    values = c("1_stabilize"="lightgreen", "2_treatment"="#DD8888", "3_recovery"="#4488DD"),
    labels = c(
      "1_stabilize" = "0-20yr (stabilize)", 
      "2_treatment" = "20-30yr (treatment)", 
      "3_recovery"  = "30-50yr (recovery)"
    ),
    name = "Phase"
  ) +
  labs(title="Herbicide effect on direct mortality",
    x="Time (years)", y = "Ratio Effect/Control - Shannon Diveristy") +
  facet_grid(Biome~stress_level, scale="free")



# 1. Agréger les données pour "Pool All Biomes"
df_relative_roll <- df_relative %>%
  arrange(Biome, stress_level, VarOutput, period, stressor, Time) %>%
  group_by(Biome, stress_level, VarOutput, period, stressor) %>%
  mutate(
    across(
      c(mean_relative, median_relative, qinf95_relative, qsup95_relative),
      # Remplacement de rollmean par rollapply
      ~ rollapply(.x, width = 30, FUN = mean, na.rm = TRUE, fill = NA, align = "center"),
      .names = "{.col}_roll"
    )
  ) %>%
  ungroup()

df_pool <- df_relative |>
  filter(stressor != "NO") |>
  # On groupe par temps, niveau de stress et stresseur (en ignorant le Biome)
  group_by(Time, stress_level, stressor, VarOutput, period) |>
  summarise(
    # On calcule la médiane inter-biomes (vous pouvez utiliser mean() si préféré)
    median_relative = median(median_relative, na.rm = TRUE),
    qinf95_relative = median(qinf95_relative, na.rm = TRUE),
    qsup95_relative = median(qsup95_relative, na.rm = TRUE),
    .groups = "drop"
  ) |>
  # Optionnel : Forcer l'ordre d'affichage des stresseurs comme sur l'image
  mutate(stressor = factor(stressor, levels = c("biomass", "SEbiomass", "survival", "sterility", "seednumber", "establishment")))

df_pool_roll <- df_pool %>%
  filter(VarOutput=="shannon") %>%
  arrange(stress_level, stressor, VarOutput, Time, period) %>%
  group_by(stress_level, stressor, VarOutput, period) %>%
  mutate(
    across(
      c(median_relative, qinf95_relative, qsup95_relative),
      # Remplacement de rollmean par rollapply
      ~ rollapply(.x, width = 5, FUN = mean, na.rm = TRUE, fill = NA, align = "center"),
      .names = "{.col}_roll"
    )
  ) %>%
  filter(!is.na(median_relative_roll), !is.na(qinf95_relative_roll), !is.na(qsup95_relative_roll)) %>%
  ungroup()
# 2. Création du graphique
ggplot(data = df_pool_roll, 
       aes(x = Time/30, y = median_relative_roll, ymin = qinf95_relative_roll, ymax = qsup95_relative_roll, 
           fill = period)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(color = "grey80", fill = NA) # Pour bien délimiter les panneaux
  ) +
  geom_ribbon(alpha = 0.5, color = NA) + # color=NA pour ne pas avoir de bordure sur le ruban
  geom_line(linewidth = 0.5) +
  # scale_fill_manual(
  #   values = c("z00" = "#555555", "z01" = "#78C850", "z02" = "#E0C048", 
  #              "z03" = "#D07040", "z04" = "#C04040", "z05" = "#8858C8"),
  #   name = "Stressor level:"
  # ) +
  # scale_color_manual(
  #   values = c("z00" = "#555555", "z01" = "#78C850", "z02" = "#E0C048",
  #              "z03" = "#D07040", "z04" = "#C04040", "z05" = "#8858C8"),
  #   name = "Stressor level:"
  # ) +
  scale_fill_manual(
    values = c("1_stabilize"="lightgreen", "2_treatment"="#DD8888", "3_recovery"="#4488DD"),
    labels = c(
      "1_stabilize" = "0-20yr (stabilize)", 
      "2_treatment" = "20-30yr (treatment)", 
      "3_recovery"  = "30-50yr (recovery)"
    ),
    name = "Phase"
  ) +
  labs(
    title = "Pool All Biomes",
    x = "Time (years)", 
    y = "Ratio Effect/Control Shannon Diversity"
  ) +
  facet_grid(stress_level ~ stressor) 
  # scale_y_continuous(breaks = seq(0.8, 1.1, by = 0.1), limits = c(0.78, 1.12))

df_pool_roll <- df_pool %>%
  filter(VarOutput=="NPFT") %>%
  arrange(stress_level, stressor, Time) %>%
  group_by(stress_level, stressor) %>%
  mutate(
    across(
      c(median_relative, qinf95_relative, qsup95_relative),
      # Remplacement de rollmean par rollapply
      ~ rollapply(.x, width = 30, FUN = mean, na.rm = TRUE, fill = NA, align = "center"),
      .names = "{.col}_roll"
    )
  ) %>%
  filter(!is.na(median_relative_roll), !is.na(qinf95_relative_roll), !is.na(qsup95_relative_roll)) %>%
  ungroup()
# 2. Création du graphique
ggplot(data = df_pool_roll, 
       aes(x = Time/30, y = median_relative_roll, ymin = qinf95_relative_roll, ymax = qsup95_relative_roll, 
           fill = period)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(color = "grey80", fill = NA) # Pour bien délimiter les panneaux
  ) +
  geom_ribbon(alpha = 0.5, color = NA) + # color=NA pour ne pas avoir de bordure sur le ruban
  geom_line(linewidth = 0.5) +
  scale_fill_manual(
    values = c("1_stabilize"="lightgreen", "2_treatment"="#DD8888", "3_recovery"="#4488DD"),
    labels = c(
      "1_stabilize" = "0-20yr (stabilize)", 
      "2_treatment" = "20-30yr (treatment)", 
      "3_recovery"  = "30-50yr (recovery)"
    ),
    name = "Phase"
  ) +
  labs(
    title = "Pool All Biomes",
    x = "Time (years)", 
    y = "Ratio Effect/Control Number PFTs"
  ) +
  facet_grid(stress_level ~ stressor) 

df_pool_roll <- df_pool %>%
  filter(VarOutput=="totMass") %>%
  arrange(stress_level, stressor, Time) %>%
  group_by(stress_level, stressor) %>%
  mutate(
    across(
      c(median_relative, qinf95_relative, qsup95_relative),
      # Remplacement de rollmean par rollapply
      ~ rollapply(.x, width = 30, FUN = mean, na.rm = TRUE, fill = NA, align = "center"),
      .names = "{.col}_roll"
    )
  ) %>%
  filter(!is.na(median_relative_roll), !is.na(qinf95_relative_roll), !is.na(qsup95_relative_roll)) %>%
  ungroup()
# 2. Création du graphique
ggplot(data = df_pool_roll, 
       aes(x = Time/30, y = median_relative_roll, ymin = qinf95_relative_roll, ymax = qsup95_relative_roll, 
           fill = period)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(color = "grey80", fill = NA) # Pour bien délimiter les panneaux
  ) +
  geom_ribbon(alpha = 0.5, color = NA) + # color=NA pour ne pas avoir de bordure sur le ruban
  geom_line(linewidth = 0.5) +
  scale_fill_manual(
    values = c("1_stabilize"="lightgreen", "2_treatment"="#DD8888", "3_recovery"="#4488DD"),
    labels = c(
      "1_stabilize" = "0-20yr (stabilize)", 
      "2_treatment" = "20-30yr (treatment)", 
      "3_recovery"  = "30-50yr (recovery)"
    ),
    name = "Phase"
  ) +
  labs(
    title = "Pool All Biomes",
    x = "Time (years)", 
    y = "Ratio Effect/Control Total Biomass"
  ) +
  facet_grid(stress_level ~ stressor) 

df_pool_roll <- df_pool %>%
  filter(VarOutput=="NInd") %>%
  arrange(stress_level, stressor, Time) %>%
  group_by(stress_level, stressor) %>%
  mutate(
    across(
      c(median_relative, qinf95_relative, qsup95_relative),
      # Remplacement de rollmean par rollapply
      ~ rollapply(.x, width = 30, FUN = mean, na.rm = TRUE, fill = NA, align = "center"),
      .names = "{.col}_roll"
    )
  ) %>%
  filter(!is.na(median_relative_roll), !is.na(qinf95_relative_roll), !is.na(qsup95_relative_roll)) %>%
  ungroup()
# 2. Création du graphique
ggplot(data = df_pool_roll, 
       aes(x = Time/30, y = median_relative_roll, ymin = qinf95_relative_roll, ymax = qsup95_relative_roll, 
           fill = period)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(color = "grey80", fill = NA) # Pour bien délimiter les panneaux
  ) +
  geom_ribbon(alpha = 0.5, color = NA) + # color=NA pour ne pas avoir de bordure sur le ruban
  geom_line(linewidth = 0.5) +
  scale_fill_manual(
    values = c("1_stabilize"="lightgreen", "2_treatment"="#DD8888", "3_recovery"="#4488DD"),
    labels = c(
      "1_stabilize" = "0-20yr (stabilize)", 
      "2_treatment" = "20-30yr (treatment)", 
      "3_recovery"  = "30-50yr (recovery)"
    ),
    name = "Phase"
  ) +
  labs(
    title = "Pool All Biomes",
    x = "Time (years)", 
    y = "Ratio Effect/Control Total Number of Individuals"
  ) +
  facet_grid(stress_level ~ stressor) 

################################################################################
#
# Pt
#
#
Pt_summary <- read_csv("output-EFSA_r20all/Pt_summary_noPFT.csv.gz")
Pt_summary <- Pt_summary |>
  mutate(
    Time_years = Time / 30,
    period = case_when(
      Time_years >= 0  & Time_years < 20 ~ "1_stabilize",
      Time_years >= 20 & Time_years < 30 ~ "2_treatment",
      Time_years >= 30 & Time_years <= 50 ~ "3_recovery",
      TRUE ~ "other"
    )
  )

# Application à Grd_summary
Pt_summary <- Pt_summary %>%
  rename(stress_type_num = stress_type) %>%
  mutate(
    code_str = sprintf("%03d", as.integer(stress_type_num)),
    # Extraire chaque chiffre
    c_hundreds = substr(code_str, 1, 1),
    c_tens     = substr(code_str, 2, 2),
    c_units    = substr(code_str, 3, 3),
    # Remplacer par le nom d'origine
    be2th_decoded        = be2th_map[c_hundreds],
    stress_level_decoded = stress_level_map[c_tens],
    stressor             = stressor_map[c_units]
  ) %>%
  select(-code_str, -c_hundreds, -c_tens, -c_units)

head(Pt_summary)

df_control <- Pt_summary %>%
  filter(stressor == "NO") %>%
  select(Biome, Time, stress_level, VarOutput, period,
         mean_control = mean, median_control = median)

# 2. Joindre avec le tableau complet et calculer le ratio
df_relative <- Pt_summary %>%
  left_join(df_control, by = c("Biome", "Time", "stress_level", "VarOutput", "period")) %>%
  mutate(
    mean_relative = mean / mean_control,
    median_relative = median / median_control,
    qinf95_relative = qinf95 / median_control,
    qsup95_relative = qsup95 / median_control,
  )

df_sel <- df_relative |>
  filter(VarOutput=="Inds", stressor=="survival") |>
  filter(!is.na(median_relative), !is.infinite(median_relative))
ggplot(data=df_sel,
       aes(x=Time/30, y=median_relative, ymin=qinf95_relative, ymax=qsup95_relative,
           fill=period)) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  geom_ribbon(alpha=0.7) +
  geom_line() +
  scale_fill_manual(
    values = c("1_stabilize"="lightgreen", "2_treatment"="#DD8888", "3_recovery"="#4488DD"),
    labels = c(
      "1_stabilize" = "0-20yr (stabilize)", 
      "2_treatment" = "20-30yr (treatment)", 
      "3_recovery"  = "30-50yr (recovery)"
    ),
    name = "Phase"
  ) +
  labs(title="Herbicide effect on direct mortality",
       x="Time (years)", y = "Ratio Effect/Control - Number of individuals") +
  facet_grid(Biome~stress_level, scale="free")

df_relative_roll <- df_relative %>%
  arrange(Biome, stress_level, VarOutput, period, stressor, Time) %>%
  group_by(Biome, stress_level, VarOutput, period, stressor) %>%
  mutate(
    across(
      c(mean_relative, median_relative, qinf95_relative, qsup95_relative),
      # Remplacement de rollmean par rollapply
      ~ rollapply(.x, width = 5, FUN = mean, na.rm = TRUE, fill = NA, align = "center"),
      .names = "{.col}_roll"
    )
  ) %>%
  ungroup()

df_sel <- df_relative_roll |>
  filter(VarOutput=="Inds", stressor=="survival") # |>
  # filter(!is.na(median_relative_roll), !is.infinite(median_relative_roll))
ggplot(data=df_sel,
       aes(x=Time/30, y=median_relative_roll,
           ymin=qinf95_relative_roll, ymax=qsup95_relative_roll,
           fill=period)) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  geom_ribbon(alpha=0.7) +
  geom_line() +
  scale_fill_manual(
    values = c("1_stabilize"="lightgreen", "2_treatment"="#DD8888", "3_recovery"="#4488DD"),
    labels = c(
      "1_stabilize" = "0-20yr (stabilize)", 
      "2_treatment" = "20-30yr (treatment)", 
      "3_recovery"  = "30-50yr (recovery)"
    ),
    name = "Phase"
  ) +
  labs(title="Herbicide effect on direct mortality",
       x="Time (years)", y = "Ratio Effect/Control - Number of individuals") +
  facet_grid(Biome~stress_level, scale="free")

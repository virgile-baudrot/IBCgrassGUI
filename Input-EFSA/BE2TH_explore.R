library(dplyr)
library(ggplot2)
library(tidygraph)
library(readr)
library(tidyverse)
library(reshape2)

df_BE2TH <- readr::read_csv("Input-EFSA/BE2TH_cleaned.csv")




# Extract just the biome columns mentioned in the plots
biome_cols <- c("ALP_I1", "ATL_I1", "BLS_I1", "BOR_I1", 
                "CON_I1", "MAC_I1", "MED_I1", "PAN_I1", "STE_I1")

# Subset the dataframe to only these columns and ensure they are numeric
df_biomes <- df_BE2TH %>%
  select(all_of(biome_cols)) %>%
  mutate(across(everything(), as.numeric))

###### BAR PLOT
pft_cols <- c("Category", "LMR", "MaxMass..mg.", "mSeed..mg.", "Dist..m.", 
              "Gmax", "Memo..weeks.", "SLA.comb.", "Palat", "Clonal", 
              "meanSpacerLength..cm.", "sdSpacerLength..cm.", "Resshare")

# 2. Préparer les données
bar_data <- df_BE2TH %>%
  # Créer une colonne PFT_ID qui fusionne les valeurs de toutes les colonnes de traits
  unite("PFT_ID", all_of(pft_cols), sep = "_", remove = FALSE) %>%
  select(PFT_ID, all_of(biome_cols)) %>%
  # Convertir les biomes en format long
  pivot_longer(cols = all_of(biome_cols), names_to = "Biome", values_to = "Presence") %>%
  # Garder uniquement les présences
  filter(Presence == 1) %>%
  group_by(Biome) %>%
  summarise(
    Number_species = n(),                 # Nombre total d'espèces (lignes)
    Number_PFT = n_distinct(PFT_ID)       # Nombre de combinaisons uniques de traits (PFTs)
  ) %>%
  mutate(Biome = str_remove(Biome, "_I1")) %>%
  pivot_longer(cols = c("Number_species", "Number_PFT"), 
               names_to = "Number", values_to = "Count")

# 3. Générer le graphique
ggplot(bar_data, aes(x = Biome, y = Count, fill = Number)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("Number_PFT" = "#d9d9d9", "Number_species" = "#969696")) +
  labs(x = "Biomes",
       y = "Number of Species or PFT") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank())


ç###### CORRELATION
# Calculate correlation matrix
cor_matrix <- cor(df_biomes, use = "pairwise.complete.obs")

# Convert to long format for ggplot
cor_long <- melt(cor_matrix)

# Plot
ggplot(cor_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 2)), size = 3) +
  scale_fill_gradient2(low = "#67a9cf", mid = "white", high = "#ef8a62", 
                       midpoint = 0, limit = c(-1, 1), name = "Correlation") +
  labs(title = "Correlation between biomes",
       x = "Biomes",
       y = "Biomes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5))


###### CO-OCCURENCE
# Calculate co-occurrence matrix using matrix multiplication (crossprod)
co_matrix <- crossprod(as.matrix(df_biomes))

# Convert to long format for ggplot
co_long <- melt(co_matrix)

# Plot
ggplot(co_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = value), size = 3) +
  scale_fill_gradient(low = "white", high = "#805ad5", 
                      name = "Number of\n Species") +
  labs(title = "Co-occurrence of Species",
       x = "Biomes",
       y = "Biomes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5))

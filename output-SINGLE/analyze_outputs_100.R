library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)

################################################################################
#
# SINGLE SPECIES - SEVERAL SCENARIOS
#
################################################################################
# 2. Define scenarios to map the numeric IDs
scenarios_name <- list(
  "Control (No Herbicide)"              = 101,
  "Biomass Reduction"                   = 102, # biomass
  "Seedling Biomass Reduction"          = 103, # SEbiomass
  "Increased Mortality (Survival)"      = 104, # survival
  "Reduced Establishment Probability"   = 105, # establishment
  "Reduced Fertility (Seed Sterility)"  = 106, # sterility
  "Reduced Seed Number"                 = 107  # seednumber
)

lookup_table <- tibble(
  scenario_name = names(scenarios_name),
  Scenario = as.character(unlist(scenarios_name)) 
)

# 3. Define the exact order for the legend to match the color gradient
legend_order <- c(
  "Control (No Herbicide)",              # Green
  "Biomass Reduction",                   # Yellow
  "Seedling Biomass Reduction",          # Light Orange
  "Reduced Establishment Probability",   # Dark Orange
  "Increased Mortality (Survival)",      # Red
  "Reduced Fertility (Seed Sterility)",  # Dark Red
  "Reduced Seed Number"                  # Burgundy
)
# 6. Color palette: Green for the control, warm colors for the stressors
custom_colors <- c(
  "Control (No Herbicide)"              = "#1a9850", # Green
  "Biomass Reduction"                   = "#fee08b", # Very light yellow/orange
  "Seedling Biomass Reduction"          = "#fdae61", # Light orange
  "Reduced Establishment Probability"   = "#f46d43", # Dark orange
  "Increased Mortality (Survival)"      = "#d73027", # Red
  "Reduced Fertility (Seed Sterility)"  = "#a50026", # Dark red
  "Reduced Seed Number"                 = "#67001f"  # Burgundy
)

################################################################################
# Pt data
################################################################################
pt_long_data <- readRDS("output-SINGLE/pt_long_data_100.rds")


# 4. Join and convert to factor with the specific gradient order
pt_long_data <- pt_long_data %>%
  left_join(lookup_table, by = "Scenario") %>%
  mutate(scenario_name = factor(scenario_name, levels = legend_order))

# 5. Rename output variables (OutputVariable) for the facets
var_labels <- c(
  "Inds"      = "Number of Individuals (Inds)",
  "seedlings" = "Number of Seedlings",
  "seeds"     = "Seed Bank (Seeds)",
  "repromass" = "Reproductive Biomass (mg)",
  "cover"     = "Cover / Area of Influence",
  "shootmass" = "Above-ground Shoot Biomass (mg)"
)

# 7. Plot 1: Raw Values (Control on top)
plot_raw <- pt_long_data %>%
  # Sort the data so the Control is plotted last (on top of everything else)
  arrange(scenario_name == "Control (No Herbicide)") %>%
  ggplot(aes(x = time/30, y = Value, color = scenario_name)) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  facet_wrap(~ OutputVariable, scales = "free_y", ncol = 2, labeller = labeller(OutputVariable = var_labels)) +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  labs(
    title = "Temporal Dynamics of Model Outputs",
    # subtitle = "Comparison of raw values between the control and herbicide stress scenarios",
    x = "Time (Years)",
    y = "Output Value",
    color = "Scenarios"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

print(plot_raw)

# 8. Plot 2: Ratio to Control
pt_long_data_dropNA <- pt_long_data[!is.na(pt_long_data$Ratio_to_Control), ]

plot_ratio <- pt_long_data_dropNA %>%
  # Same trick for the ratio plot if you want the Control line (which will be perfectly flat at 1.0) on top
  arrange(scenario_name == "Control (No Herbicide)") %>%
  ggplot(aes(x = time/30, y = Ratio_to_Control, color = scenario_name)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", alpha = 0.5) +
  geom_line(linewidth = 0.6, alpha = 0.5) +
  facet_wrap(~ OutputVariable, scales = "free_y", ncol = 2, labeller = labeller(OutputVariable = var_labels)) +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  labs(
    title = "Relative Effects of Herbicide Over Time",
    # subtitle = "Ratio (Scenario Value / Control Value) for each variable",
    x = "Time (Years)",
    y = "Ratio (Treatment / Control)",
    color = "Scenarios"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

print(plot_ratio)

################################################################################
# Grd data
################################################################################
# 1. Load the data
grd_long_data <- readRDS("output-SINGLE/grd_long_data_100.rds")

# 4. Join and convert to factor with the specific gradient order
grd_long_data <- grd_long_data %>%
  left_join(lookup_table, by = "Scenario") %>%
  mutate(scenario_name = factor(scenario_name, levels = legend_order))

# 5. Rename output variables (OutputVariable) for the facets based on IBC-grass documentation
var_labels_grd <- c(
  "totMass"         = "Total Biomass (mg)",
  "NInd"            = "Number of Individuals",
  "abovemass"       = "Above-ground Biomass (mg)",
  "belowmass"       = "Below-ground Biomass (mg)",
  "mean_ares"       = "Mean Above-ground Resources",
  "mean_bres"       = "Mean Below-ground Resources",
  "shannon"         = "Shannon Diversity Index",
  "meanShannon"     = "Mean Shannon Diversity",
  "NPFT"            = "PFT Richness (Count)",
  "meanNPFT"        = "Mean PFT Richness",
  "Cutted"          = "Harvested Biomass (Cut)",
  "NNonClonal"      = "Number of Non-Clonal Plants",
  "NClonal"         = "Number of Clonal Plants",
  "mean_generation" = "Mean Generation Age",
  "mean_genet_size" = "Mean Genet Size",
  "NGenets"         = "Number of Genets"
)


# 7. Plot 1: Raw Values (Control on top)
plot_raw <- grd_long_data %>%
  # Sort the data so the Control is plotted last (on top of everything else)
  arrange(scenario_name == "Control (No Herbicide)") %>%
  ggplot(aes(x = time/30, y = Value, color = scenario_name)) +
  geom_line(linewidth = 0.5, alpha = 0.8) +
  facet_wrap(~ OutputVariable, scales = "free_y", ncol = 4, labeller = labeller(OutputVariable = var_labels_grd)) +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  labs(
    title = "Community-Level Dynamics (Grid Outputs)",
    # subtitle = "Raw values of spatial, resource, and biodiversity metrics over time",
    x = "Time (Years)",
    y = "Output Value",
    color = "Scenarios"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

print(plot_raw)


grd_long_data_dropNA <- grd_long_data[!is.na(grd_long_data$Ratio_to_Control), ]

plot_ratio <- grd_long_data_dropNA %>%
  arrange(scenario_name == "Control (No Herbicide)") %>%
  ggplot(aes(x = time/30, y = Ratio_to_Control, color = scenario_name)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", alpha = 0.5) +
  geom_line(linewidth = 0.5, alpha = 0.8) +
  facet_wrap(~ OutputVariable, scales = "free_y", ncol = 4, labeller = labeller(OutputVariable = var_labels_grd)) +
  scale_color_manual(values = custom_colors) +
  theme_minimal() +
  labs(
    title = "Relative Herbicide Effects on Community Structure",
    # subtitle = "Ratio (Scenario Value / Control Value) for grid variables",
    x = "Time (Years)",
    y = "Ratio (Treatment / Control)",
    color = "Scenarios"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

print(plot_ratio)




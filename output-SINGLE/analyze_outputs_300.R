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
  "Control (No Herbicide)"              = "control",
  "Increased Mortality (Survival)"      = "treatment" # survival
)

lookup_table <- tibble(
  scenario_name = names(scenarios_name),
  Scenario_Type = as.character(unlist(scenarios_name)) 
)

# 3. Define the exact order for the legend to match the color gradient
legend_order <- c(
  "Control (No Herbicide)",              # Green
  "Increased Mortality (Survival)"      # Red
)
# 6. Color palette: Green for the control, warm colors for the stressors
custom_colors <- c(
  "control"              = "#1a9850", # Green
  "treatment"      = "#d73027" # Red
)

################################################################################
# Pt data
################################################################################
pt_long_data <- readRDS("output-SINGLE/pt_long_data_300.rds")

# 4. Join and convert to factor with the specific gradient order
pt_long_data <- pt_long_data %>%
  left_join(lookup_table, by = "Scenario_Type") %>%
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

pt_long_data_raw <- pt_long_data |>
  dplyr::group_by(Scenario_Type, time, PFT, OutputVariable) |>
  dplyr::summarise(
    mean = mean(Value),
    median = quantile(Value, 0.5),
    qinf95 = quantile(Value, 0.025),
    qsup95 = quantile(Value, 0.975)
  )

# pt_long_data_density <- pt_long_data |>
#   dplyr::mutate(
#     cat_time = dplyr::case_when(
#       time <= 20 ~ "pre-treatment",
#       time <= 30 ~ "treatment",
#       TRUE       ~ "post-treatment"
#     )
#   )
# ggplot(pt_long_data_density, aes(Value, fill=Scenario_Type)) +
#   scale_x_log10() +
#   scale_fill_manual(values = custom_colors) +
#   facet_grid(cat_time ~ OutputVariable, scales = "free") + 
#   geom_density(alpha=0.5)

# 7. Plot 1: Raw Values (Control on top)
plot_raw <- pt_long_data_raw %>%
  # Sort the data so the Control is plotted last (on top of everything else)
  # arrange(scenario_name == "Control (No Herbicide)") %>%
  ggplot(aes(x = time/30, y = mean, ymin=qinf95, ymax=qsup95, color = Scenario_Type, fill=Scenario_Type)) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  geom_ribbon(linewidth = 0.05, alpha = 0.5) +
  facet_wrap(~ OutputVariable, scales = "free_y", ncol = 2, labeller = labeller(OutputVariable = var_labels)) +
  scale_color_manual(values = custom_colors, name="Mean") +
  scale_fill_manual(values = custom_colors, name="95% Interval") +
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


pt_long_data_dropNA <- pt_long_data[!is.na(pt_long_data$Ratio_to_Control), ] |>
  dplyr::group_by(Scenario_Type, time, PFT, OutputVariable) |>
  dplyr::summarise(
    mean = mean(Ratio_to_Control),
    median = median(Ratio_to_Control),
    qinf95 = quantile(Ratio_to_Control, 0.025),
    qsup95 = quantile(Ratio_to_Control, 0.975)
  )
plot_ratio <- pt_long_data_dropNA %>%
  # Same trick for the ratio plot if you want the Control line (which will be perfectly flat at 1.0) on top
  # arrange(scenario_name == "Control (No Herbicide)") %>%
  ggplot(aes(x = time/30, y = mean, ymin=qinf95, ymax=qsup95, color=Scenario_Type, fill=Scenario_Type)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", alpha = 0.5) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  geom_ribbon(linewidth = 0.05, alpha = 0.5) +
  facet_wrap(~ OutputVariable, scales = "free_y", ncol = 2, labeller = labeller(OutputVariable = var_labels)) +
  scale_color_manual(values = custom_colors, name="Mean") +
  scale_fill_manual(values = custom_colors, name="95% Interval") +
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
grd_long_data <- readRDS("output-SINGLE/grd_long_data_300.rds")

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
grd_long_data_raw <- grd_long_data |>
  dplyr::group_by(Scenario_Type, time, OutputVariable) |>
  dplyr::summarise(
    mean = mean(Ratio_to_Control),
    median = median(Ratio_to_Control, na.rm=TRUE),
    qinf95 = quantile(Ratio_to_Control, 0.025, na.rm=TRUE),
    qsup95 = quantile(Ratio_to_Control, 0.975, na.rm=TRUE)
  )
plot_raw <- grd_long_data_raw %>%
  # Sort the data so the Control is plotted last (on top of everything else)
  # arrange(scenario_name == "Control (No Herbicide)") %>%
  ggplot(aes(x = time/30, y = mean, ymin=qinf95, ymax=qsup95, color=Scenario_Type, fill=Scenario_Type)) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  geom_ribbon(linewidth = 0.05, alpha = 0.5) +
  facet_wrap(~ OutputVariable, scales = "free_y", ncol = 4, labeller = labeller(OutputVariable = var_labels_grd)) +
  scale_color_manual(values = custom_colors, name="Mean") +
  scale_fill_manual(values = custom_colors, name="95% Interval") +
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


grd_long_data_dropNA <- grd_long_data[!is.na(grd_long_data$Ratio_to_Control), ] |>
  dplyr::group_by(Scenario_Type, time, OutputVariable) |>
  dplyr::summarise(
    mean = mean(Ratio_to_Control),
    median = median(Ratio_to_Control, na.rm=TRUE),
    qinf95 = quantile(Ratio_to_Control, 0.025, na.rm=TRUE),
    qsup95 = quantile(Ratio_to_Control, 0.975, na.rm=TRUE)
  )

plot_ratio <- grd_long_data_dropNA %>%
  # arrange(scenario_name == "Control (No Herbicide)") %>%
  ggplot(aes(x = time/30, y = mean, ymin=qinf95, ymax=qsup95, color=Scenario_Type, fill=Scenario_Type)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", alpha = 0.5) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  geom_ribbon(linewidth = 0.05, alpha = 0.5) +
  facet_wrap(~ OutputVariable, scales = "free_y", ncol = 4, labeller = labeller(OutputVariable = var_labels_grd)) +
  scale_color_manual(values = custom_colors, name="Mean") +
  scale_fill_manual(values = custom_colors, name="95% Interval") +
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




# §4.1.1 Single species, one run per stressor: Figures 4.1, 4.2 (Pt_) and 4.4, 4.5 (Grd_)
#   Rscript output-SINGLE/analyze_outputs_100.R [in_dir] [fig_dir]
#   in_dir: outputs/4.1_single_species/stressors (run_4.1_single_species.R stressors)
# The 7 runs share one seed, so each treatment is compared with the control run of the same
# seed: the ratio is exactly 1 before the first herbicide year.

source("Input-EFSA/04_output_analysis.R")
args <- commandArgs(trailingOnly = TRUE)
in_dir <- if (length(args) >= 1) args[1] else "outputs/4.1_single_species/stressors"
fig_dir <- if (length(args) >= 2) args[2] else "EFSA_REPORT/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
prov <- read_provenance(in_dir)

# Scenario names, legend order and colours (control green, stressors warm colours)
scenario_names <- c(
  NO = "Control (No Herbicide)",
  biomass = "Biomass Reduction",
  SEbiomass = "Seedling Biomass Reduction",
  establishment = "Reduced Establishment Probability",
  survival = "Increased Mortality (Survival)",
  sterility = "Reduced Fertility (Seed Sterility)",
  seednumber = "Reduced Seed Number"
)
custom_colors <- setNames(c("#1a9850", "#fee08b", "#fdae61", "#f46d43", "#d73027", "#a50026", "#67001f"),
                          scenario_names)

var_labels <- c(
  Inds = "Number of Individuals (Inds)",
  seedlings = "Number of Seedlings",
  seeds = "Seed Bank (Seeds)",
  repromass = "Reproductive Biomass (mg)",
  cover = "Cover / Area of Influence",
  shootmass = "Above-ground Shoot Biomass (mg)"
)
var_labels_grd <- c(
  totMass = "Total Biomass (mg)",
  NInd = "Number of Individuals",
  abovemass = "Above-ground Biomass (mg)",
  belowmass = "Below-ground Biomass (mg)",
  mean_ares = "Mean Above-ground Resources",
  mean_bres = "Mean Below-ground Resources",
  shannon = "Shannon Diversity Index",
  NPFT = "PFT Richness (Count)",
  Cutted = "Harvested Biomass (Cut)",
  mean_generation = "Mean Generation Age",
  mean_genet_size = "Mean Genet Size"
)

# long format, with the paired control value and the ratio to it
to_long <- function(dt, vars) {
  dt[, (vars) := lapply(.SD, as.numeric), .SDcols = vars]
  long <- melt(dt, id.vars = c("stressor", "Time"), measure.vars = vars,
               variable.name = "OutputVariable", value.name = "Value", variable.factor = FALSE)
  ctl <- long[stressor == "NO", .(Time, OutputVariable, Control_Value = Value)]
  long <- ctl[long, on = .(Time, OutputVariable)]
  long[, Ratio_to_Control := fifelse(Control_Value > 0, Value / Control_Value, NA_real_)]
  long[, scenario_name := factor(scenario_names[stressor], levels = scenario_names)]
  rbind(long[stressor != "NO"], long[stressor == "NO"]) # control drawn last (on top)
}

plot_lines <- function(long, y, labels, ncol, title, ylab, ref = NULL) {
  ggplot(long[!is.na(get(y))], aes(x = Time / 30, y = .data[[y]], color = scenario_name)) +
    phase_layer(prov) +
    (if (!is.null(ref)) geom_hline(yintercept = ref, linetype = "dashed", color = "black", alpha = 0.5)) +
    geom_line(linewidth = 0.5, alpha = 0.8) +
    facet_wrap(~ OutputVariable, scales = "free_y", ncol = ncol,
               labeller = labeller(OutputVariable = labels)) +
    scale_color_manual(values = custom_colors) +
    theme_minimal() +
    labs(title = title, x = "Time (Years)", y = ylab, color = "Scenarios") +
    theme(strip.text = element_text(face = "bold", size = 9),
          legend.position = "bottom", legend.title = element_text(face = "bold"))
}

# ---------------- Pt data (Figures 4.1, 4.2) ----------------
pt <- fread(file.path(in_dir, "Pt.csv.gz"))[, c("seed", "rep") := NULL][, community := 1L]
pt <- complete_pt(pt, c("community", "stressor"), "community")[, community := NULL]
pt_long <- to_long(pt, names(var_labels))

save_fig(plot_lines(pt_long, "Value", var_labels, 2, "Temporal Dynamics of Model Outputs", "Output Value"),
         file.path(fig_dir, "fig_4.1_stressors_Pt_raw.png"), 16, 16)
save_fig(plot_lines(pt_long, "Ratio_to_Control", var_labels, 2, "Relative Effects of Herbicide Over Time",
                    "Ratio (Treatment / Control)", ref = 1),
         file.path(fig_dir, "fig_4.1_stressors_Pt_ratio.png"), 16, 16)

# ---------------- Grd data (Figures 4.4, 4.5) ----------------
# meanShannon, meanNPFT, NNonClonal, NClonal and NGenets repeat one end-of-run value (R§5.2)
grd <- fread(file.path(in_dir, "Grd.csv.gz"))
grd_long <- to_long(grd, names(var_labels_grd))

save_fig(plot_lines(grd_long, "Value", var_labels_grd, 4, "Community-Level Dynamics (Grid Outputs)", "Output Value"),
         file.path(fig_dir, "fig_4.1_stressors_Grd_raw.png"), 20, 16)
save_fig(plot_lines(grd_long, "Ratio_to_Control", var_labels_grd, 4,
                    "Relative Herbicide Effects on Community Structure", "Ratio (Treatment / Control)", ref = 1),
         file.path(fig_dir, "fig_4.1_stressors_Grd_ratio.png"), 20, 16)

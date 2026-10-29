# §4.1.2 Single species, stochasticity of controls and treatment: Figures 4.6 (Pt_) and 4.7 (Grd_)
#   Rscript output-SINGLE/analyze_outputs_300.R [in_dir] [fig_dir]
#   in_dir: outputs/4.1_single_species/variability (run_4.1_single_species.R variability)
# n replicates of the control and of the survival treatment; replicate i of both uses seed i.
# Raw figures: mean and 95 % interval (2.5 and 97.5 percentiles) over replicates.
# Ratio figures:
#   treatment (red): treatment run i / control run i (paired, exactly 1 before exposure)
#   control (green): control run i / mean of all control runs, i.e. the variability between
#   runs due to the stochasticity of the model alone

source("Input-EFSA/04_output_analysis.R")
args <- commandArgs(trailingOnly = TRUE)
in_dir <- if (length(args) >= 1) args[1] else "outputs/4.1_single_species/variability"
fig_dir <- if (length(args) >= 2) args[2] else "EFSA_REPORT/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
prov <- read_provenance(in_dir)

custom_colors <- c(control = "#1a9850", treatment = "#d73027")

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
  Cutted = "Harvested Biomass (Cut)",
  mean_generation = "Mean Generation Age",
  mean_genet_size = "Mean Genet Size"
) # meanShannon, meanNPFT, NNonClonal, NClonal, NGenets: one end-of-run value (R§5.2)

to_long <- function(dt, vars) {
  dt[, (vars) := lapply(.SD, as.numeric), .SDcols = vars]
  long <- melt(dt, id.vars = c("stressor", "rep", "Time"), measure.vars = vars,
               variable.name = "OutputVariable", value.name = "Value", variable.factor = FALSE)
  long[, Scenario_Type := fifelse(stressor == "NO", "control", "treatment")]
  ctl <- long[Scenario_Type == "control", .(rep, Time, OutputVariable, Control_Value = Value)]
  long <- ctl[long, on = .(rep, Time, OutputVariable)]
  long[, Control_Mean := mean(Control_Value), by = .(Time, OutputVariable)]
  long[, Ratio := fifelse(Scenario_Type == "control",
                          fifelse(Control_Mean > 0, Value / Control_Mean, NA_real_),
                          fifelse(Control_Value > 0, Value / Control_Value, NA_real_))]
  long
}

summarise_reps <- function(long, y) {
  long[!is.na(get(y)), .(mean = mean(get(y)),
                         qinf95 = quantile(get(y), 0.025, names = FALSE),
                         qsup95 = quantile(get(y), 0.975, names = FALSE)),
       by = .(Scenario_Type, Time, OutputVariable)]
}

plot_band <- function(s, labels, ncol, title, ylab, ref = NULL) {
  ggplot(s, aes(x = Time / 30, y = mean, ymin = qinf95, ymax = qsup95,
                color = Scenario_Type, fill = Scenario_Type)) +
    phase_layer(prov) +
    (if (!is.null(ref)) geom_hline(yintercept = ref, linetype = "dashed", color = "black", alpha = 0.5)) +
    geom_ribbon(linewidth = 0.05, alpha = 0.4) +
    geom_line(linewidth = 0.6, alpha = 0.8) +
    facet_wrap(~ OutputVariable, scales = "free_y", ncol = ncol,
               labeller = labeller(OutputVariable = labels)) +
    scale_color_manual(values = custom_colors, name = "Mean") +
    scale_fill_manual(values = custom_colors, name = "95% Interval") +
    theme_minimal() +
    labs(title = title, x = "Time (Years)", y = ylab) +
    theme(strip.text = element_text(face = "bold", size = 9),
          legend.position = "bottom", legend.title = element_text(face = "bold"))
}

# ---------------- Pt data (Figure 4.6) ----------------
pt <- fread(file.path(in_dir, "Pt.csv.gz"))[, seed := NULL]
pt <- complete_pt(pt, c("stressor", "rep"), "rep") # zeros when the species is absent
pt_long <- to_long(pt, names(var_labels))
save_fig(plot_band(summarise_reps(pt_long, "Value"), var_labels, 2,
                   "Temporal Dynamics of Model Outputs", "Output Value"),
         file.path(fig_dir, "fig_4.1_variability_Pt_raw.png"), 16, 16)
save_fig(plot_band(summarise_reps(pt_long, "Ratio"), var_labels, 2,
                   "Relative Effects of Herbicide Over Time", "Ratio to control", ref = 1),
         file.path(fig_dir, "fig_4.1_variability_Pt_ratio.png"), 16, 16)

# ---------------- Grd data (Figure 4.7) ----------------
grd <- fread(file.path(in_dir, "Grd.csv.gz"))
grd_long <- to_long(grd, names(var_labels_grd))
save_fig(plot_band(summarise_reps(grd_long, "Value"), var_labels_grd, 3,
                   "Community-Level Dynamics (Grid Outputs)", "Output Value"),
         file.path(fig_dir, "fig_4.1_variability_Grd_raw.png"), 18, 16)
save_fig(plot_band(summarise_reps(grd_long, "Ratio"), var_labels_grd, 3,
                   "Relative Herbicide Effects on Community Structure", "Ratio to control", ref = 1),
         file.path(fig_dir, "fig_4.1_variability_Grd_ratio.png"), 18, 16)

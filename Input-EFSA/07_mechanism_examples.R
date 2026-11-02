# §4.4 Weekly dynamics of one community: one species and the whole community, all stressors
#
#   Rscript Input-EFSA/07_mechanism_examples.R [in_root] [fig_dir]
#   defaults: outputs/main-2026-10-27/4.4_stressor_effects EFSA_REPORT/figures
# For each case below, two figures (level z05 for the 6 stressors, one column each):
#   fig_4.4_mechanism_<biome>_<MC>_<species>.png  every Pt_ output of the species
#   fig_4.4_mechanism_<biome>_<MC>_community.png  every Grd_ output that varies in time
# Each panel shows the control (black), the null run (grey) and the treatment (red).
# Weekly values (30 weeks per year), years first_year to last_year.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
in_root <- if (length(args) >= 1) args[1] else "outputs/main-2026-10-27/4.4_stressor_effects"
fig_dir <- if (length(args) >= 2) args[2] else "EFSA_REPORT/figures"

cases <- list(
  list(biome = "ALP_I1", MC = 3, species = c("Veronica_polita")),
  list(biome = "ALP_I1", MC = 8, species = c("Scleranthus_annuus", "Agrostis_gigantea")),
  # annual and perennial with nearly the same effect, both abundant in the control (MaxAge)
  list(biome = "MAC_I1", MC = 7, species = c("Lepidium_didymum", "Polypogon_viridis"))
)
level <- "z05"
first_year <- 18
last_year <- 36

STRESSORS <- c("biomass", "SEbiomass", "survival", "establishment", "sterility", "seednumber")
run_labels <- c(NO = "control", null = "null run", treatment = "treatment (z05)")
run_cols <- setNames(c("#1a1a1a", "#9e9e9e", "#d73027"), run_labels)
PT_LABELS <- c(Inds = "Plants (Inds)", seedlings = "Seedlings", seeds = "Seed bank (seeds)",
               repromass = "Reproductive\nmass (mg)", cover = "Cover", shootmass = "Shoot mass (mg)")
GRD_LABELS <- c(NInd = "Plants (NInd)", totMass = "Total mass (mg)", abovemass = "Above-ground\nmass (mg)",
                belowmass = "Below-ground\nmass (mg)", mean_ares = "Mean above-ground\nresource",
                mean_bres = "Mean below-ground\nresource", shannon = "Shannon (Grd)", NPFT = "NPFT (Grd)",
                Cutted = "Cut mass (mg)", mean_generation = "Mean generation", mean_genet_size = "Mean genet size")

# one copy of the control and null run per stressor column, next to the treatment
by_stressor <- function(dt) {
  base <- dt[stressor %in% c("NO", "null")]
  rbind(dt[level == lev & stressor %in% STRESSORS][, `:=`(run = run_labels[["treatment"]], column = stressor)],
        rbindlist(lapply(STRESSORS, function(s) copy(base)[, `:=`(run = run_labels[stressor], column = s)])))
}

plot_runs <- function(long, labels, title, subtitle) {
  long[, `:=`(variable = factor(variable, names(labels), labels), column = factor(column, STRESSORS),
              run = factor(run, run_labels))]
  ggplot(long, aes(Time / 30, value, colour = run)) +
    theme_minimal(base_size = 8) +
    annotate("rect", xmin = 20, xmax = 30, ymin = -Inf, ymax = Inf, fill = "#d73027", alpha = 0.07) +
    geom_line(data = long[run != run_labels[["treatment"]]], linewidth = 0.3) +
    geom_line(data = long[run == run_labels[["treatment"]]], linewidth = 0.35) +
    facet_grid(variable ~ column, scales = "free_y", switch = "y") +
    scale_colour_manual(values = run_cols, name = NULL) +
    scale_x_continuous(breaks = seq(first_year, last_year, 4)) +
    labs(x = "Time (years)", y = NULL, title = title, subtitle = subtitle) +
    theme(legend.position = "top", panel.grid.minor = element_blank(), strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 1), strip.text.x = element_text(face = "bold"),
          panel.spacing = unit(0.5, "lines"))
}

lev <- level
for (cs in cases) {
  d <- file.path(in_root, cs$biome, sprintf("MC%02d", cs$MC))
  com <- fread(file.path(d, "community.csv"))[order(-u)][, rank := .I]
  keep <- function(dt) dt[Time > first_year * 30 & Time <= last_year * 30 &
                            (level == lev | stressor %in% c("NO", "null"))]
  pt <- keep(fread(file.path(d, "Pt.csv.gz")))
  grd <- keep(fread(file.path(d, "Grd.csv.gz")))
  sub <- sprintf("Control (black), null run (grey) and treatment at %s (red); weekly values; shaded: herbicide years 21-30 (one pulse in week 1 of each year)", lev)

  for (sp in cs$species) {
    info <- com[Species == sp]
    x <- pt[PFT == sp]
    # weeks without a row: the species had no plant and no seed
    grid <- unique(pt[, .(level, stressor, Time)])
    x <- x[grid, on = .(level, stressor, Time)]
    for (v in names(PT_LABELS)) set(x, j = v, value = fifelse(is.na(x[[v]]), 0, as.numeric(x[[v]])))
    long <- melt(by_stressor(x), id.vars = c("run", "column", "Time"), measure.vars = names(PT_LABELS),
                 variable.name = "variable", variable.factor = FALSE)
    title <- sprintf("%s community MC%02d: %s (rank %d of 20, MaxAge %d, effect %.3f at %s), all stressors",
                     cs$biome, cs$MC, gsub("_", " ", sp), info$rank, info$MaxAge, info[[lev]], lev)
    ggsave(file.path(fig_dir, sprintf("fig_4.4_mechanism_%s_MC%02d_%s.png", cs$biome, cs$MC, sp)),
           plot_runs(long, PT_LABELS, title, sub), width = 30, height = 20, units = "cm", dpi = 200, bg = "white")
  }

  grd[, (names(GRD_LABELS)) := lapply(.SD, as.numeric), .SDcols = names(GRD_LABELS)]
  vars <- names(GRD_LABELS)[vapply(names(GRD_LABELS), function(v) uniqueN(grd[[v]]) > 1, logical(1))]
  long <- melt(by_stressor(grd), id.vars = c("run", "column", "Time"), measure.vars = vars,
               variable.name = "variable", variable.factor = FALSE)
  title <- sprintf("%s community MC%02d: community outputs (Grd_), all stressors", cs$biome, cs$MC)
  ggsave(file.path(fig_dir, sprintf("fig_4.4_mechanism_%s_MC%02d_community.png", cs$biome, cs$MC)),
         plot_runs(long, GRD_LABELS[vars], title, sub), width = 30, height = 4 + 2.4 * length(vars),
         units = "cm", dpi = 200, bg = "white")
  message(cs$biome, " MC", cs$MC, ": Grd variables plotted: ", paste(vars, collapse = ", "))
}

# Figures of the report §4.4 (A4 page width), from the tables written by 04, 06, 08 and 09
#
#   Rscript Input-EFSA/16_report_4.4_figures.R [in_root] [fig_dir]
#   defaults: outputs/main-2026-10-27/4.4_stressor_effects EFSA_REPORT/figures/4.4
# Run 04_output_analysis.R, 06_composition_shift.R (z05 and z09) and 09_population_at_risk.R
# first. The time axis starts at year 15: before year 21 the treatments are identical to their
# control (same seed), so the ratios are exactly 1 and the first years add nothing.
#   fig_4.4_<endpoint>_survival_by_biome.png   paired ratio, survival, one row per biome
#   fig_4.4_<endpoint>_pooled.png              paired ratio, 90 communities, all stressors
#   fig_4.4_par_<endpoint>_pooled.png          PFT at risk (PaR) over time (NInd, richness)
#   fig_4.4_par_NInd_z01_z04_z09.png           the same, 3 levels and 3 stressors (larger SETAC fig2)
#   fig_4.4_par_vs_level.png                   PFT at risk against the SSD level
#   fig_4.4_composition_time_<level>.png       sensitive / other species, year by year
#   fig_4.4_species_sterility.png              change of each species against its effect (SETAC fig3)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
in_root <- if (length(args) >= 1) args[1] else "outputs/main-2026-10-27/4.4_stressor_effects"
fig_dir <- if (length(args) >= 2) args[2] else "EFSA_REPORT/figures/4.4"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
prov <- jsonlite::read_json(file.path(in_root, "BOR_I1", "MC03", "provenance.json"))
herb_years <- (prov$Tinit + 1):(prov$Tinit + prov$HerbDuration)
x_from <- 15
base <- 9 # pt
STRESSORS <- c("biomass", "SEbiomass", "survival", "establishment", "sterility", "seednumber")
RUN_ORDER <- c("null", STRESSORS)
run_lab <- c(null = "null run\n(no herbicide)", setNames(STRESSORS, STRESSORS))

save_fig <- function(p, name, width, height) {
  ggsave(file.path(fig_dir, name), p, width = width, height = height, units = "cm", dpi = 300, bg = "white")
  message("written: ", name)
}
# herbicide years shaded, x axis from year x_from
herb_layer <- annotate("rect", xmin = min(herb_years) - 0.5, xmax = max(herb_years) + 0.5, ymin = -Inf, ymax = Inf,
                       fill = "#d73027", alpha = 0.10)
x_years <- scale_x_continuous(breaks = c(20, 30, 40, 50), limits = c(x_from, prov$Tmax), expand = c(0, 0.5))
theme_rep <- theme_minimal(base_size = base) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "top", plot.caption = element_text(hjust = 0, size = base - 1.5),
        plot.caption.position = "plot")

# ---- 1. Paired ratios treatment / control, with the null band ----
ENDPOINTS <- c(NInd = "number of individuals", shootmass = "above-ground biomass", richness = "number of PFTs",
               shannon = "Shannon diversity")
null_cap <- paste("Black line: median over the communities; orange: 95 % interval of treatment / control.",
                  "\nGrey: 95 % interval of null run / control (no herbicide, new random stream from year 21).",
                  "\nShaded: herbicide years. Week 20 of each year.")
eff_plot <- function(s, null, rows, cols, ep) {
  ggplot(s, aes(year)) +
    theme_rep +
    herb_layer +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30", linewidth = 0.3) +
    geom_ribbon(data = null, aes(ymin = q025, ymax = q975), fill = "grey50", alpha = 0.35) +
    geom_ribbon(aes(ymin = q025, ymax = q975), fill = "#C75314", alpha = 0.35) +
    geom_line(aes(y = median), linewidth = 0.35) +
    x_years +
    facet_grid(as.formula(paste(rows, "~", cols))) +
    labs(x = "Time (years)", y = paste("Treatment / control:", ENDPOINTS[[ep]]), caption = null_cap) +
    theme(axis.text.x = element_text(size = base - 2.5), axis.text.y = element_text(size = base - 2),
          strip.text.y = element_text(size = base - 1.5), panel.spacing = unit(0.3, "lines"))
}
read_eff <- function(f) {
  e <- fread(file.path(in_root, f))[, year := Time %/% 30 + 1][year >= x_from]
  e
}
bb <- read_eff("effects_by_biome.csv")
pp <- read_eff("effects_pooled.csv")
for (ep in names(ENDPOINTS)) {
  s <- bb[endpoint == ep & stressor == "survival"]
  save_fig(eff_plot(s, bb[endpoint == ep & stressor == "null", -c("level")], "biome", "level", ep),
           sprintf("fig_4.4_%s_survival_by_biome.png", ep), 24, 20)
  s <- pp[endpoint == ep & stressor %in% STRESSORS][, stressor := factor(stressor, STRESSORS)]
  save_fig(eff_plot(s, pp[endpoint == ep & stressor == "null", -c("level", "stressor")], "level", "stressor", ep),
           sprintf("fig_4.4_%s_pooled.png", ep), 17, 22)
}

# ---- 2. PFT at risk (PaR, after the Population at Risk of Baudrot et al. 2021) over time ----
GROUPS <- c("5 % PaR", "10 % PaR", "20 % PaR", "50 % PaR", "Rest (50 % least sensitive)", "Community")
group_cols <- setNames(c("#67000d", "#cb181d", "#fb6a4a", "#fcae91", "#4575b4", "#1a1a1a"), GROUPS)
PAR_EP <- c(NInd = "number of individuals", richness = "number of PFTs present", shootmass = "above-ground biomass")
par_cap <- paste("90 communities of 20 PFTs, all biomes pooled. x % PaR (PFT at risk): the x % most sensitive PFTs of each",
                 "community\n(5 %: 1 PFT, 10 %: 2, 20 %: 4, 50 %: 10). Ratio of the sums over the communities;",
                 "band: 95 % bootstrap interval\n(communities resampled). Shaded: herbicide years. Week 20 of each year.")
pr <- fread(file.path(in_root, "population_at_risk_ratio.csv"))[year >= x_from]
LEVELS <- sort(setdiff(unique(pr$level), "null"))
null_rows <- function(s) rbindlist(lapply(LEVELS, function(l) copy(s[stressor == "null"])[, level := l]))
pr <- rbind(pr[stressor != "null"], null_rows(pr))
pr[, `:=`(stressor = factor(stressor, RUN_ORDER), group = factor(group, GROUPS))]
par_plot <- function(s, ep, ncol_legend = 6) {
  ggplot(s, aes(year, colour = group, fill = group)) +
    theme_rep +
    herb_layer +
    geom_hline(yintercept = 1, colour = "grey40", linewidth = 0.3) +
    geom_ribbon(aes(ymin = lo, ymax = hi), colour = NA, alpha = 0.15) +
    geom_line(aes(y = ratio), linewidth = 0.4) +
    x_years +
    scale_colour_manual(values = group_cols, name = NULL) +
    scale_fill_manual(values = group_cols, name = NULL) +
    guides(colour = guide_legend(ncol = ncol_legend), fill = guide_legend(ncol = ncol_legend)) +
    labs(x = "Time (years)", y = paste("Treatment / control:", PAR_EP[[ep]]), caption = par_cap) +
    theme(axis.text.x = element_text(size = base - 2.5), axis.text.y = element_text(size = base - 2),
          panel.spacing = unit(0.3, "lines"), legend.text = element_text(size = base - 1))
}
for (ep in c("NInd", "richness"))
  save_fig(par_plot(pr[endpoint == ep], ep) + facet_grid(level ~ stressor, labeller = labeller(stressor = run_lab)),
           sprintf("fig_4.4_par_%s_pooled.png", ep), 17, 23)

# larger version of the SETAC fig2: 3 stressors, 3 levels
s <- pr[endpoint == "NInd" & stressor %in% c("biomass", "survival", "establishment") & level %in% c("z01", "z04", "z09")]
save_fig(par_plot(droplevels(s), "NInd", 3) + facet_grid(level ~ stressor), "fig_4.4_par_NInd_z01_z04_z09.png", 16, 13)

# ---- 3. PFT at risk against the SSD level (herbicide years pooled) ----
# P4 of the review of Ch. 3-4 (make_proposed_figures.R), table P4_par_vs_level.csv
p4 <- fread("EFSA_REPORT/revision/drafts/REVIEW_Ch3-4_20261030_figs/P4_par_vs_level.csv")
G4 <- c("5 % PaR (1 species)", "20 % PaR (4 species)", "Rest (10 least sensitive)", "Community (20 species)") # names in the csv
G4_lab <- setNames(c("5 % PaR (1 PFT)", "20 % PaR (4 PFTs)", "Rest (10 least sensitive PFTs)", "Community (20 PFTs)"), G4)
gcol <- setNames(c("#b2182b", "#e8743b", "#2c6fbb", "#1a1a1a"), G4)
null4 <- p4[stressor == "null" & group == G4[1]]
p4 <- p4[stressor %in% STRESSORS][, `:=`(z = as.numeric(sub("z0", "", level)) / 10,
                                         stressor = factor(stressor, STRESSORS), group = factor(group, G4))]
p <- ggplot(p4, aes(z, r, colour = group, fill = group)) +
  theme_rep +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = null4$lo, ymax = null4$hi, fill = "#e6e6e6") +
  geom_hline(yintercept = 1, colour = "grey40", linewidth = 0.3) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.6) + geom_point(size = 1) +
  facet_wrap(~ stressor, nrow = 2) +
  scale_colour_manual(values = gcol, labels = G4_lab, name = NULL) + scale_fill_manual(values = gcol, labels = G4_lab, name = NULL) +
  scale_x_continuous(breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) +
  guides(colour = guide_legend(nrow = 2), fill = guide_legend(nrow = 2)) +
  labs(x = "SSD level z (95th percentile of the PFT effects)",
       y = "Treatment / control: number of individuals\n(week 20, herbicide years 21-30)",
       caption = paste("90 communities pooled; ratio of the sums over the communities, 95 % bootstrap interval",
                       "(500 resamples).\nGrey: 95 % interval of the null run for the 5 % PaR, the reference for a zero effect.")) +
  theme(panel.spacing = unit(0.6, "lines"))
save_fig(p, "fig_4.4_par_vs_level.png", 16, 12)

# ---- 4. Sensitive and other species, year by year ----
CGROUPS <- c("Sensitive species", "Other species", "Community total")
ccols <- setNames(c("#d73027", "#4575b4", "#1a1a1a"), CGROUPS)
mean_ci <- function(v) list(mean = mean(v), lo = mean(v) - 1.96 * sd(v) / sqrt(length(v)),
                            hi = mean(v) + 1.96 * sd(v) / sqrt(length(v)))
for (lev in c("z05", "z09")) {
  f <- file.path(in_root, sprintf("composition_shift_%s_yearly.csv", lev))
  if (!file.exists(f)) next
  y <- fread(f)[year >= x_from, mean_ci(r), by = .(year, stressor, group)]
  y[, `:=`(stressor = factor(stressor, RUN_ORDER), group = factor(group, CGROUPS))]
  p <- ggplot(y, aes(year, mean, colour = group, fill = group)) +
    theme_rep +
    herb_layer +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_ribbon(aes(ymin = lo, ymax = hi), colour = NA, alpha = 0.2) +
    geom_line(linewidth = 0.45) +
    x_years +
    facet_wrap(~ stressor, nrow = 2, labeller = as_labeller(run_lab)) +
    scale_colour_manual(values = ccols, name = NULL) + scale_fill_manual(values = ccols, name = NULL) +
    labs(x = "Time (years)", y = "Change in individuals vs control\n(% of the group's own individuals)",
         caption = sprintf(paste("90 communities; sensitive species: effect >= 0.2 at %s. Mean and 95 %% confidence",
                                 "interval of the mean over the communities.\nShaded: herbicide years. Week 20 of each year."), lev))
  save_fig(p, sprintf("fig_4.4_composition_time_%s.png", lev), 17, 11)
}

# ---- 5. Shoot mass of each species against its effect, sterility, z05 (SETAC fig3 without the purple bands) ----
cache <- file.path(in_root, "species_herbyears_z05.csv")
stopifnot(file.exists(cache)) # written by 11_setac2027_figures.R
x <- fread(cache)
x <- x[x[stressor == "NO", .(biome, MC, PFT, C = shootmass)], on = .(biome, MC, PFT)][stressor != "NO"]
x[, d := 100 * (shootmass - C) / sum(C), by = .(biome, MC, stressor)]
x[, rank := frank(-C, ties.method = "first"), by = .(biome, MC, stressor)]
RANKS_E <- c("1-3", "4-10", "11-20")
s3 <- x[stressor == "sterility"][, dom := cut(rank, c(0, 3, 10, 20), RANKS_E)]
green_cols <- setNames(c("#00441b", "#41ab5d", "#a1d99b"), RANKS_E)
x_breaks <- c(0.01, 0.05, 0.1, 0.2, 0.5, 1)
p <- ggplot(s3[order(-rank)], aes(effect, d, colour = dom)) +
  theme_rep +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_point(size = 0.7, alpha = 0.5) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.9, span = 0.5, key_glyph = "path") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.5, linetype = "22", key_glyph = "path") +
  scale_colour_manual(values = green_cols, name = "Rank of the species in the control (1 = largest shoot mass)") +
  scale_x_sqrt(breaks = x_breaks, labels = format(x_breaks, drop0trailing = TRUE), expand = c(0, 0.01)) +
  coord_cartesian(ylim = quantile(s3$d, c(0.002, 0.998))) +
  guides(colour = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1,
                               override.aes = list(linewidth = 0.9, alpha = 1))) +
  labs(x = "Effect of the species at z05 (square-root scale)",
       y = "Change in shoot mass of the species vs control\n(% of the shoot mass of the control community)",
       caption = paste("Sterility, z05. One point per species and community (90 communities x 20 species);",
                       "\nshoot mass in week 20 summed over the herbicide years 21-30. Solid: local mean (loess);",
                       "\ndashed: linear regression, by rank class. y axis cut at the 0.2 and 99.8 percentiles."))
save_fig(p, "fig_4.4_species_sterility.png", 16, 12)

# ---- 6. Table of effect classes (TABLE 4.1), as in Table 12 of the previous version of the report ----
# Herbicide years 21-30, mean of the yearly values, 90 communities pooled. Reference: null run.
#   mid: median over the communities (5 % PaR: ratio of the sums); lo: 2.5th percentile (5 % PaR: bootstrap)
#   Class 5: mid below the lower bound of the null run
#   Class 4: mid at least 0.05 below the mid of the null run
#   Class 3: lo at least 0.04 below the lower bound of the null run
#   Class 2: lo 0.02 to 0.04 below the lower bound of the null run
#   Class 1: otherwise
# The biomes are not classified one by one: with 10 communities, the classes of a biome follow the noise.
cls_herb <- function(s) s[year %in% herb_years, .(mid = mean(mid), lo = mean(lo)), by = .(endpoint, stressor, level)]
ce <- fread(file.path(in_root, "effects_pooled.csv"))[, year := Time %/% 30 + 1][endpoint %in% names(ENDPOINTS)]
ce <- cls_herb(ce[, .(endpoint, stressor, level, year, mid = median, lo = q025)])
cp <- fread(file.path(in_root, "population_at_risk_ratio.csv"))[endpoint == "NInd" & group == "5 % PaR"]
cp <- cls_herb(cp[, .(endpoint = "PaR5_NInd", stressor, level, year, mid = ratio, lo)])
ca <- rbind(ce, cp)
ca <- ca[stressor != "null"][ca[stressor == "null", .(endpoint, m0 = mid, lo0 = lo)], on = "endpoint"]
ca[, class := fifelse(mid < lo0, 5L, fifelse(m0 - mid >= 0.05, 4L,
                fifelse(lo0 - lo >= 0.04, 3L, fifelse(lo0 - lo >= 0.02, 2L, 1L))))]
tab <- dcast(ca, endpoint + stressor ~ level, value.var = "class")
tab[, `:=`(endpoint = factor(endpoint, c("shannon", "richness", "shootmass", "NInd", "PaR5_NInd")),
           stressor = factor(stressor, c("biomass", "SEbiomass", "establishment", "seednumber", "sterility", "survival")))]
setorder(tab, endpoint, stressor)
fwrite(tab, file.path(fig_dir, "table_4.4_effect_classes.csv"))
fwrite(ca[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 4) else v)],
       file.path(fig_dir, "table_4.4_effect_classes_values.csv"))
message("written: table_4.4_effect_classes.csv")
print(tab)

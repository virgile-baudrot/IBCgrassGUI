# §4.4 Response of each species against its effect, for several Pt_ outputs
#
#   Rscript Input-EFSA/08_species_response.R [in_root] [fig_dir] [level]
#   defaults: outputs/main-2026-10-27/4.4_stressor_effects EFSA_REPORT/figures z05
# One point per species and community (90 x 20 = 1 800 per panel), one panel per run type
# (null run and the 6 stressors at the level). For variable v:
#   x: effect of the species at the level (its SSD draw)
#   y: (v of the species in the treatment - v of the species in the control)
#      / (v summed over all species of the control community) x 100,
#      week 20 of each herbicide year, then mean over the herbicide years.
# y is a share of the community, so the 20 species of a community add up to the change of
# the community total of v. Figures: fig_4.4_composition_species_<level>[_<variable>].png
# (no suffix for Inds, as before).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
in_root <- if (length(args) >= 1) args[1] else "outputs/main-2026-10-27/4.4_stressor_effects"
fig_dir <- if (length(args) >= 2) args[2] else "EFSA_REPORT/figures"
lev <- if (length(args) >= 3) args[3] else "z05"
threshold <- 0.2
week <- 20
RUN_ORDER <- c("null", "biomass", "SEbiomass", "survival", "establishment", "sterility", "seednumber")
VARS <- c(Inds = "individuals", cover = "cover", shootmass = "shoot mass")

dirs <- dirname(list.files(in_root, "^Pt\\.csv\\.gz$", recursive = TRUE, full.names = TRUE))
prov <- jsonlite::read_json(file.path(dirs[1], "provenance.json"))
herb_years <- (prov$Tinit + 1):(prov$Tinit + prov$HerbDuration)

sp <- rbindlist(lapply(dirs, function(d) {
  com <- fread(file.path(d, "community.csv"))
  pt <- fread(file.path(d, "Pt.csv.gz"), select = c("level", "stressor", "Time", "PFT", names(VARS)))
  pt <- pt[Time %% 30 == week & ceiling(Time / 30) %in% herb_years & (level == lev | stressor %in% c("NO", "null"))]
  full <- CJ(stressor = unique(pt$stressor), Time = unique(pt$Time), PFT = com$Species)
  pt <- pt[, lapply(.SD, sum), by = .(stressor, Time, PFT), .SDcols = names(VARS)][full, on = .(stressor, Time, PFT)]
  for (v in names(VARS)) set(pt, j = v, value = fifelse(is.na(pt[[v]]), 0, as.numeric(pt[[v]])))
  long <- melt(pt, id.vars = c("stressor", "Time", "PFT"), variable.name = "variable", variable.factor = FALSE)
  ctl <- long[stressor == "NO", .(Time, PFT, variable, ctl = value)]
  long <- ctl[long[stressor != "NO"], on = .(Time, PFT, variable)]
  long[, total_ctl := sum(ctl), by = .(stressor, Time, variable)]
  long[, .(change = mean(100 * (value - ctl) / total_ctl)), by = .(stressor, PFT, variable)][
    data.table(PFT = com$Species, effect = com[[lev]]), on = "PFT"][
    , `:=`(biome = basename(dirname(d)), MC = basename(d))]
}))
sp[, stressor := factor(stressor, RUN_ORDER)]
fwrite(sp, file.path(in_root, sprintf("species_change_vs_effect_%s.csv", lev)))

for (v in names(VARS)) {
  x <- sp[variable == v]
  p <- ggplot(x, aes(effect, change)) +
    theme_minimal(base_size = 9) +
    annotate("rect", xmin = threshold, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#d73027", alpha = 0.06) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_point(size = 0.6, alpha = 0.35, colour = "grey30") +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "#d73027", linewidth = 0.6, span = 0.5) +
    scale_x_log10(breaks = c(0.001, 0.01, 0.1, 0.2, 0.5, 1), labels = c("0.001", "0.01", "0.1", "0.2", "0.5", "1")) +
    coord_cartesian(ylim = quantile(x$change, c(0.002, 0.998))) +
    facet_wrap(~ stressor, nrow = 2, labeller = as_labeller(c(null = "null run (no herbicide; x: effect it would have)",
                                                               setNames(RUN_ORDER[-1], RUN_ORDER[-1])))) +
    labs(x = sprintf("Species effect at %s (log scale; shaded: sensitive, >= %.1f)", lev, threshold),
         y = sprintf("Change in %s of the species vs control\n(%% of the %s of the control community)", VARS[[v]], VARS[[v]]),
         title = sprintf("%s: %s of each species of the 90 communities, mean over the herbicide years", lev, VARS[[v]]),
         subtitle = "One point per species and community. Red line: local mean (loess)") +
    theme(panel.grid.minor = element_blank())
  suffix <- if (v == "Inds") "" else paste0("_", v)
  ggsave(file.path(fig_dir, sprintf("fig_4.4_composition_species_%s%s.png", lev, suffix)), p,
         width = 24, height = 14, units = "cm", dpi = 200, bg = "white")
}

# group sums (sensitive / other / total) per variable, to compare with the individuals
grp <- sp[, .(change = sum(change)), by = .(biome, MC, stressor, variable, sens = effect >= threshold)]
grp <- rbind(grp[, .(biome, MC, stressor, variable, group = fifelse(sens, "sensitive", "other"), change)],
             grp[, .(change = sum(change), group = "total"), by = .(biome, MC, stressor, variable)])
print(dcast(grp[, .(m = round(mean(change), 2)), by = .(variable, stressor, group)],
            variable + stressor ~ group, value.var = "m"))

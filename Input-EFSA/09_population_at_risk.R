# §4.4 Effects on the x % Population at Risk (PaR), all biomes pooled
#
#   Rscript Input-EFSA/09_population_at_risk.R [in_root] [fig_dir]
#   defaults: outputs/main-2026-10-27/4.4_stressor_effects EFSA_REPORT/figures
# Adapted from the Average Population at Risk of Baudrot et al. (2021, Ecotox. Environ.
# Saf. 207:111215): the x % Population at Risk is the x % of the population that suffers
# the most, and its effect is averaged over that part only. Here the population is the
# species of a community: species are ranked by their sensitivity u (the effect at every
# level is proportional to the level, so the rank is the same at all levels) and the x %
# PaR is the x % most sensitive species, as for the HC5 of a species sensitivity
# distribution. With 20 species per community: PaR5 = 1 species, PaR10 = 2, PaR20 = 4,
# PaR50 = 10. The groups are nested; "Rest" is the 10 least sensitive species.
# The most sensitive species is often rare (median 0.6 individuals per community in the
# control, week 20), so a treatment / control ratio community by community is mostly 0/0
# or 1/0. The ratio is taken over all at-risk individuals of the 90 communities, i.e.
# sum(treatment) / sum(control), with a 95 % interval from a bootstrap of the communities
# (the same resamples for every group, run and year). The null run (no herbicide, new
# random stream from year Tinit+1) is shown as a first column: its ratio is what a zero
# effect gives.
# Endpoints (week 20): number of individuals, shoot mass and number of species present in
# the group (richness).
# Figures:
#   fig_4.4_par_<endpoint>_pooled.png           treatment / control
#   fig_4.4_par_<endpoint>_absolute_pooled.png  treatment and control, mean per community

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
in_root <- if (length(args) >= 1) args[1] else "outputs/main-2026-10-27/4.4_stressor_effects"
fig_dir <- if (length(args) >= 2) args[2] else "EFSA_REPORT/figures"
week <- 20
n_boot <- 1000
set.seed(2021)
RUN_ORDER <- c("null", "biomass", "SEbiomass", "survival", "establishment", "sterility", "seednumber")
PAR <- c(PaR5 = 0.05, PaR10 = 0.10, PaR20 = 0.20, PaR50 = 0.50)
GROUPS <- c("5 % PaR", "10 % PaR", "20 % PaR", "50 % PaR", "Rest (50 % least sensitive)", "Community")
group_cols <- setNames(c("#67000d", "#cb181d", "#fb6a4a", "#fcae91", "#4575b4", "#1a1a1a"), GROUPS)
ENDPOINTS <- c(NInd = "number of individuals", shootmass = "above-ground biomass (mg)",
               richness = "number of species present")

dirs <- dirname(list.files(in_root, "^Pt\\.csv\\.gz$", recursive = TRUE, full.names = TRUE))
prov <- jsonlite::read_json(file.path(dirs[1], "provenance.json"))
herb_years <- (prov$Tinit + 1):(prov$Tinit + prov$HerbDuration)

# endpoints per group, community, run and year; absent species counted as 0.
# sp: individuals of each species summed over the herbicide years (for the curve x -> R_x)
res <- lapply(dirs, function(d) {
  com <- fread(file.path(d, "community.csv"))
  n_sp <- nrow(com)
  rank <- setNames(frank(-com$u, ties.method = "first"), com$Species)
  pt <- fread(file.path(d, "Pt.csv.gz"), select = c("level", "stressor", "Time", "PFT", "Inds", "shootmass"))
  pt <- pt[Time %% 30 == week, .(Inds = sum(Inds), shootmass = sum(shootmass)), by = .(level, stressor, Time, PFT)]
  full <- CJ(run = unique(paste(pt$level, pt$stressor)), Time = unique(pt$Time), PFT = com$Species)
  full[, c("level", "stressor") := tstrsplit(run, " ")][, run := NULL]
  pt <- pt[full, on = .(level, stressor, Time, PFT)]
  pt[is.na(Inds), `:=`(Inds = 0, shootmass = 0)]
  pt[, rk := rank[PFT]]
  members <- c(lapply(PAR, function(p) seq_len(round(p * n_sp))), list((n_sp / 2 + 1):n_sp, seq_len(n_sp)))
  out <- rbindlist(Map(function(g, m) pt[rk %in% m, .(NInd = sum(Inds), shootmass = sum(shootmass),
                                                      richness = sum(Inds > 0)),
                                         by = .(level, stressor, Time)][, group := g],
                       GROUPS, members))
  sp <- pt[ceiling(Time / 30) %in% herb_years, .(N = sum(Inds)), by = .(level, stressor, PFT)]
  id <- list(biome = basename(dirname(d)), MC = as.integer(sub("MC", "", basename(d))))
  list(grp = out[, names(id) := id], sp = sp[, names(id) := id])
})
grp <- rbindlist(lapply(res, `[[`, "grp"))
sp <- rbindlist(lapply(res, `[[`, "sp"))
rm(res)
grp[, year := ceiling(Time / 30)][, Time := NULL]
long <- melt(grp, measure.vars = names(ENDPOINTS), variable.name = "endpoint", variable.factor = FALSE)
ctl <- long[stressor == "NO", .(biome, MC, year, group, endpoint, control = value)]
long <- ctl[long[stressor != "NO"], on = .(biome, MC, year, group, endpoint)]
long[stressor == "null", level := "null"]
fwrite(long, file.path(in_root, "population_at_risk.csv.gz"))

# pooled ratio and bootstrap over communities: one row per community, one column per cell
coms <- unique(long[, .(biome, MC)])[, id := .I]
long <- coms[long, on = .(biome, MC)]
W <- t(rmultinom(n_boot, nrow(coms), rep(1, nrow(coms)))) # n_boot x communities
boot_ratio <- function(x) {
  x[, cell := .GRP, by = .(endpoint, level, stressor, year, group)]
  cells <- unique(x[, .(cell, endpoint, level, stressor, year, group)])
  tr <- as.matrix(dcast(x, id ~ cell, value.var = "value")[order(id)][, -1])
  ct <- as.matrix(dcast(x, id ~ cell, value.var = "control")[order(id)][, -1])
  q <- function(m, p) apply(m, 2, quantile, p, na.rm = TRUE, names = FALSE)
  bt <- W %*% tr
  bc <- W %*% ct
  n <- nrow(tr)
  cells[, `:=`(ratio = colSums(tr) / colSums(ct), lo = q(bt / bc, 0.025), hi = q(bt / bc, 0.975),
               treatment = colMeans(tr), treatment_lo = q(bt / n, 0.025), treatment_hi = q(bt / n, 0.975),
               control = colMeans(ct), control_lo = q(bc / n, 0.025), control_hi = q(bc / n, 0.975))]
}
cells <- rbindlist(lapply(split(long, by = "endpoint"), boot_ratio))
ratio <- cells[, .(endpoint, level, stressor, year, group, ratio, lo, hi)]
# mean per community and 95 % bootstrap interval of the mean, treatment and control in long format
absolute <- rbind(cells[, .(endpoint, level, stressor, year, group, run = "treatment",
                            value = treatment, lo = treatment_lo, hi = treatment_hi)],
                  cells[, .(endpoint, level, stressor, year, group, run = "control",
                            value = control, lo = control_lo, hi = control_hi)])
fwrite(ratio, file.path(in_root, "population_at_risk_ratio.csv"))

# the null run is drawn as a first column, the same in every row
LEVELS <- sort(setdiff(unique(long$level), "null"))
null_rows <- function(s) rbindlist(lapply(LEVELS, function(l) copy(s[stressor == "null"])[, level := l]))
prep <- function(s) {
  s <- rbind(s[stressor != "null"], null_rows(s))
  s[, `:=`(stressor = factor(stressor, RUN_ORDER), group = factor(group, GROUPS))]
}
ratio <- prep(ratio)
absolute <- prep(absolute)
facet_lab <- labeller(stressor = c(null = "null run\n(no herbicide)", setNames(RUN_ORDER[-1], RUN_ORDER[-1])))
note <- sprintf("%d communities of %d species, all biomes pooled. x %% PaR: the x %% most sensitive species of each community (5 %%: 1 species, 10 %%: 2, 20 %%: 4, 50 %%: 10). Week %d.",
                nrow(coms), prov$n_species, week)
# ymin: 0 on a log scale, where -Inf cannot be transformed and the window would be dropped
base_plot <- function(s, ymin = -Inf) {
  ggplot(s, aes(year, colour = group, fill = group)) +
    theme_minimal(base_size = 8) +
    annotate("rect", xmin = prov$Tinit + 0.5, xmax = prov$Tinit + prov$HerbDuration + 0.5,
             ymin = ymin, ymax = Inf, fill = "#d73027", alpha = 0.07) +
    facet_grid(level ~ stressor, labeller = facet_lab) +
    scale_colour_manual(values = group_cols, name = NULL) +
    scale_fill_manual(values = group_cols, name = NULL) +
    guides(colour = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
    theme(legend.position = "top", panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
}

for (e in names(ENDPOINTS)) {
  p <- base_plot(ratio[endpoint == e]) +
    geom_hline(yintercept = 1, colour = "grey40", linewidth = 0.3) +
    geom_ribbon(aes(ymin = lo, ymax = hi), colour = NA, alpha = 0.15) +
    geom_line(aes(y = ratio), linewidth = 0.4) +
    labs(x = "Time (years)", y = paste("Treatment / control:", ENDPOINTS[[e]]),
         title = sprintf("Population at Risk: treatment / control, %s", ENDPOINTS[[e]]),
         subtitle = "Shaded: herbicide years. Ratio of the sums over the communities; band: 95 % bootstrap interval (communities resampled)",
         caption = note)
  ggsave(file.path(fig_dir, sprintf("fig_4.4_par_%s_pooled.png", e)), p,
         width = 24, height = 30, units = "cm", dpi = 200, bg = "white")

  a <- absolute[endpoint == e][, run := factor(run, c("treatment", "control"))]
  p <- base_plot(a, ymin = if (e != "richness") 0 else -Inf) +
    geom_ribbon(aes(ymin = lo, ymax = hi, group = interaction(group, run)), colour = NA, alpha = 0.15) +
    geom_line(aes(y = value, linetype = run, linewidth = run)) +
    scale_linetype_manual(values = c(treatment = "solid", control = "22"), name = NULL) +
    scale_linewidth_manual(values = c(treatment = 0.45, control = 0.3), name = NULL) +
    (if (e != "richness") scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale()))) +
    labs(x = "Time (years)", y = paste("Mean per community:", ENDPOINTS[[e]], if (e != "richness") "(log scale)"),
         title = sprintf("Population at Risk: absolute values, %s", ENDPOINTS[[e]]),
         subtitle = "Shaded: herbicide years. Solid: treatment (or null run); dashed: control. Mean over the communities; band: 95 % bootstrap interval of the mean",
         caption = note)
  ggsave(file.path(fig_dir, sprintf("fig_4.4_par_%s_absolute_pooled.png", e)), p,
         width = 24, height = 30, units = "cm", dpi = 200, bg = "white")
  message("written: fig_4.4_par_", e)
}

# ratio over the herbicide years, z05 and z09
print(dcast(ratio[endpoint == "NInd" & year %in% herb_years & level %in% c("z05", "z09"),
                  .(min_ratio = round(min(ratio), 2)), by = .(level, stressor, group)],
            level + stressor ~ group, value.var = "min_ratio"))

# ---- the whole curve x -> R_x (Shorrocks) ----
# Unit: a species of a community (all communities pooled). Outcome y: its individuals summed
# over the herbicide years, treatment / control; species without individual in the control
# over these years are left out (same species for every run). Unlike the PaR groups above,
# the species are ranked by their own outcome, separately in each run, as the theorem
# requires: R_x is the mean of y over the x % of the species with the lowest y. If
# R_x(herbicide) <= R_x(null run) for every x, the null run dominates at the second order.
ctl_sp <- sp[stressor == "NO", .(biome, MC, PFT, N_ctl = N)]
y <- ctl_sp[sp[stressor != "NO"], on = .(biome, MC, PFT)][N_ctl > 0]
y[, `:=`(y = N / N_ctl, run = fifelse(stressor == "null", "null", paste(level, stressor)))]
x_grid <- exp(seq(log(0.005), 0, length.out = 200))
r_curve <- function(v) {
  v <- sort(v)
  cumsum(v)[ceiling(x_grid * length(v))] / ceiling(x_grid * length(v))
}
curves <- y[, .(x = x_grid, R = r_curve(y)), by = .(run, level, stressor)]
# 95 % band of the null curve: communities resampled
null_y <- y[stressor == "null"]
null_y[, cid := .GRP, by = .(biome, MC)]
by_com <- split(null_y$y, null_y$cid)
boot <- sapply(seq_len(200), function(i) r_curve(unlist(by_com[sample(length(by_com), replace = TRUE)])))
null_band <- data.table(x = x_grid, lo = apply(boot, 1, quantile, 0.025), hi = apply(boot, 1, quantile, 0.975))
null_curve <- curves[stressor == "null", .(x, R)]
fwrite(curves, file.path(in_root, "population_at_risk_curve.csv"))

stressor_cols <- c(biomass = "#E69F00", SEbiomass = "#CC79A7", survival = "#D55E00",
                   establishment = "#009E73", sterility = "#0072B2", seednumber = "#56B4E9")
sc <- curves[stressor != "null"][, stressor := factor(stressor, names(stressor_cols))]
p <- ggplot(sc, aes(x, R)) +
  theme_minimal(base_size = 9) +
  geom_vline(xintercept = PAR, linetype = "dotted", colour = "grey60", linewidth = 0.3) +
  geom_ribbon(data = null_band, aes(x, ymin = lo, ymax = hi), inherit.aes = FALSE, fill = "grey50", alpha = 0.3) +
  geom_line(data = null_curve, colour = "grey30", linewidth = 0.5, linetype = "22") +
  geom_line(aes(colour = stressor), linewidth = 0.5) +
  facet_wrap(~ level, ncol = 3) +
  scale_x_log10(breaks = c(0.01, PAR, 1), labels = scales::label_percent()) +
  scale_colour_manual(values = stressor_cols, name = NULL) +
  labs(x = "x: share of the species with the lowest outcome (log scale)",
       y = expression(R[x] * ": mean treatment / control over these species"),
       title = expression("The whole curve " * x %->% R[x] * ": individuals, herbicide years pooled"),
       subtitle = "Grey dashed: null run (no herbicide) and 95 % interval (communities resampled). A curve below the null run for every x: the null run is preferred by every risk-averse decision maker",
       caption = sprintf("%d species of %d communities (species without individual in the control over years %d-%d left out). Dotted: x = 5, 10, 20, 50 %%. Week %d.",
                         uniqueN(y[, .(biome, MC, PFT)]), nrow(coms), min(herb_years), max(herb_years), week)) +
  theme(legend.position = "top", panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
        plot.subtitle = element_text(size = 7.5))
ggsave(file.path(fig_dir, "fig_4.4_par_curve_NInd.png"), p, width = 24, height = 20, units = "cm", dpi = 200, bg = "white")

# dominance check: share of x where the herbicide curve lies at or below the null run
dom <- curves[stressor != "null"][null_curve, on = "x", R_null := i.R][, .(below = mean(R <= R_null + 1e-12)), by = .(level, stressor)]
print(dcast(dom, level ~ stressor, value.var = "below"))

# the same, as the difference to the null run: dominance reads as a curve below 0 for every x
sc[null_curve, on = "x", delta := R - i.R]
band <- null_band[null_curve, on = "x"][, `:=`(lo = lo - R, hi = hi - R)]
p <- ggplot(sc, aes(x, delta)) +
  theme_minimal(base_size = 9) +
  geom_vline(xintercept = PAR, linetype = "dotted", colour = "grey60", linewidth = 0.3) +
  geom_ribbon(data = band, aes(x, ymin = lo, ymax = hi), inherit.aes = FALSE, fill = "grey50", alpha = 0.3) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.3) +
  geom_line(aes(colour = stressor), linewidth = 0.5) +
  facet_wrap(~ level, ncol = 3) +
  scale_x_log10(breaks = c(0.01, PAR, 1), labels = scales::label_percent()) +
  scale_colour_manual(values = stressor_cols, name = NULL) +
  labs(x = "x: share of the species with the lowest outcome (log scale)",
       y = expression(R[x] * "(herbicide) - " * R[x] * "(null run)"),
       title = expression("The whole curve " * x %->% R[x] * ", difference to the null run: individuals, herbicide years pooled"),
       subtitle = "Grey: 95 % interval of the null run (communities resampled). A curve below 0 for every x: the null run dominates at the second order",
       caption = sprintf("Below x of about 8 %% the species are absent in every run (R_x = 0). Dotted: x = 5, 10, 20, 50 %%. Week %d.", week)) +
  theme(legend.position = "top", panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
        plot.subtitle = element_text(size = 7.5))
ggsave(file.path(fig_dir, "fig_4.4_par_curve_NInd_difference.png"), p, width = 24, height = 20, units = "cm", dpi = 200, bg = "white")

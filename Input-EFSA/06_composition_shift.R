# §4.4 The herbicide changes the composition of the community, not its size
#
#   Rscript Input-EFSA/06_composition_shift.R [in_root] [fig_dir] [level]
#   defaults: outputs/main-2026-10-27/4.4_stressor_effects EFSA_REPORT/figures z05
# Species of each community are split into "sensitive" (effect >= 0.2 at the level) and
# "other" species. For each group and for the whole community, the relative change of the
# number of individuals (week 20) between a run and the control is
#   r = 100 * (N_run - N_control) / ((N_run + N_control) / 2)
# i.e. a % of the group's own individuals, with a symmetric denominator (no bias between
# run and control). The total is then a weighted mean of the two groups and lies between
# them. The null run (no herbicide, new random stream from year Tinit+1) is shown as a
# reference column: its r is what a zero effect gives. It is not subtracted community by
# community: the null run is one more random realisation, and subtracting it adds its own
# noise and chance offset to every stressor (checked: +5 % for the sensitive species under
# biomass and SEbiomass). Means over communities with 95 % confidence intervals of the
# mean (normal approximation).
# Figures:
#   fig_4.4_composition_decomposition_<level>.png           herbicide years pooled
#   fig_4.4_composition_decomposition_<level>_by_biome.png  the same, one panel per biome
#   fig_4.4_composition_time_<level>.png                    year by year

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
GROUPS <- c("Sensitive species", "Other species", "Community total")
group_cols <- setNames(c("#d73027", "#4575b4", "#1a1a1a"), GROUPS)

dirs <- dirname(list.files(in_root, "^Pt\\.csv\\.gz$", recursive = TRUE, full.names = TRUE))
prov <- jsonlite::read_json(file.path(dirs[1], "provenance.json"))
herb_years <- (prov$Tinit + 1):(prov$Tinit + prov$HerbDuration)

# individuals per species, run and year (week 20); absent species counted as 0
species <- rbindlist(lapply(dirs, function(d) {
  com <- fread(file.path(d, "community.csv"))
  pt <- fread(file.path(d, "Pt.csv.gz"), select = c("biome", "MC", "level", "stressor", "Time", "PFT", "Inds"))
  pt <- pt[Time %% 30 == week & (level == lev | stressor %in% c("NO", "null"))]
  full <- CJ(stressor = unique(pt$stressor), Time = unique(pt$Time), PFT = com$Species)
  pt <- pt[, .(Inds = sum(Inds)), by = .(stressor, Time, PFT)][full, on = .(stressor, Time, PFT)]
  pt[is.na(Inds), Inds := 0]
  pt[data.table(PFT = com$Species, effect = com[[lev]]), effect := i.effect, on = "PFT"]
  pt[, `:=`(biome = basename(dirname(d)), MC = as.integer(sub("MC", "", basename(d))), year = ceiling(Time / 30))]
}))
ctl <- species[stressor == "NO", .(biome, MC, year, PFT, Inds_ctl = Inds)]
species <- ctl[species[stressor != "NO"], on = .(biome, MC, year, PFT)]
species[, total_ctl := sum(Inds_ctl), by = .(biome, MC, year, stressor)]
species[, group := fifelse(effect >= threshold, GROUPS[1], GROUPS[2])]

# individuals per group (and total), community, run and year
grp <- species[, .(N = sum(Inds), C = sum(Inds_ctl)), by = .(biome, MC, year, stressor, group)]
grp <- rbind(grp, grp[, .(N = sum(N), C = sum(C), group = GROUPS[3]), by = .(biome, MC, year, stressor)])
grp[, group := factor(group, GROUPS)]
rel <- function(n, c) fifelse(n + c > 0, 100 * (n - c) / ((n + c) / 2), 0)
mean_ci <- function(v) list(mean = mean(v), lo = mean(v) - 1.96 * sd(v) / sqrt(length(v)),
                            hi = mean(v) + 1.96 * sd(v) / sqrt(length(v)), n = length(v))
# herbicide years pooled: ratio of sums over the years, per community
pooled <- grp[year %in% herb_years, .(r = rel(sum(N), sum(C))), by = .(biome, MC, stressor, group)]
yearly <- grp[, .(r = rel(N, C)), by = .(biome, MC, year, stressor, group)]
pooled[, stressor := factor(stressor, RUN_ORDER)]
yearly[, stressor := factor(stressor, RUN_ORDER)]
fwrite(pooled, file.path(in_root, sprintf("composition_shift_%s.csv", lev)))
fwrite(yearly, file.path(in_root, sprintf("composition_shift_%s_yearly.csv", lev)))

n_com <- uniqueN(pooled[, .(biome, MC)])
share_sens <- species[stressor == "null" & year %in% herb_years,
                      .(s = sum(Inds_ctl[group == GROUPS[1]]) / sum(Inds_ctl)), by = .(biome, MC)][, mean(s)]
note <- sprintf("%d communities; sensitive: effect >= %.1f at %s (%.1f species per community, %.0f %% of the control individuals). Week %d. Mean and 95 %% CI of the mean.",
                n_com, threshold, lev, species[stressor == "null" & year == 1, sum(group == GROUPS[1]) / n_com],
                100 * share_sens, week)
ylab <- "Change in individuals vs control\n(% of the group's own individuals)"

point_plot <- function(s, base_size, point_size) {
  ggplot(s, aes(stressor, mean, colour = group)) +
    theme_minimal(base_size = base_size) +
    annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.5) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.6),
                    size = point_size, linewidth = 0.6) +
    scale_colour_manual(values = group_cols, name = NULL) +
    scale_x_discrete(labels = c(null = "null run\n(no herbicide)")) +
    theme(legend.position = "top", panel.grid.major.x = element_blank())
}

# ---- A: herbicide years pooled ----
a <- pooled[, mean_ci(r), by = .(stressor, group)]
pa <- point_plot(a, 9, 0.35) +
  labs(x = NULL, y = ylab,
       title = sprintf("%s: change in individuals of the sensitive species, the other species and the whole community", lev),
       subtitle = sprintf("Herbicide years %d-%d pooled. Shaded: null run, the reference for a zero effect", min(herb_years), max(herb_years)),
       caption = note) +
  theme(plot.title = element_text(size = 9.5))
ggsave(file.path(fig_dir, sprintf("fig_4.4_composition_decomposition_%s.png", lev)), pa,
       width = 20, height = 11, units = "cm", dpi = 200, bg = "white")

# ---- A by biome: 10 communities per biome ----
ab <- pooled[, mean_ci(r), by = .(biome, stressor, group)]
pab <- point_plot(ab, 8, 0.2) +
  facet_wrap(~ biome, ncol = 3) +
  labs(x = NULL, y = ylab,
       title = sprintf("%s: change in individuals of the sensitive and other species, by biome", lev),
       subtitle = sprintf("Herbicide years %d-%d pooled; 10 communities per biome, mean and 95 %% CI of the mean",
                          min(herb_years), max(herb_years)),
       caption = sprintf("Sensitive: effect >= %.1f at %s. Shaded: null run (no herbicide). Week %d.", threshold, lev, week)) +
  theme(strip.text = element_text(face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, sprintf("fig_4.4_composition_decomposition_%s_by_biome.png", lev)), pab,
       width = 24, height = 22, units = "cm", dpi = 200, bg = "white")
fwrite(dcast(ab, biome + stressor ~ group, value.var = "mean"),
       file.path(in_root, sprintf("composition_decomposition_%s_by_biome.csv", lev)))

# ---- B: year by year ----
b <- yearly[, mean_ci(r), by = .(year, stressor, group)]
pb <- ggplot(b, aes(year, mean, colour = group, fill = group)) +
  theme_minimal(base_size = 9) +
  annotate("rect", xmin = prov$Tinit + 0.5, xmax = prov$Tinit + prov$HerbDuration + 0.5, ymin = -Inf, ymax = Inf,
           fill = "#d73027", alpha = 0.07) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_ribbon(aes(ymin = lo, ymax = hi), colour = NA, alpha = 0.2) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ stressor, nrow = 2, labeller = as_labeller(c(null = "null run (no herbicide)",
                                                             setNames(RUN_ORDER[-1], RUN_ORDER[-1])))) +
  scale_colour_manual(values = group_cols, name = NULL) +
  scale_fill_manual(values = group_cols, name = NULL) +
  labs(x = "Time (years)", y = ylab,
       title = sprintf("%s: change in individuals of the sensitive and other species, year by year", lev),
       subtitle = "Shaded: herbicide years. Each year: relative change of the group vs control (week 20)", caption = note) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, sprintf("fig_4.4_composition_time_%s.png", lev)), pb,
       width = 24, height = 14, units = "cm", dpi = 200, bg = "white")

# The figure of each species against its effect is made by 08_species_response.R

print(dcast(a, stressor ~ group, value.var = c("mean", "lo", "hi"))[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 1) else v)])
# check: the community total lies between the two groups, community by community
raw <- dcast(pooled, biome + MC + stressor ~ group, value.var = "r")
setnames(raw, GROUPS, c("s", "o", "t"))
raw <- raw[!is.na(s) & !is.na(o)] # communities without sensitive species have one group only
message("total between the groups: ", raw[, sum(t >= pmin(s, o) - 1e-9 & t <= pmax(s, o) + 1e-9)], " / ", nrow(raw))
message(note)

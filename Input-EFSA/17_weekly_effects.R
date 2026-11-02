# §4.4 Week-by-week effects around the herbicide application (exploration)
#
#   Rscript Input-EFSA/17_weekly_effects.R [in_root] [fig_dir]
#   defaults: outputs/main-2026-10-27/4.4_stressor_effects EFSA_REPORT/figures/4.4/weekly
# The model writes the Pt_ file every week (30 weeks per year); scripts 04 and 09 keep only
# week 20. The herbicide acts once per year, in week week_start (1): each effect factor of a
# plant is used once and then reset to 0 (CTDPlant.cpp). This script reads every week of the
# years around the exposure to see how large each effect is just after the application and
# how much of it is left in week 20.
# Groups of each community: 5 % PaR (the most sensitive PFT), 20 % PaR (4 PFTs), rest (10
# least sensitive PFTs) and the whole community. Ratio treatment / control of the sums over
# the 90 communities, 95 % bootstrap interval over the communities; the null run (no
# herbicide, new random stream from year 21) is the reference for a zero effect.
# Outputs:
#   weekly_effects.csv.gz                       ratio per week, group, endpoint, run
#   fig_4.4_weekly_year21_<endpoint>.png        first herbicide year, weeks 1-30
#   fig_4.4_weekly_years20-32_PaR5.png          weekly ratio of the 5 % PaR, years 20-32, z09
#   fig_4.4_weekly_profile.png                  mean ratio by week over the herbicide years 21-30

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
in_root <- if (length(args) >= 1) args[1] else "outputs/main-2026-10-27/4.4_stressor_effects"
fig_dir <- if (length(args) >= 2) args[2] else "EFSA_REPORT/figures/4.4/weekly"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(2021)
n_boot <- 200
years_keep <- 20:32
LEVELS <- c("z05", "z09")
STRESSORS <- c("biomass", "SEbiomass", "survival", "establishment", "sterility", "seednumber")
GROUPS <- c("5 % PaR (1 PFT)", "20 % PaR (4 PFTs)", "Rest (10 least sensitive PFTs)", "Community (20 PFTs)")
gcol <- setNames(c("#b2182b", "#e8743b", "#2c6fbb", "#1a1a1a"), GROUPS)

dirs <- dirname(list.files(in_root, "^Pt\\.csv\\.gz$", recursive = TRUE, full.names = TRUE))
prov <- jsonlite::read_json(file.path(dirs[1], "provenance.json"))
t0 <- (min(years_keep) - 1) * 30
t1 <- max(years_keep) * 30

grp <- rbindlist(lapply(dirs, function(d) {
  com <- fread(file.path(d, "community.csv"))
  n_sp <- nrow(com)
  rank <- setNames(frank(-com$u, ties.method = "first"), com$Species)
  pt <- fread(file.path(d, "Pt.csv.gz"), select = c("level", "stressor", "Time", "PFT", "Inds", "shootmass"))
  pt <- pt[Time > t0 & Time <= t1 & (stressor %in% c("NO", "null") | level %in% LEVELS)]
  pt <- pt[, .(Inds = sum(Inds), shootmass = sum(shootmass)), by = .(level, stressor, Time, PFT)]
  pt[, rk := rank[PFT]]
  members <- list(1, 1:4, (n_sp / 2 + 1):n_sp, seq_len(n_sp))
  # absent PFTs contribute 0 to the sums, so no completion is needed for sums
  out <- rbindlist(Map(function(g, m) pt[rk %in% m, .(NInd = sum(Inds), shootmass = sum(shootmass)),
                                         by = .(level, stressor, Time)][, group := g], GROUPS, members))
  out[, `:=`(biome = basename(dirname(d)), MC = basename(d))]
  out
}))
# complete the missing (community, run, week, group) cells with 0
runs <- unique(grp[, .(level, stressor)])
full <- CJ(com = unique(paste(grp$biome, grp$MC)), r = seq_len(nrow(runs)), Time = (t0 + 1):t1, group = GROUPS)
full <- cbind(full, runs[full$r])[, r := NULL][, c("biome", "MC") := tstrsplit(com, " ")][, com := NULL]
grp <- grp[full, on = .(biome, MC, level, stressor, Time, group)]
for (v in c("NInd", "shootmass")) set(grp, which(is.na(grp[[v]])), v, 0)
long <- melt(grp, measure.vars = c("NInd", "shootmass"), variable.name = "endpoint", variable.factor = FALSE)
ctl <- long[stressor == "NO", .(biome, MC, Time, group, endpoint, control = value)]
long <- ctl[long[stressor != "NO"], on = .(biome, MC, Time, group, endpoint)]
long[stressor == "null", level := "null"]

# ratio of sums and bootstrap over communities (same resamples for every cell)
coms <- unique(long[, .(biome, MC)])[, id := .I]
long <- coms[long, on = .(biome, MC)]
W <- matrix(0L, n_boot, nrow(coms))
for (b in seq_len(n_boot)) W[b, ] <- tabulate(sample(nrow(coms), replace = TRUE), nrow(coms))
res <- long[, {
  tr <- numeric(nrow(coms)); ct <- numeric(nrow(coms))
  tr[id] <- value; ct[id] <- control
  bt <- as.vector(W %*% tr) / as.vector(W %*% ct)
  .(ratio = sum(tr) / sum(ct), lo = quantile(bt, 0.025, na.rm = TRUE), hi = quantile(bt, 0.975, na.rm = TRUE))
}, by = .(endpoint, level, stressor, group, Time)]
res[, `:=`(year = (Time - 1) %/% 30 + 1, week = (Time - 1) %% 30 + 1)]
fwrite(res, file.path(fig_dir, "weekly_effects.csv.gz"))

theme_rep <- theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"), legend.position = "top",
        plot.caption = element_text(hjust = 0, size = 7.5), plot.caption.position = "plot")
EP <- c(NInd = "number of individuals", shootmass = "above-ground biomass")
cap <- paste("Ratio of the sums over the 90 communities; band: 95 % bootstrap interval (communities resampled).",
             "\nGrey: null run (no herbicide), same groups. Dotted line: week 20, the week used in the report.")

# 1. first herbicide year, weeks 1-30, levels z05 and z09
for (ep in names(EP)) {
  s <- res[endpoint == ep & year == 21 & stressor != "null"][, stressor := factor(stressor, STRESSORS)]
  n <- res[endpoint == ep & year == 21 & stressor == "null", -c("stressor", "level")]
  s[, group := factor(group, GROUPS)]; n[, group := factor(group, GROUPS)]
  p <- ggplot(s, aes(week)) + theme_rep +
    geom_hline(yintercept = 1, colour = "grey40", linewidth = 0.3) +
    geom_vline(xintercept = 20, linetype = "dotted") +
    geom_ribbon(data = n, aes(ymin = lo, ymax = hi), fill = "grey60", alpha = 0.35) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = level), alpha = 0.2) +
    geom_line(aes(y = ratio, colour = level), linewidth = 0.5) +
    facet_grid(group ~ stressor, labeller = labeller(group = label_wrap_gen(14))) +
    scale_colour_manual(values = c(z05 = "#e08214", z09 = "#7f3b08"), name = "Level") +
    scale_fill_manual(values = c(z05 = "#e08214", z09 = "#7f3b08"), name = "Level") +
    scale_x_continuous(breaks = c(1, 10, 20, 30)) +
    labs(x = "Week of year 21 (first herbicide year, application in week 1)",
         y = paste("Treatment / control:", EP[[ep]]), caption = cap) +
    theme(axis.text = element_text(size = 6.5), strip.text.y = element_text(size = 7))
  ggsave(file.path(fig_dir, sprintf("fig_4.4_weekly_year21_%s.png", ep)), p, width = 18, height = 14,
         units = "cm", dpi = 300, bg = "white")
}

# 2. 5 % PaR, z09, every week of years 20-32
s <- res[group == GROUPS[1] & level %in% c("z09", "null")]
s[, run := fifelse(stressor == "null", "null run", "z09")]
s[, stressor_f := factor(fifelse(stressor == "null", NA_character_, stressor), STRESSORS)]
n <- s[stressor == "null", -c("stressor_f")]
s <- s[stressor != "null"]
p <- ggplot(s, aes((Time - 1) / 30 + 1)) + theme_rep +
  annotate("rect", xmin = 21, xmax = 31, ymin = -Inf, ymax = Inf, fill = "#d73027", alpha = 0.08) +
  geom_hline(yintercept = 1, colour = "grey40", linewidth = 0.3) +
  geom_ribbon(data = n, aes(ymin = lo, ymax = hi), fill = "grey60", alpha = 0.4) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#b2182b", alpha = 0.2) +
  geom_line(aes(y = ratio), colour = "#b2182b", linewidth = 0.35) +
  geom_point(data = s[week == 20], aes(y = ratio), size = 0.8) +
  facet_grid(endpoint ~ stressor_f, scales = "free_y", labeller = labeller(endpoint = EP)) +
  scale_x_continuous(breaks = seq(20, 32, 2)) +
  labs(x = "Time (years; 30 weeks per year)", y = "Treatment / control, 5 % PaR (most sensitive PFT), z09",
       caption = paste("Weekly values. Black points: week 20 of each year (the values used in the report).",
                       "Shaded: herbicide years 21-30.\nGrey: null run. Band: 95 % bootstrap interval.")) +
  theme(axis.text = element_text(size = 6.5))
ggsave(file.path(fig_dir, "fig_4.4_weekly_years20-32_PaR5.png"), p, width = 18, height = 10, units = "cm",
       dpi = 300, bg = "white")

# 3. mean weekly profile over the herbicide years 21-30 (ratio of the summed values)
herb <- (prov$Tinit + 1):(prov$Tinit + prov$HerbDuration)
pr <- long[((Time - 1) %/% 30 + 1) %in% herb, .(ratio = sum(value) / sum(control)),
           by = .(endpoint, level, stressor, group, week = (Time - 1) %% 30 + 1)]
pr[, group := factor(group, GROUPS)]
n <- pr[stressor == "null", -c("stressor", "level")]
s <- pr[stressor != "null"][, stressor := factor(stressor, STRESSORS)]
p <- ggplot(s, aes(week, ratio, colour = group)) + theme_rep +
  geom_hline(yintercept = 1, colour = "grey40", linewidth = 0.3) +
  geom_vline(xintercept = 20, linetype = "dotted") +
  geom_line(data = n, aes(group = group), colour = "grey60", linewidth = 0.4) +
  geom_line(linewidth = 0.5) +
  facet_grid(endpoint + level ~ stressor, scales = "free_y",
             labeller = labeller(endpoint = c(NInd = "individuals", shootmass = "biomass"))) +
  scale_colour_manual(values = gcol, name = NULL) +
  scale_x_continuous(breaks = c(1, 10, 20, 30)) +
  guides(colour = guide_legend(nrow = 2)) +
  labs(x = "Week of the year (application in week 1)", y = "Treatment / control, herbicide years 21-30 pooled",
       caption = "Ratio of the sums over the 90 communities and the 10 herbicide years. Grey: null run (4 groups).") +
  theme(axis.text = element_text(size = 6.5), strip.text.y = element_text(size = 6.5))
ggsave(file.path(fig_dir, "fig_4.4_weekly_profile.png"), p, width = 18, height = 16, units = "cm", dpi = 300,
       bg = "white")

# numbers for the text: week 2-5 and week 20, year 21 and herbicide years pooled
fwrite(pr, file.path(fig_dir, "weekly_profile.csv"))
options(width = 200)
tab <- pr[group %in% GROUPS[c(1, 4)] & week %in% c(1, 2, 3, 5, 10, 20, 30)]
print(dcast(tab, endpoint + group + stressor + level ~ week, value.var = "ratio"), digits = 3)
message("written in ", fig_dir)

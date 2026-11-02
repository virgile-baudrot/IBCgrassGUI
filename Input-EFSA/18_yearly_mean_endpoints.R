# §4.4 Endpoints as the mean over the 30 weeks of each year (instead of week 20)
#
#   Rscript Input-EFSA/18_yearly_mean_endpoints.R [in_root] [fig_dir]
#   defaults: outputs/main-2026-10-27/4.4_stressor_effects EFSA_REPORT/figures/4.4/yearly_mean
# Scripts 04, 09 and 16 take the endpoints in week 20 of each year. The weekly outputs
# (script 17) show that the effects of survival and biomass are largest just after the
# application in week 1 and partly gone in week 20. Here each endpoint is computed every
# week from the living plants of the Pt_ file (as in 04_output_analysis.R) and then averaged
# over the 30 weeks of the year. The ratios and figures are the same as in script 16:
#   fig_4.4_<endpoint>_pooled_yearmean.png   paired ratio treatment / control, 90 communities,
#                                            median and 95 % interval, null run in grey
#   fig_4.4_par_NInd_pooled_yearmean.png     PFT at risk (PaR), ratio of the sums over the
#                                            communities, 95 % bootstrap interval
#   effects_pooled_yearmean.csv, population_at_risk_ratio_yearmean.csv

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
in_root <- if (length(args) >= 1) args[1] else "outputs/main-2026-10-27/4.4_stressor_effects"
fig_dir <- if (length(args) >= 2) args[2] else "EFSA_REPORT/figures/4.4/yearly_mean"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(2021)
n_boot <- 1000
x_from <- 15
STRESSORS <- c("biomass", "SEbiomass", "survival", "establishment", "sterility", "seednumber")
RUN_ORDER <- c("null", STRESSORS)
GROUPS <- c("5 % PaR", "10 % PaR", "20 % PaR", "50 % PaR", "Rest (50 % least sensitive)", "Community")

dirs <- dirname(list.files(in_root, "^Pt\\.csv\\.gz$", recursive = TRUE, full.names = TRUE))
prov <- jsonlite::read_json(file.path(dirs[1], "provenance.json"))
herb_years <- (prov$Tinit + 1):(prov$Tinit + prov$HerbDuration)

# ---- 1. weekly endpoints, then mean over the weeks of each year ----
res <- lapply(dirs, function(d) {
  com <- fread(file.path(d, "community.csv"))
  n_sp <- nrow(com)
  rank <- setNames(frank(-com$u, ties.method = "first"), com$Species)
  pt <- fread(file.path(d, "Pt.csv.gz"), select = c("level", "stressor", "Time", "PFT", "Inds", "shootmass"))
  pt <- pt[Time > (x_from - 1) * 30]
  pt <- pt[, .(Inds = sum(Inds), shootmass = sum(shootmass)), by = .(level, stressor, Time, PFT)]
  # species without plant or seed have no row: complete with zeros
  full <- CJ(run = unique(paste(pt$level, pt$stressor)), Time = unique(pt$Time), PFT = com$Species)
  full[, c("level", "stressor") := tstrsplit(run, " ")][, run := NULL]
  pt <- pt[full, on = .(level, stressor, Time, PFT)]
  pt[is.na(Inds), `:=`(Inds = 0, shootmass = 0)]
  pt[, year := (Time - 1) %/% 30 + 1]
  ep <- pt[, .(NInd = as.numeric(sum(Inds)), shootmass = sum(shootmass), richness = as.numeric(sum(Inds > 0)),
               shannon = { p <- Inds[Inds > 0] / sum(Inds); -sum(p * log(p)) }),
           by = .(level, stressor, year, Time)]
  ep <- ep[, lapply(.SD, mean), by = .(level, stressor, year), .SDcols = c("NInd", "shootmass", "richness", "shannon")]
  pt[, rk := rank[PFT]]
  members <- list(1, 1:2, 1:4, 1:10, (n_sp / 2 + 1):n_sp, seq_len(n_sp))
  par <- rbindlist(Map(function(g, m) pt[rk %in% m, .(NInd = sum(Inds)), by = .(level, stressor, Time, year)][
    , .(NInd = mean(NInd)), by = .(level, stressor, year)][, group := g], GROUPS, members))
  id <- list(biome = basename(dirname(d)), MC = basename(d))
  list(ep = ep[, names(id) := id], par = par[, names(id) := id])
})
ep <- rbindlist(lapply(res, `[[`, "ep"))
par <- rbindlist(lapply(res, `[[`, "par"))
rm(res)

# ---- 2. paired ratios, median and 95 % interval over the 90 communities ----
long <- melt(ep, id.vars = c("biome", "MC", "level", "stressor", "year"), variable.name = "endpoint",
             variable.factor = FALSE)
ctl <- long[stressor == "NO", .(biome, MC, year, endpoint, control = value)]
eff <- ctl[long[stressor != "NO"], on = .(biome, MC, year, endpoint)]
eff[, effect := fifelse(control > 0, value / control, NA_real_)]
pp <- eff[, .(median = median(effect, na.rm = TRUE), q025 = quantile(effect, 0.025, na.rm = TRUE, names = FALSE),
              q975 = quantile(effect, 0.975, na.rm = TRUE, names = FALSE), n = sum(!is.na(effect))),
          by = .(endpoint, stressor, level, year)]
fwrite(pp, file.path(fig_dir, "effects_pooled_yearmean.csv"))

# ---- 3. PFT at risk: ratio of the sums, bootstrap over the communities ----
ctl <- par[stressor == "NO", .(biome, MC, year, group, control = NInd)]
par <- ctl[par[stressor != "NO"], on = .(biome, MC, year, group)]
par[stressor == "null", level := "null"]
coms <- unique(par[, .(biome, MC)])[, id := .I]
par <- coms[par, on = .(biome, MC)]
W <- matrix(0L, n_boot, nrow(coms))
for (b in seq_len(n_boot)) W[b, ] <- tabulate(sample(nrow(coms), replace = TRUE), nrow(coms))
pr <- par[, {
  tr <- numeric(nrow(coms)); ct <- numeric(nrow(coms))
  tr[id] <- NInd; ct[id] <- control
  bt <- as.vector(W %*% tr) / as.vector(W %*% ct)
  .(ratio = sum(tr) / sum(ct), lo = quantile(bt, 0.025, na.rm = TRUE, names = FALSE),
    hi = quantile(bt, 0.975, na.rm = TRUE, names = FALSE))
}, by = .(level, stressor, group, year)]
fwrite(pr, file.path(fig_dir, "population_at_risk_ratio_yearmean.csv"))

# ---- 4. figures (same layout as script 16) ----
base <- 9
run_lab <- c(null = "null run\n(no herbicide)", setNames(STRESSORS, STRESSORS))
save_fig <- function(p, name, width, height) {
  ggsave(file.path(fig_dir, name), p, width = width, height = height, units = "cm", dpi = 300, bg = "white")
  message("written: ", name)
}
herb_layer <- annotate("rect", xmin = min(herb_years) - 0.5, xmax = max(herb_years) + 0.5, ymin = -Inf, ymax = Inf,
                       fill = "#d73027", alpha = 0.10)
x_years <- scale_x_continuous(breaks = c(20, 30, 40, 50), limits = c(x_from, prov$Tmax), expand = c(0, 0.5))
theme_rep <- theme_minimal(base_size = base) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "top", plot.caption = element_text(hjust = 0, size = base - 1.5),
        plot.caption.position = "plot")
ENDPOINTS <- c(shootmass = "above-ground biomass", richness = "number of PFTs", shannon = "Shannon diversity",
               NInd = "number of individuals")
null_cap <- paste("Endpoint: mean over the 30 weeks of each year. Black line: median over the communities;",
                  "\norange: 95 % interval of treatment / control. Grey: 95 % interval of null run / control",
                  "\n(no herbicide, new random stream from year 21). Shaded: herbicide years.")
for (ep_n in names(ENDPOINTS)) {
  s <- pp[endpoint == ep_n & stressor %in% STRESSORS][, stressor := factor(stressor, STRESSORS)]
  null <- pp[endpoint == ep_n & stressor == "null", -c("level", "stressor")]
  p <- ggplot(s, aes(year)) + theme_rep + herb_layer +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30", linewidth = 0.3) +
    geom_ribbon(data = null, aes(ymin = q025, ymax = q975), fill = "grey50", alpha = 0.35) +
    geom_ribbon(aes(ymin = q025, ymax = q975), fill = "#C75314", alpha = 0.35) +
    geom_line(aes(y = median), linewidth = 0.35) +
    x_years + facet_grid(level ~ stressor) +
    labs(x = "Time (years)", y = paste("Treatment / control:", ENDPOINTS[[ep_n]], "(yearly mean)"), caption = null_cap) +
    theme(axis.text.x = element_text(size = base - 2.5), axis.text.y = element_text(size = base - 2),
          strip.text.y = element_text(size = base - 1.5), panel.spacing = unit(0.3, "lines"))
  save_fig(p, sprintf("fig_4.4_%s_pooled_yearmean.png", ep_n), 17, 22)
}

group_cols <- setNames(c("#67000d", "#cb181d", "#fb6a4a", "#fcae91", "#4575b4", "#1a1a1a"), GROUPS)
par_cap <- paste("Endpoint: mean over the 30 weeks of each year. 90 communities of 20 PFTs, all biomes pooled.",
                 "\nx % PaR (PFT at risk): the x % most sensitive PFTs of each community (5 %: 1 PFT, 10 %: 2, 20 %: 4,",
                 "50 %: 10).\nRatio of the sums over the communities; band: 95 % bootstrap interval. Shaded: herbicide years.")
LEVELS <- sort(setdiff(unique(pr$level), "null"))
q <- pr[year >= x_from]
q <- rbind(q[stressor != "null"], rbindlist(lapply(LEVELS, function(l) copy(q[stressor == "null"])[, level := l])))
q[, `:=`(stressor = factor(stressor, RUN_ORDER), group = factor(group, GROUPS))]
p <- ggplot(q, aes(year, colour = group, fill = group)) + theme_rep + herb_layer +
  geom_hline(yintercept = 1, colour = "grey40", linewidth = 0.3) +
  geom_ribbon(aes(ymin = lo, ymax = hi), colour = NA, alpha = 0.15) +
  geom_line(aes(y = ratio), linewidth = 0.4) +
  x_years +
  scale_colour_manual(values = group_cols, name = NULL) + scale_fill_manual(values = group_cols, name = NULL) +
  guides(colour = guide_legend(ncol = 6), fill = guide_legend(ncol = 6)) +
  facet_grid(level ~ stressor, labeller = labeller(stressor = run_lab)) +
  labs(x = "Time (years)", y = "Treatment / control: number of individuals (yearly mean)", caption = par_cap) +
  theme(axis.text.x = element_text(size = base - 2.5), axis.text.y = element_text(size = base - 2),
        panel.spacing = unit(0.3, "lines"), legend.text = element_text(size = base - 1))
save_fig(p, "fig_4.4_par_NInd_pooled_yearmean.png", 17, 23)

# ---- 5. comparison with week 20, herbicide years 21-30 ----
options(width = 200)
w20 <- fread(file.path(in_root, "effects_pooled.csv"))[, year := Time %/% 30 + 1][year %in% herb_years]
cmp <- merge(pp[year %in% herb_years, .(yearmean = mean(median), lo_yearmean = mean(q025)), by = .(endpoint, stressor, level)],
             w20[, .(week20 = mean(median), lo_week20 = mean(q025)), by = .(endpoint, stressor, level)],
             by = c("endpoint", "stressor", "level"))
print(dcast(cmp[level %in% c("z00", "z05", "z09")], endpoint + stressor ~ level, value.var = c("week20", "yearmean")),
      digits = 3)
p20 <- fread(file.path(in_root, "population_at_risk_ratio.csv"))[year %in% herb_years & group %in% GROUPS[c(1, 5, 6)]]
cp <- merge(pr[year %in% herb_years & group %in% GROUPS[c(1, 5, 6)], .(yearmean = mean(ratio)), by = .(group, stressor, level)],
            p20[, .(week20 = mean(ratio)), by = .(group, stressor, level)], by = c("group", "stressor", "level"))
print(dcast(cp[level %in% c("null", "z05", "z09")], group + stressor ~ level, value.var = c("week20", "yearmean")),
      digits = 3)
fwrite(cmp, file.path(fig_dir, "compare_week20_yearmean.csv"))
fwrite(cp, file.path(fig_dir, "compare_week20_yearmean_par.csv"))

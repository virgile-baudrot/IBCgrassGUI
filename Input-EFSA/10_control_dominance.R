# Control runs: dominance of a few species and number of species that persist, by the
# number of species at the start (report section "Analysis of Control", old Fig. 4.4-4.11)
#
#   Rscript Input-EFSA/10_control_dominance.R [in_4.3] [in_4.4] [fig_dir] [example_biome]
#   defaults: outputs/control-2026-10-28/4.3_control_design
#             outputs/main-2026-10-27/4.4_stressor_effects
#             EFSA_REPORT/figures/4.2
#             BOR_I1
# in_4.3: outputs of run_4.3_control_design.R, designs full, s50, s25 and s10 (nested random
#   subsets of the biome pool, 5 replicates, 173 x 173 cells, 50 years, no herbicide).
# in_4.4: the controls (level z00, stressor NO) of the §4.4 communities of 20 species give the
#   fifth initial size. Optional: skipped if the folder is missing.
# All values are read at week 20 of each year (before seed release, as in §4.4). The dominance
# of a species is its mean number of individuals at week 20 of the last 10 years.

source(file.path(Sys.getenv("IBC_ROOT", "."), "Input-EFSA", "04_output_analysis.R"))
suppressPackageStartupMessages(library(patchwork))

args <- commandArgs(trailingOnly = TRUE)
arg <- function(i, default) if (length(args) >= i) args[i] else default
in43 <- arg(1, "outputs/control-2026-10-28/4.3_control_design")
in44 <- arg(2, "outputs/main-2026-10-27/4.4_stressor_effects")
fig_dir <- arg(3, "EFSA_REPORT/figures/4.2")
example <- arg(4, "BOR_I1")
week <- 20
last_years <- 10
n_top <- 10 # dominants of the full pool highlighted in the subset rankings
dir.create(file.path(fig_dir, "by_biome"), recursive = TRUE, showWarnings = FALSE)

DESIGNS <- c("full", "s50", "s25", "r20", "s10")

# ------------------------------------------------------------------
# Reading
# ------------------------------------------------------------------

read_design <- function(dir) {
  pt <- fread(file.path(dir, "Pt.csv.gz"))[, seed := NULL]
  pt <- pt[Time %% 30 == week]
  # communities.csv lists the species of each replicate, also those never present in Pt_
  com <- fread(file.path(dir, "communities.csv"))[, .(rep, PFT = Species)]
  com[, `:=`(biome = pt$biome[1], design = pt$design[1])]
  times <- unique(pt[, .(biome, design, rep, Time)])
  full <- merge(times, com, by = c("biome", "design", "rep"), allow.cartesian = TRUE)
  pt <- pt[full, on = c("biome", "design", "rep", "Time", "PFT")]
  for (v in intersect(PT_VARS, names(pt))) set(pt, which(is.na(pt[[v]])), v, 0)
  pt
}

dirs <- dirname(list.files(in43, "^Pt\\.csv\\.gz$", recursive = TRUE, full.names = TRUE))
if (length(dirs) == 0) stop("no Pt.csv.gz under ", in43)
pt <- rbindlist(lapply(dirs, read_design), use.names = TRUE)
prov <- read_provenance(dirs[1])
tmax <- prov$Tmax

# §4.4 controls: 20 random species, 10 communities per biome (MC used as replicate)
if (dir.exists(in44)) {
  f44 <- list.files(in44, "^Pt\\.csv\\.gz$", recursive = TRUE, full.names = TRUE)
  pt44 <- rbindlist(lapply(f44, function(f) {
    x <- fread(cmd = sprintf("zcat %s | awk -F, 'NR == 1 || ($5 == \"z00\" && $6 == \"NO\" && $7 %% 30 == %d)'",
                             shQuote(f), week))
    com <- fread(file.path(dirname(f), "community.csv"))[, .(PFT = Species)]
    full <- CJ(Time = unique(x$Time), PFT = com$PFT)
    x <- x[full, on = c("Time", "PFT")]
    for (v in intersect(PT_VARS, names(x))) set(x, which(is.na(x[[v]])), v, 0)
    prov44 <- read_provenance(dirname(f))
    x[, `:=`(biome = prov44$biome, design = "r20", rep = prov44$MC)]
    x[, intersect(names(pt), names(x)), with = FALSE]
  }), use.names = TRUE)
  pt <- rbind(pt, pt44, fill = TRUE)
}

n_start <- pt[, .(n0 = uniqueN(PFT)), by = .(biome, design, rep)]
size_label <- function(design, n0) fifelse(design == "full", paste0("full pool (", n0, ")"),
                                           paste0(n0, " species"))
pt[, design := factor(design, levels = intersect(DESIGNS, unique(design)))]

# ------------------------------------------------------------------
# Summaries
# ------------------------------------------------------------------

ep <- pt[, .(richness = sum(Inds > 0), NInd = sum(Inds), shootmass = sum(shootmass)),
         by = .(biome, design, rep, Time)]
ep <- n_start[ep, on = c("biome", "design", "rep")]
ep[, year := Time %/% 30 + 1]

# dominance of each species in each run: mean over the last years
last <- pt[Time > (tmax - last_years) * 30,
           .(Inds = mean(Inds), shootmass = mean(shootmass), seedlings = mean(seedlings), seeds = mean(seeds)),
           by = .(biome, design, rep, PFT)]
last[, share := Inds / sum(Inds), by = .(biome, design, rep)]
setorder(last, biome, design, rep, -Inds)
last[, rank := seq_len(.N), by = .(biome, design, rep)]

# top n_top species of the full pool (mean over the replicates)
top_full <- last[design == "full", .(Inds = mean(Inds)), by = .(biome, PFT)][order(biome, -Inds)]
top_full <- top_full[, .SD[seq_len(min(n_top, .N))], by = biome]

# one row per run: species present at the end (mean over the last years of the number of
# species with living plants), species holding 50 % and 80 % of the plants
end_tab <- last[, .(n50 = which(cumsum(share) >= 0.5)[1],
                    n80 = which(cumsum(share) >= 0.8)[1],
                    top_share = max(share)), by = .(biome, design, rep)]
end_tab <- n_start[end_tab, on = c("biome", "design", "rep")]
end_tab <- ep[year > tmax - last_years, .(present = mean(richness)), by = .(biome, design, rep)][
  end_tab, on = c("biome", "design", "rep")]
# share of the plants that belong to a top species of the full pool, in the subsets
end_tab <- merge(end_tab,
                 last[top_full[, .(biome, PFT)], on = c("biome", "PFT"), nomatch = 0L][
                   , .(n_top_in = .N, top_full_share = sum(share)), by = .(biome, design, rep)],
                 by = c("biome", "design", "rep"), all.x = TRUE)
fwrite(end_tab[order(biome, design, rep)], file.path(fig_dir, "control_end_of_run.csv"))
fwrite(end_tab[, .(n0 = round(mean(n0)), present_mean = mean(present), present_min = min(present),
                   present_max = max(present), n80_mean = mean(n80), top_share_mean = mean(top_share),
                   n_runs = .N), by = .(design)][order(design)],
       file.path(fig_dir, "control_end_of_run_by_size.csv"))

# ------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------

design_colours <- c(full = "#08306b", s50 = "#2171b5", s25 = "#6baed6", r20 = "#d94801", s10 = "#9ecae1")
design_names <- c(full = "full pool", s50 = "50 species", s25 = "25 species", r20 = "20 species (§4.4)",
                  s10 = "10 species")

# number of individuals per species and year, one panel per replicate (old Fig. 4.5, 4.6, 4.9, 4.10)
# alphabetic = TRUE: species in alphabetic order from the top instead of by abundance
plot_heatmap <- function(x, title, alphabetic = FALSE) {
  ord <- if (alphabetic) sort(unique(x$PFT), decreasing = TRUE) else
    x[Time > (tmax - last_years) * 30, .(m = mean(Inds)), by = PFT][order(m), PFT]
  x <- copy(x)[, PFT := factor(PFT, levels = ord)]
  x[, Inds_plot := fifelse(Inds > 0, Inds, NA_real_)]
  n_sp <- length(ord)
  ggplot(x, aes(x = Time %/% 30 + 1, y = PFT, fill = Inds_plot)) +
    theme_minimal(base_size = 8) +
    geom_tile() +
    scale_fill_viridis_c(trans = "log10", na.value = "grey88", name = "Individuals\n(log scale)") +
    scale_x_continuous(breaks = seq(0, tmax, 10), expand = c(0, 0)) +
    facet_grid(~ rep, labeller = labeller(rep = function(r) paste("replicate", r))) +
    labs(x = "Time (years)", y = NULL, title = title,
         caption = paste("Week 20 of each year. Grey: no living plant. Species",
                         if (alphabetic) "in alphabetic order." else "ordered by their mean abundance in the last 10 years.")) +
    theme(axis.text.y = element_text(size = if (n_sp > 60) 3.5 else if (n_sp > 30) 5 else 7),
          panel.grid = element_blank(), panel.spacing = unit(2, "mm"))
}

# number of individuals per species at the end of the run, one panel per replicate
# (old Fig. 4.7, 4.8, 4.11). Bars of the dominants of the full pool are dark.
plot_ranking <- function(x, top, title) {
  x <- copy(x)[, top := PFT %in% top]
  x[, key := paste(rep, PFT, sep = "__")]
  x[, key := factor(key, levels = x[order(rep, -Inds), key])]
  ggplot(x, aes(x = key, y = Inds, fill = top)) +
    theme_minimal(base_size = 8) +
    geom_col() +
    scale_x_discrete(labels = function(k) sub("^[0-9]+__", "", k)) +
    scale_fill_manual(values = c(`TRUE` = "#08306b", `FALSE` = "#9ecae1"), breaks = c("TRUE", "FALSE"),
                      labels = c(`TRUE` = paste("among the", n_top, "dominants of the full pool"),
                                 `FALSE` = "other species"), name = NULL) +
    facet_wrap(~ rep, ncol = 1, scales = "free_x", labeller = labeller(rep = function(r) paste("replicate", r))) +
    labs(x = NULL, y = "Mean number of individuals, last 10 years", title = title) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                     size = if (uniqueN(x$PFT) > 60) 3.5 else 5),
          legend.position = "top", panel.grid.major.x = element_blank())
}

# full pool, the four Pt_ variables, mean and replicates (old Fig. 4.4)
plot_ranking_vars <- function(x, title) {
  ord <- x[, .(m = mean(Inds)), by = PFT][order(-m), PFT]
  vars <- c(Inds = "Individuals", shootmass = "Shoot mass (mg)", seedlings = "Seedlings", seeds = "Seeds")
  fill <- c(Inds = "steelblue", shootmass = "#F54927", seedlings = "#5AA33B", seeds = "#953BA3")
  long <- melt(x, id.vars = c("rep", "PFT"), measure.vars = names(vars), variable.factor = FALSE)
  long[, PFT := factor(PFT, levels = ord)]
  ps <- lapply(names(vars), function(v) {
    d <- long[variable == v]
    ggplot(d, aes(x = PFT, y = value)) +
      theme_minimal(base_size = 8) +
      stat_summary(fun = mean, geom = "col", fill = fill[[v]], alpha = 0.8) +
      geom_point(size = 0.4, alpha = 0.6) +
      labs(x = NULL, y = vars[[v]]) +
      theme(axis.text.x = if (v == "seeds") element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3.5) else
        element_blank(), panel.grid.major.x = element_blank())
  })
  wrap_plots(ps, ncol = 1) +
    plot_annotation(title = title,
                    caption = "Mean over the last 10 years (week 20). Bars: mean of the 5 replicates; points: replicates. Species in the order of the number of individuals.")
}

# species present over time, all biomes, by initial number of species
# (legend without the reference to §4.4 for the 20 species)
plot_richness_time <- function(e) {
  s <- e[, .(mean = mean(richness), lo = min(richness), hi = max(richness)), by = .(biome, design, year)]
  names_time <- replace(design_names, "r20", "20 species")
  ggplot(s, aes(x = year, colour = design, fill = design)) +
    theme_minimal(base_size = 9) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(aes(y = mean), linewidth = 0.5) +
    scale_colour_manual(values = design_colours, breaks = DESIGNS, labels = names_time, name = "At the start") +
    scale_fill_manual(values = design_colours, breaks = DESIGNS, labels = names_time, name = "At the start") +
    scale_y_continuous(limits = c(0, NA)) +
    facet_wrap(~ biome, ncol = 3) +
    labs(x = "Time (years)", y = "Number of species with living plants",
         caption = "Week 20 of each year. Line: mean of the replicates; band: minimum to maximum.") +
    theme(legend.position = "top")
}

# species present at the end against species at the start
plot_richness_end <- function(t) {
  ggplot(t, aes(x = n0, y = present)) +
    theme_minimal(base_size = 9) +
    geom_line(data = data.table(n0 = c(10, max(t$n0)), present = c(10, max(t$n0))),
              linetype = "dashed", colour = "grey50") +
    geom_jitter(aes(colour = design), width = 0.01, height = 0.15, size = 1.2, alpha = 0.7) +
    stat_summary(aes(group = design), fun = mean, geom = "point", shape = 3, size = 3) +
    scale_x_log10(breaks = c(10, 20, 25, 50, 100)) +
    scale_colour_manual(values = design_colours, breaks = DESIGNS, labels = design_names, name = "At the start") +
    facet_wrap(~ biome, ncol = 3) +
    labs(x = "Number of species at the start (log scale)",
         y = "Species with living plants, mean of the last 10 years",
         caption = "Points: replicates; crosses: mean. Dashed line: all species persist.") +
    theme(legend.position = "top")
}

# number of individuals over time, all biomes, by initial number of species
# (legend without the reference to §4.4 for the 20 species)
plot_total_time <- function(e) {
  s <- e[, .(mean = mean(NInd), lo = min(NInd), hi = max(NInd)), by = .(biome, design, year)]
  names_time <- replace(design_names, "r20", "20 species")
  ggplot(s, aes(x = year, colour = design, fill = design)) +
    theme_minimal(base_size = 9) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(aes(y = mean), linewidth = 0.5) +
    scale_colour_manual(values = design_colours, breaks = DESIGNS, labels = names_time, name = "At the start") +
    scale_fill_manual(values = design_colours, breaks = DESIGNS, labels = names_time, name = "At the start") +
    facet_wrap(~ biome, ncol = 3) +
    labs(x = "Time (years)", y = "Number of individuals (173 x 173 cells)",
         caption = "Week 20 of each year. Line: mean of the replicates; band: minimum to maximum.") +
    theme(legend.position = "top")
}

# ------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------

save_fig(plot_richness_time(ep), file.path(fig_dir, "fig_4.2_richness_time_all_biomes.png"), 18, 18)
save_fig(plot_richness_end(end_tab), file.path(fig_dir, "fig_4.2_richness_end_vs_start.png"), 18, 18)
save_fig(plot_total_time(ep), file.path(fig_dir, "fig_4.2_NInd_time_all_biomes.png"), 18, 18)

for (b in unique(pt$biome)) {
  sub <- if (b == example) fig_dir else file.path(fig_dir, "by_biome")
  tag <- file.path(sub, paste0("fig_4.2_", b))
  for (d in intersect(c("full", "s50", "s25", "s10"), as.character(unique(pt[biome == b, design])))) {
    n0 <- n_start[biome == b & design == d, n0[1]]
    lab <- paste0(b, ", ", size_label(d, n0), " at the start, 173 x 173 cells, controls")
    save_fig(plot_heatmap(pt[biome == b & design == d], lab), paste0(tag, "_heatmap_", d, ".png"),
             24, if (n0 > 60) 22 else if (n0 > 30) 16 else 11)
    save_fig(plot_ranking(last[biome == b & design == d], top_full[biome == b, PFT], lab),
             paste0(tag, "_ranking_", d, ".png"), 18, 22)
  }
  # 20 species of §4.4 (10 communities), example biome only, to illustrate
  if (b == example && "r20" %in% pt[biome == b, design]) {
    n_r20 <- uniqueN(pt[biome == b & design == "r20", PFT])
    lab <- paste0(b, ", 20 species at the start (controls of §4.4), 173 x 173 cells")
    h <- if (n_r20 > 60) 22 else if (n_r20 > 30) 16 else 11
    save_fig(plot_heatmap(pt[biome == b & design == "r20"], lab), paste0(tag, "_heatmap_r20.png"), 40, h)
    save_fig(plot_heatmap(pt[biome == b & design == "r20"], lab, alphabetic = TRUE),
             paste0(tag, "_heatmap_r20_alphabetic.png"), 40, h)
  }
  if ("full" %in% pt[biome == b, design])
    save_fig(plot_ranking_vars(last[biome == b & design == "full"],
                               paste0(b, ", full pool at the start, 173 x 173 cells, controls")),
             paste0(tag, "_ranking_variables_full.png"), 24, 18)
}

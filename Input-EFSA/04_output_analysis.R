# Analysis and figures of the simulations of §4.3 and §4.4 (and shared helpers for §4.1)
#
#   Rscript Input-EFSA/04_output_analysis.R [out_4.3] [out_4.4] [fig_dir]
#   defaults: outputs/4.3_control_design outputs/4.4_stressor_effects EFSA_REPORT/figures
# Reads the Pt.csv.gz / Grd.csv.gz written by run_4.3_control_design.R and
# run_4.4_stressor_effects.R. The file can also be sourced for its functions
# (output-SINGLE/analyze_outputs_*.R do so).
#
# Endpoints are computed from Pt_ (living plants per species), not from the shannon and
# NPFT columns of Grd_, which count litter and seeds (R§5.2). Pt_ has no row for a species
# without plant or seed: these rows are added with zeros before any sum.
# Effects are paired: each treatment run is divided by the control run of the same
# community (same species, same seed), then summarised over communities by the median and
# the 2.5 and 97.5 percentiles. "Pooled" figures pool the communities of all biomes.
# After the first herbicide year a treatment run follows another random stream than its
# control, so its ratio also contains the stochasticity of the model. The null run of each
# community (no herbicide, new random stream from year Tinit+1) gives that part alone:
# its interval (grey) is drawn behind every effect and is the reference for a zero effect.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

STRESSOR_ORDER <- c("NO", "biomass", "SEbiomass", "survival", "establishment", "sterility", "seednumber")
PT_VARS <- c("Inds", "seedlings", "seeds", "repromass", "cover", "shootmass")
# Grd_ columns that hold one end-of-run value on every row (R§5.2): never plotted
GRD_CONSTANT <- c("meanShannon", "meanNPFT", "NNonClonal", "NClonal", "NGenets")
ENDPOINT_LABELS <- c(NInd = "Number of individuals", shootmass = "Above-ground biomass (mg)",
                     richness = "Species richness", shannon = "Shannon diversity",
                     bray_curtis = "Bray-Curtis distance to control")

# ------------------------------------------------------------------
# Reading and endpoints
# ------------------------------------------------------------------

read_provenance <- function(dir) jsonlite::read_json(file.path(dir, "provenance.json"))

# Add the missing Pt_ rows (species without plant or seed) with zeros.
# run_cols identify a run, community_cols the set of species shared by its runs.
complete_pt <- function(pt, run_cols, community_cols = run_cols) {
  species <- unique(pt[, c(community_cols, "PFT"), with = FALSE])
  times <- unique(pt[, c(run_cols, "Time"), with = FALSE])
  full <- merge(times, species, by = community_cols, allow.cartesian = TRUE)
  pt <- pt[full, on = c(run_cols, "Time", "PFT")]
  for (v in intersect(PT_VARS, names(pt))) set(pt, which(is.na(pt[[v]])), v, 0)
  pt
}

# Community endpoints of each run and time step, from living plants (Inds)
pt_endpoints <- function(pt, run_cols) {
  pt[, .(NInd = as.numeric(sum(Inds)),
         shootmass = as.numeric(sum(shootmass)),
         richness = as.numeric(sum(Inds > 0)),
         shannon = {
           p <- Inds[Inds > 0] / sum(Inds)
           -sum(p * log(p))
         }), by = c(run_cols, "Time")]
}

# Bray-Curtis distance (on Inds per species) between each treatment run and its control
paired_bray_curtis <- function(pt, community_cols) {
  ctl <- pt[stressor == "NO", c(community_cols, "Time", "PFT", "Inds"), with = FALSE]
  setnames(ctl, "Inds", "Inds_control")
  x <- ctl[pt[stressor != "NO"], on = c(community_cols, "Time", "PFT")]
  x[, .(value = sum(abs(Inds - Inds_control)) / sum(Inds + Inds_control)),
    by = c(community_cols, "level", "stressor", "Time")]
}

# Ratio treatment / paired control for each endpoint (NA when the control is 0)
paired_effects <- function(ep, community_cols, endpoints = c("NInd", "shootmass", "richness", "shannon")) {
  long <- melt(ep, id.vars = c(community_cols, "level", "stressor", "Time"),
               measure.vars = endpoints, variable.name = "endpoint", variable.factor = FALSE)
  ctl <- long[stressor == "NO", c(community_cols, "Time", "endpoint", "value"), with = FALSE]
  setnames(ctl, "value", "control")
  eff <- ctl[long[stressor != "NO"], on = c(community_cols, "Time", "endpoint")]
  eff[, effect := fifelse(control > 0, value / control, NA_real_)]
  eff
}

summarise_effects <- function(eff, by) {
  eff[, .(median = median(effect, na.rm = TRUE),
          q025 = quantile(effect, 0.025, na.rm = TRUE, names = FALSE),
          q975 = quantile(effect, 0.975, na.rm = TRUE, names = FALSE),
          n = sum(!is.na(effect))), by = c(by, "Time")]
}

# All communities of a §4.4 output folder. week: within-year week kept (20 = before seed
# release, R§5.2); NA keeps every week (memory: about 1 M Pt rows per community).
read_44 <- function(root, week = 20) {
  dirs <- dirname(list.files(root, "^Pt\\.csv\\.gz$", recursive = TRUE, full.names = TRUE))
  if (length(dirs) == 0) stop("no Pt.csv.gz under ", root)
  run_cols <- c("biome", "MC", "level", "stressor")
  res <- lapply(dirs, function(d) {
    pt <- fread(file.path(d, "Pt.csv.gz"))[, seed := NULL]
    if (!is.na(week)) pt <- pt[Time %% 30 == week]
    pt <- complete_pt(pt, run_cols, c("biome", "MC"))
    list(ep = pt_endpoints(pt, run_cols), bc = paired_bray_curtis(pt, c("biome", "MC")))
  })
  ep <- rbindlist(lapply(res, `[[`, "ep"))
  eff <- paired_effects(ep, c("biome", "MC"))
  bc <- rbindlist(lapply(res, `[[`, "bc"))[, `:=`(endpoint = "bray_curtis", effect = value)]
  list(ep = ep, eff = rbind(eff, bc, fill = TRUE), prov = read_provenance(dirs[1]))
}

# One design of §4.3 (all weeks)
read_43 <- function(dir) {
  pt <- fread(file.path(dir, "Pt.csv.gz"))[, seed := NULL]
  pt <- complete_pt(pt, c("biome", "design", "rep"))
  list(pt = pt, grd = fread(file.path(dir, "Grd.csv.gz")), prov = read_provenance(dir))
}

# ------------------------------------------------------------------
# Plot helpers
# ------------------------------------------------------------------

# Herbicide years Tinit+1 ... Tinit+HerbDuration, shown on every time plot
phase_layer <- function(prov) {
  list(annotate("rect", xmin = prov$Tinit, xmax = prov$Tinit + prov$HerbDuration,
                ymin = -Inf, ymax = Inf, fill = "#d73027", alpha = 0.10),
       geom_vline(xintercept = c(prov$Tinit, prov$Tinit + prov$HerbDuration),
                  linetype = "dotted", colour = "grey40", linewidth = 0.3))
}

save_fig <- function(p, file, width = 16, height = 12) {
  ggsave(file, p, width = width, height = height, units = "cm", dpi = 300, bg = "white")
  message("written: ", file)
}

# ------------------------------------------------------------------
# §4.3 control runs
# ------------------------------------------------------------------

# Number of individuals per species over time, one panel per replicate (Fig. 4.12, 4.13, 4.16, 4.17)
plot_species_heatmap <- function(pt) {
  ggplot(pt) +
    theme_minimal() +
    scale_fill_gradientn(colors = c("red", "blue", "cyan", "green"),
                         values = scales::rescale(c(0, 0.001, 0.5, 1)), name = "# Inds") +
    labs(x = "Time (years)", y = "Species") +
    geom_tile(aes(x = Time / 30, y = PFT, fill = Inds)) +
    theme(axis.text.y = element_text(size = 4)) +
    facet_grid(~ rep, labeller = label_both)
}

# Species ranking: sum over the run of individuals, shoot mass, seedlings or seeds
# (Fig. 4.11, 4.14, 4.15, 4.18)
species_sums <- function(pt) {
  pt[, .(individuals = sum(Inds), shootmass = sum(shootmass), seedlings = sum(seedlings),
         seeds = sum(seeds)), by = .(PFT, rep)]
}

plot_species_ranking <- function(sums, var = "individuals", top = 30) {
  fill <- c(individuals = "steelblue", shootmass = "#F54927", seedlings = "#5AA33B", seeds = "#953BA3")
  ord <- sums[, .(tot = sum(get(var))), by = PFT][order(-tot)][seq_len(min(top, .N)), PFT]
  ggplot(sums[PFT %in% ord], aes(x = factor(PFT, levels = ord), y = .data[[var]])) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) +
    labs(x = "Species", y = paste("Sum of", var, "over the run")) +
    geom_col(fill = fill[[var]]) +
    facet_grid(rep ~ ., labeller = label_both)
}

# Community endpoints of the control runs, one line per replicate, one panel per design
plot_control_dynamics <- function(ep, ep_name) {
  ggplot(ep, aes(x = Time / 30, y = .data[[ep_name]], group = rep)) +
    theme_minimal() +
    labs(x = "Time (years)", y = ENDPOINT_LABELS[[ep_name]]) +
    geom_line(alpha = 0.4) +
    facet_wrap(~ design, scales = "free_y")
}

# ------------------------------------------------------------------
# §4.4 effects
# ------------------------------------------------------------------

# Median and 95 % interval over communities of the paired effect.
# null: summary of the null runs, without the facet columns that it must repeat over
plot_effect <- function(s, ep_name, prov, rows, cols, title = NULL, null = NULL) {
  ref <- if (ep_name == "bray_curtis") 0 else 1
  ylab <- if (ep_name == "bray_curtis") ENDPOINT_LABELS[[ep_name]] else
    paste("Treatment / control:", tolower(ENDPOINT_LABELS[[ep_name]]))
  ggplot(s[endpoint == ep_name], aes(x = Time / 30)) +
    theme_minimal() +
    phase_layer(prov) +
    geom_hline(yintercept = ref, linetype = "dashed", colour = "grey30") +
    (if (!is.null(null)) geom_ribbon(data = null[endpoint == ep_name], aes(ymin = q025, ymax = q975),
                                     fill = "grey50", alpha = 0.35)) +
    (if (!is.null(null)) geom_line(data = null[endpoint == ep_name], aes(y = median), colour = "grey45",
                                   linewidth = 0.4)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), fill = "#C75314", alpha = 0.35) +
    geom_line(aes(y = median), linewidth = 0.4) +
    labs(x = "Time (years)", y = ylab, title = title,
         caption = if (!is.null(null)) "Grey: null run / control (no herbicide, new random stream from the first herbicide year)") +
    facet_grid(as.formula(paste(rows, "~", cols)))
}

# Effect at the end of the herbicide years and at the end of the run (pooled over biomes)
effect_table <- function(eff, prov) {
  yrs <- c(end_exposure = prov$Tinit + prov$HerbDuration, end_run = prov$Tmax)
  x <- eff[ceiling(Time / 30) %in% yrs]
  x[, period := names(yrs)[match(ceiling(Time / 30), yrs)]]
  dcast(summarise_effects(x, c("endpoint", "stressor", "level", "period"))[, Time := NULL],
        endpoint + stressor + level ~ period, value.var = c("median", "q025", "q975"))
}

# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  arg <- function(i, default) if (length(args) >= i) args[i] else default
  out43 <- arg(1, "outputs/4.3_control_design")
  out44 <- arg(2, "outputs/4.4_stressor_effects")
  fig_dir <- arg(3, "EFSA_REPORT/figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  # §4.3: one set of figures per biome and design
  for (bdir in list.dirs(out43, recursive = FALSE)) {
    designs <- intersect(c("full", "nodom", "r20"), basename(list.dirs(bdir, recursive = FALSE)))
    eps <- list()
    for (d in designs) {
      x <- read_43(file.path(bdir, d))
      tag <- paste0("fig_4.3_", basename(bdir), "_", d)
      save_fig(plot_species_heatmap(x$pt), paste0(fig_dir, "/", tag, "_heatmap.png"), 24, 14)
      sums <- species_sums(x$pt)
      for (v in c("individuals", "shootmass", "seedlings", "seeds"))
        save_fig(plot_species_ranking(sums, v), paste0(fig_dir, "/", tag, "_ranking_", v, ".png"), 16, 16)
      eps[[d]] <- pt_endpoints(x$pt, c("design", "rep"))
    }
    if (length(eps)) {
      ep <- rbindlist(eps)[, design := factor(design, levels = designs)]
      for (e in c("NInd", "shootmass", "richness", "shannon"))
        save_fig(plot_control_dynamics(ep, e),
                 paste0(fig_dir, "/fig_4.3_", basename(bdir), "_control_", e, ".png"), 16, 7)
    }
  }

  # §4.4: per biome (one stressor) and pooled over biomes (all stressors)
  if (dir.exists(out44)) {
    x <- read_44(out44)
    null_eff <- x$eff[stressor == "null"]
    eff <- x$eff[stressor != "null"][, stressor := factor(stressor, levels = STRESSOR_ORDER[-1])]
    by_biome <- summarise_effects(eff, c("endpoint", "biome", "stressor", "level"))
    pooled <- summarise_effects(eff, c("endpoint", "stressor", "level"))
    null_by_biome <- summarise_effects(null_eff, c("endpoint", "biome"))
    null_pooled <- summarise_effects(null_eff, "endpoint")
    # the null rows are written with stressor "null" and level "z00"
    fwrite(rbind(by_biome, null_by_biome[, `:=`(stressor = "null", level = "z00")], use.names = TRUE),
           file.path(out44, "effects_by_biome.csv"))
    fwrite(rbind(pooled, null_pooled[, `:=`(stressor = "null", level = "z00")], use.names = TRUE),
           file.path(out44, "effects_pooled.csv"))
    fwrite(rbind(effect_table(eff, x$prov), effect_table(null_eff, x$prov)),
           file.path(out44, "effects_table.csv"))
    null_by_biome[, c("stressor", "level") := NULL]
    null_pooled[, c("stressor", "level") := NULL]
    # figure sizes were set for 5 levels (z01-z05): one level = one column (by biome) or row (pooled)
    n_lev <- uniqueN(eff$level)
    for (e in names(ENDPOINT_LABELS)) {
      for (s in c("survival", "sterility"))
        save_fig(plot_effect(by_biome[stressor == s], e, x$prov, "biome", "level", paste("Stressor:", s),
                             null = null_by_biome),
                 paste0(fig_dir, "/fig_4.4_", e, "_", s, "_by_biome.png"), 3 + 2.6 * n_lev, 20)
      save_fig(plot_effect(pooled, e, x$prov, "level", "stressor", "All biomes pooled", null = null_pooled),
               paste0(fig_dir, "/fig_4.4_", e, "_pooled.png"), 18, 3 + 2.2 * n_lev)
    }
  }
}

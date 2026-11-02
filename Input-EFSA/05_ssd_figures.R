# §4.4 Species sensitivity of each community: SSD curves and exact species values
#
#   Rscript Input-EFSA/05_ssd_figures.R [in_root] [fig_dir] [biome ...]
#   with biomes: only the figures gathering the communities of these biomes
#   defaults: outputs/main-2026-10-27/4.4_stressor_effects EFSA_REPORT/figures/4.4_ssd
# Reads the community.csv written by run_4.4_stressor_effects.R (one per biome and MC).
# One figure per community, all of them in one PDF, and one table of all values:
#   left:  the 5 SSDs (lognormal, sdlog 1, 95th percentile = z) and the 20 species on each.
#          A species has the same normal draw u at every level, so it sits at the same
#          height (pnorm(u)) on the 5 curves. Effects are capped at 0.999 (open circles).
#   right: the exact effect of each species at each level (the value used in the runs).
# Per biome, all communities together (log and linear x axis): one panel per level, the
# theoretical SSD and the empirical distribution of the 20 species of each community.

source(file.path(Sys.getenv("IBC_ROOT", "."), "Input-EFSA", "02_simulation_setting.R"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
in_root <- if (length(args) >= 1) args[1] else "outputs/main-2026-10-27/4.4_stressor_effects"
fig_dir <- if (length(args) >= 2) args[2] else "EFSA_REPORT/figures/4.4_ssd"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

sdlog <- 1
cap <- 0.999
# one blue hue, light to dark with z (sequential)
z_colours <- setNames(colorRampPalette(c("#9ecae1", "#08306b"))(length(Z_LEVELS)), names(Z_LEVELS))

files <- list.files(in_root, "^community\\.csv$", recursive = TRUE, full.names = TRUE)
if (length(files) == 0) stop("no community.csv under ", in_root)
all <- rbindlist(lapply(files, function(f) {
  x <- fread(f)
  prov <- jsonlite::read_json(file.path(dirname(f), "provenance.json"))
  x[, `:=`(biome = prov$biome, MC = prov$MC, seed = prov$seed)]
}))
setorder(all, biome, MC, -u)
all[, rank := seq_len(.N), by = .(biome, MC)] # 1 = most sensitive species
long <- melt(all, measure.vars = names(Z_LEVELS), variable.name = "level", value.name = "effect",
             variable.factor = FALSE)
long[, `:=`(capped = effect >= cap, cum = pnorm(u))]
fwrite(all[, c("biome", "MC", "seed", "rank", "Species", "MaxAge", "u", names(Z_LEVELS)), with = FALSE],
       file.path(fig_dir, "ssd_species_values.csv"))

# theoretical SSD curves: P(effect <= x) = plnorm(x, mu(z), sdlog)
grid_x <- 10^seq(-4, 0, length.out = 400)
curves <- rbindlist(lapply(names(Z_LEVELS), function(l)
  data.table(level = l, x = grid_x, cum = plnorm(grid_x, mu_lognorm(Z_LEVELS[[l]], sdlog), sdlog))))

plot_community <- function(b, mc) {
  sp <- long[biome == b & MC == mc]
  tab <- all[biome == b & MC == mc]
  x_min <- min(1e-3, 10^floor(log10(min(sp$effect))))

  ssd <- ggplot() +
    theme_minimal(base_size = 9) +
    geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_line(data = curves[x >= x_min], aes(x, cum, colour = level), linewidth = 0.6) +
    geom_point(data = sp[capped == FALSE], aes(effect, cum, colour = level), size = 1.6) +
    geom_point(data = sp[capped == TRUE], aes(effect, cum, colour = level), size = 1.8, shape = 1) +
    geom_text_repel(data = sp[level == "z05"], aes(effect, cum, label = rank), size = 2.2,
                    colour = "grey20", nudge_x = 0.12, direction = "y", min.segment.length = 0,
                    segment.colour = "grey70", segment.size = 0.2, max.overlaps = Inf, seed = 1) +
    scale_x_log10(limits = c(x_min, 1.6), breaks = 10^(-4:0), labels = function(v) format(v, scientific = FALSE, drop0trailing = TRUE)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_colour_manual(values = z_colours, labels = paste0(names(Z_LEVELS), " (95th pct = ", Z_LEVELS, ")"),
                        name = NULL) +
    labs(x = "Species effect (log scale)", y = "Fraction of species with a smaller effect",
         subtitle = "SSD curves and species (numbers: rank, 1 = most sensitive; open circle: capped at 0.999)") +
    theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
    guides(colour = guide_legend(nrow = 3))

  cells <- melt(tab, id.vars = c("rank", "Species", "MaxAge"), measure.vars = c("u", names(Z_LEVELS)),
                variable.name = "column", value.name = "value", variable.factor = FALSE)
  cells[, label := fifelse(column == "u", sprintf("%.3f", value),
                           fifelse(value >= cap, sprintf("%.4f*", value), sprintf("%.4f", value)))]
  cells[column != "u", shade := log10(value)]
  cells[, column := factor(column, levels = c("u", names(Z_LEVELS)))]
  row_lab <- tab[, sprintf("%2d  %s (MaxAge %d)", rank, gsub("_", " ", Species), MaxAge)]
  cells[, row := factor(rank, levels = rev(tab$rank), labels = rev(row_lab))]

  values <- ggplot(cells, aes(column, row)) +
    theme_minimal(base_size = 9) +
    geom_tile(aes(fill = shade), colour = "white", linewidth = 0.6) +
    geom_text(aes(label = label), size = 2.4, colour = "grey10") +
    scale_fill_gradient(low = "#f7fbff", high = "#9ecae1", na.value = "grey95", guide = "none") +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL, subtitle = "Exact species effects used in the runs (* capped at 0.999)") +
    theme(panel.grid = element_blank(), axis.text.y = element_text(hjust = 0, size = 7),
          axis.text.x = element_text(face = "bold"))

  (ssd | values) + plot_layout(widths = c(1, 1.6)) +
    plot_annotation(title = sprintf("%s, community MC%02d (seed %d): species sensitivity distributions",
                                    b, mc, tab$seed[1]),
                    caption = "Lognormal SSD, sdlog = 1, 95th percentile = z; effect = exp(mu(z) + sdlog * u), same u for every z")
}

# All communities of one biome: theoretical SSD (coloured) and, for each community, the
# empirical distribution of its 20 species (grey steps; point i at the i-th smallest effect)
plot_biome <- function(b, log_x = TRUE) {
  sp <- long[biome == b][order(MC, level, effect)]
  sp[, ecdf := seq_len(.N) / .N, by = .(MC, level)]
  n_mc <- uniqueN(sp$MC)
  x_min <- if (log_x) min(1e-3, 10^floor(log10(min(sp$effect)))) else 0
  cv <- curves[x >= max(x_min, 1e-4)]
  if (!log_x) cv <- rbind(data.table(level = names(Z_LEVELS), x = 0, cum = 0), cv)
  ggplot(sp, aes(effect, ecdf)) +
    theme_minimal(base_size = 9) +
    geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_step(aes(group = MC), colour = "grey55", linewidth = 0.25, alpha = 0.7) +
    geom_point(data = sp[capped == FALSE], colour = "grey35", size = 0.7) +
    geom_point(data = sp[capped == TRUE], colour = "grey10", size = 1.6, shape = 1) +
    geom_line(data = cv, aes(x, cum, colour = level), linewidth = 0.9) +
    (if (log_x) scale_x_log10(limits = c(x_min, 1), breaks = 10^(-4:0),
                              labels = function(v) format(v, scientific = FALSE, drop0trailing = TRUE))
     else scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25))) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_colour_manual(values = z_colours, guide = "none") +
    facet_wrap(~ level, nrow = 2, labeller = as_labeller(setNames(
      paste0(names(Z_LEVELS), ": 95th percentile = ", Z_LEVELS), names(Z_LEVELS)))) +
    labs(x = paste0("Species effect", if (log_x) " (log scale)"), y = "Fraction of species with a smaller effect",
         title = sprintf("%s, %d communities of 20 species: species sensitivity distributions", b, n_mc),
         subtitle = "Coloured: theoretical SSD (lognormal, sdlog = 1). Grey: the 20 species of each community (one step line per community; open circle: capped at 0.999)") +
    theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
          panel.spacing.x = unit(0.8, "lines"))
}

biomes_only <- if (length(args) >= 3) args[-(1:2)] else NULL
for (b in if (is.null(biomes_only)) unique(all$biome) else biomes_only) {
  for (log_x in c(TRUE, FALSE))
    ggsave(file.path(fig_dir, sprintf("fig_4.4_ssd_%s_all_MC_%s.png", b, if (log_x) "log" else "linear")),
           plot_biome(b, log_x), width = 32, height = 18, units = "cm", dpi = 200, bg = "white")
}
if (!is.null(biomes_only)) quit(save = "no")

comms <- unique(all[, .(biome, MC)])
pdf_file <- file.path(fig_dir, "fig_4.4_ssd_all_communities.pdf")
pdf(pdf_file, width = 34 / 2.54, height = 15 / 2.54)
for (i in seq_len(nrow(comms))) {
  p <- plot_community(comms$biome[i], comms$MC[i])
  print(p)
  ggsave(file.path(fig_dir, sprintf("fig_4.4_ssd_%s_MC%02d.png", comms$biome[i], comms$MC[i])), p,
         width = 34, height = 15, units = "cm", dpi = 200, bg = "white")
}
invisible(dev.off())
message(nrow(comms), " communities: ", fig_dir, " (PNG per community, ", basename(pdf_file),
        ", ssd_species_values.csv)")

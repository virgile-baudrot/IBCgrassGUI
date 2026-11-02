# §4.3.1 Ranking of species in the control runs of one biome (report Fig. 4.11, one figure
# per variable: individuals, shoot mass, seedlings, seeds)
#
#   Rscript Input-EFSA/12_species_ranking.R [in_dir] [fig_dir] [top]
#   defaults: outputs/control-2026-10-28/4.3_control_design/BOR_I1/full
#             EFSA_REPORT/figures/4.3
#             30
# in_dir: one design of run_4.3_control_design.R (5 replicates, 173 x 173 cells, 50 years,
#   no herbicide). Each variable is summed over all the weeks of the run, as in the July
#   figure. Species are ordered by the sum over the 5 replicates (pooled), then the bars are
#   drawn for each replicate (one row per replicate). Only the top species are shown.

source(file.path(Sys.getenv("IBC_ROOT", "."), "Input-EFSA", "04_output_analysis.R"))

args <- commandArgs(trailingOnly = TRUE)
arg <- function(i, default) if (length(args) >= i) args[i] else default
in_dir <- arg(1, "outputs/control-2026-10-28/4.3_control_design/BOR_I1/full")
fig_dir <- arg(2, "EFSA_REPORT/figures/4.3")
top <- as.integer(arg(3, 30))
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pt <- fread(file.path(in_dir, "Pt.csv.gz"))
biome <- pt$biome[1]
design <- pt$design[1]
sums <- pt[, .(individuals = sum(Inds), shootmass = sum(shootmass), seedlings = sum(seedlings),
               seeds = sum(seeds)), by = .(PFT, rep)]
# species never present in Pt_ of a replicate get 0
sums <- sums[CJ(PFT = unique(sums$PFT), rep = unique(sums$rep)), on = c("PFT", "rep")]
for (v in c("individuals", "shootmass", "seedlings", "seeds")) set(sums, which(is.na(sums[[v]])), v, 0)

fills <- c(individuals = "steelblue", shootmass = "#F54927", seedlings = "#5AA33B", seeds = "#953BA3")
ylabs <- c(individuals = "Number of individuals", shootmass = "Shoot mass (mg)",
           seedlings = "Number of seedlings", seeds = "Number of seeds")

plot_ranking <- function(var) {
  ord <- sums[, .(tot = sum(get(var))), by = PFT][order(-tot)][seq_len(min(top, .N)), PFT]
  x <- sums[PFT %in% ord][, PFT := factor(PFT, levels = ord, labels = gsub("_", " ", ord))]
  ggplot(x, aes(x = PFT, y = .data[[var]])) +
    theme_minimal(base_size = 9) +
    geom_col(fill = fills[[var]], width = 0.8) +
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                       n.breaks = 3) +
    facet_grid(rep ~ ., labeller = labeller(rep = function(r) paste("Sim.", r))) +
    labs(x = NULL, y = paste0(ylabs[[var]], "\n(sum over the run)")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7, face = "italic"),
          axis.text.y = element_text(size = 6),
          strip.text.y = element_text(angle = 0, face = "bold"),
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          panel.spacing.y = unit(1, "mm"),
          plot.margin = margin(4, 4, 2, 22, "mm"))
}

for (v in names(fills))
  save_fig(plot_ranking(v), file.path(fig_dir, sprintf("fig_4.3_%s_%s_ranking_%s.png", biome, design, v)),
           27, 9)

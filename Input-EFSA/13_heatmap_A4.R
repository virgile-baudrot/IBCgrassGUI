# §4.3.1 Number of individuals per species and year in the control runs of one biome, sized
# for one A4 page (16 x 23 cm), the 5 replicates (same content as fig_4.2_*_heatmap_full)
#
#   Rscript Input-EFSA/13_heatmap_A4.R [in_dir] [fig_dir] [reps]
#   defaults: outputs/control-2026-10-28/4.3_control_design/BOR_I1/full
#             EFSA_REPORT/figures/4.3
#             1,2,3,4,5

source(file.path(Sys.getenv("IBC_ROOT", "."), "Input-EFSA", "04_output_analysis.R"))

args <- commandArgs(trailingOnly = TRUE)
arg <- function(i, default) if (length(args) >= i) args[i] else default
in_dir <- arg(1, "outputs/control-2026-10-28/4.3_control_design/BOR_I1/full")
fig_dir <- arg(2, "EFSA_REPORT/figures/4.3")
reps <- as.integer(strsplit(arg(3, "1,2,3,4,5"), ",")[[1]])
week <- 20
last_years <- 10
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

prov <- read_provenance(in_dir)
tmax <- prov$Tmax
pt <- fread(file.path(in_dir, "Pt.csv.gz"))[Time %% 30 == week & rep %in% reps]
# species of the pool never present in Pt_ get 0
com <- fread(file.path(in_dir, "communities.csv"))[rep %in% reps, .(rep, PFT = Species)]
full <- merge(unique(pt[, .(rep, Time)]), com, by = "rep", allow.cartesian = TRUE)
pt <- pt[full, on = c("rep", "Time", "PFT")]
pt[is.na(Inds), Inds := 0]

ord <- pt[Time > (tmax - last_years) * 30, .(m = mean(Inds)), by = PFT][order(m), PFT]
pt[, PFT := factor(PFT, levels = ord, labels = gsub("_", " ", ord))]
pt[, Inds_plot := fifelse(Inds > 0, Inds, NA_real_)]

p <- ggplot(pt, aes(x = Time %/% 30 + 1, y = PFT, fill = Inds_plot)) +
  theme_minimal(base_size = 7) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log10", na.value = "grey88", name = "Individuals (log scale)") +
  scale_x_continuous(breaks = seq(10, tmax - 10, 10), expand = c(0, 0)) +
  facet_grid(~ rep, labeller = labeller(rep = function(r) paste("Sim.", r))) +
  labs(x = "Time (years)", y = NULL) +
  theme(axis.text.y = element_text(size = 5.5, face = "italic"),
        axis.text.x = element_text(size = 5),
        strip.text = element_text(face = "bold"),
        panel.grid = element_blank(), panel.spacing = unit(1, "mm"),
        legend.position = "bottom", legend.key.height = unit(2.5, "mm"),
        legend.key.width = unit(12, "mm"), legend.title = element_text(vjust = 0.8))

save_fig(p, file.path(fig_dir, sprintf("fig_4.3_%s_%s_heatmap_A4.png", pt$biome[1], pt$design[1])), 16, 23)

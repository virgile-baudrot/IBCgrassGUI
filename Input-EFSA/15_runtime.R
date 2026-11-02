# Computing time of one IBC-grass run (50 years, controls) against the grid size and the
# number of species (Chapter 1, cost of the simulations)
#
#   Rscript Input-EFSA/15_runtime.R [grid_dir] [control_dir] [fig_dir]
#   defaults: outputs/grid-2026-10-30/4.3_control_design/BOR_I1   (grid series, no trampling)
#             outputs/control-2026-10-28/4.3_control_design       (173 x 173, trampling on)
#             EFSA_REPORT/figures/1
# Times are those written in runs.csv: wall time of one run when 5 runs share one Scaleway
# Serverless job of 6 vCPU (one run per core, AMD EPYC). A power law t = a * cells^b is fitted
# on all the grids (10 x 10 to 350 x 350, no trampling); the table extrapolates it up to 1000 x 1000 cells.
# The figure shows the grid series only (no trampling, up to 350 x 350, the largest grid run).

source(file.path(Sys.getenv("IBC_ROOT", "."), "Input-EFSA", "04_output_analysis.R"))
suppressPackageStartupMessages(library(patchwork))

args <- commandArgs(trailingOnly = TRUE)
arg <- function(i, default) if (length(args) >= i) args[i] else default
grid_dir <- arg(1, "outputs/grid-2026-10-30/4.3_control_design/BOR_I1")
ctl_dir <- arg(2, "outputs/control-2026-10-28/4.3_control_design")
fig_dir <- arg(3, "EFSA_REPORT/figures/1")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

read_runs <- function(f) {
  r <- fread(f)[ok == TRUE]
  com <- fread(file.path(dirname(f), "communities.csv"))
  r[, n_species := com[, .N, by = rep][r, on = "rep", N]]
  r[, .(biome, design, rep, CellNum, tramp, n_species, minutes = seconds / 60)]
}
grid <- rbindlist(lapply(list.files(grid_dir, "^runs\\.csv$", recursive = TRUE, full.names = TRUE), read_runs))
ctl <- rbindlist(lapply(list.files(ctl_dir, "^runs\\.csv$", recursive = TRUE, full.names = TRUE), read_runs))
bor173 <- ctl[biome == "BOR_I1" & design == "full"]

# power law on the area, all the grids
fit <- lm(log(minutes) ~ log(CellNum^2), data = grid)
b <- coef(fit)[[2]]
pred <- data.table(CellNum = sort(unique(c(seq(10, 1000, 5), 173))))
pred[, minutes := exp(predict(fit, newdata = pred))]
tab <- unique(pred[CellNum %in% c(100, 173, 250, 350, 500, 700, 1000)])[order(CellNum)]
max_side <- max(grid$CellNum)
tab[, hours := round(minutes / 60, 1)]
fwrite(rbind(grid[, .(minutes = median(minutes), n = .N, source = "measured"), by = .(CellNum, tramp)],
             tab[, .(CellNum, tramp = 0, minutes, n = NA, source = "power law")], fill = TRUE),
       file.path(fig_dir, "runtime_grid_size.csv"))

breaks <- c(10, 25, 50, 100, 173, 250, 350)
time_breaks <- c(1 / 60, 0.1, 1, 10, 60, 180, 600, 60 * 24, 60 * 72)
time_labels <- c("1 s", "6 s", "1 min", "10 min", "1 h", "3 h", "10 h", "1 day", "3 days")
pa <- ggplot(grid, aes(x = CellNum, y = minutes)) +
  geom_line(data = pred[CellNum <= max_side], linetype = "dashed", colour = "grey40") +
  geom_point(size = 1.4, alpha = 0.7) +
  scale_x_log10(breaks = breaks) +
  scale_y_log10(breaks = time_breaks, labels = time_labels) +
  labs(x = "Grid side (cells, log scale)", y = "Time of one run (log scale)",
       title = sprintf("Boreal pool (88 species), time ~ area^%.1f", b)) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

pb <- ggplot(ctl, aes(x = n_species, y = minutes)) +
  geom_point(size = 1.2, alpha = 0.5, colour = "#2c6fbb") +
  expand_limits(y = 0, x = 0) +
  labs(x = "Number of species at the start", y = "Time of one run (minutes)",
       title = "B. 173 x 173 cells, trampling on, 9 biomes") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 9, face = "bold"))

p <- pa +
  plot_annotation(caption = paste0("50 years, controls, no trampling. Wall time of one run, 5 runs in parallel ",
                                   "on a Scaleway job of 6 vCPU.\nDashed: power law fitted on all the grids."),
                  theme = theme(plot.caption = element_text(size = 6.5)))
save_fig(p, file.path(fig_dir, "fig_1_runtime.png"), 13, 10)
print(tab)

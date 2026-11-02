# §4.3.2 Effect of the grid size on dominance: Boreal pool, controls, trampling off,
# grids of 10 x 10 to 173 x 173 cells (run_4.3_control_design.R, designs g010-g173)
#
#   Rscript Input-EFSA/14_grid_size.R [in_dir] [ref_dir] [fig_dir]
#   defaults: outputs/grid-2026-10-30/4.3_control_design/BOR_I1
#             outputs/control-2026-10-28/4.3_control_design/BOR_I1/full
#             EFSA_REPORT/figures/4.3
# ref_dir: the full pool at 173 x 173 with trampling (tramp = 0.1, as in §4.4), drawn in grey
#   for comparison (other seeds). Skipped if missing.
# All values at week 20 of each year; dominance = mean number of plants over the last 10 years.

source(file.path(Sys.getenv("IBC_ROOT", "."), "Input-EFSA", "04_output_analysis.R"))
suppressPackageStartupMessages(library(patchwork))

args <- commandArgs(trailingOnly = TRUE)
arg <- function(i, default) if (length(args) >= i) args[i] else default
in_dir <- arg(1, "outputs/grid-2026-10-30/4.3_control_design/BOR_I1")
ref_dir <- arg(2, "outputs/control-2026-10-28/4.3_control_design/BOR_I1/full")
fig_dir <- arg(3, "EFSA_REPORT/figures/4.3")
week <- 20
last_years <- 10
n_top <- 25
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

read_size <- function(dir, label) {
  prov <- read_provenance(dir)
  pt <- fread(file.path(dir, "Pt.csv.gz"))[Time %% 30 == week]
  com <- fread(file.path(dir, "communities.csv"))[, .(rep, PFT = Species)]
  full <- merge(unique(pt[, .(rep, Time)]), com, by = "rep", allow.cartesian = TRUE)
  pt <- pt[full, on = c("rep", "Time", "PFT")]
  pt[is.na(Inds), Inds := 0]
  pt[, `:=`(size = label, cells = prov$CellNum^2, tmax = prov$Tmax)]
  pt[, .(size, cells, tmax, rep, Time, PFT, Inds)]
}

dirs <- sort(list.dirs(in_dir, recursive = FALSE))
dirs <- dirs[file.exists(file.path(dirs, "Pt.csv.gz"))]
sizes <- sapply(dirs, function(d) read_provenance(d)$CellNum)
labels <- paste0(sizes, " x ", sizes)
pt <- rbindlist(Map(read_size, dirs, labels))
has_ref <- file.exists(file.path(ref_dir, "Pt.csv.gz"))
ref_label <- "173 x 173, trampling on"
if (has_ref) pt <- rbind(pt, read_size(ref_dir, ref_label))
size_levels <- c(labels, if (has_ref) ref_label)
pt[, size := factor(size, levels = size_levels)]
tmax <- pt$tmax[1]

# ---------------- summaries ----------------
last <- pt[Time > (tmax - last_years) * 30, .(Inds = mean(Inds)), by = .(size, cells, rep, PFT)]
last[, share := if (sum(Inds) > 0) Inds / sum(Inds) else 0, by = .(size, rep)]
setorder(last, size, rep, -Inds)
run <- last[, .(present = sum(Inds > 0),
                n80 = which(cumsum(share) >= 0.8)[1],
                top_share = max(share),
                density = sum(Inds) / cells[1]), by = .(size, rep)]
# mean number of species with living plants in a year, last 10 years
present <- pt[Time > (tmax - last_years) * 30, .(n = sum(Inds > 0)), by = .(size, rep, Time)][
  , .(present_year = mean(n)), by = .(size, rep)]
run <- present[run, on = c("size", "rep")]
# similarity of the 5 simulations of one size: mean Bray-Curtis dissimilarity between pairs
bc <- last[, {
  m <- dcast(.SD, rep ~ PFT, value.var = "Inds", fill = 0)[, rep := NULL]
  m <- as.matrix(m)
  pr <- combn(nrow(m), 2)
  .(bc = mean(apply(pr, 2, function(k) sum(abs(m[k[1], ] - m[k[2], ])) / sum(m[k[1], ] + m[k[2], ]))))
}, by = size]
fwrite(run[bc, on = "size"][order(size, rep)], file.path(fig_dir, "grid_size_BOR_I1.csv"))

# ---------------- figure 1: summary against grid size ----------------
# x: side of the grid (log scale), trampling off only (the 173 x 173 reference with trampling is left out)
sides <- last[, .(side = sqrt(cells[1])), by = size]
run_off <- sides[run, on = "size"][size != ref_label]
bc_off <- sides[bc, on = "size"][size != ref_label]
breaks <- sort(unique(run_off$side))
point_panel <- function(d, y, ylab) {
  s <- d[, .(m = mean(get(y))), by = side]
  ggplot(d, aes(x = side, y = .data[[y]])) +
    geom_line(data = s, aes(y = m), linewidth = 0.4) +
    geom_point(size = 1.3, alpha = 0.7) +
    scale_x_log10(breaks = breaks) +
    expand_limits(y = 0) +
    labs(x = "Grid side (cells, log scale)", y = ylab) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank())
}
p1 <- point_panel(run_off, "present_year", "Species with living plants") +
  point_panel(run_off, "n80", "Species holding 80 % of the plants") +
  point_panel(run_off, "top_share", "Share of the first species") +
  point_panel(bc_off, "bc", "Bray-Curtis between simulations") +
  plot_layout(nrow = 1) +
  plot_annotation(caption = paste0("Boreal pool (88 species), controls, trampling off. Mean over the last 10 years, ",
                                   "week 20. Points: simulations; line through the means."),
                  theme = theme(plot.caption = element_text(size = 7)))
save_fig(p1, file.path(fig_dir, "fig_4.3_BOR_I1_grid_size_summary.png"), 24, 9)

# ---------------- figure 2: share of each species by grid size ----------------
sh <- last[, .(share = mean(share), n_rep = sum(Inds > 0)), by = .(size, PFT)]
ord <- sh[size == labels[length(labels)]][order(-share), PFT]
top <- unique(c(ord[seq_len(n_top)],
                sh[share >= 0.02, PFT])) # also species that dominate only at another size
sh <- sh[PFT %in% top][, PFT := factor(PFT, levels = rev(top), labels = rev(gsub("_", " ", top)))]
sh[, lab := fifelse(share >= 0.005, sprintf("%.0f", 100 * share), "")]
p2 <- ggplot(sh, aes(x = size, y = PFT, fill = 100 * share)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = lab, colour = share > 0.08), size = 2.2) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"), guide = "none") +
  scale_fill_viridis_c(option = "mako", direction = -1, limits = c(0, NA), name = "% of the plants") +
  labs(x = "Grid (cells)", y = NULL,
       caption = paste0("Mean share of the plants over the last 10 years (week 20) and the 5 simulations. ",
                        "Species ordered by their share at ", labels[length(labels)], ".")) +
  theme_minimal(base_size = 8) +
  theme(axis.text.y = element_text(face = "italic"), axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank(), legend.position = "bottom", legend.key.height = unit(2.5, "mm"),
        plot.caption = element_text(size = 6.5))
save_fig(p2, file.path(fig_dir, "fig_4.3_BOR_I1_grid_size_species.png"), 16, 16)

# ---------------- figure 3: number of species over time ----------------
rt <- pt[, .(n = sum(Inds > 0)), by = .(size, rep, Time)]
p3 <- ggplot(rt, aes(x = Time / 30, y = n, group = rep)) +
  geom_line(linewidth = 0.3, alpha = 0.7) +
  facet_wrap(~ size, nrow = if (length(size_levels) > 6) 2 else 1) +
  labs(x = "Time (years)", y = "Species with living plants (week 20)") +
  theme_minimal(base_size = 9)
save_fig(p3, file.path(fig_dir, "fig_4.3_BOR_I1_grid_size_richness_time.png"), 24,
         if (length(size_levels) > 6) 12 else 8)

# §4.3 Analysis of control simulations (choice of the community design)
#
#   Rscript Input-EFSA/run_4.3_control_design.R [biome] [design] [out_dir] [ncores]
#   biome  : default BOR_I1 (report §4.3)
#   design : full  - the whole species pool of the biome           (§4.3.1, §4.3.2)
#            nodom - the pool without the 20 dominant species       (§4.3.3)
#                    (dominants ranked from the "full" runs, which must exist)
#            r20   - 20 species drawn at random for each replicate  (§4.3.4)
#            all   - full, then nodom, then r20
#            s50, s25, s10 - 50, 25 or 10 species drawn at random for each replicate,
#                    nested: replicate i of s10 is in s25, which is in s50 (one random
#                    order of the pool per biome and replicate, its first n species)
#            nested - full, then s50, s25 and s10 (revised §4.2/§4.3, Oct 2026)
#
# Controls only (no herbicide code, HerbEffectType = 0), 173 x 173 cells, 50 years,
# n_rep replicates per design, each with its own seed (logged in runs.csv).
# The 50 x 50 grid of the first version of §4.3 is not used: with tramp = 0.1 it has no
# trampling at all (R§5.3), so it is not comparable with the 173 x 173 runs of §4.4.

source(file.path(Sys.getenv("IBC_ROOT", "."), "Input-EFSA", "02_simulation_setting.R"))

args <- commandArgs(trailingOnly = TRUE)
biome <- if (length(args) >= 1) args[1] else "BOR_I1"
design <- if (length(args) >= 2) args[2] else "all"
out_root <- if (length(args) >= 3) args[3] else "outputs/4.3_control_design"
ncores <- as.integer(if (length(args) >= 4) args[4] else parallel::detectCores())
stopifnot(biome %in% BIOMES, design %in% c("full", "nodom", "r20", "all", "s50", "s25", "s10", "nested"))

# ---------------- design parameters ----------------
p <- default_params(CellNum = 173, Tmax = 50, Tinit = 20, HerbDuration = 10, cut = 0)
n_rep <- 5
n_species <- 20 # r20 design
n_dominant <- 20 # nodom design
rank_years <- 10 # dominance = mean number of individuals at week 20 of the last 10 years
exclude_crops <- FALSE # decision D6
for (k in c("CellNum", "Tmax")) { # test runs only
  v <- Sys.getenv(paste0("TEST_", toupper(k)))
  if (v != "") p[[k]] <- as.numeric(v)
}
# ---------------------------------------------------

pool <- read_species_pool(biome, exclude_crops = exclude_crops)
design_id <- c(full = 1L, nodom = 2L, r20 = 3L, s50 = 4L, s25 = 5L, s10 = 6L)
n_nested <- c(s50 = 50L, s25 = 25L, s10 = 10L)

# rank species by their mean number of individuals at week 20 of the last rank_years years
rank_dominance <- function(pt, Tmax) {
  last <- pt[Time %% 30 == 20 & Time > (Tmax - rank_years) * 30]
  # species absent from a run have no row in Pt_: average over all runs and years, zeros included
  n_obs <- uniqueN(last[, .(rep, Time)])
  dom <- last[, .(mean_Inds = sum(Inds) / n_obs), by = PFT][order(-mean_Inds)]
  dom[, rank := .I]
  dom
}

run_design <- function(d) {
  out_dir <- file.path(out_root, biome, d)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  species <- pool
  if (d == "nodom") {
    f <- file.path(out_root, biome, "full", "dominance.csv")
    if (!file.exists(f)) stop("run the 'full' design first: ", f, " is missing")
    dom <- fread(f)
    species <- pool[!pool$Species_model %in% dom[rank <= n_dominant, PFT], ]
  }
  seeds <- 300000L + 1000L * match(biome, BIOMES) + 100L * design_id[[d]] + seq_len(n_rep)
  communities <- lapply(seq_len(n_rep), function(i) {
    if (d == "r20") {
      set.seed(seeds[i])
      return(species[sort(sample(nrow(species), n_species)), ])
    }
    if (d %in% names(n_nested)) {
      # same order for every nested design of this biome and replicate
      set.seed(400000L + 1000L * match(biome, BIOMES) + i)
      return(species[sort(sample(nrow(species))[seq_len(n_nested[[d]])]), ])
    }
    species
  })
  jobs <- lapply(seq_len(n_rep), function(i) list(
    run_dir = file.path(out_dir, "runs", sprintf("rep%02d", i)),
    fieldedge = build_Fieldedge(communities[[i]], sens = 0),
    herbfact = build_HerbFact(p$HerbDuration, "NO"),
    p = p, seed = seeds[i], herb = 0, scenario = design_id[[d]], MC = i))
  cat(biome, d, ":", n_rep, "runs of", sapply(communities, nrow)[1], "species\n")

  status <- run_many(jobs, ncores)
  status[, `:=`(biome = biome, design = d, rep = seq_len(n_rep))]
  fwrite(status, file.path(out_dir, "runs.csv"))
  fwrite(rbindlist(lapply(seq_len(n_rep), function(i)
    data.table(rep = i, Species = communities[[i]]$Species_model, MaxAge = communities[[i]]$MaxAge))),
    file.path(out_dir, "communities.csv"))
  write_provenance(out_dir, p, list(script = "run_4.3_control_design.R", biome = biome, design = d,
                                    seeds = seeds, n_species = if (d %in% names(n_nested)) n_nested[[d]] else n_species,
                                    n_dominant = n_dominant,
                                    exclude_crops = exclude_crops))
  if (!all(status$ok)) stop(sum(!status$ok), " run(s) failed in ", out_dir)

  out <- collect_outputs(status[, .(run_dir, biome, design, rep, seed)])
  fwrite(out$Pt, file.path(out_dir, "Pt.csv.gz"))
  fwrite(out$Grd, file.path(out_dir, "Grd.csv.gz"))
  if (d == "full") fwrite(rank_dominance(out$Pt, p$Tmax), file.path(out_dir, "dominance.csv"))
  unlink(file.path(out_dir, "runs"), recursive = TRUE)
}

designs <- switch(design, all = c("full", "nodom", "r20"), nested = c("full", "s50", "s25", "s10"), design)
for (d in designs) run_design(d)

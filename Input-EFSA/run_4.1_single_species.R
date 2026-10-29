# §4.1 Effect and variability of the outputs with a single species (single PFT)
#
#   Rscript Input-EFSA/run_4.1_single_species.R [part] [out_dir] [ncores]
#   part : stressors   - §4.1.1: 900 x 900 cells (81 m2), 1 control + 1 run per stressor
#          variability - §4.1.2: 173 x 173 cells (3 m2), n_rep controls and n_rep survival runs
#          all         - both
#
# Species: the first species of BE2TH_cleaned.csv (Abutilon theophrasti, an annual).
# Effect: TXT mode, sens = 1 and HerbFact = `effect` in the herbicide years for the stressor
# of the run, so every individual gets the effect `effect`.
# Herbicide years: Tinit+1 ... Tinit+HerbDuration (21-30), recovery to Tmax (40).
# Seeds: in "stressors" all 7 runs share one seed; in "variability" replicate i of the
# control and replicate i of the treatment share seed i (identical up to year 20), and the
# spread across replicates shows the stochasticity of the model.

source(file.path(Sys.getenv("IBC_ROOT", "."), "Input-EFSA", "02_simulation_setting.R"))

args <- commandArgs(trailingOnly = TRUE)
part <- if (length(args) >= 1) args[1] else "all"
out_root <- if (length(args) >= 2) args[2] else "outputs/4.1_single_species"
ncores <- as.integer(if (length(args) >= 3) args[3] else parallel::detectCores())
stopifnot(part %in% c("stressors", "variability", "all"))

# ---------------- design parameters ----------------
species_row <- 1 # Abutilon theophrasti
effect <- 0.9 # effect of the stressor on every individual, each herbicide year
p <- default_params(Tmax = 40, Tinit = 20, HerbDuration = 10, cut = 1)
cellnum <- c(stressors = 900, variability = 173)
n_rep <- 30 # variability part
variability_stressor <- "survival"
for (k in c("Tmax")) { # test runs only
  v <- Sys.getenv(paste0("TEST_", toupper(k)))
  if (v != "") p[[k]] <- as.numeric(v)
}
if (Sys.getenv("TEST_CELLNUM") != "") cellnum[] <- as.numeric(Sys.getenv("TEST_CELLNUM"))
# ---------------------------------------------------

species <- read_species_pool()[species_row, ]
cat("species:", species$Species, "MaxAge", species$MaxAge, "\n")
fieldedge <- build_Fieldedge(species, sens = 1)

run_part <- function(pt) {
  out_dir <- file.path(out_root, pt)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pp <- modifyList(p, list(CellNum = cellnum[[pt]]))
  if (pt == "stressors") {
    design <- data.table(stressor = c("NO", STRESSORS), rep = 1L, seed = 410001L)
  } else {
    design <- CJ(stressor = c("NO", variability_stressor), rep = seq_len(n_rep))
    design[, seed := 420000L + rep]
  }
  design[, herb := as.integer(stressor != "NO")]
  design[, run_dir := file.path(out_dir, "runs", sprintf("%s_rep%02d", stressor, rep))]
  jobs <- lapply(seq_len(nrow(design)), function(i) with(design[i], list(
    run_dir = run_dir, fieldedge = fieldedge,
    herbfact = build_HerbFact(pp$HerbDuration, stressor, level = effect),
    p = pp, seed = seed, herb = herb, scenario = match(stressor, c("NO", STRESSORS)) - 1, MC = rep)))

  status <- run_many(jobs, ncores)
  status[, `:=`(stressor = design$stressor, rep = design$rep)]
  fwrite(status, file.path(out_dir, "runs.csv"))
  write_provenance(out_dir, pp, list(script = "run_4.1_single_species.R", part = pt,
                                     species = species$Species, effect = effect))
  if (!all(status$ok)) stop(sum(!status$ok), " run(s) failed in ", out_dir)

  out <- collect_outputs(design[, .(run_dir, stressor, rep, seed)])
  fwrite(out$Pt, file.path(out_dir, "Pt.csv.gz"))
  fwrite(out$Grd, file.path(out_dir, "Grd.csv.gz"))
  unlink(file.path(out_dir, "runs"), recursive = TRUE)
}

for (pt in if (part == "all") c("stressors", "variability") else part) run_part(pt)

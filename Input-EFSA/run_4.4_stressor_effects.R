# §4.4 Stressor effects on 20-species communities, 9 BE2TH biomes
#
# One call = one community (biome, MC): 1 control + 1 null run + 9 SSD levels x 6 stressors
# = 56 runs.
#   Rscript Input-EFSA/run_4.4_stressor_effects.R BOR_I1 1 [out_dir] [ncores]
# (or environment variables BIOME, MC, OUT_DIR, NCORES, as in the Scaleway job)
# The full design is 9 biomes x nMC communities, i.e. 9 * nMC calls (see scaleway/).
#
# Extra levels only (environment variable LEVELS, e.g. LEVELS=z06,z07,z08,z09): the
# treatments of these levels only, without control and null run (6 runs per level). The
# community, seed and u are those of the full design, so the runs are paired with the
# control of an earlier submission. The 26 Oct 2026 submission (main-2026-10-26) had z01-z05;
# z06-z09 were run on 27 Oct 2026 and joined to it with:
#   Rscript Input-EFSA/run_4.4_stressor_effects.R merge <old_root> <new_root> <out_root>
#
# Design
# - Community: 20 species drawn at random from the biome pool, the same for all runs.
# - Seed: one per community (seed = 1000 * biome index + MC). All runs use the same
#   RNG seed, so control and treatments are identical up to the first herbicide year.
#   From the first herbicide event on, the treatment draws different random numbers
#   (an affected plant dies or grows differently), so its random stream and that of the
#   control diverge: after year Tinit, treatment / control also contains the
#   stochasticity of the model, even for a negligible effect.
# - Null run (level z00, stressor "null"): a control (no herbicide) whose RNG is seeded
#   again (reseed = seed + 500000) at the start of year Tinit+1. It shares the state of
#   the control at the end of year Tinit and then follows another random stream, as a
#   treatment does. null / control is the ratio expected for a zero effect: its interval
#   over communities is the reference against which the treatment effects are read.
# - SSD: one standard normal draw u per species; the species effect at level z is
#   exp(mu(z) + sdlog * u), capped at 0.999 (95th percentile of the SSD = z). A species
#   has the same effect for the 6 stressors (one sensitivity per species).
# - Effect: TXT mode, HerbFact = 1 in the herbicide years for the stressor of the run,
#   so the effect of species i is sens[i] every herbicide year.
# - Herbicide years: Tinit+1 ... Tinit+HerbDuration (21-30), then recovery to Tmax.

source(file.path(Sys.getenv("IBC_ROOT", "."), "Input-EFSA", "02_simulation_setting.R"))

args <- commandArgs(trailingOnly = TRUE)

# levels run by this call: all, or those of the environment variable LEVELS
run_levels <- function() {
  v <- Sys.getenv("LEVELS")
  lv <- if (v == "") names(Z_LEVELS) else strsplit(v, ",")[[1]]
  stopifnot(all(lv %in% names(Z_LEVELS)))
  lv
}

# Join extra levels (run with LEVELS) to an earlier submission, community by community:
#   Rscript Input-EFSA/run_4.4_stressor_effects.R merge <old_root> <new_root> <out_root>
# Checks before writing: same species and same u, the new runs use the same seed, and every new
# treatment equals the old control (Grd_ and Pt_) up to the end of year Tinit (pairing).
if (length(args) >= 1 && args[1] == "merge") {
  stopifnot(length(args) == 4)
  old_root <- args[2]; new_root <- args[3]; out_root <- args[4]
  dirs <- sub(paste0("^", new_root, "/"), "", dirname(list.files(new_root, "^runs\\.csv$", recursive = TRUE,
                                                                 full.names = TRUE)))
  for (d in dirs) {
    od <- file.path(old_root, d); nd <- file.path(new_root, d); md <- file.path(out_root, d)
    oc <- fread(file.path(od, "community.csv")); nc <- fread(file.path(nd, "community.csv"))
    stopifnot(identical(oc$Species, nc$Species), identical(oc$u, nc$u),
              isTRUE(all.equal(oc, nc[, names(oc), with = FALSE], tolerance = 1e-12)))
    oruns <- fread(file.path(od, "runs.csv")); nruns <- fread(file.path(nd, "runs.csv"))
    stopifnot(all(nruns$ok), uniqueN(c(oruns$seed, nruns$seed)) == 1,
              !any(nruns$level %in% oruns$level))
    tinit_weeks <- nruns$Tinit[1] * 30
    dir.create(md, recursive = TRUE, showWarnings = FALSE)
    for (f in c("Grd", "Pt")) {
      o <- fread(file.path(od, paste0(f, ".csv.gz"))); n <- fread(file.path(nd, paste0(f, ".csv.gz")))
      # Grd_ columns that repeat an end-of-run value on every row are not compared
      skip <- c("biome", "MC", "seed", "reseed", "level", "stressor",
                "meanShannon", "meanNPFT", "NNonClonal", "NClonal", "NGenets")
      cols <- setdiff(intersect(names(o), names(n)), skip)
      ctrl <- o[stressor == "NO" & Time <= tinit_weeks, ..cols]
      key <- intersect(c("Time", "PFT"), cols)
      for (r in split(n[Time <= tinit_weeks], by = c("level", "stressor"))) {
        a <- setorderv(copy(ctrl), key); b <- setorderv(r[, ..cols], key)
        if (!isTRUE(all.equal(a, b, check.attributes = FALSE)))
          stop(d, " ", f, ": ", r$level[1], " ", r$stressor[1], " differs from the control before Tinit")
      }
      fwrite(rbind(o, n, use.names = TRUE), file.path(md, paste0(f, ".csv.gz")))
    }
    fwrite(nc, file.path(md, "community.csv"))
    fwrite(rbind(oruns, nruns, use.names = TRUE), file.path(md, "runs.csv"))
    file.copy(file.path(od, "provenance.json"), file.path(md, "provenance.json"), overwrite = TRUE)
    file.copy(file.path(nd, "provenance.json"), file.path(md, "provenance_z06-z09.json"), overwrite = TRUE)
    cat(d, "merged\n")
  }
  quit(status = 0)
}

# Integrity check of the whole design (step S4.3):
#   Rscript Input-EFSA/run_4.4_stressor_effects.R check [out_dir] [nMC]
# Every (biome, MC) must have its 56 successful runs (6 per level with LEVELS) of Tmax*30 weeks
# and its Pt/Grd files.
if (length(args) >= 1 && args[1] == "check") {
  out_root <- if (length(args) >= 2) args[2] else "outputs/4.4_stressor_effects"
  nMC <- if (length(args) >= 3) as.integer(args[3]) else 10L
  levels_run <- run_levels()
  n_expected <- 6L * length(levels_run) + if (setequal(levels_run, names(Z_LEVELS))) 2L else 0L
  expected <- CJ(biome = BIOMES, MC = seq_len(nMC))
  res <- expected[, {
    d <- file.path(out_root, biome, sprintf("MC%02d", MC))
    f <- file.path(d, "runs.csv")
    r <- if (file.exists(f)) fread(f) else NULL
    .(n_runs = if (is.null(r)) 0L else nrow(r),
      n_ok = if (is.null(r)) 0L else sum(r$ok),
      weeks = if (is.null(r)) NA_character_ else paste(unique(r$n_weeks), collapse = "/"),
      outputs = file.exists(file.path(d, "Pt.csv.gz")) && file.exists(file.path(d, "Grd.csv.gz")))
  }, by = .(biome, MC)]
  res[, complete := n_runs == n_expected & n_ok == n_expected & outputs]
  fwrite(res, file.path(out_root, "check.csv"))
  print(res[, .(complete = sum(complete), missing_or_failed = sum(!complete)), by = biome])
  if (any(!res$complete)) print(res[complete == FALSE])
  quit(status = as.integer(any(!res$complete)))
}

biome <- if (length(args) >= 1) args[1] else Sys.getenv("BIOME")
MC <- as.integer(if (length(args) >= 2) args[2] else Sys.getenv("MC"))
out_root <- if (length(args) >= 3) args[3] else Sys.getenv("OUT_DIR", "outputs/4.4_stressor_effects")
ncores <- as.integer(if (length(args) >= 4) args[4] else Sys.getenv("NCORES", parallel::detectCores()))
stopifnot(biome %in% BIOMES, !is.na(MC), MC >= 1)

# ---------------- design parameters (report §3.1.3 and §4.4) ----------------
p <- default_params(CellNum = 173, Tmax = 50, Tinit = 20, HerbDuration = 10, cut = 0)
n_species <- 20
sdlog <- 1
exclude_crops <- FALSE # decision D6
# test runs only: e.g. TEST_CELLNUM=50 TEST_TMAX=6 TEST_TINIT=2 TEST_HERBDURATION=2
for (k in c("CellNum", "Tmax", "Tinit", "HerbDuration")) {
  v <- Sys.getenv(paste0("TEST_", toupper(k)))
  if (v != "") p[[k]] <- as.numeric(v)
}
# -----------------------------------------------------------------------------

seed <- 1000L * match(biome, BIOMES) + MC
out_dir <- file.path(out_root, biome, sprintf("MC%02d", MC))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("biome", biome, "MC", MC, "seed", seed, "->", out_dir, "\n")

# community and species sensitivities
set.seed(seed)
pool <- read_species_pool(biome, exclude_crops = exclude_crops)
community <- pool[sort(sample(nrow(pool), n_species)), ]
u <- rnorm(n_species)
sens <- sapply(Z_LEVELS, function(z) ssd_effect(u, z, sdlog))
fwrite(data.table(Species = community$Species_model, MaxAge = community$MaxAge, u = u, sens),
       file.path(out_dir, "community.csv"))

# the 56 runs (or the treatments of LEVELS)
reseed <- seed + 500000L
levels_run <- run_levels()
runs <- if (setequal(levels_run, names(Z_LEVELS)))
  list(list(level = "z00", stressor = "NO", sens = 0, herb = 0, scenario = 0, reseed = 0),
       list(level = "z00", stressor = "null", sens = 0, herb = 0, scenario = 1, reseed = reseed)) else list()
for (i in match(levels_run, names(Z_LEVELS))) for (j in seq_along(STRESSORS)) {
  runs[[length(runs) + 1]] <- list(level = names(Z_LEVELS)[i], stressor = STRESSORS[j],
                                   sens = sens[, i], herb = 1, scenario = 10 * i + j, reseed = 0)
}
jobs <- lapply(runs, function(r) list(
  run_dir = file.path(out_dir, "runs", paste(r$level, r$stressor, sep = "_")),
  fieldedge = build_Fieldedge(community, sens = r$sens, stressor = r$stressor),
  herbfact = build_HerbFact(p$HerbDuration, r$stressor, level = 1),
  p = p, seed = seed, herb = r$herb, scenario = r$scenario, MC = MC, reseed = r$reseed))

status <- run_many(jobs, ncores)
status[, `:=`(biome = biome, level = sapply(runs, `[[`, "level"), stressor = sapply(runs, `[[`, "stressor"))]
fwrite(status, file.path(out_dir, "runs.csv"))
write_provenance(out_dir, p, list(script = "run_4.4_stressor_effects.R", biome = biome, MC = MC,
                                  seed = seed, reseed_null = reseed, levels = levels_run,
                                  n_species = n_species, sdlog = sdlog,
                                  exclude_crops = exclude_crops))

if (!all(status$ok)) {
  print(status[ok == FALSE, .(level, stressor, exit_status, n_weeks)])
  stop(sum(!status$ok), " run(s) failed for ", biome, " MC ", MC, ". Run directories kept for inspection.")
}

out <- collect_outputs(status[, .(run_dir, biome, MC, seed, reseed, level, stressor)])
fwrite(out$Pt, file.path(out_dir, "Pt.csv.gz"))
fwrite(out$Grd, file.path(out_dir, "Grd.csv.gz"))
unlink(file.path(out_dir, "runs"), recursive = TRUE)
cat("done:", nrow(status), "runs,", round(sum(status$seconds) / 60), "CPU minutes\n")

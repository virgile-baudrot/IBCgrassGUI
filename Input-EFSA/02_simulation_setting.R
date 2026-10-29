# Shared functions for the EFSA simulations of the report (Ch. 4)
#   run_4.1_single_species.R   -> §4.1
#   run_4.3_control_design.R   -> §4.3
#   run_4.4_stressor_effects.R -> §4.4
# Every run has its own directory, a logged RNG seed (19th model argument) and is
# checked after the run (exit status 0 and Tmax*30 weeks in Grd_). A failed run is
# logged and kept as failed: communities are never re-drawn after a failure.

library(data.table)

BIOMES <- c("STE_I1", "BLS_I1", "ALP_I1", "ATL_I1", "BOR_I1", "CON_I1", "PAN_I1", "MED_I1", "MAC_I1")
STRESSORS <- c("biomass", "SEbiomass", "survival", "establishment", "sterility", "seednumber")
# z06-z09 were added on 27 Oct 2026 (run separately, see LEVELS in run_4.4_stressor_effects.R):
# the community and the SSD draw u of a (biome, MC) do not depend on the number of levels.
Z_LEVELS <- c(z01 = 0.1, z02 = 0.2, z03 = 0.3, z04 = 0.4, z05 = 0.5,
              z06 = 0.6, z07 = 0.7, z08 = 0.8, z09 = 0.9)

# root of the repository (scripts are run from it, or from the Docker image /opt/IBCgrass)
ibc_root <- function() normalizePath(Sys.getenv("IBC_ROOT", "."))
ibc_bin <- function() file.path(ibc_root(), "Model-files", "IBCgrassGUI")

# git commit of the code used (IBC_COMMIT is set in the Docker image, where .git is absent)
ibc_commit <- function() {
  commit <- Sys.getenv("IBC_COMMIT", "")
  if (commit == "") {
    commit <- tryCatch(system("git rev-parse --short HEAD", intern = TRUE), error = function(e) "unknown")
    dirty <- tryCatch(length(system("git status --porcelain -- Model-files Input-EFSA", intern = TRUE)) > 0,
                      error = function(e) FALSE)
    if (dirty) commit <- paste0(commit, "-dirty")
  }
  commit
}

# ------------------------------------------------------------------
# Species pool and Fieldedge.txt
# ------------------------------------------------------------------

read_species_pool <- function(biome = NULL, exclude_crops = FALSE,
                              file = file.path(ibc_root(), "Input-EFSA", "BE2TH_cleaned.csv")) {
  df <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
  if (!is.null(biome)) df <- df[df[[biome]] == 1, ]
  if (exclude_crops) df <- df[df$Category != "Crop", ]
  df$Species_model <- clean_species_name(df$Species)
  stopifnot(!anyDuplicated(df$Species_model)) # the model silently drops duplicated names
  df
}

clean_species_name <- function(x) {
  x <- gsub(" ", "_", x)
  x <- gsub("\\.", "", x)
  iconv(x, "latin1", "ASCII//TRANSLIT", sub = "")
}

# Fieldedge.txt for the species in df.
# sens: effect of each species (0-1). In TXT mode (EffectModel = 0) the effect of species i
#   in herbicide year t is HerbFact[t, stressor] * sens[i].
# stressor: the EC50 column of this stressor is the EC50 giving the same effect in DR mode
#   (EffectModel = 2, rate 1, slope 4); the other EC50 are 0 ("no effect" in the C++ code).
#   The EC50 columns are not used in TXT mode.
build_Fieldedge <- function(df, sens = 1, stressor = "NO") {
  n <- nrow(df)
  sens <- rep_len(sens, n)
  ec50 <- function(s) if (s == stressor) build_dose_response(sens) else 0
  data.frame(
    ID = 1:n,
    Species = df$Species_model,
    MaxAge = df$MaxAge, # 1 annual, 2 biennial, 100 perennial (see 00_life_span.R)
    AllocSeed = 0.05,
    LMR = df$LMR,
    m0 = df$mSeed..mg.,
    MaxMass = df$MaxMass..mg.,
    mSeed = df$mSeed..mg.,
    Dist = df$Dist..m.,
    pEstab = 0.5,
    Gmax = df$Gmax,
    SLA = df$SLA.comb.,
    palat = 0.5, # not used without grazing (graz = 0), must be > 0
    memo = df$Memo..weeks.,
    RAR = 1,
    growth = 0.25,
    mThres = 0.2,
    clonal = df$Clonal,
    propSex = df$Clonal, # read but not used by the model
    meanSpacerLength = df$meanSpacerLength..cm.,
    sdSpacerLength = ifelse(df$sdSpacerLength..cm. == 17.5, 12.5, df$sdSpacerLength..cm.),
    Resshare = df$Resshare,
    AllocSpacer = ifelse(df$Clonal == 1, 0.05, 0),
    mSpacer = ifelse(df$Clonal == 1, 70, 0),
    sens = sens,
    allocroot = 1,
    allocshoot = 1,
    EC50_biomass = ec50("biomass"), slope_biomass = 4,
    EC50_SEbiomass = ec50("SEbiomass"), slope_SEbiomass = 4,
    EC50_survival = ec50("survival"), slope_survival = 4,
    EC50_establishment = ec50("establishment"), slope_establishment = 4,
    EC50_sterility = ec50("sterility"), slope_sterility = 4,
    EC50_seednumber = ec50("seednumber"), slope_seednumber = 4,
    FlowerWeek = 16, # default in publications
    DispWeek = 20, # FlowerWeek + 4
    GermPeriod = 1
  )
}

# ------------------------------------------------------------------
# Effects: species sensitivity distribution (SSD) and HerbFact.txt
# ------------------------------------------------------------------

# meanlog of a lognormal whose 95th percentile is z
mu_lognorm <- function(z, sdlog = 1) {
  log(z) - qnorm(0.95) * sdlog
}

# Species effect from the SSD. u is one standard normal draw per species, drawn once per
# community: the same u for every z level keeps the ranking of species sensitivity across
# levels (z is a gradient within a community). Effects of 1 have no finite EC50: cap below 1.
ssd_effect <- function(u, z, sdlog = 1, cap = 0.999) {
  pmin(exp(mu_lognorm(z, sdlog) + sdlog * u), cap)
}

ec50_foo <- function(dose_response, rate, slope) {
  (rate^slope / dose_response - rate^slope)^(1 / slope)
}

# EC50 giving the species effect at application rate 1 and slope 4 (EffectModel = 2).
# Effect 0 -> EC50 = 0, which means "no effect" in the C++ model.
build_dose_response <- function(effect) {
  rxp <- ec50_foo(effect, 1, 4)
  rxp[effect == 0] <- 0
  rxp
}

# HerbFact.txt: one row per herbicide year, one column per stressor (C++ column order).
# The column of the stressor gets `level` every year, the others 0. stressor "NO": all 0.
build_HerbFact <- function(HerbDuration, stressor = "NO", level = 1) {
  col <- function(s) rep(if (s == stressor) level else 0, HerbDuration)
  data.frame(
    Biomass = col("biomass"),
    Mortality = col("survival"),
    SeedlingBiomass = col("SEbiomass"),
    Establishment = col("establishment"),
    SeedSterility = col("sterility"),
    SeedNumber = col("seednumber")
  )
}

# ------------------------------------------------------------------
# Running the model
# ------------------------------------------------------------------

# Model parameters (command-line arguments 1-16). Scripts override what they need.
default_params <- function(...) {
  p <- list(
    ModelVersion = 3,
    CellNum = 173, # 173 x 173 cells of 1 cm2 = 3 m2
    Tmax = 50, # years
    Tinit = 20, # years before the first herbicide year
    SeedInput = 10, # seeds per species and year from outside
    belowres = 90,
    abres = 100,
    abampl = 0,
    tramp = 0.1, # weekly trampling (see R§5.3)
    graz = 0,
    cut = 0,
    week_start = 1, # week of herbicide effect in each herbicide year
    HerbDuration = 10, # herbicide years: Tinit+1 ... Tinit+HerbDuration
    EffectModel = 0 # 0 = TXT (HerbFact.txt x sens)
  )
  modifyList(p, list(...))
}

# One model run in its own directory. Returns one row (data.table) describing the run.
# herb: 0 = control (no herbicide code at all), 1 = treatment
# reseed: if not 0, the model draws a new random stream from the start of year Tinit+1
# (20th model argument; null runs of §4.4)
run_ibc <- function(run_dir, fieldedge, herbfact, p, seed, herb = 1, scenario = 1, MC = 1,
                    reseed = 0, guard = 3600) {
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  # the model appends to its output files: remove anything left by an earlier attempt
  unlink(list.files(run_dir, "^(Pt|Grd|Ind)_.*\\.txt$", full.names = TRUE))
  write.table(fieldedge, file.path(run_dir, "Fieldedge.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(herbfact, file.path(run_dir, "HerbFact.txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
  args <- c(p$ModelVersion, p$CellNum, p$Tmax, p$Tinit, "Fieldedge.txt", p$SeedInput,
            p$belowres, p$abres, p$abampl, p$tramp, p$graz, p$cut, p$week_start,
            p$HerbDuration, herb, p$EffectModel, scenario, MC, seed, if (reseed != 0) reseed)
  cmd <- sprintf("cd %s && timeout %d %s %s > model.log 2>&1",
                 shQuote(run_dir), guard, shQuote(ibc_bin()), paste(args, collapse = " "))
  elapsed <- system.time(status <- system(cmd))[["elapsed"]]
  grd <- list.files(run_dir, "^Grd_.*\\.txt$", full.names = TRUE)
  pt <- list.files(run_dir, "^Pt_.*\\.txt$", full.names = TRUE)
  n_weeks <- if (length(grd) == 1) length(readLines(grd)) - 1L else NA_integer_
  ok <- status == 0 && length(grd) == 1 && length(pt) == 1 && identical(n_weeks, as.integer(p$Tmax * 30))
  row <- data.table(run_dir = run_dir, seed = seed, reseed = reseed, herb = herb, scenario = scenario, MC = MC,
                    exit_status = status, n_weeks = n_weeks, seconds = round(elapsed),
                    ok = ok, as.data.table(p))
  fwrite(row, file.path(run_dir, "run.csv"))
  row
}

# Run a list of runs (each a list of arguments of run_ibc) on ncores cores.
run_many <- function(runs, ncores = parallel::detectCores()) {
  res <- parallel::mclapply(runs, function(r) do.call(run_ibc, r),
                            mc.cores = ncores, mc.preschedule = FALSE)
  failed <- vapply(res, inherits, logical(1), "try-error")
  if (any(failed)) stop("R error in ", sum(failed), " run(s): ", res[failed][[1]])
  rbindlist(res)
}

# Read the Pt_ and Grd_ outputs of successful runs and add the run descriptors.
# meta: data.table with one row per run, a column run_dir and the descriptor columns.
collect_outputs <- function(meta) {
  read_one <- function(prefix, d) {
    f <- list.files(d, paste0("^", prefix, "_.*\\.txt$"), full.names = TRUE)
    fread(f, sep = "\t")
  }
  desc <- setdiff(names(meta), "run_dir")
  pt <- rbindlist(lapply(seq_len(nrow(meta)), function(i)
    cbind(meta[i, ..desc], read_one("Pt", meta$run_dir[i]))))
  grd <- rbindlist(lapply(seq_len(nrow(meta)), function(i)
    cbind(meta[i, ..desc], read_one("Grd", meta$run_dir[i]))))
  list(Pt = pt, Grd = grd)
}

# Write a session record next to the outputs (code version, parameters, R session).
write_provenance <- function(out_dir, p, extra = list()) {
  info <- c(list(date = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                 commit = ibc_commit(),
                 model_md5 = unname(tools::md5sum(ibc_bin())),
                 R = R.version.string),
            p, extra)
  jsonlite::write_json(info, file.path(out_dir, "provenance.json"), auto_unbox = TRUE, pretty = TRUE)
}

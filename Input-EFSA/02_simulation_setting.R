# LIBRARIES
#library(foreach)
#library(doParallel)
#library(labeling)

# Function to build application rate vector for IBCgrass
# application_rate: numeric, application rate value
# duration: integer, number of weeks the application rate is applied
build_AppRate <- function(application_rate, herbicide_duration){
    data.frame(rep(application_rate, herbicide_duration))
}

# Function to save simulation configuration to an RDS file
save_configuration <- function(filepath, ModelVersion, CellNum, Tmax, InitDuration, NamePftFile,
                               SeedInput, belowres, abres, abampl, tramp, graz, cut,
                               week_start, HerbDuration, HerbEffectType,
                               EffectModel, scenario, MC, nMC){
    config <- list(
        ModelVersion = ModelVersion,
        CellNum = CellNum,
        Tmax = Tmax,
        InitDuration = InitDuration,
        NamePftFile = NamePftFile,
        SeedInput = SeedInput,
        belowres = belowres,
        abres = abres,
        abampl = abampl,
        tramp = tramp,
        graz = graz,
        cut = cut,
        week_start = week_start,
        HerbDuration = HerbDuration,
        HerbEffectType = HerbEffectType,
        EffectModel = EffectModel,
        scenario = scenario,
        MC = MC,
        nMC = nMC
    )
    dfconifg = as.data.frame(t(unlist(config)))
    saveRDS(config, file = paste0(filepath, ".rds"))
    write.csv(dfconifg, file = paste0(filepath, ".csv"), row.names=FALSE)
}


ec50_foo = function(dose_response, rate, slope){
    ec50 = (rate^slope/dose_response - rate^slope)^(1/slope)
    return(ec50)
}
# Species Sensitivity Distribution (SSD): one effect level (0-1) per species
build_effect_dist <- function(n, stressor_type="low"){
    if(stressor_type=="z00"){
      out_dist = rep(0, n)
    }
    if(stressor_type=="high"){
      out_dist = runif(n, 0.25, 0.50)
    }
    if(stressor_type=="low"){
      out_dist = runif(n, 0.01, 0.25)
    }
    if(stressor_type=="medium"){
      out_dist = pmax(rnorm(n, 0.25, 0.05), 0)
    }
    if(stressor_type=="z01"){
      out_dist = rlnorm(n, meanlog = mu_lognorm(0.1), sdlog = 1)
    }
    if(stressor_type=="z02"){
      out_dist = rlnorm(n, meanlog = mu_lognorm(0.2), sdlog = 1)
    }
    if(stressor_type=="z03"){
      out_dist = rlnorm(n, meanlog = mu_lognorm(0.3), sdlog = 1)
    }
    if(stressor_type=="z04"){
      out_dist = rlnorm(n, meanlog = mu_lognorm(0.4), sdlog = 1)
    }
    if(stressor_type=="z05"){
      out_dist = rlnorm(n, meanlog = mu_lognorm(0.5), sdlog = 1)
    }
    # an effect of 1 (100%) has no finite EC50, cap just below
    return(pmin(out_dist, 0.999))
}
# EC50 giving the species effect at application rate 1 and slope 4 (EffectModel=2)
build_dose_response <- function(effect){
    rxp = ec50_foo(effect,1,4)
    # no effect: EC50 = 0 means "no effect" in the C++ model
    rxp[effect == 0] = 0
  return(rxp)
}

build_id_scenario <- function(be2th, stress_type, stressor){
    ids_100 = 100*switch(be2th,
                    "STE_I1"=1,
                    "BLS_I1"=2,
                    "ALP_I1"=3,
                    "ATL_I1"=4,
                    "BOR_I1"=5,
                    "CON_I1"=6,
                    "PAN_I1"=7,
                    "MED_I1"=8,
                    "MAC_I1"=9
    )
    ids_10 = 10*switch(stress_type,
                   "low"=1,
                   "high"=2,
                   "medium"=3,
                   "z01" = 4,
                   "z02" = 5,
                   "z03" = 6,
                   "z04" = 7,
                   "z05" = 8,
                   "z00" = 9,
    )
    ids_1 = switch(stressor,
        "NO"=0,
        "biomass"=1, 
        "SEbiomass"=2, 
        "survival"=3, 
        "establishment"=4, 
        "sterility"=5,
        "seednumber"=6)
    return(ids_100+ids_10+ids_1)
}

mu_lognorm <- function(z, sdlog=1){
  log(z) - 1.645 * sdlog
}

# Function to build HerbFact data frame for EFSA simulations
# The effect of species i in year t is HerbFact[t] * sens[i] (C++ EffectModel=0).
# The species effect (SSD) is in sens (see build_Fieldedge_BE2TH), so the yearly
# factor is 1: same effect every year, as in the dose-response mode with rate=1.
build_HerbFact <- function(HerbDuration, stressor, stressor_type){

    level_vec = rep(ifelse(stressor_type=="z00", 0, 1), HerbDuration)
    no_vec = rep(0, HerbDuration)
    
    data.frame(
        Biomass=if("biomass" %in% stressor){level_vec}else{no_vec},
        Mortality=if("survival" %in% stressor){level_vec}else{no_vec},
        SeedlingBiomass=if("SEbiomass" %in% stressor){level_vec}else{no_vec},
        Establishment=if("establishment" %in% stressor){level_vec}else{no_vec},
        SeedSterility=if("sterility" %in% stressor){level_vec}else{no_vec},
        SeedNumber=if("seednumber" %in% stressor){level_vec}else{no_vec}
    )
}

# Function to build FieldEdge data frame for EFSA simulations using BE2TH database
# df: data frame with plant traits
# stressor: character vector with stressor names
build_Fieldedge_BE2TH = function(dfBE2TH, group = "", stressor="NO", stressor_type="low"){

    df = dfBE2TH[dfBE2TH[[group]]==1,]
    n = nrow(df)
    # one effect per species from the SSD, used both as sens (txt, EffectModel=0)
    # and as EC50 (dose-response, EffectModel=2) so that both modes give the same effect
    effect = build_effect_dist(n, stressor_type)
    rxp = build_dose_response(effect)
    
    SPECIES = gsub(" ", "_", df$Species)
    SPECIES = gsub("\\.", "", SPECIES)
    SPECIES = iconv(SPECIES, "latin1", "ASCII//TRANSLIT", sub="")

    df <- data.frame(
        ID=1:n,
        Species=SPECIES,
        MaxAge=df$MaxAge, # 1 annual, 2 biennial, 100 perennial (see 00_life_span.R)
        AllocSeed=0.05,
        LMR=df$LMR,
        m0=df$mSeed..mg.,
        MaxMass=df$MaxMass..mg.,
        mSeed=df$mSeed..mg.,
        Dist=df$Dist..m.,
        pEstab=0.5,
        Gmax=df$Gmax,
        SLA=df$SLA.comb.,
        # palat=df$Palat,
        palat=0.5, # cannot be 0.0,  exclude after
        memo=df$Memo..weeks.,
        RAR=1,
        growth=0.25,
        mThres=0.2,
        clonal=df$Clonal,
        propSex=df$Clonal, # no error
        meanSpacerLength=df$meanSpacerLength..cm.,
        sdSpacerLength=ifelse(df$sdSpacerLength..cm. == 17.5, 12.5,  df$sdSpacerLength..cm.),
        Resshare=df$Resshare,
        AllocSpacer=ifelse(df$Clonal==1, 0.05 , 0),
        mSpacer=ifelse(df$Clonal==1, 70, 0),
        sens=effect, # species effect in txt mode (multiplies HerbFact.txt)
        allocroot=1,
        allocshoot=1,
        EC50_biomass = if("biomass" %in% stressor){rxp}else{0},
        slope_biomass = 4,
        EC50_SEbiomass = if("SEbiomass" %in% stressor){rxp}else{0},
        slope_SEbiomass = 4,
        EC50_survival = if("survival" %in% stressor){rxp}else{0},
        slope_survival = 4,
        EC50_establishment = if("establishment" %in% stressor){rxp}else{0},
        slope_establishment = 4,
        EC50_sterility = if("sterility" %in% stressor){rxp}else{0},
        slope_sterility = 4,
        EC50_seednumber = if("seednumber" %in% stressor){rxp}else{0},
        slope_seednumber = 4,
        FlowerWeek=16, # DEFAULT IN PUBLICATIONS
        DispWeek=20, # IS = FlowerWeek+4
        GermPeriod=1
        # Overwintering=1
    )
    return(df)
}
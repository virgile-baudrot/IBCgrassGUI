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
build_dose_response <- function(n, stressor_type="low"){
    if(stressor_type=="high"){
        out_dist = runif(n, 0.25, 0.50)
    }
    if(stressor_type=="low"){
        out_dist = runif(n, 0.01, 0.25)
    }
    if(stressor_type=="medium"){
        out_dist = pmax(rnorm(n, 0.25, 0.05), 0)
    }
    rxp = ec50_foo(out_dist,1,4)
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
                   "medium"=3
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

# Function to build HerbFact data frame for EFSA simulations
build_HerbFact <- function(HerbDuration, stressor, stressor_type){
    if(stressor_type=="low") level_vec = runif(HerbDuration, 0.01,0.25)
    if(stressor_type=="medium") level_vec = pmax(rnorm(HerbDuration, 0.25, 0.05), 0)
    if(stressor_type=="high") level_vec = runif(HerbDuration, 0.25,0.50)
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
    rxp = build_dose_response(n, stressor_type)
    
    SPECIES = gsub(" ", "_", df$Species)
    SPECIES = gsub("\\.", "", SPECIES)
    SPECIES = iconv(SPECIES, "latin1", "ASCII//TRANSLIT", sub="")

    df <- data.frame(
        ID=1:n,
        Species=SPECIES,
        MaxAge=100,
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
        sens=1, # MUST BE 1 for sensitivity !!! SPECIFICALLY FOR TXT !!!
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
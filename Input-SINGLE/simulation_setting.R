# simulation_setting.R

# Function to build application rate vector for IBCgrass
build_AppRate <- function(application_rate, herbicide_duration){
  data.frame(rep(application_rate, herbicide_duration))
}

# Function to build HerbFact.txt with a fixed effect value for a specific trait
build_HerbFact_fixed <- function(HerbDuration, effect_name, effect_value) {
  level_vec <- rep(effect_value, HerbDuration)
  no_vec <- rep(0, HerbDuration)
  
  data.frame(
    Biomass = if(effect_name == "biomass") level_vec else no_vec,
    Mortality = if(effect_name == "survival") level_vec else no_vec,
    SeedlingBiomass = if(effect_name == "SEbiomass") level_vec else no_vec,
    Establishment = if(effect_name == "establishment") level_vec else no_vec,
    SeedSterility = if(effect_name == "sterility") level_vec else no_vec,
    SeedNumber = if(effect_name == "seednumber") level_vec else no_vec
  )
}

# Function to build the correctly formatted Fieldedge.txt from the raw CSV
build_Fieldedge_Single <- function(df) {
  n <- nrow(df)
  rxp <- 0 # Set EC50 to 0 because we apply effects via HerbFact.txt (TXT way)
  
  # Clean species names
  SPECIES = gsub(" ", "_", df$Species)
  SPECIES = gsub("\\.", "", SPECIES)
  SPECIES = iconv(SPECIES, "latin1", "ASCII//TRANSLIT", sub="")
  
  data.frame(
    ID=1:n, Species=SPECIES, MaxAge=100, AllocSeed=0.05, LMR=df$LMR,
    m0=df$mSeed..mg., MaxMass=df$MaxMass..mg., mSeed=df$mSeed..mg.,
    Dist=df$Dist..m., pEstab=0.5, Gmax=df$Gmax, SLA=df$SLA.comb.,
    palat=0.5, memo=df$Memo..weeks., RAR=1, growth=0.25, mThres=0.2,
    clonal=df$Clonal, propSex=df$Clonal, meanSpacerLength=df$meanSpacerLength..cm.,
    sdSpacerLength=ifelse(df$sdSpacerLength..cm. == 17.5, 12.5,  df$sdSpacerLength..cm.),
    Resshare=df$Resshare, AllocSpacer=ifelse(df$Clonal==1, 0.05 , 0),
    mSpacer=ifelse(df$Clonal==1, 70, 0), sens=1, allocroot=1, allocshoot=1,
    EC50_biomass = rxp, slope_biomass = 4,
    EC50_SEbiomass = rxp, slope_SEbiomass = 4,
    EC50_survival = rxp, slope_survival = 4,
    EC50_establishment = rxp, slope_establishment = 4,
    EC50_sterility = rxp, slope_sterility = 4,
    EC50_seednumber = rxp, slope_seednumber = 4,
    FlowerWeek=16, DispWeek=20, GermPeriod=1
  )
}
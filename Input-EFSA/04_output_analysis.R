library(ggplot2)
library(dplyr)
library(stringr)
PATH_OUTPUT = "output-EFSA/"

sim_repo <- list.files(PATH_OUTPUT, full.names=FALSE)

Pt_files <- lapply(sim_repo, function(r){
    files = list.files(paste0(PATH_OUTPUT, r, "/Model-files/"), full.names=TRUE)
    files[grepl("Pt_", files)]
})
names(Pt_files) <- sim_repo

Grd_files <- lapply(sim_repo, function(r){
    files = list.files(paste0(PATH_OUTPUT, r, "/Model-files/"), full.names=TRUE)
    files[grepl("Grd_", files)]
})
names(Grd_files) <- sim_repo

build_DF <- function(ls_files){
    df = lapply(seq_along(ls_files), function(id){
        d = ls_files[[id]]
        combined_df <- lapply(seq_along(d), function(i) {
            df <- read.table(d[[i]], sep = "\t", header = TRUE)
            df$simulations = i
            return(df)
        }) %>%
            dplyr::bind_rows()
        combined_df$modality = names(Grd_files)[[id]]
        return(combined_df)
    }) %>%
        dplyr::bind_rows()
    splitter = str_split(df$modality, "_", simplify = TRUE)
    df$community = splitter[,1]
    df$stressor = splitter[,4]
    df$stressor_level = splitter[,3]
    df$stressor = factor(
        df$stressor,
        levels = c("NO", "biomass", "SEbiomass", "survival", "sterility", "seednumber","establishment"))
    return(df)
}


############# GRD
Grd_DF <- build_DF(Grd_files)

ggplot() +
    theme_minimal() +
    labs(x="Time (weeks)", y="Number of Individuals") + 
    scale_x_continuous(limits = c(448,500)) +
    scale_color_manual(values = c("#115522", "#553322", "#aa5522")) +
    geom_line(data = Grd_DF,
              aes(x = Time, y = NInd, group = simulations, color = stressor_level),
              alpha = 0.2  ) + 
    facet_grid(community ~ stressor)

ggplot() +
    theme_minimal() +
    labs(x="Time (weeks)", y="Total Biomass") + 
    scale_x_continuous(limits = c(0,500)) +
    scale_color_manual(values = c("#115522", "#553322", "#aa5522")) +
    geom_line(data = Grd_DF,
              aes(x = Time, y = totMass, group = simulations, color = stressor_level),
              alpha = 0.2  ) + 
    facet_grid(community ~ stressor)

ggplot() +
    theme_minimal() +
    labs(x="Time (weeks)", y="Shannon Diversity") + 
    scale_x_continuous(limits = c(0,500)) +
    scale_color_manual(values = c("#115522", "#553322", "#aa5522")) +
    geom_line(data = Grd_DF,
              aes(x = Time, y = shannon, group = simulations, color = stressor_level),
              alpha = 0.2  ) + 
    facet_grid(community ~ stressor)

############# PT
Pt_DF <- build_DF(Pt_files)
PFTs = unique(Pt_DF$PFT) 
Pt_DF_ = Pt_DF[Pt_DF$PFT == PFTs[[1]], ]

ggplot() +
    theme_minimal() +
    labs(x="Time (weeks)", y="seeds") + 
    scale_x_continuous(limits = c(448,500)) +
    scale_color_manual(values = c("#115522", "#553322", "#aa5522")) +
    geom_line(data = Pt_DF_,
              aes(x = Time, y = seedlings, group = simulations, color = stressor_level),
              alpha = 0.2  ) + 
    facet_grid(community ~ stressor)
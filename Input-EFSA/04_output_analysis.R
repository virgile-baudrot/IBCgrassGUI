build_DF <- function(ls_files){
  sim_repo <- basename(names(ls_files))
  df = lapply(seq_along(ls_files), function(id){
    d = ls_files[[id]]
    combined_df <- lapply(seq_along(d), function(i) {
      df <- read.table(d[[i]], sep = "\t", header = TRUE)
      df$simulations = i
      return(df)
    }) %>%
      dplyr::bind_rows()
    combined_df$modality = sim_repo[[id]]
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

# -------------------
# POPULATION CONTROL
# -------------------
p_popC_ind <- function(Grd_DF, time_burn=NA){
  if(!is.na(time_burn)){
    Grd_DF_burnin <- Grd_DF[Grd_DF$Time > time_burn*30,]
    p <- ggplot() +
      theme_minimal() +
      labs(x="Time (years)", y="Number of Individuals") + 
      scale_x_continuous(limits=c(time_burn, NA)) +
      geom_line(data = Grd_DF_burnin,
                aes(x = Time / 30, y = NInd, group = simulations),
                alpha = 0.2  ) + 
      facet_grid(~ community)
  } else{
    p <- ggplot() +
      theme_minimal() +
      labs(x="Time (years)", y="Number of Individuals") + 
      geom_line(data = Grd_DF,
                aes(x = Time / 30, y = NInd, group = simulations),
                alpha = 0.2  ) + 
      facet_grid(~ community)
  }
  return(p)
}

p_popC_mass <- function(Grd_DF, time_burn=NA){
  if(!is.na(time_burn)){
    Grd_DF_burnin <- Grd_DF[Grd_DF$Time > time_burn*30,]
    p <- ggplot() +
      theme_minimal() +
      labs(x="Time (years)", y="Total Mass (g)") + 
      scale_x_continuous(limits=c(time_burn, NA)) +
      geom_line(data = Grd_DF_burnin,
                aes(x = Time / 30, y = totMass, group = simulations),
                alpha = 0.2  ) + 
      facet_grid(~ community)
  } else{
    p <- ggplot() +
      theme_minimal() +
      labs(x="Time (years)", y="Total Mass (g)") + 
      geom_line(data = Grd_DF,
                aes(x = Time / 30, y = totMass, group = simulations),
                alpha = 0.2  ) + 
      facet_grid(~ community)
  }
  return(p)
}

p_popC_pft <- function(Grd_DF, time_burn=NA){
  if(!is.na(time_burn)){
    Grd_DF_burnin <- Grd_DF[Grd_DF$Time > time_burn*30,]
    p <- ggplot() +
      theme_minimal() +
      labs(x="Time (years)", y="Number of PFT") + 
      scale_x_continuous(limits=c(time_burn, NA)) +
      geom_line(data = Grd_DF_burnin,
                aes(x = Time / 30, y = NPFT, group = simulations),
                alpha = 0.2  ) + 
      facet_wrap(~ community, scale="free")
  } else{
    p <- ggplot() +
      theme_minimal() +
      labs(x="Time (years)", y="Number of PFT") + 
      geom_line(data = Grd_DF,
                aes(x = Time / 30, y = NPFT, group = simulations),
                alpha = 0.2  ) + 
      facet_wrap(~ community, scale="free")
  }
  return(p)
}

p_popC_div <- function(Grd_DF, time_burn=NA){
  if(!is.na(time_burn)){
    Grd_DF_burnin <- Grd_DF[Grd_DF$Time > time_burn*30,]
    p <- ggplot() +
      theme_minimal() +
      labs(x="Time (years)", y="Shannon Diversity") + 
      scale_x_continuous(limits=c(time_burn, NA)) +
      geom_line(data = Grd_DF_burnin,
                aes(x = Time / 30, y = shannon, group = simulations),
                alpha = 0.2  ) + 
      facet_grid(~ community)
  } else{
    p <- ggplot() +
      theme_minimal() +
      labs(x="Time (years)", y="Shannon Diversity") + 
      geom_line(data = Grd_DF,
                aes(x = Time / 30, y = shannon, group = simulations),
                alpha = 0.2  ) + 
      facet_grid(~ community)
  }
  return(p)
}

# -------------------
# INDIVIDUALS
# -------------------
p_ind_ind <- function(Pt_DF_reduce){
  mx_inds = max(Pt_DF_reduce$Inds)
  Pt_DF_reduce |>
    ggplot() +
    theme_minimal() +
    scale_fill_gradientn(
      colors = c("red", "blue", "cyan", "green"),
      values = scales::rescale(c(0, 0.001, 0.5, 1)),
      name = "# Inds"
    ) +
    labs(x="Time (years)", y="Individuals/PFT") + 
    geom_tile(aes(
      x = Time / 30, y = PFT, fill=Inds)) +
    facet_grid(~ simulations)
}


df_ind_sum <- function(Pt_DF_reduce){
  Pt_DF_sum = Pt_DF_reduce %>%
    dplyr::group_by(PFT, simulations) |>
    dplyr::summarise(
      sum_perPFT = sum(Inds),
      sum_shootmass = sum(shootmass),
      sum_seedling = sum(seedlings),
      sum_seed = sum(seeds)) |>
    dplyr::arrange(sum_perPFT) 
  return(Pt_DF_sum)
}

p_indC <- function(Pt_DF_sum, level){
  if(level=="individuals"){
    p <- ggplot(data = Pt_DF_sum) +
      theme_minimal() +
      labs(x="PFT (Species)", y="Sum of individuals / PFT") +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      ) +
      geom_bar(
        stat = "identity", fill = "steelblue",
        aes(x = reorder(PFT, -sum_perPFT), y = sum_perPFT)
      ) +
      facet_grid(simulations~.)
  }
  if(level=="shootmass"){
    p <- ggplot(data = Pt_DF_sum) +
      theme_minimal() +
      labs(x="PFT (Species)", y="Sum of shootmass / PFT") +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      ) +
      geom_bar(
        stat = "identity", fill = "#F54927",
        aes(x = reorder(PFT, -sum_shootmass), y = sum_shootmass)
      ) +
      facet_grid(simulations~.)
  }
  if(level=="seedling"){
    p <- ggplot(data = Pt_DF_sum) +
      theme_minimal() +
      labs(x="PFT (Species)", y="Sum of sum_seedling / PFT") +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      ) +
      geom_bar(
        stat = "identity", fill = "#5AA33B",
        aes(x = reorder(PFT, -sum_seedling), y = sum_seedling)
      ) +
      facet_grid(simulations~.)
  }
  if(level=="seeds"){
    p <- ggplot(data = Pt_DF_sum) +
      theme_minimal() +
      labs(x="PFT (Species)", y="Sum of seeds / PFT") +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      ) +
      geom_bar(
        stat = "identity", fill = "#953BA3",
        aes(x = reorder(PFT, -sum_seed), y = sum_seed)
      ) +
      facet_grid(simulations~.)
  }
  return(p)
} 

# -------------------
# POPULATION EFFECT
# -------------------
p_popE_ind <- function(Grd_DF, time_burn=NA){
  if(!is.na(time_burn)){
    Grd_DF_burnin <- Grd_DF[Grd_DF$Time > time_burn*30,]
    p <- ggplot() +
      theme_minimal() +
      labs(x="Time (years)", y="Number of Individuals") + 
      scale_x_continuous(limits=c(time_burn, NA)) +
      geom_line(data = Grd_DF_burnin,
                aes(x = Time / 30, y = NInd, group = simulations),
                alpha = 0.2  ) + 
      facet_grid(stressor ~ stressor_level)
  } else{
    p <- ggplot() +
      theme_minimal() +
      labs(x="Time (years)", y="Number of Individuals") + 
      geom_line(data = Grd_DF,
                aes(x = Time / 30, y = NInd, group = simulations),
                alpha = 0.2  ) + 
      facet_grid(stressor ~ stressor_level)
  }
  return(p)
}

p_popE_mass <- function(Grd_DF, time_burn=NA){
  if(!is.na(time_burn)){
    Grd_DF_burnin <- Grd_DF[Grd_DF$Time > time_burn*30,]
    p <- ggplot() +
      theme_minimal() +
      labs(x="Time (years)", y="Total Mass (g)") + 
      scale_x_continuous(limits=c(time_burn, NA)) +
      geom_line(data = Grd_DF_burnin,
                aes(x = Time / 30, y = totMass, group = simulations),
                alpha = 0.2  ) + 
      facet_grid(stressor ~ stressor_level)
  } else{
    p <- ggplot() +
      theme_minimal() +
      labs(x="Time (years)", y="Total Mass (g)") + 
      geom_line(data = Grd_DF,
                aes(x = Time / 30, y = totMass, group = simulations),
                alpha = 0.2  ) + 
      facet_grid(stressor ~ stressor_level)
  }
  return(p)
}

p_popE_pft <- function(Grd_DF, time_burn=NA){
  if(!is.na(time_burn)){
    Grd_DF_burnin <- Grd_DF[Grd_DF$Time > time_burn*30,]
    p <- ggplot() +
      theme_minimal() +
      labs(x="Time (years)", y="Number of PFT") + 
      scale_x_continuous(limits=c(time_burn, NA)) +
      geom_line(data = Grd_DF_burnin,
                aes(x = Time / 30, y = NPFT, group = simulations),
                alpha = 0.2  ) + 
      facet_wrap(stressor ~ stressor_level, scale="free")
  } else{
    p <- ggplot() +
      theme_minimal() +
      labs(x="Time (years)", y="Number of PFT") + 
      geom_line(data = Grd_DF,
                aes(x = Time / 30, y = NPFT, group = simulations),
                alpha = 0.2  ) + 
      facet_wrap(stressor ~ stressor_level, scale="free")
  }
  return(p)
}

p_popE_div <- function(Grd_DF, time_burn=NA){
  if(!is.na(time_burn)){
    Grd_DF_burnin <- Grd_DF[Grd_DF$Time > time_burn*30,]
    p <- ggplot() +
      theme_minimal() +
      labs(x="Time (years)", y="Shannon Diversity") + 
      scale_x_continuous(limits=c(time_burn, NA)) +
      geom_line(data = Grd_DF_burnin,
                aes(x = Time / 30, y = shannon, group = simulations),
                alpha = 0.2  ) + 
      facet_grid(stressor ~ stressor_level)
  } else{
    p <- ggplot() +
      theme_minimal() +
      labs(x="Time (years)", y="Shannon Diversity") + 
      geom_line(data = Grd_DF,
                aes(x = Time / 30, y = shannon, group = simulations),
                alpha = 0.2  ) + 
      facet_grid(stressor ~ stressor_level)
  }
  return(p)
}

df_effect_mean <- function(Grd_DF, time_burn=35){
  Grd_DF_burnin = Grd_DF[Grd_DF$Time > time_burn*30,]
  Grd_DF_burnin_sum = Grd_DF_burnin |>
    dplyr::group_by(Time, stressor_level, stressor) |>
    dplyr::summarise(
      mean_NInd = mean(NInd),
      mean_NPFT = mean(NPFT),
      mean_totMass = mean(totMass)
    )
  return(Grd_DF_burnin_sum)
}

# ----------------------
# INDIVIDUALS
# ----------------------
p_indE_ind <- function(Pt_DF, stress, stress_level){
  Pt_DF_reduce <- Pt_DF[Pt_DF$stressor==stress & Pt_DF$stressor_level==stress_level, ]
  
  mx_inds = max(Pt_DF_reduce$Inds)
  Pt_DF_reduce |>
    ggplot() +
    theme_minimal() +
    labs(title=paste("Stressor:", stress, ", level:", stress_level)) +
    scale_fill_gradientn(
      colors = c("red", "blue", "cyan", "green"),
      values = scales::rescale(c(0, 0.001, 0.5, 1)),
      name = "# Inds"
    ) +
    labs(x="Time (years)", y="Individuals/PFT") + 
    geom_tile(aes(
      x = Time / 30, y = PFT, fill=Inds)) +
    facet_grid(~ simulations)
  
}

p_ind_mean <- function(Grd_DF_burnin_sum, item){
  if(item=="Individuals"){
    p <- Grd_DF_burnin_sum |>
      ggplot() +
      theme_minimal() +
      labs(x="Time (years)", y="Mean Number of Individuals") + 
      scale_color_manual(
        values=c("#64BA2B", "#A7BA2B", "#BA912B", "#BA622B", "#BA392B")) +
      scale_x_continuous(limits=c(35, NA)) +
      geom_line(aes(x = Time / 30, y = mean_NInd, color=stressor_level),
                alpha=0.5) + 
      facet_grid(~ stressor)
  }
  if(item=="Mass") {
    p <- Grd_DF_burnin_sum |>
      ggplot() +
      theme_minimal() +
      labs(x="Time (years)", y="Mean Total Mass") + 
      scale_color_manual(
        values=c("#64BA2B", "#A7BA2B", "#BA912B", "#BA622B", "#BA392B")) +
      scale_x_continuous(limits=c(35, NA)) +
      geom_line(aes(x = Time / 30, y = mean_totMass, color=stressor_level),
                alpha=0.5) + 
      facet_grid(~ stressor)
  }
  if(item=="PFT"){
    p <- Grd_DF_burnin_sum |>
      ggplot() +
      theme_minimal() +
      labs(x="Time (years)", y="Mean Number of PFT") + 
      scale_color_manual(
        values=c("#64BA2B", "#A7BA2B", "#BA912B", "#BA622B", "#BA392B")) +
      scale_x_continuous(limits=c(35, NA)) +
      geom_line(aes(x = Time / 30, y = mean_NPFT, color=stressor_level),
                alpha=0.5) + 
      facet_grid(~ stressor)
  }
  return(p)
}

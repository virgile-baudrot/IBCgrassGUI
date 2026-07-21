build_DF <- function(ls_files){
  sim_repo <- basename(names(ls_files))
  df = lapply(seq_along(ls_files), function(id){
    d = ls_files[[id]]
    combined_df <- lapply(seq_along(d), function(i) {
      df <- tryCatch({
        read.table(d[[i]], sep = "\t", header = TRUE)
      }, error = function(e) {
        warning(paste("Erreur de lecture sur le fichier :", d[[i]]))
        return(NULL)
      })
      
      if (!is.null(df)) {
        df$simulations = i
        return(df)
      } else {
        return(NULL)
      }
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


# ----------------------
# POP
# ----------------------

pop_effect_prepare <- function(Grd_DF, QUANT=0.5, ratio=TRUE, rollmean=FALSE, rollmean_k=30){
  Grd_DF_z00 = Grd_DF |>
    dplyr::filter(stressor_level == "z00") |>
    group_by(Time, community, stressor) |>
    dplyr::mutate(
      NInd_control_median = median(NInd, na.rm=TRUE),
      totMass_control_median = median(totMass, na.rm=TRUE),
      abovemass_control_median = median(abovemass, na.rm=TRUE),
      NPFT_control_median = median(NPFT, na.rm=TRUE),
      shannon_control_median = median(shannon, na.rm=TRUE)
    ) |>
    dplyr::mutate(
      stressor_level=NULL,
      NInd=NULL,
      totMass=NULL,
      abovemass=NULL,
      NPFT=NULL,
      shannon=NULL
    )
  
  Grd_DF_merge = dplyr::left_join(Grd_DF, Grd_DF_z00, by=c("Time", "community", "stressor"))
  
  if(ratio){
    Grd_DF_diff = Grd_DF_merge |>
      dplyr::mutate(
        Nind_diff = NInd / NInd_control_median,
        totMass_diff= totMass / totMass_control_median,
        abovemass_diff= abovemass / abovemass_control_median,
        NPFT_diff = NPFT / NPFT_control_median,
        shannon_diff = shannon / shannon_control_median
      )
  } else{
    Grd_DF_diff = Grd_DF_merge |>
      dplyr::mutate(
        Nind_diff = NInd - NInd_control_median,
        totMass_diff= totMass - totMass_control_median,
        abovemass_diff= abovemass - abovemass_control_median,
        NPFT_diff = NPFT - NPFT_control_median,
        shannon_diff = shannon - shannon_control_median
      )
  }
  
  if(length(QUANT) == 1){
    Grd_DF_diff_sum <- Grd_DF_diff |>
      dplyr::group_by(Time, community, stressor, stressor_level) |>
      dplyr::summarise(
        q_NInd_diff = quantile(Nind_diff, QUANT, na.rm=TRUE),
        q_totMass_diff = quantile(totMass_diff, QUANT, na.rm=TRUE),
        q_abovemass_diff = quantile(abovemass_diff, QUANT, na.rm=TRUE),
        q_NPFT_diff = quantile(NPFT_diff, QUANT, na.rm=TRUE),
        q_shannon_diff = quantile(shannon_diff, QUANT, na.rm=TRUE)
      )
  } else{
    Grd_DF_diff_sum <- Grd_DF_diff |>
      dplyr::group_by(Time, community, stressor, stressor_level) |>
      dplyr::summarise(
        q_NInd_diff = quantile(Nind_diff, QUANT[1], na.rm=TRUE) - quantile(Nind_diff, QUANT[2], na.rm=TRUE),
        q_totMass_diff = quantile(totMass_diff, QUANT[1], na.rm=TRUE) - quantile(totMass_diff, QUANT[2], na.rm=TRUE),
        q_abovemass_diff = quantile(abovemass_diff, QUANT[1], na.rm=TRUE) - quantile(abovemass_diff, QUANT[2], na.rm=TRUE),
        q_NPFT_diff = quantile(NPFT_diff, QUANT[1], na.rm=TRUE) - quantile(NPFT_diff, QUANT[2], na.rm=TRUE),
        q_shannon_diff = quantile(shannon_diff, QUANT[1], na.rm=TRUE) - quantile(shannon_diff, QUANT[2], na.rm=TRUE)
      )
  }
  
  if(rollmean){
    Grd_DF_diff_sum <- Grd_DF_diff_sum |>
      dplyr::group_by(community, stressor, stressor_level) |>
      dplyr::mutate(
        q_NInd_diff = zoo::rollmean(q_NInd_diff, k = rollmean_k, fill = NA, align = "right"),
        q_totMass_diff = zoo::rollmean(q_totMass_diff, k = rollmean_k, fill = NA, align = "right"),
        q_abovemass_diff = zoo::rollmean(q_abovemass_diff, k = rollmean_k, fill = NA, align = "right"),
        q_NPFT_diff = zoo::rollmean(q_NPFT_diff, k = rollmean_k, fill = NA, align = "right"),
        q_shannon_diff = zoo::rollmean(q_shannon_diff, k = rollmean_k, fill = NA, align = "right")
      )
  }
  
  return(Grd_DF_diff_sum)
}

plot_eff_001 <- function(Grd_DF_diff_sum){
  ggplot(data=Grd_DF_diff_sum) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(x = "Time (years)", y = "Ratio Effect/Control Number of Individuals") +
    scale_color_manual(
      name="Stressor level:",
      values=c("black", "#44C714", "#C7C714", "#C75314", "#C71414", "#7314C7")
    )+
    geom_line(
      aes(x=Time/30, y = q_NInd_diff, color=stressor_level),
      alpha=0.7) +
    facet_grid(stressor~.)
}

plot_eff_002 <- function(Grd_DF_diff_sum){
  ggplot(data=Grd_DF_diff_sum) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(x = "Time (years)", y = "Ratio Effect/Control Total Mass (g)") +
    #scale_x_continuous(limits=c(5,NA)) +
    scale_color_manual(
      name="Stressor level:",
      values=c("black", "#44C714", "#C7C714", "#C75314", "#C71414", "#7314C7")
    )+
    geom_line(
      aes(x=Time/30, y = q_totMass_diff, color=stressor_level),
      alpha=0.7) +
    facet_grid(stressor~.)
}

plot_eff_003 <- function(Grd_DF_diff_sum){
  ggplot(data=Grd_DF_diff_sum) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(x = "Time (years)", y = "Ratio Effect/Control Above Mass (g)") +
    scale_color_manual(
      name="Stressor level:",
      values=c("black", "#44C714", "#C7C714", "#C75314", "#C71414", "#7314C7")
    )+
    geom_line(
      aes(x=Time/30, y = q_abovemass_diff, color=stressor_level),
      alpha=0.7) +
    facet_grid(stressor~.)
}

plot_eff_004 <- function(Grd_DF_diff_sum){
  ggplot(data=Grd_DF_diff_sum) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(x = "Time (years)", y = "Ratio Effect/Control Above Mass (g)") +
    scale_color_manual(
      name="Stressor level:",
      values=c("black", "#44C714", "#C7C714", "#C75314", "#C71414", "#7314C7")
    )+
    geom_line(
      aes(x=Time/30, y = q_abovemass_diff, color=stressor_level),
      alpha=0.7) +
    facet_grid(stressor~stressor_level)
}

plot_eff_005 <- function(Grd_DF_diff_sum){
  ggplot(data=Grd_DF_diff_sum) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(x = "Time (years)", y = "Ratio Effect/Control Number PTF") +
    scale_color_manual(
      name="Stressor level:",
      values=c("black", "#44C714", "#C7C714", "#C75314", "#C71414", "#7314C7")
    )+
    geom_line(
      aes(x=Time/30, y = q_NPFT_diff, color=stressor_level),
      alpha=0.7) +
    facet_grid(stressor~.)
}

plot_eff_006 <- function(Grd_DF_diff_sum){
  ggplot(data=Grd_DF_diff_sum) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(x = "Time (years)", y = "Ratio Effect/Control Shannon diveristy") +
    scale_color_manual(
      name="Stressor level:",
      values=c("black", "#44C714", "#C7C714", "#C75314", "#C71414", "#7314C7")
    )+
    geom_line(
      aes(x=Time/30, y = q_shannon_diff, color=stressor_level),
      alpha=0.7) +
    facet_grid(stressor~.)
}

plot_eff_006b <- function(Grd_DF_diff_sum){
  ggplot(data=Grd_DF_diff_sum) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(x = "Time (years)", y = "Ratio Effect/Control Shannon diveristy") +
    scale_color_manual(
      name="Stressor level:",
      values=c("black", "#44C714", "#C7C714", "#C75314", "#C71414", "#7314C7")
    )+
    geom_line(
      aes(x=Time/30, y = q_shannon_diff, color=stressor_level),
      alpha=0.7) +
    facet_grid(stressor~stressor_level)
}

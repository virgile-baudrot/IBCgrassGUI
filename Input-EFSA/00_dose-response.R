library(ggplot2)

dose_response <- function(rate, EC50, slope){
    response = rate^slope/(EC50^slope+rate^slope)
    return(response)
}

ec50_foo = function(dose_response, rate, slope){
    ec50 = (rate^slope/dose_response - rate^slope)^(1/slope)
    return(ec50)
}

mu_lognorm <- function(z, sdlog=1){
  log(z) - 1.645 * sdlog
}



# PLOTTING DOSE_RESPONSE
application_rate = seq(0, 2, by=0.01)
grids_par = expand.grid(
    EC50 = seq(0, 4.1, by=0.5),
    slope = 4
)

ls <- apply(grids_par, 1, function(p){
    data.frame(
        rate = application_rate,
        EC50 = p["EC50"],
        slope = p["slope"],
        response = dose_response(application_rate, p["EC50"], p["slope"])
    )
})
do.call("rbind", ls) -> dfDR


ggplot(dfDR, aes(x = rate, y = response, color = as.factor(EC50))) +
  geom_line() +
  theme_minimal() +
  labs(title = "Dose-Response Curve, slope=4", x = "Application Rate", y = "Response") +
  geom_vline(xintercept = 1, linetype=2, color="red") +
  scale_color_brewer(name="EC50", palette = "RdYlGn")

# PLOTTING RESPONSE WITH FIXED APPLICATION RATE
application_rate_fixed = 1
grids_par = expand.grid(
    EC50 = seq(0, 4, by=0.01),
    slope = 4
)
ls <- apply(grids_par, 1, function(p){
    data.frame(
        rate = application_rate_fixed,
        EC50 = p["EC50"],
        slope = p["slope"],
        response = dose_response(application_rate_fixed, p["EC50"], p["slope"])
    )
})
dfEC50 <- do.call("rbind", ls)

ggplot(dfEC50, aes(x = EC50, y = response, group=slope,)) +
    theme_minimal() +
    labs(title = "Response ~ EC50 - Fixed Application Rate=1 & Slope=4",
        x = "EC50",  y = "Response" ) +
    geom_line( color="#115566")


effect_dist = lapply(seq(0.1,0.5, by=0.1), function(p){
  data.frame(
    p=p,
    # rdist=rlnorm(1e5, meanlog = mu_lognorm(p), sdlog = 1)
    pdist=seq(0.01,0.99, by=0.01),
    qdist=qlnorm(seq(0.01,0.99, by=0.01), meanlog = mu_lognorm(p), sdlog = 1)
  )
})
df_effect_dist = do.call("rbind", effect_dist)

ggplot(data=df_effect_dist) +
  theme_minimal() +
  labs(title="Species Sensitivity Distribution",
       x="Species Distribution", y="Effect Level") +
  scale_color_manual(
    name="Quantile 95",
    values=c("#114411", "#115599", "#995599", "#995511", "#BB1111")
  ) +
  geom_line(aes(x=pdist, y=qdist, color=as.factor(p)), alpha=0.7) +
  geom_hline(
    yintercept = c(0.1, 0.2, 0.3, 0.4, 0.5),
    color=c("#114411", "#115599", "#995599", "#995511", "#BB1111"),
    linetype=4
  )



x = ec50_foo(out_dist,1,4)
y = dose_response(1, x, 4)
plot(hist(y))
dfpt = data.frame(
    ec50=x,
    dr=y
)

ggplot() +
    theme_minimal() +
    labs(title = "Response ~ EC50 - with Fixed Application Rate = 1",
         x = "EC50",
         y = "Response" ) +
    scale_color_manual(
        name="slope",
        values=c("#115566", "#665511", "#116655")) +
    geom_line(
        data=dfEC50,
        aes(x = EC50, y = response, group=slope, color=as.factor(slope))) + 
    geom_point(data=dfpt,
        aes(x=ec50, y=dr), alpha=0.1, size=4)
#     
ggplot() +
    theme_minimal() +
    labs(x = "EC50") +
    geom_histogram(data=dfpt, aes(ec50))

ggplot() +
    theme_minimal() +
    labs(x = "Response") +
    geom_histogram(data=dfpt, aes(dr))


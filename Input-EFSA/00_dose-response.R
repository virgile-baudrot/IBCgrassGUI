library(ggplot2)

dose_response <- function(rate, EC50, slope){
    response = rate^slope/(EC50^slope+rate^slope)
    return(response)
}

ec50_foo = function(dose_response, rate, slope){
    ec50 = (rate^slope/dose_response - rate^slope)^(1/slope)
    return(ec50)
}

# PLOTTING DOSE_RESPONSE
application_rate = seq(0, 2, by=0.01)
grids_par = expand.grid(
    EC50 = seq(0, 4.1, by=0.5),
    slope = c(1,4,7)
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
ggplot(dfDR, aes(x = rate, y = response, group = EC50, color = EC50)) +
    theme_minimal() +
    labs(
        title = "Dose-Response Curve",
        x = "Application Rate",
        y = "Response"
    ) +
    geom_line() +
    facet_grid(
        . ~ slope,
        labeller = labeller(
            slope = function(s) paste("slope:", s)
        )
    )

# PLOTTING RESPONSE WITH FIXED APPLICATION RATE
application_rate_fixed = 1
grids_par = expand.grid(
    EC50 = seq(0, 4, by=0.01),
    slope = c(1,4,7)
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


ggplot(dfEC50, aes(x = EC50, y = response, group=slope, color=as.factor(slope))) +
    theme_minimal() +
    labs(title = "Response ~ EC50 - with Fixed Application Rate = 1",
        x = "EC50",
        y = "Response" ) +
    scale_color_manual(
        name="slope",
        values=c("#115566", "#665511", "#116655")) +
    geom_line()


y = dose_response(1, runif(1e3,1.5,4), 4)
plot(hist(y))

out_dist = runif(1e4,0.25,0.50)
out_dist = runif(1e4,0.01,0.25)
out_dist = rnorm(1e4,0.25,0.05)

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


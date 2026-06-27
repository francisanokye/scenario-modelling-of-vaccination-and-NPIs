library(macpan2)
library(tidyverse);theme_set(theme_bw())
library(shellpipes)

loadEnvironments()

NL_POP <- 510550
upper_plot_time <- 300
burn_in_start  <- omicron_start - lubridate::days(150)
upper_date <- burn_in_start + upper_plot_time

housedat <- data.frame(date = as.Date("2022-01-11")
	, value = NL_POP * 0.209
	, matrix = "cum_booster"
)

# simulate trajectories from posterior
sims = (rdsRead("calibrate")
	|> mp_trajectory_sd(conf.int = TRUE, back_transform = TRUE)
#	|> filter(time >= offset0)
	|> mutate(date = burn_in_start + (time - 1))
)

serodat <- rdsRead("seroprevdata.rds") 
fitsero <- rdsRead("fitsero.rds")

print(head(serodat))

print(head(fitsero))

trajs <- (ggplot(sims, aes(date, value))
	+ geom_line(linewidth = 1,color="red")
	+ geom_point(data=serodat)
	+ geom_point(data=fitsero)
	+ geom_point(data = housedat,color="blue")
	+ geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4)
	+ facet_wrap(~matrix, scales = "free")
)

print(trajs)

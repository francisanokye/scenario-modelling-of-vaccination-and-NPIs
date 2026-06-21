library(macpan2)
library(tidyverse);theme_set(theme_bw())
library(shellpipes)

loadEnvironments()

upper_plot_time <- 300
upper_date <- day0 + upper_plot_time

housedat <- data.frame(date = as.Date("2022-01-11")
	, value = 100000
	, matrix = "cum_booster"
)

# simulate trajectories from posterior
sims = (rdsRead("calibrate")
	|> mp_trajectory_sd(conf.int = TRUE, back_transform = TRUE)
#	|> filter(time >= offset0)
	|> mutate(date = day0 + time)
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

library(macpan2)
library(tidyverse);theme_set(theme_bw())
library(shellpipes)

loadEnvironments()

upper_date <- burn_in_start + upper_plot_time

housedat <- data.frame(date = as.Date("2022-01-11")
	, value = NL_POP * 0.209
	, matrix = "cum_booster"
)

# simulate trajectories from posterior
sims = (rdsRead("calibrate")
	|> mp_trajectory_sd(conf.int = TRUE, back_transform = TRUE)
	|> mutate(date = burn_in_start + (time - 1))
)

serodat <- rdsRead("seroprevdata.rds") 
fitsero <- rdsRead("fitsero.rds")

print(head(serodat))

print(head(fitsero))

trajs <- (ggplot(sims, aes(date, value))
	+ geom_line(linewidth = 2,color="red")
	+ geom_point(data=serodat)
	+ geom_point(data=fitsero)
	+ geom_point(data = housedat,color="blue", size = 5)
	+ geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4)
	+ facet_wrap(~matrix, scales = "free")
	+ theme(
	  axis.text.x        = element_text(size = 25, hjust = 0.5),
	  axis.text.y        = element_text(size = 25),
	  axis.title.y       = element_text(size = 25),
	  legend.position    = "bottom",
	  legend.justification = "center",
	  legend.direction   = "vertical",
	  strip.text         = element_text(size = 25),
	  strip.background   = element_blank(),
	  legend.text        = element_text(size = 25),
	  legend.background  = element_blank(),
	  panel.grid.major.x = element_blank(),
	  plot.background    = element_blank(),
	  plot.margin        = margin(12, 35, 12, 12)
	)
)

png("../figures/calibration_plots.png", width = 5000, height = 3500, res = 300, bg = "white",type = "cairo")
print(trajs)
dev.off()



library(macpan2)
library(tidyverse);theme_set(theme_bw())
library(shellpipes)
library(ggthemes)

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
	  axis.text.x        = element_text(size = 25, hjust = 1, angle = 45),
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

# png("../figures/calibration_plots.png", width = 5000, height = 3500, res = 300, bg = "white",type = "cairo")
# print(trajs)
# dev.off()

simm <- (sims 
  |> filter(date >= omicron_start, matrix == "newR")
)

# Make sure intervention_removed is a Date
intervention_date <- as.Date(intervention_removed)

p <- (ggplot() 
  + geom_point(data = fitsero, aes(x = date, y = value, color = "serosurveillance"), size = 3.5, show.legend = TRUE) 
  + geom_ribbon(data = simm, aes(x = date, ymin = conf.low, ymax = conf.high),fill = "red",alpha = 0.40,show.legend = FALSE) 
  + geom_line(data = simm,aes(x = date,y = value,color = "model"),linewidth = 1.5,show.legend = TRUE) 
  + geom_vline(xintercept = intervention_date,colour = "black",linetype = "dashed",linewidth = 1.3,show.legend = FALSE) 
  + scale_colour_manual(
    name = NULL,
    values = c("serosurveillance" = "black","model" = "red"),
    breaks = c("serosurveillance","model"),
    labels = c("serosurveillance","model"),
    guide = guide_legend(
      override.aes = list(
        linetype = c(0, 1),
        shape = c(16, NA),
        linewidth = c(0, 1.2)
      )
    )
  ) 
  + labs(y = "Number of infections",title = "Infection incidence estimated from serosurveillance and model fit") 
  + scale_x_date(limits = c(omicron_start, analysis_end),date_breaks = "2 weeks", date_labels = "%d %b %Y",expand = c(0.01, 0)) 
  + theme_clean() 
  + theme(
    axis.text.x = element_text(size = 22,angle = 45,hjust = 1),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 30,color = "black"),
    plot.subtitle = element_blank(),
    legend.position = c(0.75, 0.45),
    legend.justification = c(0.5, 0.5),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 25),
    legend.background = element_blank(),
    plot.title = element_text(size = 30,hjust = 0.5,face = "plain"),
    plot.background = element_blank(),
    plot.tag = element_text(size = 35,face = "bold"),
    plot.tag.position = c(0.01, 0.98)
    )
)

png("../figures/SI_Fig1.png",width = 5000,height = 2500,res = 300,bg = "white",type = "cairo")
print(p)
dev.off()
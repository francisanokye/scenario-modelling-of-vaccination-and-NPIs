library(macpan2)
library(tidyverse);theme_set(theme_bw())
library(shellpipes)


offset0 <- 150
upper_plot_time <- 300

# simulate trajectories from posterior
sims = (rdsRead("calibrate")
        |> mp_trajectory_sd(conf.int = TRUE, back_transform = TRUE)
        |> dplyr::filter(time >= offset0)
)


# function to filter plotting time
lim_dat = \(x) dplyr::filter(x, time <= upper_plot_time)

# plot observed and fitted incidence
trajs <- (ggplot()
          + geom_point(aes(time, value), data = lim_dat(rdsRead("seroprevdata")))
          + geom_point(aes(time, value), data = lim_dat(rdsRead("fitsero")))
          + geom_line(aes(time, value), linewidth = 1, data = lim_dat(sims),color="red")
          + geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high), alpha = 0.4, data = lim_dat(sims))
          + facet_wrap(~matrix, scales = "free")
          + scale_y_continuous(name = "")
          + scale_x_continuous(limits = c(offset0, upper_plot_time))
          + theme_bw()
)

print(trajs)

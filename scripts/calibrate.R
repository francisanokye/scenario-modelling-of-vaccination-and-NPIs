library(macpan2)
library(shellpipes)
rpcall("calibrate.Rout calibrate.R timevar_spec.rds seroprevdata.rds fitsero.rds params.rda")
library(conflicted)
library(tidyverse)
library(dplyr)
library(ggthemes)
library(broom.mixed)

# check macpan2 version
if (packageVersion("macpan2") < "2.5.0") {
  stop(
    "please install a new version of macpan2\n"
    , "https://canmod.github.io/macpan2/#installation"
  )
}

# set seed and options
set.seed(2025)
options(macpan2_log_dir = ".")
loadEnvironments()

spline_beta = TRUE

# load input data
timevar_spec <- rdsRead("timevar_spec.rds")
seroprevdata <- rdsRead("seroprevdata.rds")
fitserodata <- rdsRead("fitsero.rds")
time_steps = max(seroprevdata$time)
upper_plot_time = 300

## convert cumulative to weekly incidence
#fitserodata = (seroprevdata 
#	|> mutate(value = timevar_spec$default$N * (value - dplyr::lag(value)) / 7)
#	|> mutate(matrix = "inc")
#)

# construct spline basis for beta(t)
if (spline_beta) {
  basis_cols = 11
  basis_rows = time_steps
  X = splines::ns(1:basis_rows
                  , basis_cols
                  , intercept = TRUE
                  , Boundary.knots = c(offset0, basis_rows)
  )
  timevar_spec = mp_tmb_insert_glm_timevar(timevar_spec
                                           , parameter_name = "beta"
                                           , design_matrix = X
                                           , timevar_coef = rep(0, basis_cols)
                                           , link_function = mp_log
  )
}

# prior constructor
get_prior = function(trans) function(rng) {
  mp_normal(
    (trans(rng[1]) + trans(rng[2])) / 2
    , log((trans(rng[2]) - trans(rng[1])) / (2 * 1.96))
  )
}

# define priors
priors = list(log_beta = mp_normal(log(0.25), log(1))
              , time_var_beta = mp_normal(0, log(1))
              , log_gamma_a = get_prior(log)(prior_range$gamma_a)
              , log_gamma_i = get_prior(log)(prior_range$gamma_i)
              , logit_kappa2 = get_prior(qlogis)(prior_range$kappa2)
              , logit_kappa3 = get_prior(qlogis)(prior_range$kappa3)
              , log_sigma = get_prior(log)(prior_range$sigma)
)

# fit model using seroprevalence data
calibrator = mp_tmb_calibrator(
  spec = timevar_spec |> mp_rk4()
  , data = (fitserodata##  seroprevdata                   ## change here?!? 
            |> select(-date) 
            |> dplyr::filter(matrix == "newR") 
            |> mutate(matrix = "log_newR", value = log(value))
  )
  , time = mp_sim_bounds(1, time_steps)
  , traj = list(log_newR = mp_normal(sd = mp_fit(0.1)))
  , par = priors
  , outputs = c("log_beta_thing", "log_inc", "log_newR","logit_serop_total")
)

# optimize fit
mp_optimize(calibrator)

# simulate trajectories from posterior
sims = (calibrator
        |> mp_trajectory_sd(conf.int = TRUE, back_transform = TRUE)
        |> dplyr::filter(time >= offset0)
)


# function to filter plotting time
lim_dat = \(x) dplyr::filter(x, time <= upper_plot_time)

# plot observed and fitted incidence
A = (ggplot()
      + geom_point(aes(time, value), data = lim_dat(seroprevdata))
      + geom_point(aes(time, value), data = lim_dat(fitserodata))
      + geom_line(aes(time, value), linewidth = 1, data = lim_dat(sims),color="red")
      + geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high), alpha = 0.4, data = lim_dat(sims))
      + facet_wrap(~matrix, scales = "free", ncol = 1)
      + scale_y_continuous(name = "")
      + scale_x_continuous(limits = c(offset0, upper_plot_time))
      + theme_bw()
)

png("../figures/A.png", width = 3000, height = 2500, res = 300, bg = "white", type = "cairo")
A
dev.off()



# plot full compartment trajectories
B = (calibrator
      |> mp_optimized_spec("modified")
      |> mp_simulator(time_steps, outputs = mp_state_vars(timevar_spec))
      |> mp_trajectory()
      |> lim_dat()
      |> ggplot()
      + aes(time, value)
      + geom_line()
      + facet_wrap(~matrix, scales = "free")
      + theme_bw()
)

png("../figures/B.png", width = 3000, height = 2500, res = 300, bg = "white", type = "cairo")
B
dev.off()

# prepare vaccination data for plotting
remade_vax_data = (list(
  double_vac = data.frame(
    time = timevar_spec$integers$double_vac_changepoints
    , value = timevar_spec$default$double_vac_values
  )
  , booster_shot = data.frame(
    time = timevar_spec$integers$booster_vac_changepoints
    , value = timevar_spec$default$booster_vac_values
  )
) 
|> bind_rows(.id = "matrix")
)

# plot vaccination coverage over time
C = (calibrator
      |> mp_optimized_spec("modified")
      |> mp_simulator(time_steps, outputs = c("double_vac","booster_shot"))
      |> mp_trajectory()
      |> lim_dat()
      |> ggplot()
      + aes(time, value)
      + geom_line()
      + geom_point(data = remade_vax_data)
      + facet_wrap(~matrix, scales = "free")
      + theme_bw()
)

png("../figures/C.png", width = 3000, height = 2500, res = 300, bg = "white", type = "cairo")
C
dev.off()

# extract fitted coefficients
model_estimates = mp_tmb_coef(calibrator, conf.int = TRUE) |> dplyr::select(-term, -col, -type) |> mutate(mat = abbreviate(mat, minlength = 15))
print(model_estimates, digits = 2)

# save fitted calibrator
rdsSave(calibrator)

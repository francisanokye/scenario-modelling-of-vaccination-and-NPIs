library(macpan2)
library(shellpipes)
library(conflicted)
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggthemes)
library(forcats)
library(broom.mixed)
# 
# check macpan2 version
if (packageVersion("macpan2") < "2.5.0") {
  stop(
    "please install a new version of macpan2\n"
    , "https://canmod.github.io/macpan2/#installation"
  )
}
# 
# # set seed and options
set.seed(2025)
options(scipen = 9999)
options(macpan2_log_dir = ".")
loadEnvironments()
conflicted::conflicts_prefer(dplyr::filter)

spline_beta = TRUE

# load input data
timevar_spec <- rdsRead("timevar_spec.rds")
seroprevdata <- rdsRead("seroprevdata.rds")
fitserodata <- rdsRead("fitsero.rds")
time_steps = max(seroprevdata$time)
upper_plot_time = 300
sim_horizon = 530
## Are we calling date0 from other places? We need to be consistent with the dates
date0 = as.Date("2021-12-15")

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

## Do we really need this? Yes

# prior constructor
get_prior = function(trans) function(rng) {
  mp_norm(
    (trans(rng[1]) + trans(rng[2])) / 2
    , log((trans(rng[2]) - trans(rng[1])) / (2 * 1.96))
  )
}

# define priors
priors = list(log_beta = mp_norm(log(0.25), log(1))
              , time_var_beta = mp_norm(0, log(1))
              , log_gamma_a = get_prior(log)(prior_range$gamma_a)
              , log_gamma_i = get_prior(log)(prior_range$gamma_i)
              , logit_kappa2 = get_prior(qlogis)(prior_range$kappa2)
              , logit_kappa3 = get_prior(qlogis)(prior_range$kappa3)
              , log_sigma = get_prior(log)(prior_range$sigma)
)

# fit model using seroprevalence data
calibrator = mp_tmb_calibrator(
  spec = timevar_spec |> mp_rk4()
  , data = (fitserodata
            |> select(-date) 
            |> dplyr::filter(matrix == "newR") 
            |> mutate(matrix = "log_newR", value = log(value))
  )
  , time = mp_sim_bounds(1, time_steps)
  , traj = list(log_newR = mp_norm(sd = mp_fit(0.1)))
  , par = priors
  , outputs = c("log_beta_thing", "log_inc", "log_newR","logit_serop_total", "cum_double_vac", "cum_booster")
)

# optimize fit
mp_optimize(calibrator)


# extract fitted coefficients
model_estimates = mp_tmb_coef(calibrator, conf.int = TRUE) |> 
  dplyr::select(-term, -col, -type) |> 
  mutate(mat = abbreviate(mat, minlength = 15))

print(model_estimates, digits = 2)

rdsSave(calibrator)


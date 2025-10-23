library(macpan2)
#library(conflicted)
library(tidyverse)
library(dplyr)
library(ggthemes)
library(broom.mixed)
library(zoo)
library(ggplot2)
library(shellpipes)
rpcall("calibrate.Rout calibrate.R timevar_spec.rds seroprevdata.rds fitsero.rds params.rda")

loadEnvironments()

# check macpan2 version
if (packageVersion("macpan2") < "2.5.0") {
  stop(
    "please install a new version of macpan2\n"
    , "https://canmod.github.io/macpan2/#installation"
  )
}

# ---------------------------------------
# Calibrate one scenario spec
# ---------------------------------------

## Why do you create functions if you only use them ones?
## Do not create functions unless you are using more than 3 times
calibrate_spec <- function(spec_file = "timevar_spec",
                           sero_file = "seroprevdata",
                           fitsero_file = "fitsero",
                           spline_beta = TRUE,
                           upper_plot_time = 300) {

  # ---------------------------------------
  # Load data and setup
  # ---------------------------------------
  timevar_spec <- rdsRead(spec_file)
  seroprevdata <- rdsRead(sero_file)
  fitserodata  <- rdsRead(fitsero_file)
  time_steps   <- max(seroprevdata$time)
  #tag          <- gsub("\\.rds$", "", basename(spec_file))
  set.seed(2025); options(macpan2_log_dir = ".")

  # ---------------------------------------
  # Construct spline basis for beta(t)
  # ---------------------------------------
  if (isTRUE(spline_beta)) {
    basis_cols <- 5L
    X <- splines::ns(1:time_steps, basis_cols, intercept = TRUE,
                     Boundary.knots = c(offset0, time_steps))
    timevar_spec <- mp_tmb_insert_glm_timevar(timevar_spec,
                                              parameter_name = "beta",
                                              design_matrix = X,
                                              timevar_coef = rep(0, basis_cols),
                                              link_function = mp_log)
  }

  # ---------------------------------------
  # Priors
  # ---------------------------------------

## Too much hard coding

  get_prior <- function(trans) function(rng) {
    mp_normal((trans(rng[1]) + trans(rng[2])) / 2,
              log((trans(rng[2]) - trans(rng[1])) / (2 * 1.96)))
  }

  priors <- list(
    log_beta      = mp_normal(log(0.25), log(1)),
    time_var_beta = mp_normal(0, log(1)),
    log_gamma_a   = get_prior(log)(prior_range$gamma_a),
    log_gamma_i   = get_prior(log)(prior_range$gamma_i),
    logit_kappa2  = get_prior(qlogis)(prior_range$kappa2),
    logit_kappa3  = get_prior(qlogis)(prior_range$kappa3),
    log_sigma     = get_prior(log)(prior_range$sigma)
  )

  # ---------------------------------------
  # Prepare seroprevalence data
  # ---------------------------------------
  sero_for_fit <- seroprevdata |>
    select(-date) |>
    dplyr::filter(matrix == "serop") |>
    mutate(matrix = "logit_serop", value = qlogis(value))

  fitserodata <- seroprevdata |>
    mutate(value = timevar_spec$default$N * (value - dplyr::lag(value)) / 7) |>
    mutate(matrix = "inc")

  # ---------------------------------------
  # Calibrate and simulate
  # ---------------------------------------
  calibrator <- mp_tmb_calibrator(
    spec   = timevar_spec |> mp_rk4(),
    data   = sero_for_fit,
    time   = mp_sim_bounds(1, time_steps),
    traj   = list(logit_serop = mp_normal(sd = mp_fit(1))),
    par    = priors,
    outputs = c("log_beta_thing", "log_inc", "logit_serop")
  )
  mp_optimize(calibrator)

  sims <- calibrator |>
    mp_trajectory_sd(conf.int = TRUE, back_transform = TRUE) |>
    dplyr::filter(time >= offset0)

  # ---------------------------------------
  # uncomment below to troubleshoot and show model fit 
  # ---------------------------------------
  lim_dat <- \(x) dplyr::filter(x, time <= upper_plot_time)

  # Observed vs fitted incidence
  # print(ggplot() +
  #        geom_point(aes(time, value), data = lim_dat(seroprevdata)) +
  #        geom_point(aes(time, value), data = lim_dat(fitserodata)) +
  #        geom_line(aes(time, value), linewidth = 1, data = lim_dat(sims)) +
  #        geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high),alpha = 0.4, data = lim_dat(sims)) +
  #        facet_wrap(~matrix, scales = "free", ncol = 1) +
  #        scale_y_continuous(name = "") +
  #        scale_x_continuous(limits = c(offset0, upper_plot_time)) +
  #        theme_bw())

  # uncomment codeblock below to show full compartment trajectories
  # print(calibrator |>
  #        mp_optimized_spec("modified") |>
  #        mp_simulator(time_steps, outputs = mp_state_vars(timevar_spec)) |>
  #        mp_trajectory() |>
  #        lim_dat() |>
  #        ggplot(aes(time, value)) +
  #        geom_line() +
  #        facet_wrap(~matrix, scales = "free") +
  #        theme_bw())


  return(calibrator)
}

# ---------------------------------------
# run calibration and save
# ---------------------------------------

calibrator  <- calibrate_spec()

rdsSave(calibrator)


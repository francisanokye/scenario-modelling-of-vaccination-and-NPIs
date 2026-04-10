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

fig_dir   <- "../figures"
# diag_dir  <- file.path(fig_dir, "diagnostics_actual")
# posterior_dir <- file.path(fig_dir, "posterior")
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  invisible(path)
}

ensure_dir(fig_dir)
# ensure_dir(diag_dir)
# ensure_dir(posterior_dir)

# load input data
timevar_spec <- rdsRead("timevar_spec.rds")
seroprevdata <- rdsRead("seroprevdata.rds")
fitserodata <- rdsRead("fitsero.rds")
time_steps = max(seroprevdata$time)
upper_plot_time = 300
sim_horizon = 500
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

# simulate trajectories from posterior
sims = (calibrator
        |> mp_trajectory_sd(conf.int = TRUE, back_transform = TRUE)
        |> dplyr::filter(time >= offset0)
)


# function to filter plotting time
lim_dat = \(x) dplyr::filter(x, time <= upper_plot_time)

# plot observed and fitted incidence
trajs <- (ggplot()
          + geom_point(aes(time, value), data = lim_dat(seroprevdata))
          + geom_point(aes(time, value), data = lim_dat(fitserodata))
          + geom_line(aes(time, value), linewidth = 1, data = lim_dat(sims),color="red")
          + geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high), alpha = 0.4, data = lim_dat(sims))
          + facet_wrap(~matrix, scales = "free", ncol = 1)
          + scale_y_continuous(name = "")
          + scale_x_continuous(limits = c(offset0, upper_plot_time))
          + theme_bw()
)

png("../figures/fitted_trajs.png", width = 4500, height = 2500, res = 300, bg = "white", type = "cairo")
print(trajs)
dev.off()

# plot full compartment trajectories
print(calibrator
      |> mp_optimized_spec("modified")
      |> mp_simulator(time_steps, outputs = c(mp_state_vars(timevar_spec)), "cum_booster")
      |> mp_trajectory()
      |> lim_dat()
      |> ggplot()
      + aes(time, value)
      + geom_line()
      + facet_wrap(~matrix, scales = "free")
      + theme_bw()
)


# extract fitted coefficients
model_estimates = mp_tmb_coef(calibrator, conf.int = TRUE) |> dplyr::select(-term, -col, -type) |> mutate(mat = abbreviate(mat, minlength = 15))
print(model_estimates, digits = 2)

############### Posteriror plot##########################################

keep_params <- c("gamma_a","gamma_i","kappa2","kappa3","sigma")

lab_map <- c(
  gamma_a = "gamma[a]",
  gamma_i = "gamma[i]",
  kappa3  = "kappa[3]",
  sigma   = "sigma",
  kappa2  = "kappa[2]"
)

df_keep <- model_estimates %>%
  dplyr::filter(mat %in% keep_params) %>%
  dplyr::mutate(
    param = forcats::fct_relevel(mat, keep_params)  
  )

posterior_plt <- ggplot(df_keep, aes(x = estimate, y = param)) +
  geom_vline(xintercept = 0, linetype = 3, linewidth = 0.6) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, linewidth = 0.7) +
  geom_point(size = 2.4) +
  scale_y_discrete(labels = function(x) parse(text = lab_map[x])) +
  labs(
    title = "Posterior summaries of fitted model parameters",
    x = "Estimate (95% CI)",
    y = NULL
  ) +
  theme_clean()+
  theme(
    axis.text.x       = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.text.y       = element_text(size = 12),
    axis.title.x      = element_text(size = 14),
    axis.title.y      = element_text(size = 14 ),
    plot.title        = element_text(size = 16, hjust = 0.5, face = "plain"),
    strip.text        = element_text(size = 16, face = "plain", hjust = 0.5),
    legend.position   = "top",
    legend.direction  = "horizontal",
    legend.box        = "horizontal",
    legend.spacing.x  = unit(1.2, "cm"),
    legend.spacing.y  = unit(0.5, "cm"),
    legend.background = element_blank(),
    legend.text       = element_text(size = 12, face = "plain"),
    panel.spacing     = unit(1.5, "cm"),
    panel.grid.major.y = element_line(),
    #panel.grid.major.y = element_blank(),
    
    plot.background   = element_blank()
  )

# fn <- file.path(posterior_dir, "posterior_plot.png")
# png(fn, width = 3500, height = 1500, res = 300, bg = "white", type = "cairo")
# posterior_plt
# dev.off()

# ----------------------- settings ------------------------------------------
n_draw <- 40000L

# plotmath labels for non-time-var params
lab_map <- c(
  beta    = "beta",
  gamma_a = "gamma[a]",
  gamma_i = "gamma[i]",
  sigma   = "sigma",
  kappa2  = "kappa[2]",
  kappa3  = "kappa[3]",
  dstr_prms_sd_nR = "sigma[nR]"
)

# classify transforms
log_params   <- c("beta","gamma_a","gamma_i","sigma","dstr_prms_sd_nR")
logit_params <- c("kappa2","kappa3")
id_params    <- c("time_var_beta")  # already on fitted scale (Normal prior)

# ----------------------- match your prior constructor -----------------------
prior_from_range <- function(rng, trans) {
  mu <- (trans(rng[1]) + trans(rng[2])) / 2
  log_sd <- log((trans(rng[2]) - trans(rng[1])) / (2 * 1.96))
  sd <- exp(log_sd)
  list(mu = mu, sd = sd)
}
prior_from_mp_norm <- function(mean, log_sd) list(mu = mean, sd = exp(log_sd))

# posterior approx on transformed scale (delta method)
post_on_log <- function(est, se) {
  est <- pmax(est, 1e-12)
  list(mu = log(est), sd = abs(se / est))
}
post_on_logit <- function(est, se) {
  p <- pmin(pmax(est, 1e-8), 1 - 1e-8)
  list(mu = qlogis(p), sd = abs(se / (p * (1 - p))))
}

# ----------------------- prior specs (edit as needed) -----------------------
prior_spec <- tibble::tibble(
  mat = c("beta", "time_var_beta", "dstr_prms_sd_nR"),
  prior_mean  = c(log(0.25), 0, log(0.1)),
  prior_logsd = c(log(1),    log(1), log(1))
)

# model_estimates must have: mat, row, estimate, std.error
me <- model_estimates %>%
  dplyr::select(mat, row, estimate, std.error) %>%
  left_join(prior_spec, by = "mat")

# ----------------------- draw generator ------------------------------------

############################################################

keep_params <- c("gamma_a","gamma_i","kappa2","kappa3","sigma")
set.seed(2025)
n_draw <- 40000L

# safe unicode labels (no parsing)
lab_unicode <- c(
  gamma_a = "γₐ",
  gamma_i = "γᵢ",
  kappa2  = "κ₂",
  kappa3  = "κ₃",
  sigma   = "σ"
)

# prior constructor (matches your model)
prior_from_range <- function(rng, trans) {
  mu <- (trans(rng[1]) + trans(rng[2])) / 2
  log_sd <- log((trans(rng[2]) - trans(rng[1])) / (2 * 1.96))
  sd <- exp(log_sd)
  list(mu = mu, sd = sd)
}

# posterior approx on transformed scale (delta method)
post_on_log <- function(est, se) {
  est <- pmax(est, 1e-12)
  list(mu = log(est), sd = abs(se / est))
}
post_on_logit <- function(est, se) {
  p <- pmin(pmax(est, 1e-8), 1 - 1e-8)
  list(mu = qlogis(p), sd = abs(se / (p * (1 - p))))
}

# keep only the 5 fitted params
me <- model_estimates %>%
  dplyr::filter(mat %in% keep_params) %>%
  dplyr::select(mat, estimate, std.error)

# draw generator (prior + posterior on natural scale)
make_draws_one <- function(mat, est, se) {
  
  if (mat %in% c("gamma_a","gamma_i","sigma")) {
    pr <- prior_from_range(prior_range[[mat]], log)
    po <- post_on_log(est, se)
    
    prior_val <- exp(rnorm(n_draw, pr$mu, pr$sd))
    post_val  <- exp(rnorm(n_draw, po$mu, po$sd))
    
  } else { # kappa2, kappa3
    pr <- prior_from_range(prior_range[[mat]], qlogis)
    po <- post_on_logit(est, se)
    
    prior_val <- plogis(rnorm(n_draw, pr$mu, pr$sd))
    post_val  <- plogis(rnorm(n_draw, po$mu, po$sd))
  }
  
  tibble::tibble(
    mat   = mat,
    value = c(prior_val, post_val),
    dist  = rep(c("Prior","Posterior"), each = n_draw)
  )
}

draws_df <- purrr::pmap_dfr(
  list(me$mat, me$estimate, me$std.error),
  make_draws_one
) %>%
  mutate(
    dist = factor(dist, levels = c("Prior","Posterior")),
    mat  = factor(mat, levels = keep_params),
    facet_lab = factor(lab_unicode[as.character(mat)],
                       levels = unname(lab_unicode[keep_params]))
  )

# plot
posterior <- ggplot(draws_df, aes(x = value, fill = dist, color = dist)) +
  geom_density(alpha = 0.35, linewidth = 0.8, adjust = 1.1) +
  facet_wrap(~ facet_lab, scales = "free", ncol = 5) +
  scale_fill_manual(values = c(Prior = "grey30", Posterior = "orange")) +
  scale_color_manual(values = c(Prior = "grey30", Posterior = "orange")) +
  labs(
    title = "Prior and posterior marginal distributions",
    x = NULL, y = NULL,
    fill = NULL, color = NULL
  ) +
  theme_clean()+
  theme(
    axis.text.x       = element_text(size = 20, angle = 0, hjust = 0.5),
    axis.text.y       = element_text(size = 20),
    axis.title.x      = element_text(size = 20),
    axis.title.y      = element_text(size = 20 ),
    plot.title        = element_text(size = 24, hjust = 0.5, face = "plain"),
    strip.text        = element_text(size = 20, face = "plain", hjust = 0.5),
    legend.position   = "top",
    legend.direction  = "horizontal",
    legend.box        = "horizontal",
    legend.spacing.x  = unit(1.2, "cm"),
    legend.spacing.y  = unit(0.5, "cm"),
    legend.background = element_blank(),
    legend.text       = element_text(size = 20, face = "plain"),
    panel.spacing     = unit(1.5, "cm"),
    panel.grid.major.y = element_line(),
    plot.background   = element_blank()
  )

# fn <- file.path(posterior_dir, "prior_posterior.png")
# png(fn, width = 7500, height = 1500, res = 300, bg = "white", type = "cairo")
# posterior
# dev.off()

actual_diagnostic = (calibrator
                     |> mp_optimized_spec("modified")
                     |> mp_simulator(sim_horizon, outputs = c(mp_state_vars(timevar_spec), "total_inc", "inc"))
                     |> mp_trajectory()
                     |> dplyr::mutate(date = date0 + (time - offset0))
                     |> ggplot() +
                       # shaded window: time 0 to offset0
                       geom_rect(
                         data = data.frame(
                           xmin = date0 - offset0,
                           xmax = date0
                         ),
                         aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                         inherit.aes = FALSE,
                         fill = "grey30", alpha = 0.25
                       ) +
                       geom_vline(
                         xintercept = date0 + (231 - offset0),
                         color = "black",
                         linetype = 4,
                         linewidth = 1.0
                       ) +
                       geom_vline(
                         xintercept = date0 + (300 - offset0),
                         colour = "purple",
                         linetype = 4,
                         linewidth = 1.5
                       ) +
                       aes(date, value) +
                       geom_line() +
                       labs(
                         title = "Diagnostic compartment trajectories (Actual)"
                       ) +
                       scale_x_date(expand = c(0, 0), date_breaks = "month", date_labels = "%b %d %Y") +
                       facet_wrap(~matrix, scales = "free") +
                       theme_clean() +
                       theme(
                         axis.text.x       = element_text(size = 12,  hjust = 1, angle = 45, margin = margin(t = 8)),
                         axis.text.y       = element_text(size = 18),
                         axis.title.x      = element_text(size = 18),
                         axis.title.y      = element_text(size = 18),
                         plot.title        = element_text(size = 24, hjust = 0.5, face = "plain"),
                         strip.text        = element_text(size = 22, face = "plain", hjust = 0.5),
                         strip.background  = element_blank(),
                         legend.position   = "top",
                         legend.direction  = "horizontal",
                         legend.box        = "horizontal",
                         legend.spacing.x  = unit(1.2, "cm"),
                         legend.spacing.y  = unit(0.5, "cm"),
                         legend.background = element_blank(),
                         legend.text       = element_text(size = 22, face = "plain"),
                         panel.spacing     = unit(1.5, "cm"),
                         plot.background   = element_blank()
                       )
)

# fn <- file.path(diag_dir,paste0("diagnostic_", "actual",".png"))
# png(fn, width = 7500, height = 5000, res = 300, bg = "white", type = "cairo")
# actual_diagnostic
# dev.off()


# save fitted calibrator
rdsSave(calibrator)



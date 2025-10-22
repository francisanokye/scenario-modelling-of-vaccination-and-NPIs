# counterfactuals via policy-neutral transmission scaling with robust injection/fallback
# all comments are lowercase and concise

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(ggthemes); library(splines)
  library(shellpipes); library(macpan2)
  library(ggh4x); library(scales); library(tibble);
  library(tidyverse); library(broom.mixed); library(zoo)
})

set.seed(2025)
try(loadEnvironments(), silent = TRUE)

date0    <- as.Date("2021-12-15")
date_end <- as.Date("2022-05-22")
offset0  <- 150
upper_plot_time <- 300
start_date <- date0
basis_cols <- 5
tau_days <- 2.0

# panels: vaccination scales
vax_scales <- c("25%" = 0.25, "50%" = 0.50, "75%" = 0.75)
vax_levels <- names(vax_scales)
vax_fac    <- factor(vax_levels, levels = vax_levels)
 
als_cols <- c("als2_only", "als4_only", "no_als_only")
als_fac  <- factor(als_cols, levels = als_cols)
 
#' # ------------------------------- inputs ------------------------------------
calibrator   <- readRDS("calibrator_vax_scale_100.rds")
seroprevdata <- readRDS("seroprevdata.rds")
time_steps   <- max(seroprevdata$time)

#' # freeze the calibrated spec
spec_fixed <- mp_optimized_spec(calibrator, "modified")
if (is.null(spec_fixed$default$double_vac_values_base))
  spec_fixed$default$double_vac_values_base  <- spec_fixed$default$double_vac_values
if (is.null(spec_fixed$default$booster_vac_values_base))
  spec_fixed$default$booster_vac_values_base <- spec_fixed$default$booster_vac_values
#' 
#' N_pop <- as.numeric(spec_fixed$default$N)  # population for scaling to counts
#' 
#' # ----------------------- extract β̂(t) and policies ------------------------
# beta_hat(t) from calibrator and extract beta
beta_df <- calibrator |>
  mp_trajectory_sd(conf.int = TRUE, back_transform = TRUE) |>
  dplyr::filter(matrix == "beta_thing", time >= offset0, time <= upper_plot_time) |>
  mutate(date = seq(start_date, by = "day", length.out = length(value)),
         beta_hat = value) |>
  select(date, beta_hat) |>
  filter(date >= date0, date <= date_end)
 
base_df  <- tibble(date = beta_df$date, t = as.numeric(date - start_date) + 1)
 
# ---- simulation helpers ---------------------------------------------------
set_vax <- function(cf, scale) {
  cf$default$double_vac_values  <- scale * cf$default$double_vac_values_base
  cf$default$booster_vac_values <- scale * cf$default$booster_vac_values_base
  cf
}

simulate_series <- function(spec, scale) {
  cf <- spec %>% set_vax(scale)
  cf %>%
    mp_simulator(time_steps, outputs = "inc") %>%
    mp_trajectory() %>%
    filter(matrix == "inc") %>%
    mutate(date = date0 + (time - offset0)) %>%
    filter(date >= date0, date <= date_end) %>%
    transmute(date, inc = value)
}

# ---- baseline (observed @ 100%) -------------------------------------------
baseline <- purrr::imap_dfr(vax_scales, function(scl, vlbl) {
  simulate_series(spec_fixed, scl) |>
    dplyr::mutate(
      vax = factor(vlbl, levels = levels(vax_fac)),
      series = "observed"
    ) |>
    mutate(inc = as.integer(inc)) |>
    tidyr::crossing(tibble(als_scen = als_fac)) |>   # replicate to all als columns
    dplyr::select(date, inc, vax, als_scen, series)
}) |>
  arrange(als_scen)

print(baseline)
####################################### Counterfactual Scenarios ####################################### 

in_int  <- function(d, a, b) d >= as.Date(a) & d <= as.Date(b)
 
policy_df <- base_df %>%
  mutate(
    ALS = case_when(
      in_int(date,"2021-12-15","2021-12-24") ~ "ALS2",
      in_int(date,"2021-12-25","2022-01-03") ~ "ALS3",
      in_int(date,"2022-01-04","2022-02-07") ~ "ALS4",
      in_int(date,"2022-02-08","2022-03-14") ~ "ALS3",
      in_int(date,"2022-03-15","2022-05-22") ~ "NoALS",
      TRUE ~ "ALS2"
    ),
    school = if_else(in_int(date,"2021-12-20","2022-01-25"), "closed", "open")
  )
 
# multipliers (rebased so no-als = 1)
R_ALS4 <- 1.52; R_ALS3 <- 2.12; R_NoALS <- 2.32
R_open <- 2.03; R_closed <- 1.44
R_ALS2 <- sqrt(R_ALS3 * R_NoALS)

m_als    <- c(ALS4 = R_ALS4/R_NoALS, ALS3 = R_ALS3/R_NoALS, ALS2 = R_ALS2/R_NoALS, NoALS = 1.00)
m_school <- c(open = 1.00, closed = R_closed/R_open)

# smooth categorical timelines to c¹ series
smooth_step <- function(t, t0, tau = 2.0) plogis((as.numeric(t - t0))/tau)

make_smooth_series <- function(date, lev, map, tau = tau_days) {
  lev <- as.character(lev)
  seg_idx    <- c(1L, which(lev[-1L] != lev[-length(lev)]) + 1L)
  seg_dates  <- date[seg_idx]
  seg_levels <- lev[seg_idx]
  m <- rep(unname(map[seg_levels[1L]]), length(date))
  if (length(seg_idx) >= 2L) {
    for (j in 2:length(seg_idx)) {
      t0 <- seg_dates[j]
      m_target <- unname(map[seg_levels[j]])
      s <- smooth_step(date, t0, tau)
      m <- m*(1 - s) + m_target*s
    }
  }
  m
}

# recover policy-neutral beta0_hat
df_est <- policy_df %>%
  left_join(beta_df, by = "date") %>%
  mutate(
    m_ALS_hat    = make_smooth_series(date, ALS, map = m_als, tau = tau_days),
    m_school_hat = make_smooth_series(date, school, map = m_school, tau = tau_days),
    m_obs_hat    = m_ALS_hat * m_school_hat,
    beta0_hat    = beta_hat / pmax(m_obs_hat, 1e-8),
    beta_recon   = beta0_hat * m_obs_hat
  )

stopifnot(!any(is.na(df_est$beta0_hat)))

# ------------------- target curves for 3 counterfactuals -------------------
# counterfactual betas (schools as observed; ALS constant throughout) ----
k_ALS <- c(ALS2 = unname(m_als["ALS2"]),
           ALS4 = unname(m_als["ALS4"]),
           NoALS = unname(m_als["NoALS"]))  

df_est_cf <- df_est |>
    mutate(
    date = date,
    beta_cf_ALS2  = beta0_hat * (k_ALS["ALS2"]  * m_school_hat),
    beta_cf_ALS4  = beta0_hat * (k_ALS["ALS4"]  * m_school_hat),
    beta_cf_NoALS = beta0_hat * (k_ALS["NoALS"] * m_school_hat)
  )

############# -------------------- fit time_var_beta for each scenario ------------------
# use the same basis size p as in the calibrated spec
p <- length(spec_fixed$default$time_var_beta)
nT <- nrow(df_est)
tgrid <- if (nT > 1) seq.int(from = offset0, to = time_steps, length.out = nT) else offset0

# note: splines::ns() has no 'degree' argument
X <- splines::ns(tgrid, df = p, intercept = TRUE, Boundary.knots = c(offset0, time_steps))

log_beta0 <- as.numeric(spec_fixed$default$log_beta)

fit_timevar_coef <- function(beta_cf, log_beta0, X, lambda_frac = 1e-6) {
  y   <- log(pmax(beta_cf, 1e-12)) - log_beta0
  XtX <- crossprod(X); Xty <- crossprod(X, y)
  lam <- lambda_frac * mean(diag(XtX))
  as.numeric(solve(XtX + lam * diag(ncol(X)), Xty))
}

thetaALS2  <- fit_timevar_coef(df_est_cf$beta_cf_ALS2,  log_beta0, X)
thetaALS4  <- fit_timevar_coef(df_est_cf$beta_cf_ALS4,  log_beta0, X)
thetaNoALS <- fit_timevar_coef(df_est_cf$beta_cf_NoALS, log_beta0, X)

# print(thetaALS2)
# print(thetaALS4)
# print(thetaNoALS)
# print(spec_fixed$default$time_var_beta)


## ---- helper: ensure vax vectors numeric and set a scale -------------------
# --- helpers ---------------------------------------------------------------

# sanitize_beta_path <- function(spec) {
#   # keep only time_var_beta as the beta driver; zero out other beta-related pieces if present
#   if (!is.null(spec$default$linear_pred_beta)) spec$default$linear_pred_beta <- rep(0, length(spec$default$linear_pred_beta))
#   if (!is.null(spec$default$beta_t))           spec$default$beta_t           <- NULL
#   spec
# }

set_vax_numeric <- function(spec) {
  for (nm in c("double_vac_values_base","booster_vac_values_base",
               "double_vac_values","booster_vac_values")) {
    if (!is.null(spec$default[[nm]])) spec$default[[nm]] <- as.numeric(spec$default[[nm]])
  }
  spec
}

set_vax <- function(spec, scale) {
  spec <- set_vax_numeric(spec)
  spec$default$double_vac_values  <- scale * spec$default$double_vac_values_base
  spec$default$booster_vac_values <- scale * spec$default$booster_vac_values_base
  spec
}

simulate_series <- function(spec, scale, time_steps, date0, date_end, offset0) {
  spec <- set_vax(spec, scale)
  spec %>%
    mp_simulator(time_steps, outputs = "inc") %>%
    mp_trajectory() %>%
    dplyr::filter(matrix == "inc") #%>%
    #dplyr::mutate(date = date0 + (time - offset0)) %>%
    #dplyr::filter(date >= date0, date <= date_end) %>%
    #dplyr::transmute(date, inc = value)
}

# --- core loop (rows = vax; cols = ALS scenarios) --------------------------

simulate_counterfactual_grid <- function(spec_fixed, thetas, vax_scales,
                                         time_steps, date0, date_end, offset0) {
  
  dplyr::bind_rows(lapply(names(thetas), function(scen_name) {
    theta <- thetas[[scen_name]]
    
    # clone, inject β(t), neutralize policy multipliers (no double counting), sanitize beta path
    spec_cf <- spec_fixed
    spec_cf$default$time_var_beta <- as.numeric(theta)
    if (!is.null(spec_cf$default$m_ALS))    spec_cf$default$m_ALS    <- rep(1, time_steps)
    if (!is.null(spec_cf$default$m_school)) spec_cf$default$m_school <- rep(1, time_steps)
    #spec_cf <- sanitize_beta_path(spec_cf)
    
    dplyr::bind_rows(lapply(names(vax_scales), function(vlbl) {
      scl <- vax_scales[[vlbl]]
      simulate_series(spec_cf, scl, time_steps, date0, date_end, offset0) #%>%
        # dplyr::mutate(
        #   # match your baseline presentation: integer incidence
        #   inc      = as.integer(round(pmax(inc, 0))),
        #   vax      = factor(vlbl, levels = names(vax_scales)),
        #   als_scen = factor(scen_name, levels = c("als2_only", "als4_only", "no_als_only")),
        #   series   = "counterfactual"
        # )
    }))
  }))
}

# --- usage -----------------------------------------------------------------

thetas <- list(
  "als2_only" = thetaALS2,
  "als4_only" = thetaALS4,
  "no_als_only" = thetaNoALS
)
vax_scales <- c("75%"=.75,"50%"=.50,"25%"=.25)

grid_sims <- simulate_counterfactual_grid(
  spec_fixed = spec_fixed,
  thetas     = thetas,
  vax_scales = vax_scales,
  time_steps = time_steps,
  date0 = date0, date_end = date_end, offset0 = offset0
)

## ---- optional: round to whole numbers for presentation --------------------
grid_sims <- grid_sims #|> dplyr::mutate(inc = round(inc))

print(grid_sims)

quit()


plot_df <- dplyr::bind_rows(
  grid_sims, baseline
)

## ---- plot: 3×3 panels with axes on every panel ----------------------------
p_3x3 <- ggplot(plot_df, aes(x = date, y = inc, group = series, color = series)) +
  geom_line(linewidth = 1.1) +
  ggh4x::facet_grid2(rows = vars(vax), cols = vars(als_scen), axes = "all") +
  scale_x_date(limits = c(date0, date_end), expand = c(0,0),
               date_breaks = "2 weeks", date_labels = "%b %d") +
  labs(title = "counterfactual daily infections by vaccination level and ALS policy",
       x = NULL, y = "estimated daily infections") +
  theme_clean()+
  theme(
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(size = 22, hjust = 0.5, face = "plain"),
    strip.text   = element_text(size = 20, face = "plain"),
    legend.position = "bottom",
    panel.spacing   = unit(0.7, "lines"),
    plot.background = element_blank()
  ) 

print(p_3x3)

png("../figures/observed_vs_counterfactuals.png",width = 5200, height = 4200, res = 300, bg = "white", type = "cairo")
print(p_3x3)
dev.off()

# quick check: reconstructed beta(t) from spec_cf
# recon_beta <- function(spec, n) {
#   # internal predictor: log_beta + X %*% time_var_beta
#   log_beta0 <- as.numeric(spec$default$log_beta)
#   theta     <- as.numeric(spec$default$time_var_beta)
#   p         <- length(theta)
#   tgrid     <- offset0 + seq_len(n) - 1L
#   X         <- splines::ns(tgrid, df = p, intercept = TRUE, Boundary.knots = c(offset0, time_steps))
#   as.numeric(exp(log_beta0 + as.vector(X %*% theta)))
# }
# plot(recon_beta(spec_cf, as.integer(difftime(date_end, date0, units="days")) + 1), type="l")





quit()

########################
set_vax <- function(spec, scale) {
  # coerce possible non-numeric storage to numeric
  for (nm in c("double_vac_values_base","booster_vac_values_base",
               "double_vac_values","booster_vac_values")) {
    if (!is.null(spec$default[[nm]])) {
      spec$default[[nm]] <- as.numeric(spec$default[[nm]])
    }
  }
  spec$default$double_vac_values  <- scale * spec$default$double_vac_values_base
  spec$default$booster_vac_values <- scale * spec$default$booster_vac_values_base
  spec
}

simulate_series <- function(spec, scale) {
  spec %>%
    set_vax(scale) %>%
    mp_simulator(time_steps, outputs = "inc") %>%
    mp_trajectory() %>%
    dplyr::filter(matrix == "inc") %>%
    dplyr::mutate(date = date0 + (time - offset0)) %>%
    dplyr::filter(date >= date0, date <= date_end) %>%
    dplyr::transmute(date, inc = value)
}

# build a spec that uses those coefficients and neutralizes policy multipliers
coef_timevar_beta <- c(thetaALS2,thetaALS4,thetaNoALS)
spec_cf <- spec_fixed
spec_cf$default$time_var_beta <- as.numeric(thetaALS2)
if (!is.null(spec_cf$default$m_ALS)) spec_cf$default$m_ALS <- rep(1, time_steps)
if (!is.null(spec_cf$default$m_school)) spec_cf$default$m_school <- rep(1, time_steps)

counterfactual <- purrr::imap_dfr(vax_scales, function(scl, vlbl) {
  simulate_series(spec_cf, scl) |>
    dplyr::mutate(
      vax = factor(vlbl, levels = levels(vax_fac)),
      series = "counterfactual"
    ) |>
    mutate(inc = as.integer(inc)) |>
    tidyr::crossing(tibble(als_scen = als_fac)) |>   # replicate to all als columns
    dplyr::select(date, inc, vax, als_scen, series)
}) |>
  arrange(als_scen)


sim_counterfactual <- purrr::imap_dfr(vax_scales, function(scl, vlbl) {
  simulate_series(spec_cf, scl) %>%
    dplyr::mutate(
      vax    = factor(vlbl, levels = levels(vax_fac)),
      series = "counterfactual"
    ) %>%
    tidyr::crossing(tibble(als_scen = als_fac)) %>%   # replicate across als columns if needed
    dplyr::select(date, inc, vax, als_scen, series)
}) %>%
  dplyr::arrange(als_scen)

print(sim_counterfactual)

quit()







# --- helpers: make vax numeric-safe and simulate at a given vax scale -------
set_vax <- function(spec, scale) {
  # coerce possible non-numeric storage to numeric
  for (nm in c("double_vac_values_base","booster_vac_values_base",
               "double_vac_values","booster_vac_values")) {
    if (!is.null(spec$default[[nm]])) {
      spec$default[[nm]] <- as.numeric(spec$default[[nm]])
    }
  }
  spec$default$double_vac_values  <- scale * spec$default$double_vac_values_base
  spec$default$booster_vac_values <- scale * spec$default$booster_vac_values_base
  spec
}

simulate_series <- function(spec, scale) {
  spec %>%
    set_vax(scale) %>%
    mp_simulator(time_steps, outputs = "inc") %>%
    mp_trajectory() %>%
    dplyr::filter(matrix == "inc") %>%
    dplyr::mutate(date = date0 + (time - offset0)) %>%
    dplyr::filter(date >= date0, date <= date_end) %>%
    dplyr::transmute(date, inc = value)
}

# --- choose the counterfactual coefficient vector you want to run -----------
# e.g., use one of the precomputed vectors: theta_ALS2, theta_ALS4, theta_NoALS
coef_vec_cf <- thetaALS4  # <-- set this to the scenario you want

# build a spec that uses those coefficients and neutralizes policy multipliers
spec_cf <- spec_fixed
spec_cf$default$time_var_beta <- as.numeric(coef_vec_cf)
if (!is.null(spec_cf$default$m_ALS))    spec_cf$default$m_ALS    <- rep(1, time_steps)
if (!is.null(spec_cf$default$m_school)) spec_cf$default$m_school <- rep(1, time_steps)

# --- iterate over the vaccination scales correctly with imap_dfr -----------
# imap_dfr passes (value, name) → here: (scl numeric, vlbl "25%"/"50%"/"75%")
sim_counterfactual <- purrr::imap_dfr(vax_scales, function(scl, vlbl) {
  simulate_series(spec_cf, scl) %>%
    dplyr::mutate(
      vax    = factor(vlbl, levels = levels(vax_fac)),
      series = "counterfactual"
    ) %>%
    tidyr::crossing(tibble(als_scen = als_fac)) %>%   # replicate across als columns if needed
    dplyr::select(date, inc, vax, als_scen, series)
}) %>%
  dplyr::arrange(als_scen)

print(sim_counterfactual)


quit()
#' construct_spline_basis <- function(beta_cf) {
#'   #' construct spline (or identity) coefficients for beta(t)
#'   #'
#'   #' @description
#'   #' given a target daily transmission curve \eqn{\beta_{\mathrm{cf}}(t)}, return the
#'   #' coefficient vector (length = p) expected by the model for `matrix_coef_beta`.
#'   #' if \eqn{p = n_T} (days), the identity basis is used; otherwise a b-spline basis
#'   #' on a scaled time axis \eqn{x \in [0,1]} is fit under a log link.
#'   #'
#'   #' @param beta_cf numeric vector of target \eqn{\beta_{\mathrm{cf}}(t)} over the
#'   #' full simulation window (length \eqn{n_T}).
#'   #'
#'   #' @returns numeric vector of length \eqn{p} containing coefficients on the log scale,
#'   #' ready to assign to `matrix_coef_beta`.
#'   #'
#'   #' @details
#'   #' - computes \eqn{y = \log(\max(\beta_{\mathrm{cf}}, 1e-10))} to enforce positivity.
#'   #' - uses `coef_len` (must be defined in scope) as the intended coefficient count \eqn{p}.
#'   #' - if \eqn{p = n_T}, returns identity coefficients: \eqn{y}.
#'   #' - else builds a degree-3 b-spline design with `df = p`, `intercept = TRUE`,
#'   #'   `Boundary.knots = c(0, 1)` on \eqn{x \in [0,1]} and solves \eqn{X c \approx y}
#'   #'   using ridge, then QR, then SVD pseudoinverse as fallbacks.
#'   #'
#'   #' @section prerequisites:
#'   #' define `coef_len` before calling (e.g., `coef_len <- length(spec$all_matrices()$matrix_coef_beta)`).
#'   #'
#'   #' @importFrom splines bs
#'   #' @seealso splines::bs
#'   #'
#'   #' @examples
#'   #' \dontrun{
#'   #' # coef_len <- length(spec$all_matrices()$matrix_coef_beta)
#'   #' coefs <- construct_spline_basis(beta_cf = df_est_cf$beta_hat_ALS2)
#'   #' # then write 'coefs' back into the spec's matrix_coef_beta
#'   #' }
#'   
#'   y <- log(pmax(beta_cf, 1e-10))  # model uses log link for beta
#'   nT <- length(beta_cf)
#'   p  <- length(spec_fixed$all_matrices()$matrix_coef_beta)
#'   
#'   if (p == nT) return(as.numeric(y))  # identity basis
#'   
#'   x <- (seq_len(nT) - 1) / (nT - 1)  # scaled 0..1
#'   X <- splines::bs(x, df = p, degree = 3, intercept = TRUE, Boundary.knots = c(0, 1))
#'   
#'   XtX <- crossprod(X); Xty <- crossprod(X, y)
#'   lam <- 1e-6 * mean(diag(XtX))
#'   coef_hat <- try(solve(XtX + lam * diag(p), Xty), silent = TRUE)
#'   if (inherits(coef_hat, "try-error")) {
#'     coef_hat <- try(qr.solve(X, y), silent = TRUE)
#'     if (inherits(coef_hat, "try-error")) {
#'       sv  <- svd(X)
#'       tol <- max(dim(X)) * .Machine$double.eps * max(sv$d)
#'       dI  <- ifelse(sv$d > tol, 1 / sv$d, 0)
#'       coef_hat <- sv$v %*% (dI * (t(sv$u) %*% y))
#'     }
#'   }
#'   as.numeric(coef_hat)
#' }
#' 
#' # build target beta(t) for each constant-als scenario (already computed as df_est_cf$beta_hat_*):
#' matrix_coef_ALS2  <- construct_spline_basis(df_est_cf$beta_cf_ALS2)
#' matrix_coef_ALS4  <- construct_spline_basis(df_est_cf$beta_cf_ALS4)
#' matrix_coef_NoALS <- construct_spline_basis(df_est_cf$beta_cf_NoALS)

################### Update Spec ################### 

# spec_ALS2 <- mp_tmb_update(
#   spec_fixed,
#   phase = "during",
#   at = 1L,
#   default = list(matrix_coef_beta ~ matrix_coef_ALS2)
# )

# ALS-2 throughout @ 25%, 50% and 75%
# als2_025 <- simulate_series(spec_ALS2, 0.25) %>% mutate(series = "als2_025")
# als2_050 <- simulate_series(spec_ALS2, 0.50) %>% mutate(series = "als2_050")
# als2_075 <- simulate_series(spec_ALS2, 0.75) %>% mutate(series = "als2_075")

# ALS-4 throughout @ 25%, 50% and 75%
# als4_025 <- simulate_series(spec_ALS4, 0.25) %>% mutate(series = "als4_025")
# als4_050 <- simulate_series(spec_ALS4, 0.50) %>% mutate(series = "als4_050")
# als4_075 <- simulate_series(spec_ALS4, 0.75) %>% mutate(series = "als4_075")
# 
# # ALS-2 throughout @ 25%, 50% and 75%
# alsnoals_025 <- simulate_series(spec_NoALS, 0.25) %>% mutate(series = "alsnoals_025")
# alsnoals_050 <- simulate_series(spec_NoALS, 0.50) %>% mutate(series = "alsnoals_050")
# alsnoals_075 <- simulate_series(spec_NoALS, 0.75) %>% mutate(series = "alsnoals_075")

# print(head(baseline))
# print(head(als2_025))
# print(head(als4_025))
# print(head(alsnoals_025))

quit()

# counterfactuals via policy-neutral beta with calibrated basis (clean, end-to-end)
# all comments are in simple lowercase, short, and explicit.

# counterfactuals via policy-neutral beta with calibrated basis (clean, end-to-end)
# all comments are in simple lowercase, short, and explicit.

# counterfactuals with per-scenario time_var_beta fitting + counts output
# comments are simple and in lowercase

#' suppressPackageStartupMessages({
#'   library(dplyr); library(tidyr); library(purrr)
#'   library(ggplot2); library(ggthemes); library(splines)
#'   library(macpan2)
#' })
#' 
#' set.seed(2025)
#' 
#' # ------------------------------- config ------------------------------------
#' date0    <- as.Date("2021-12-15")
#' date_end <- as.Date("2022-05-22")
#' offset0  <- 150
#' upper_plot_time <- 300
#' start_date <- date0
#' 
#' # panels: vaccination scales
#' vax_scales <- c("25%" = 0.25, "50%" = 0.50, "75%" = 0.75)
#' vax_levels <- names(vax_scales)
#' vax_fac    <- factor(vax_levels, levels = vax_levels)
#' 
#' als_cols <- c("als2_only", "als4_only", "noals_only")
#' als_fac  <- factor(als_cols, levels = als_cols)
#' 
#' # ------------------------------- inputs ------------------------------------
#' calibrator   <- readRDS("calibrator_vax_scale_100.rds")
#' seroprevdata <- readRDS("seroprevdata.rds")
#' time_steps   <- max(seroprevdata$time)
#' 
#' # freeze the calibrated spec
#' spec_fixed <- mp_optimized_spec(calibrator, "modified")
#' if (is.null(spec_fixed$default$double_vac_values_base))
#'   spec_fixed$default$double_vac_values_base  <- spec_fixed$default$double_vac_values
#' if (is.null(spec_fixed$default$booster_vac_values_base))
#'   spec_fixed$default$booster_vac_values_base <- spec_fixed$default$booster_vac_values
#' 
#' N_pop <- as.numeric(spec_fixed$default$N)  # population for scaling to counts
#' 
#' # ----------------------- extract β̂(t) and policies ------------------------
#' 
#' beta_df <- calibrator |>
#'   mp_trajectory_sd(conf.int = TRUE, back_transform = TRUE) |>
#'   dplyr::filter(matrix == "beta_thing", time >= offset0, time <= upper_plot_time) |>
#'   mutate(date = seq(start_date, by = "day", length.out = length(value)),
#'          beta_hat = value) |>
#'   select(date, beta_hat) |>
#'   filter(date >= date0, date <= date_end)
#' 
#' 
#' # ---- simulation helpers ---------------------------------------------------
#' set_vax <- function(cf, scale) {
#'   cf$default$double_vac_values  <- scale * cf$default$double_vac_values_base
#'   cf$default$booster_vac_values <- scale * cf$default$booster_vac_values_base
#'   cf
#' }
#' 
#' simulate_series <- function(spec, scale) {
#'   cf <- spec %>% set_vax(scale)
#'   cf %>%
#'     mp_simulator(time_steps, outputs = "inc") %>%
#'     mp_trajectory() %>%
#'     filter(matrix == "inc") %>%
#'     mutate(date = date0 + (time - offset0)) %>%
#'     filter(date >= date0, date <= date_end) %>%
#'     transmute(date, inc = value)
#' }
#' 
#' # ---- baseline (observed @ 100%) -------------------------------------------
#' baseline <- purrr::imap_dfr(vax_scales, function(scl, vlbl) {
#'   simulate_series(spec_fixed, scl) |>                    
#'     dplyr::mutate(
#'       vax = factor(vlbl, levels = levels(vax_fac)),
#'       series = "observed"
#'     ) |>
#'     mutate(inc = as.integer(inc)) |>
#'     tidyr::crossing(tibble(als_scen = als_fac)) |>   # replicate to all als columns
#'     dplyr::select(date, inc, vax, als_scen, series)
#' }) |>
#'   arrange(als_scen)
#' 
#' 
#' base_df <- tibble(date = beta_df$date)
#' in_int  <- function(d, a, b) d >= as.Date(a) & d <= as.Date(b)
#' 
#' policy_df <- base_df %>%
#'   mutate(
#'     ALS = case_when(
#'       in_int(date,"2021-12-15","2021-12-24") ~ "ALS2",
#'       in_int(date,"2021-12-25","2022-01-03") ~ "ALS3",
#'       in_int(date,"2022-01-04","2022-02-07") ~ "ALS4",
#'       in_int(date,"2022-02-08","2022-03-14") ~ "ALS3",
#'       in_int(date,"2022-03-15","2022-05-22") ~ "NoALS",
#'       TRUE ~ "ALS2"
#'     ),
#'     school = if_else(in_int(date,"2021-12-20","2022-01-25"), "closed", "open")
#'   )
#' 
#' # multipliers (rebased so no-als = 1)
#' R_ALS4 <- 1.52; R_ALS3 <- 2.12; R_NoALS <- 2.32
#' R_open <- 2.03; R_closed <- 1.44
#' R_ALS2 <- sqrt(R_ALS3 * R_NoALS)
#' 
#' m_als    <- c(ALS4 = R_ALS4/R_NoALS, ALS3 = R_ALS3/R_NoALS, ALS2 = R_ALS2/R_NoALS, NoALS = 1.00)
#' m_school <- c(open = 1.00, closed = R_closed/R_open)
#' 
#' # smooth categorical timelines to c¹ series
#' smooth_step <- function(t, t0, tau = 2.0) plogis((as.numeric(t - t0))/tau)
#' make_smooth_series <- function(date, lev, map, tau = 2.0) {
#'   lev <- as.character(lev)
#'   idx <- c(1L, which(lev[-1L] != lev[-length(lev)]) + 1L)
#'   segd <- date[idx]; segl <- lev[idx]
#'   m <- rep(unname(map[segl[1L]]), length(date))
#'   if (length(idx) >= 2L) {
#'     for (j in 2:length(idx)) {
#'       s <- smooth_step(date, segd[j], tau)
#'       m <- m*(1 - s) + unname(map[segl[j]])*s
#'     }
#'   }
#'   m
#' }
#' 
#' # recover policy-neutral beta0_hat
#' df_est <- policy_df %>%
#'   left_join(beta_df, by = "date") %>%
#'   mutate(
#'     m_ALS_hat    = make_smooth_series(date, ALS,    map = m_als,    tau = 2.0),
#'     m_school_hat = make_smooth_series(date, school, map = m_school, tau = 2.0),
#'     m_obs_hat    = m_ALS_hat * m_school_hat,
#'     beta0_hat    = beta_hat / pmax(m_obs_hat, 1e-8)
#'   )
#' 
#' # ------------------- target curves for 3 counterfactuals -------------------
#' # keep school(t) as observed, make ALS constant throughout
#' k_ALS <- c(ALS2 = unname(m_als["ALS2"]),
#'            ALS4 = unname(m_als["ALS4"]),
#'            NoALS = unname(m_als["NoALS"]))
#' targets <- list(
#'   ALS2  = df_est$beta0_hat * (k_ALS["ALS2"]  * df_est$m_school_hat),
#'   ALS4  = df_est$beta0_hat * (k_ALS["ALS4"]  * df_est$m_school_hat),
#'   NoALS = df_est$beta0_hat * (k_ALS["NoALS"] * df_est$m_school_hat)
#' )
#' 
#' # -------------------- fit time_var_beta for each scenario ------------------
#' # use the same basis size p as in the calibrated spec
#' p <- length(spec_fixed$default$time_var_beta)
#' nT <- nrow(df_est)
#' tgrid <- if (nT > 1) seq.int(from = offset0, to = time_steps, length.out = nT) else offset0
#' 
#' # note: splines::ns() has no 'degree' argument
#' X <- splines::ns(tgrid, df = p, intercept = TRUE,
#'                  Boundary.knots = c(offset0, time_steps))
#' log_beta0 <- as.numeric(spec_fixed$default$log_beta)
#' 
#' fit_timevar_coef <- function(beta_cf, log_beta0, X, lambda_frac = 1e-6) {
#'   y   <- log(pmax(beta_cf, 1e-12)) - log_beta0
#'   XtX <- crossprod(X); Xty <- crossprod(X, y)
#'   lam <- lambda_frac * mean(diag(XtX))
#'   as.numeric(solve(XtX + lam * diag(ncol(X)), Xty))
#' }
#' 
#' theta_ALS2  <- fit_timevar_coef(targets$ALS2,  log_beta0, X)
#' theta_ALS4  <- fit_timevar_coef(targets$ALS4,  log_beta0, X)
#' theta_NoALS <- fit_timevar_coef(targets$NoALS, log_beta0, X)
#' 
#' print(spec_fixed$default$time_var_beta)
#' print(theta_ALS2)
#' print(theta_ALS4)
#' print(theta_NoALS)
#' 
#' 
#' print(spec_fixed %>%
#'         mp_simulator(time_steps, outputs = "inc") %>%
#'         mp_trajectory() %>%
#'         dplyr::filter(matrix == "inc") %>%
#'         dplyr::mutate(date = date0 + (time - offset0)) %>%
#'         dplyr::filter(date >= date0, date <= date_end) %>%
#'         dplyr::transmute(date, inc = as.integer(value))
#'       
#' )
#' 
#' 
#' 
#' print(spec_fixed %>%
#'   mp_simulator(time_steps, default = list(timevar_beta ~ theta_ALS2), outputs = "inc") %>%
#'   mp_trajectory() %>%
#'   dplyr::filter(matrix == "inc") %>%
#'   dplyr::mutate(date = date0 + (time - offset0)) %>%
#'   dplyr::filter(date >= date0, date <= date_end) %>%
#'   dplyr::transmute(date, inc = as.integer(value))
#' )
#' 
#' 
#' 
#' quit()
#' 
#' 
#' # --------------------- build specs (neutralize policy series) --------------
#' get_series_len <- function(spec) {
#'   d <- spec$default
#'   if (!is.null(d$time_var_beta)) return(length(d$time_var_beta))
#'   if (!is.null(d$m_ALS))         return(length(d$m_ALS))
#'   if (!is.null(d$m_school))      return(length(d$m_school))
#'   stop("cannot infer series length")
#' }
#' L_series <- get_series_len(spec_fixed)
#' onesL    <- rep(1, L_series)
#' 
#' make_spec_with_theta <- function(spec, theta, ones) {
#'   s <- spec
#'   stopifnot(length(theta) == length(s$default$time_var_beta))
#'   s$default$time_var_beta <- as.numeric(theta)
#'   if (!is.null(s$default$m_ALS))    s$default$m_ALS    <- ones
#'   if (!is.null(s$default$m_school)) s$default$m_school <- ones
#'   s
#' }
#' 
#' spec_ALS2   <- make_spec_with_theta(spec_fixed, theta_ALS2,  onesL)
#' spec_ALS4   <- make_spec_with_theta(spec_fixed, theta_ALS4,  onesL)
#' spec_NoALS  <- make_spec_with_theta(spec_fixed, theta_NoALS, onesL)
#' 
#' # ------------------------- simulation helpers ------------------------------
#' set_vax_rate <- function(spec, scale) {
#'   spec$default$double_vac_values  <- scale * spec$default$double_vac_values_base
#'   spec$default$booster_vac_values <- scale * spec$default$booster_vac_values_base
#'   spec
#' }
#' 
#' # convert flow to counts if it looks like a fraction
#' to_counts <- function(df, N) {
#'   mx <- max(df$inc, na.rm = TRUE)
#'   if (is.finite(mx) && mx < 1) df$inc <- df$inc * N
#'   df
#' }
#' 
#' simulate_series <- function(spec, scale, N) {
#'   spec |> set_vax_rate(scale) |>
#'     mp_simulator(time_steps, outputs = "inc") |>
#'     mp_trajectory() |>
#'     dplyr::filter(matrix == "inc") |>
#'     dplyr::mutate(date = date0 + (time - offset0)) |>
#'     dplyr::filter(date >= date0, date <= date_end) |>
#'     dplyr::transmute(date, inc = value) |>
#'     to_counts(N)
#' }
#' 
#' # ----------------------- baseline and counterfactuals ----------------------
#' baseline <- simulate_series(spec_fixed, 1.00, N_pop) %>% mutate(series = "observed_100")
#' 
#' cf_specs <- list("ALS-2" = spec_ALS2, "ALS-4" = spec_ALS4, "No-ALS" = spec_NoALS)
#' 
#' cf_sims <- imap_dfr(cf_specs, function(sp, als_name) {
#'   imap_dfr(vax_scales, function(scale, vlab) {
#'     simulate_series(sp, scale, N_pop) %>%
#'       mutate(vax = vlab, als = als_name, series = "counterfactual")
#'   })
#' })
#' 
#' # --------------------------- panel plot (3 x 3) ----------------------------
#' # make a strict panel ordering: rows = vaccination %, cols = als scenario
#' vax_levels <- names(vax_scales)                 # "25%","50%","75%"
#' als_levels <- c("ALS-2","ALS-4","No-ALS")
#' 
#' # replicate baseline (observed_100) into every facet so it overlays in red
#' baseline_rep <- baseline %>%
#'   tidyr::crossing(
#'     tibble(vax = factor(vax_levels, levels = vax_levels),
#'            als = factor(als_levels, levels = als_levels))
#'   ) %>%
#'   mutate(series = "observed_100")
#' 
#' # align factor levels on the cf data
#' cf_sims_f <- cf_sims %>%
#'   mutate(
#'     vax = factor(vax, levels = vax_levels),
#'     als = factor(als, levels = als_levels)
#'   )
#' 
#' # combine for plotting
#' plot_df <- bind_rows(
#'   cf_sims_f %>% mutate(series = "counterfactual"),
#'   baseline_rep
#' )
#' 
#' # human-friendly panel titles
#' plot_df <- plot_df %>%
#'   mutate(panel_title = paste0("vaccination ", vax, "  |  ", als, " throughout"))
#' 
#' # color mapping
#' cols <- c("counterfactual" = "#1f77b4",   # blue
#'           "observed_100"  = "#d62728")   # red
#' 
#' # build the figure: 3 rows (25,50,75) × 3 cols (als-2,4,no-als)
#' p_panels <- ggplot(plot_df, aes(x = date, y = inc, group = series, color = series)) +
#'   geom_line(linewidth = 0.9, alpha = 0.95) +
#'   ggh4x::facet_wrap2(~ panel_title, nrow = 3, ncol = 3, axes = "all") +  # axes on every panel
#'   scale_color_manual(values = cols, name = NULL,
#'                      breaks = c("counterfactual","observed_100"),
#'                      labels = c("counterfactual", "observed (100%)")) +
#'   scale_x_date(expand = c(0,0), date_breaks = "2 weeks", date_labels = "%b %d") +
#'   labs(
#'     title = "counterfactual incidence (blue) vs observed (red), by vaccination level and als policy",
#'     x = NULL,
#'     y = "daily infections (count)"
#'   ) +
#'   theme_clean(base_size = 12) +
#'   theme(
#'     strip.text   = element_text(size = 12),
#'     axis.text.x  = element_text(size = 9, angle = 45, hjust = 1),
#'     axis.text.y  = element_text(size = 9),
#'     axis.title.y = element_text(size = 11),
#'     plot.title   = element_text(size = 13, hjust = 0.5),
#'     legend.position = "top",
#'     panel.spacing = unit(0.8, "lines")
#'   )
#' 
#' # show on screen
#' print(p_panels)
#' 
#' png("../figures/counterfactual_trials.png",width = 5200, height = 4200, res = 300, bg = "white", type = "cairo")
#' print(p_panels)
#' dev.off()
#' 
#' quit()
#' 
#' 
#' 
#' 
#' 
#' 
#' # ---- bind the two dataframes for plotting --------------------------------
#' plot_df <- bind_rows(
#'   grid_sims,           # series == 'counterfactual' (blue)
#'   observed_from_grid   # series == 'observed'       (red)
#' )
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' simulate_with_theta <- function(spec, theta, ones, time_steps, date0, offset0, date_end) {
#'   spec2 <- make_spec_with_theta(spec, theta, ones)
#'   mp_simulator(spec2, time_steps = time_steps, outputs = "inc") |>
#'     mp_trajectory() |>
#'     dplyr::filter(matrix == "inc") |>
#'     dplyr::mutate(date = date0 + (time - offset0)) |>
#'     dplyr::filter(date >= date0, date <= date_end) |>
#'     dplyr::select(date, value) %>%
#'     rename(inc = value)
#' }
#' 
#' sim_ALS2  <- simulate_with_theta(spec_fixed, theta_ALS2,  onesL, time_steps, date0, offset0, date_end) %>% mutate(scenario = "als-2 throughout")
#' sim_ALS4  <- simulate_with_theta(spec_fixed, theta_ALS4,  onesL, time_steps, date0, offset0, date_end) %>% mutate(scenario = "als-4 throughout")
#' sim_NoALS <- simulate_with_theta(spec_fixed, theta_NoALS, onesL, time_steps, date0, offset0, date_end) %>% mutate(scenario = "no-als throughout")
#' 
#' cf_inc <- bind_rows(sim_ALS2, sim_ALS4, sim_NoALS)
#' 
#' # ---- optional: quick plot (can be removed in batch) -----------------------
#' p_cf <- ggplot(cf_inc, aes(date, inc, color = scenario)) +
#'   geom_line(linewidth = 1.1) +
#'   scale_color_manual(values = c("als-2 throughout" = "#1b9e77",
#'                                 "als-4 throughout" = "#d95f02",
#'                                 "no-als throughout" = "#7570b3")) +
#'   labs(x = NULL, y = "estimated infections", color = NULL,
#'        title = "counterfactual simulations: constant ALS throughout (schools as observed)") +
#'   theme_clean() +
#'   theme(plot.title = element_text(hjust = 0.5))
#' print(p_cf)
#' 
#' quit()
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' set_matrix_coef_beta <- function(coef_vec, neutralize_policy = TRUE) {
#'   #' set matrix_coef_beta into a macpan2 spec (safe copy–edit–write)
#'   #'
#'   #' @description
#'   #' updates `matrix_coef_beta` in a macpan2 model specification by:
#'   #' (1) reading matrices via `all_matrices()`, (2) replacing the beta-coefficient
#'   #' vector, and (3) writing matrices back using `set_matrices()` or
#'   #' `with_matrices()` (falls back to assigning `spec$matrices`). optionally
#'   #' neutralizes policy multipliers (`m_ALS`, `m_school`) to avoid double counting
#'   #' when policy effects are already embedded in the supplied beta curve.
#'   #'
#'   #' @param spec a macpan2 model spec (e.g., from `mp_optimized_spec()`).
#'   #' @param coef_vec numeric vector for `matrix_coef_beta`; length must equal
#'   #'   `NROW(spec$all_matrices()$matrix_coef_beta)`. values should be on the log scale
#'   #'   if your model links beta via `log()`.
#'   #' @param neutralize_policy logical; if `TRUE` (default), set `m_ALS` and
#'   #'   `m_school` sequences in `inputs/default` to 1 to prevent double counting.
#'   #'
#'   #' @return the modified spec, ready to simulate.
#'   #'
#'   #' @examples
#'   #' \dontrun{
#'   #' # coef_vec <- construct_spline_basis(beta_cf)  # your fitted coefficients
#'   #' spec_cf <- set_matrix_coef_beta(coef_vec, neutralize_policy = TRUE)
#'   #' # sim <- simulate_inc(spec_cf)
#'   #' }
#'   #' @export
#'   
#'   # read current matrices (copy)
#'   mats <- spec_fixed$all_matrices()
#'   if (is.null(mats$matrix_coef_beta))
#'     stop("spec_fixed$all_matrices() lacks 'matrix_coef_beta'.")
#'   
#'   # shape check + replace
#'   p <- NROW(mats$matrix_coef_beta)
#'   if (length(coef_vec) != p)
#'     stop("length(coef_vec) = ", length(coef_vec),
#'          " but NROW(matrix_coef_beta) = ", p)
#'   mats$matrix_coef_beta <- matrix(as.numeric(coef_vec), ncol = 1)
#'   
#'   # write back using public api (prefer update$matrices; else replace$matrices; else direct)
#'   upd <- spec_fixed[["update"]]
#'   rep <- spec_fixed[["replace"]]
#'   if (!is.null(upd) && is.list(upd) && is.function(upd$matrices)) {
#'     spec_fixed <- upd$matrices(mats)
#'   } else if (!is.null(rep) && is.list(rep) && is.function(rep$matrices)) {
#'     spec_fixed <- rep$matrices(mats)
#'   } else {
#'     # last resort: assign to $matrices (works on many builds)
#'     spec_fixed$matrices <- mats
#'   }
#'   
#'   spec_fixed
#' }
#' 
#' 
#' # create spec copies with scenario-specific coefficients written in
#' spec_only_ALS2  <- set_matrix_coef_beta(coef_ALS2)
#' spec_only_ALS4  <- set_matrix_coef_beta(coef_ALS4)
#' spec_only_NoALS <- set_matrix_coef_beta(coef_NoALS)
#' 
#' print(spec_only_ALS4)
#' 
#' quit()
#' 
#' # helper: simulate incidence for a given spec (keeps your date window)
#' simulate_inc <- function(spec) {
#'   spec %>%
#'     mp_simulator(time_steps, outputs = "inc") %>%
#'     mp_trajectory() %>%
#'     dplyr::filter(matrix == "inc") %>%
#'     dplyr::mutate(date = date0 + (time - offset0)) %>%
#'     dplyr::filter(date >= date0, date <= date_end) %>%
#'     dplyr::transmute(date, inc = value)
#' }
#' 
#' # simulate observed @100% (time-varying policies) for overlay if needed
#' baseline_100 <- simulate_series_observed(1.00) %>% mutate(series = "observed")
#' 
#' 
#' 
#' 
#' 
#' 
#' # write back via all_matrices(); also neutralize policy multipliers (avoid double count)
#' set_matrix_coef_beta <- function(spec, coef_vec) {
#'   mats <- spec$all_matrices()
#'   stopifnot("matrix_coef_beta" %in% names(mats))
#'   mats$matrix_coef_beta <- matrix(as.numeric(coef_vec), ncol = 1)
#'   
#'   if (!is.null(spec$set_matrices) && is.function(spec$set_matrices)) {
#'     spec <- spec$set_matrices(mats)
#'   } else if (!is.null(spec$with_matrices) && is.function(spec$with_matrices)) {
#'     spec <- spec$with_matrices(mats)
#'   } else {
#'     spec$matrices$matrix_coef_beta <- mats$matrix_coef_beta
#'   }
#'   
#'   nT <- length(coef_vec)
#'   set1 <- function(lst, key) { if (!is.null(lst[[key]])) lst[[key]] <- rep(1.0, nT); lst }
#'   if (!is.null(spec$inputs))  { spec$inputs  <- set1(spec$inputs,  "m_ALS"); spec$inputs  <- set1(spec$inputs,  "m_school") }
#'   if (!is.null(spec$default)) { spec$default <- set1(spec$default, "m_ALS"); spec$default <- set1(spec$default, "m_school") }
#'   spec
#' }
#' 
#' # build target beta vectors from your df_est_cf (school(t) kept as observed)
#' beta_cf_ALS2  <- df_est_cf$beta_hat_ALS2
#' beta_cf_ALS4  <- df_est_cf$beta_hat_ALS4
#' beta_cf_NoALS <- df_est_cf$beta_hat_NoALS
#' 
#' # fit coefficients for each constant-als scenario
#' coef_ALS2  <- fit_coef_from_beta(beta_cf_ALS2,  spec_fixed)
#' coef_ALS4  <- fit_coef_from_beta(beta_cf_ALS4,  spec_fixed)
#' coef_NoALS <- fit_coef_from_beta(beta_cf_NoALS, spec_fixed)
#' 
#' # create scenario-specific specs with updated matrix_coef_beta
#' spec_ALS2  <- set_matrix_coef_beta(spec_fixed, coef_ALS2)
#' spec_ALS4  <- set_matrix_coef_beta(spec_fixed, coef_ALS4)
#' spec_NoALS <- set_matrix_coef_beta(spec_fixed, coef_NoALS)
#' 
#' # optional: vacc scaling helper (reuses your set_vax)
#' # spec_ALS2  <- set_vax(spec_ALS2,  1.00)
#' # spec_ALS4  <- set_vax(spec_ALS4,  1.00)
#' # spec_NoALS <- set_vax(spec_NoALS, 1.00)
#' 
#' # simulate helper (unchanged)
#' simulate_inc <- function(spec) {
#'   spec %>%
#'     mp_simulator(time_steps, outputs = "inc") %>%
#'     mp_trajectory() %>%
#'     dplyr::filter(matrix == "inc") %>%
#'     dplyr::mutate(date = date0 + (time - offset0)) %>%
#'     dplyr::filter(date >= date0, date <= date_end) %>%
#'     dplyr::transmute(date, inc = value)
#' }
#' 
#' # examples:
#' # sim_ALS2  <- simulate_inc(spec_ALS2)
#' # sim_ALS4  <- simulate_inc(spec_ALS4)
#' # sim_NoALS <- simulate_inc(spec_NoALS)
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' # ----- helpers
#' set_vax <- function(cf, vax_rate) {
#'   cf$default$double_vac_values  <- vax_rate * cf$default$double_vac_values_base
#'   cf$default$booster_vac_values <- vax_rate * cf$default$booster_vac_values_base
#'   cf
#' }
#' simulate_series_observed <- function(vax_rate) {
#'   cf <- spec_fixed %>% set_vax(vax_rate)
#'   cf %>%
#'     mp_simulator(time_steps, outputs = "inc") %>%
#'     mp_trajectory() %>%
#'     filter(matrix == "inc") %>%
#'     mutate(date = date0 + (time - offset0)) %>%
#'     filter(date >= date0, date <= date_end) %>%
#'     transmute(date, inc = value)
#' }
#' 
#' df_vax_rate <- simulate_series_observed(0.25)
#' 
#' 
#' 
#' quit()
#' 
#' make_beta_cf <- function(df_est, als_level, m_als) {
#'   k <- unname(m_als[[als_level]])
#'   df_est %>% transmute(date, beta_cf = k * beta0_hat)
#' }
#' 
#' # robust setter: try to inject beta_cf; if impossible, flag fallback
#' try_set_beta_in_spec <- function(spec, dates, beta_cf_vec) {
#'   # return list(spec=..., mode="inject" | "fallback")
#'   out <- list(spec = spec, mode = "fallback")
#'   # path 1: spec$inputs$beta_t
#'   if (!is.null(spec$inputs) && !is.null(spec$inputs$beta_t)) {
#'     spec$inputs$beta_t <- beta_cf_vec
#'     if (!is.null(spec$inputs$m_ALS))    spec$inputs$m_ALS    <- rep(1.0, length(dates))
#'     if (!is.null(spec$inputs$m_school)) spec$inputs$m_school <- rep(1.0, length(dates))
#'     out$spec <- spec; out$mode <- "inject"; return(out)
#'   }
#'   # path 2: spec$inputs$beta
#'   if (!is.null(spec$inputs) && !is.null(spec$inputs$beta)) {
#'     spec$inputs$beta <- beta_cf_vec
#'     if (!is.null(spec$inputs$m_ALS))    spec$inputs$m_ALS    <- rep(1.0, length(dates))
#'     if (!is.null(spec$inputs$m_school)) spec$inputs$m_school <- rep(1.0, length(dates))
#'     out$spec <- spec; out$mode <- "inject"; return(out)
#'   }
#'   # path 3: some models keep beta under spec$default$beta_t
#'   if (!is.null(spec$default) && !is.null(spec$default$beta_t)) {
#'     spec$default$beta_t <- beta_cf_vec
#'     if (!is.null(spec$default$m_ALS))    spec$default$m_ALS    <- rep(1.0, length(dates))
#'     if (!is.null(spec$default$m_school)) spec$default$m_school <- rep(1.0, length(dates))
#'     out$spec <- spec; out$mode <- "inject"; return(out)
#'   }
#'   out
#' }
#' 
#' simulate_with_beta_or_fallback <- function(spec, beta_cf_df, vax_scale, df_est, als_level) {
#'   dates <- seq(date0, date_end, by = "day")
#'   beta_cf_df <- inner_join(tibble(date = dates), beta_cf_df, by = "date")
#'   
#'   spec_v <- set_vax(spec, vax_scale)
#'   res <- try_set_beta_in_spec(spec_v, dates, beta_cf_df$beta_cf)
#'   
#'   if (res$mode == "inject") {
#'     # true injection path
#'     res$spec %>%
#'       mp_simulator(time_steps, outputs = "inc") %>%
#'       mp_trajectory() %>%
#'       filter(matrix == "inc") %>%
#'       mutate(date = date0 + (time - offset0)) %>%
#'       filter(date >= date0, date <= date_end) %>%
#'       transmute(date, inc = value, mode = "inject")
#'   } else {
#'     # fallback: proportional incidence scaling (clearly marked)
#'     # build ratio = (k_ALS * m_school_obs) / (m_ALS_obs * m_school_obs) = k_ALS / m_ALS_obs
#'     k <- unname(m_als[[als_level]])
#'     df_ratio <- df_est %>%
#'       transmute(date, ratio = pmax(k / pmax(m_ALS_hat, 1e-8), 0))
#'     basev <- spec_v %>%
#'       mp_simulator(time_steps, outputs = "inc") %>%
#'       mp_trajectory() %>%
#'       filter(matrix == "inc") %>%
#'       mutate(date = date0 + (time - offset0)) %>%
#'       filter(date >= date0, date <= date_end) %>%
#'       transmute(date, inc = value)
#'     basev %>%
#'       left_join(df_ratio, by = "date") %>%
#'       mutate(inc = pmax(inc * ratio, 0), mode = "fallback")
#'   }
#' }
#' 
#' # ----- build constant-als beta_cf schedules
#' beta_cf_list <- list(
#'   ALS2  = make_beta_cf(df_est, "ALS2",  m_als),
#'   ALS4  = make_beta_cf(df_est, "ALS4",  m_als),
#'   NoALS = make_beta_cf(df_est, "NoALS", m_als)
#' )
#' 
#' als_levels_const <- c("ALS2","ALS4","NoALS")
#' als_labels_fac   <- c("als-2 throughout","als-4 throughout","no-als throughout")
#' 
#' # ----- simulate the 3×3 grid (preferring injection; falling back if needed)
#' grid_sims <- map_dfr(names(vax_scales), function(vlbl) {
#'   scl <- vax_scales[[vlbl]]
#'   map_dfr(als_levels_const, function(alsL) {
#'     simulate_with_beta_or_fallback(spec_fixed, beta_cf_list[[alsL]], scl, df_est, alsL) %>%
#'       mutate(
#'         vax      = factor(vlbl, levels = levels(vax_fac)),
#'         als_scen = factor(dplyr::case_match(alsL,
#'                                             "ALS2" ~ "als-2 throughout",
#'                                             "ALS4" ~ "als-4 throughout",
#'                                             "NoALS"~ "no-als throughout"),
#'                           levels = levels(als_fac)),
#'         series   = "counterfactual"
#'       )
#'   })
#' })
#' 
#' # optional: print which path was used
#' if (any(grid_sims$mode == "fallback")) {
#'   message("note: spec does not expose a writable beta(t); used proportional fallback for some/all scenarios.")
#' }
#' 
#' # observed comparator at 100% (as fitted policies)
#' baseline <- simulate_series_observed(1.00) %>% mutate(series = "observed")
#' 
#' # replicate observed into every facet for overlay
#' observed_from_grid <- grid_sims %>%
#'   select(date, vax, als_scen) %>%
#'   left_join(baseline %>% select(date, inc), by = "date") %>%
#'   mutate(series = "observed") %>%
#'   relocate(series, .after = inc)
#' 
#' plot_df <- bind_rows(grid_sims, observed_from_grid)
#' 
#' # ----- panel titles, annotations, plot
#' panel_levels <- as.vector(outer(levels(vax_fac), levels(als_fac),
#'                                 function(v, a) paste0("vaccination ", v, " || ", a)))
#' plot_df2 <- plot_df %>%
#'   mutate(panel_title = paste0("vaccination ", vax, " || ", als_scen),
#'          panel_title = factor(panel_title, levels = panel_levels))
#' 
#' cf_ann <- plot_df2 %>%
#'   filter(series == "counterfactual") %>%
#'   group_by(vax, als_scen) %>%
#'   summarize(
#'     cumulative_infections = sum(inc, na.rm = TRUE),
#'     peak_infections       = max(inc, na.rm = TRUE),
#'     peak_date             = date[which.max(inc)],
#'     .groups = "drop"
#'   ) %>%
#'   mutate(
#'     panel_title = factor(paste0("vaccination ", vax, " || ", als_scen),
#'                          levels = panel_levels),
#'     ann_lab = paste0(
#'       "peak date: ", format(peak_date, "%B %d, %Y"), "\n",
#'       "peak cases: ", scales::comma(round(peak_infections)), "\n",
#'       "total cases: ", scales::comma(round(cumulative_infections))
#'     )
#'   )
#' 
#' ann_pos <- plot_df2 %>%
#'   group_by(vax, als_scen) %>%
#'   summarize(.groups = "drop") %>%
#'   mutate(
#'     panel_title = factor(paste0("vaccination ", vax, " || ", als_scen),
#'                          levels = panel_levels),
#'     xpos = date0 + 5, ypos = 6100
#'   )
#' ann_df <- left_join(cf_ann, ann_pos, by = c("vax","als_scen","panel_title"))
#' 
#' p_3x3 <- ggplot() +
#'   geom_line(
#'     data = filter(plot_df2, series == "counterfactual"),
#'     aes(x = date, y = inc, group = 1),
#'     color = "blue", linewidth = 1.5
#'   ) +
#'   geom_line(
#'     data = filter(plot_df2, series == "observed"),
#'     aes(x = date, y = inc, group = 1),
#'     color = "red", linewidth = 1.5
#'   ) +
#'   geom_text(
#'     data = ann_df,
#'     aes(x = xpos, y = ypos, label = ann_lab),
#'     inherit.aes = FALSE,
#'     size = 5, lineheight = 0.98, hjust = 0, vjust = 1
#'   ) +
#'   ggh4x::facet_wrap2(~ panel_title, nrow = 3, ncol = 3, axes = "all") +
#'   labs(title = "counterfactual scenario (blue) vs estimated infections (red)",
#'        y = "estimated infections", x = NULL) +
#'   scale_x_date(limits = c(date0, date_end), expand = c(0,0),
#'                date_breaks = "2 weeks", date_labels = "%b %d") +
#'   theme_clean() +
#'   theme(
#'     axis.text.x  = element_text(size = 16, angle = 45, hjust = 1),
#'     axis.text.y  = element_text(size = 16),
#'     axis.title.y = element_text(size = 18),
#'     plot.title   = element_text(size = 22, hjust = 0.5),
#'     strip.text   = element_text(size = 20),
#'     legend.position = "bottom",
#'     panel.spacing   = unit(0.7, "lines"),
#'     plot.background = element_blank()
#'   )
#' 
#' png("../figures/observed_vs_counterfactuals_beta0_scaled_robust.png",
#'     width = 5200, height = 4200, res = 300, bg = "white", type = "cairo")
#' print(p_3x3)
#' dev.off()

# ----- methods note
# preferred path: inject beta_cf(t) = k_ALS * beta0(t) into the simulator while neutralizing other policy multipliers.
# fallback (only if spec exposes no writable beta path): proportional scaling of incidence by k_ALS / m_ALS_obs(t).
# the script auto-detects and reports if fallback was used.

# --------------------------- methods note (brief) ---------------------------
# we extract β̂(t) from the fitted calibrator (“beta_thing”), compute policy-neutral β0(t)
# by dividing β̂(t) by smooth als×school multipliers, then impose constant policies by
# setting β_cf(t) = k_ALS · β0(t) (als-2, als-4, no-als) and re-simulate under each
# vaccination scale (75%, 50%, 25%). observed 100% vaccination with time-varying policies
# is overlaid in every panel. annotations report per-panel total infections, peak size, and peak date.


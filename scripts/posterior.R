suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(ggthemes)
  library(shellpipes)
  library(tidyverse)
  library(macpan2)
  library(grid)
})

set.seed(2025)
options(macpan2_log_dir = ".")
loadEnvironments()

# ---- Output directory ----
fig_dir <- "../figures"

# ---- Inputs ----
calibrator <- readRDS("calibrate.rds")
sim_horizon <- 530

date0 <- as.Date("2021-12-15")
offset0 <- 150

# ---- Extract fitted coefficients ----
model_estimates <- mp_tmb_coef(calibrator, conf.int = TRUE) |>
  dplyr::select(-term, -col, -type) |>
  mutate(mat = abbreviate(mat, minlength = 15))

print(model_estimates, digits = 2)

# =====================================================================
# Posterior summary plot
# =====================================================================

keep_params <- c("gamma_a", "gamma_i", "kappa2", "kappa3", "sigma")

lab_map <- c(
  gamma_a = "gamma[a]",
  gamma_i = "gamma[i]",
  kappa2  = "kappa[2]",
  kappa3  = "kappa[3]",
  sigma   = "sigma"
)

df_keep <- model_estimates %>%
  filter(mat %in% keep_params) %>%
  filter(is.finite(estimate), is.finite(conf.low), is.finite(conf.high)) %>%
  mutate(param = forcats::fct_relevel(mat, keep_params))

posterior_plt <- ggplot(df_keep, aes(x = estimate, y = param)) +
  geom_vline(xintercept = 0, linetype = 3, linewidth = 0.6) +
  geom_errorbar(
    aes(xmin = conf.low, xmax = conf.high),
    orientation = "y",
    width = 0.2,
    linewidth = 0.7
  ) +
  geom_point(size = 2.4) +
  scale_y_discrete(labels = function(x) parse(text = lab_map[x])) +
  labs(
    title = "Posterior summaries of fitted model parameters",
    x = "Estimate (95% CI)",
    y = NULL
  ) +
  theme_clean() +
  theme(
    axis.text.x        = element_text(size = 12),
    axis.text.y        = element_text(size = 12),
    axis.title.x       = element_text(size = 14),
    axis.title.y       = element_text(size = 14),
    plot.title         = element_text(size = 16, hjust = 0.5),
    panel.grid.major.y = element_line(),
    plot.background    = element_blank()
  )

png(
  file.path(fig_dir, "posterior_plot.png"),
  width = 3500,
  height = 1500,
  res = 300,
  bg = "white",
  type = "cairo"
)
print(posterior_plt)
dev.off()

# =====================================================================
# Prior and posterior marginal distributions
# =====================================================================

set.seed(2025)
n_draw <- 40000L

lab_unicode <- c(
  gamma_a = "γₐ",
  gamma_i = "γᵢ",
  kappa2  = "κ₂",
  kappa3  = "κ₃",
  sigma   = "σ"
)

prior_from_range <- function(rng, trans) {
  mu <- (trans(rng[1]) + trans(rng[2])) / 2
  log_sd <- log((trans(rng[2]) - trans(rng[1])) / (2 * 1.96))
  sd <- exp(log_sd)
  list(mu = mu, sd = sd)
}

post_on_log <- function(est, se) {
  est <- pmax(est, 1e-12)
  list(mu = log(est), sd = abs(se / est))
}

post_on_logit <- function(est, se) {
  p <- pmin(pmax(est, 1e-8), 1 - 1e-8)
  list(mu = qlogis(p), sd = abs(se / (p * (1 - p))))
}

me <- model_estimates %>%
  filter(mat %in% keep_params) %>%
  select(mat, estimate, std.error) %>%
  filter(is.finite(estimate), is.finite(std.error), std.error > 0)

make_draws_one <- function(mat, est, se) {
  
  if (mat %in% c("gamma_a", "gamma_i", "sigma")) {
    
    pr <- prior_from_range(prior_range[[mat]], log)
    po <- post_on_log(est, se)
    
    prior_val <- exp(rnorm(n_draw, pr$mu, pr$sd))
    post_val  <- exp(rnorm(n_draw, po$mu, po$sd))
    
  } else {
    
    pr <- prior_from_range(prior_range[[mat]], qlogis)
    po <- post_on_logit(est, se)
    
    prior_val <- plogis(rnorm(n_draw, pr$mu, pr$sd))
    post_val  <- plogis(rnorm(n_draw, po$mu, po$sd))
  }
  
  tibble(
    mat   = mat,
    value = c(prior_val, post_val),
    dist  = rep(c("Prior", "Posterior"), each = n_draw)
  )
}

draws_df <- purrr::pmap_dfr(
  list(me$mat, me$estimate, me$std.error),
  make_draws_one
) %>%
  mutate(
    dist = factor(dist, levels = c("Prior", "Posterior")),
    mat  = factor(mat, levels = keep_params),
    facet_lab = factor(
      lab_unicode[as.character(mat)],
      levels = unname(lab_unicode[keep_params])
    )
  )

posterior <- ggplot(draws_df, aes(x = value, fill = dist, color = dist)) +
  geom_density(alpha = 0.35, linewidth = 0.8, adjust = 1.1) +
  facet_wrap(~ facet_lab, scales = "free", ncol = 5) +
  scale_fill_manual(values = c(Prior = "grey30", Posterior = "orange")) +
  scale_color_manual(values = c(Prior = "grey30", Posterior = "orange")) +
  labs(
    title = "Prior and posterior marginal distributions",
    x = NULL,
    y = NULL,
    fill = NULL,
    color = NULL
  ) +
  theme_clean() +
  theme(
    axis.text.x     = element_text(size = 20),
    axis.text.y     = element_text(size = 20),
    plot.title      = element_text(size = 24, hjust = 0.5),
    strip.text      = element_text(size = 20),
    legend.position = "top",
    legend.text     = element_text(size = 20),
    panel.spacing   = unit(1.5, "cm"),
    plot.background = element_blank()
  )

png(
  file.path(fig_dir, "prior_posterior.png"),
  width = 7500,
  height = 1500,
  res = 300,
  bg = "white",
  type = "cairo"
)
print(posterior)
dev.off()
suppressPackageStartupMessages({
  library(macpan2)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(ggplot2)
  library(scales)
  library(ggthemes)
  library(shellpipes)
})

loadEnvironments()

# Read time-varying model spec produced by timevar_spec.R
timevar_spec <- rdsRead("timevar_spec")

# outputs from transmission model
out_vars <- c(
  "S", "E", "A", "I", "R",
  "inc_cohort1", "inc_cohort2", "inc_cohort3", "inc",
  "serop_total", "lambda", "beta_thing"
)

traj <- timevar_spec |>
  mp_simulator(time_steps = 80, outputs = out_vars) |>
  mp_trajectory() |>
  mutate(
    matrix = as.character(matrix),
    row = as.integer(as.character(forcats::as_factor(row))),
    age = factor(row + 1, levels = 1:3, labels = age_labels)
  )

# one-sided discrete convolution
convolve_delay <- function(x, w) {
  n <- length(x)
  k <- length(w)
  out <- numeric(n)
  for (t in seq_len(n)) {
    idx <- seq_len(min(t, k))
    out[t] <- sum(x[t - idx + 1] * w[idx])
  }
  out
}

# extract infection flows by age and cohort
inc_long <- traj |>
  filter(matrix %in% c("inc_cohort1", "inc_cohort2", "inc_cohort3")) |>
  select(time, age, matrix, value) |>
  pivot_wider(names_from = matrix, values_from = value)

risk_tbl <- tibble(
  age = factor(age_labels, levels = age_labels),
  p_hosp1 = p_hosp1, p_hosp2 = p_hosp2, p_hosp3 = p_hosp3,
  p_d1    = p_d1,    p_d2    = p_d2,    p_d3    = p_d3
)

# same-day expected severe outcomes from infections
sev_instant <- inc_long |>
  left_join(risk_tbl, by = "age") |>
  mutate(
    hosp_from_inf  = p_hosp1 * inc_cohort1 + p_hosp2 * inc_cohort2 + p_hosp3 * inc_cohort3,
    death_from_inf = p_d1    * inc_cohort1 + p_d2    * inc_cohort2 + p_d3    * inc_cohort3
  )

# delayed outcomes
sev_delayed <- sev_instant |>
  group_by(age) |>
  arrange(time, .by_group = TRUE) |>
  mutate(
    hosp_adm  = convolve_delay(hosp_from_inf,  w_hosp),
    death_obs = convolve_delay(death_from_inf, w_death)
  ) |>
  ungroup()

# long form for severe outputs
sev_long <- sev_delayed |>
  select(time, age, hosp_from_inf, death_from_inf, hosp_adm, death_obs) |>
  pivot_longer(
    cols = -c(time, age),
    names_to = "matrix",
    values_to = "value"
  )

# keep transmission outputs for plotting
traj_plot <- traj |>
  filter(matrix %in% c("S", "E", "A", "I", "R", "inc", "serop_total", "lambda", "beta_thing")) |>
  select(time, age, matrix, value)

plot_df <- bind_rows(traj_plot, sev_long) |>
  mutate(
    matrix = factor(
      matrix,
      levels = c(
        "S", "E", "A", "I", "R",
        "inc",
        "hosp_from_inf", "death_from_inf",
        "hosp_adm", "death_obs",
        "serop_total", "lambda", "beta_thing"
      )
    )
  )

matrix_labeller <- c(
  S = "Susceptible",
  E = "Exposed",
  A = "Asymptomatic",
  I = "Symptomatic",
  R = "Recovered",
  inc = "Infections (Total)",
  hosp_from_inf = "Expected hospitalizations from infections",
  death_from_inf = "Expected deaths from infections",
  hosp_adm = "Hospital admissions (delayed)",
  death_obs = "Deaths (delayed)",
  serop_total = "Seroprevalence (R/N)",
  lambda = "Transmission rate (lambda)",
  beta_thing = "Transmission parameter (beta)"
)

p <- ggplot(plot_df, aes(x = time, y = value, colour = age)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(
    ~ matrix,
    scales = "free",
    labeller = labeller(matrix = matrix_labeller)
  ) +
  labs(
    title = "Age-structured SEAIR with delayed severe outcomes",
    subtitle = paste("Total pop =", label_number(scale_cut = cut_short_scale())(sum(N_age))),
    x = "Days from start",
    y = "Value"
  ) +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    plot.title = element_text(size = 25, hjust = 0.5),
    strip.text = element_text(size = 20),
    strip.background = element_rect(fill = "grey90", color = NA),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box.just = "left",
    legend.background = element_blank(),
    legend.text = element_text(size = 20),
    plot.background = element_blank(),
    plot.margin = margin(20, 20, 20, 40),
    panel.spacing.x = unit(1.2, "cm"),
    panel.grid.major.x = element_blank()
  )

print(p)


png("../figures/severe.png", width = 9500, height = 4500, res = 300, bg = "white", type = "cairo")
print(p)
dev.off()

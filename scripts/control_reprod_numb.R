suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(ggplot2)
  library(tidyr)
  library(zoo)
  library(ggthemes)
  library(patchwork)
  library(shellpipes)
  library(tidyverse)
  library(macpan2)
  library(grid)
})

theme_set(theme_bw())

set.seed(2025)
options(macpan2_log_dir = ".")
loadEnvironments()

# ---- inputs ----
calibrator     <- readRDS("calibrate.rds")
seroprevdata   <- readRDS("seroprevdata.rds")

beta_ALS4only  <- rdsRead("beta_als4")
beta_noALS     <- rdsRead("policy_neutral")

upper_plot_time <- 300
offset0         <- 150
anchor_start    <- as.Date("2021-12-15")
anchor_end      <- as.Date("2022-05-22")
cut_date        <- as.Date("2022-03-15")

# ---- model constants ----
kappa1 <- params[["kappa1"]]
kappa2 <- params[["kappa2"]]
kappa3 <- params[["kappa3"]]

gamma_i <- params[["gamma_i"]]
gamma_a <- params[["gamma_a"]]

mu   <- params[["mu"]]
zeta <- params[["zeta"]]

bracket_term <- (mu / gamma_i) + ((1 - mu) * zeta / gamma_a)

# ---- time-varying vaccination-cohort proportions ----
baseRc <- (mp_trajectory_sd(calibrator,conf.int=TRUE, back_transform=TRUE) 
	|>filter(matrix %in% c("p1", "p2", "p3","beta_thing")
		, time >= offset0
		) 
  	|> select(time, matrix, value, conf.low,conf.high) 
  	|> pivot_longer(-c(time,matrix),names_to = "est", values_to = "value")
  	|> transmute(NULL
  		, time 
		, vartype = paste0(matrix,"_",est)
		, value
	)
  	|> pivot_wider(names_from = vartype, values_from = value)
  	|> transmute(NULL
		, date = anchor_start + time -offset0  
	 	, beta = beta_thing_value
		, beta_lwr = beta_thing_conf.low
		, beta_upr = beta_thing_conf.high
		, p1 = p1_value
		, p2 = p2_value
		, p3 = p3_value
		, susceptibility = p1 * kappa1 + p2 * kappa2 + p3 * kappa3
#		, susceptibility = 0.15 * kappa1 + 0.85 * kappa2 + p3 * kappa3
		, psum = p1 + p2 + p3
		, mult_t = bracket_term * susceptibility
		, Rc = beta*mult_t
		, Rc_lwr = beta_lwr * mult_t
		, Rc_upr = beta_upr * mult_t
  		) 
	|> filter(date >= anchor_start, date <= anchor_end) 
)

print(dim(baseRc))

newRc <- function(beta_vec,scenario){
	Rcdat <- (baseRc
		|> select(date,beta,Rc,Rc_lwr,Rc_upr)
		|> mutate(beta_ratio = beta_vec/beta
			, new_beta = beta_ratio* beta
			, Rc = beta_ratio * Rc
			, Rc_lwr = beta_ratio * Rc_lwr
			, Rc_upr = beta_ratio * Rc_upr
			, scenario = scenario
		)
	)
}


noALS_Rc <- newRc(beta_noALS[-c(1:(offset0-1))],scenario="noALS")
ALS4only_Rc <- newRc(beta_ALS4only[-c(1:(offset0-1))],scenario="strict")



comboRc <- (baseRc
	|> select(date,Rc,Rc_lwr,Rc_upr)
	|> mutate(scenario = "base")
	|> bind_rows(noALS_Rc)
	|> bind_rows(ALS4only_Rc)
)

print(comboRc)

gg <- (ggplot(comboRc, aes(date,Rc,color=scenario,fill=scenario))
	+ geom_line()
	+ geom_ribbon(aes(ymin=Rc_lwr, ymax=Rc_upr),alpha=0.2)
)

print(gg)

quit()

# ============================================================

beta_actual_df <- calibrator |>
  mp_trajectory_sd(conf.int = TRUE, back_transform = TRUE) |>
  filter(
    matrix == "beta_thing",
    time >= offset0
  ) |>
  arrange(time) |>
  mutate(
    date = seq.Date(
      from = anchor_start,
      by = "day",
      length.out = n()
    )
  ) |>
  filter(date >= anchor_start, date <= anchor_end) |>
  transmute(date, beta = value, scenario = "Actual NPIs")

beta_vec_to_df <- function(beta_vec, scenario_label) {
  tibble(beta = as.numeric(beta_vec)) |>
    mutate(time = row_number()) |>
    filter(time >= offset0) |>
    mutate(
      date = seq.Date(
        from = anchor_start,
        by = "day",
        length.out = n()
      )
    ) |>
    filter(date >= anchor_start, date <= anchor_end) |>
    transmute(date, beta, scenario = scenario_label)
}

beta_noALS_df    <- beta_vec_to_df(beta_noALS, "No NPIs")
beta_ALS4only_df <- beta_vec_to_df(beta_ALS4only, "Strict NPIs")

beta_all_df <- bind_rows(
  beta_actual_df,
  beta_noALS_df,
  beta_ALS4only_df
) |>
  mutate(
    scenario = factor(
      scenario,
      levels = c("Actual NPIs", "No NPIs", "Strict NPIs")
    )
  )

p_beta <- ggplot() +
  geom_line(
    data = beta_all_df |> filter(scenario != "Actual NPIs"),
    aes(x = date, y = beta, color = scenario),
    linewidth = 2.5
  ) +
  geom_line(
    data = beta_all_df |> filter(scenario == "Actual NPIs"),
    aes(x = date, y = beta, color = scenario),
    linewidth = 2.5
  ) +
  geom_vline(
    xintercept = cut_date,
    linetype = "dashed",
    linewidth = 2
  ) +
  scale_color_manual(
    name = NULL,
    breaks = c("Actual NPIs", "No NPIs", "Strict NPIs"),
    values = c(
      "Actual NPIs" = "black",
      "No NPIs" = "red",
      "Strict NPIs" = "#00C5CD"
    )
  ) +
  scale_x_date(
    limits = c(anchor_start, anchor_end),
    date_breaks = "2 weeks",
    date_labels = "%d %b %Y",
    expand = c(0.01, 0)
  ) +
  labs(
    title = "Transmission rate",
    y = expression(beta(t)),
    x = NULL,
    tag = "B"
  ) +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 30),
    axis.title.y = element_text(size = 30, color = "black"),
    legend.position = c(0.85, 0.45),
    legend.justification = c(0.5, 0.5),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 25),
    legend.background = element_blank(),
    plot.title = element_text(size = 30, hjust = 0.5),
    plot.background = element_blank(),
    plot.tag = element_text(size = 35, face = "bold"),
    plot.tag.position = c(0.01, 0.98)
  )

# ============================================================
# panel C: weighted susceptibility
# ============================================================

p_susc <- ggplot(weighted_susc_df, aes(x = date, y = weighted_susc)) +
  geom_line(
    linewidth = 2.5,
    colour = "#1B9E77"
  ) +
  geom_vline(
    xintercept = cut_date,
    linetype = "dashed",
    linewidth = 2
  ) +
  scale_x_date(
    limits = c(anchor_start, anchor_end),
    date_breaks = "2 weeks",
    date_labels = "%d %b %Y",
    expand = c(0.01, 0)
  ) +
  labs(
    title = "Weighted susceptibility",
    y = "Weighted susceptibility",
    x = NULL,
    tag = "C"
  ) +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 30),
    axis.title.y = element_text(size = 30, color = "black"),
    plot.title = element_text(size = 30, hjust = 0.5),
    plot.background = element_blank(),
    plot.tag = element_text(size = 35, face = "bold"),
    plot.tag.position = c(0.01, 0.98)
  )

# ============================================================
# save separate Rc plot
# ============================================================

dir.create("../figures", showWarnings = FALSE, recursive = TRUE)

png(
  "../figures/control_reprod_numb.png",
  width = 5000,
  height = 1600,
  res = 300,
  bg = "white",
  type = "cairo"
)

print(p_rc)

dev.off()

# ============================================================
# combine and save all panels
# ============================================================

combined_plot <- p_rc + p_beta + p_susc +
  plot_layout(widths = c(1.5, 1, 1))

png(
  "../figures/reprod_numb_with_beta_susceptibility.png",
  width = 9500,
  height = 2000,
  res = 300,
  bg = "white",
  type = "cairo"
)

print(combined_plot)

dev.off()

# optional diagnostic check
print(
  weighted_susc_df |>
    summarise(
      min_psum = min(psum, na.rm = TRUE),
      max_psum = max(psum, na.rm = TRUE),
      min_weighted_susc = min(weighted_susc, na.rm = TRUE),
      max_weighted_susc = max(weighted_susc, na.rm = TRUE)
    )
)

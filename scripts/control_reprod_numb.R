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
base_beta      <- rdsRead("base_beta")

upper_plot_time <- 300

offset          <- 150
omicron_start   <- as.Date("2021-12-15")
burn_in_start   <- omicron_start - offset
last_date       <- as.Date("2022-05-22")

intervention_removed <- as.Date("2022-03-14")

time_omicron_start <- as.numeric(omicron_start - burn_in_start) + 1
time_last_date     <- as.numeric(last_date     - burn_in_start) + 1

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
                    , time >= time_omicron_start
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
                        , date = burn_in_start + (time - 1)
                        , beta = beta_thing_value
                        , beta_lwr = beta_thing_conf.low
                        , beta_upr = beta_thing_conf.high
                        , p1 = p1_value
                        , p2 = p2_value
                        , p3 = p3_value
                        #		, susceptibility = p1 * kappa1 + p2 * kappa2 + p3 * kappa3
                        , susceptibility = 0.15 * kappa1 + 0.85 * kappa2 + 0 * kappa3
                        , psum = p1 + p2 + p3
                        , mult_t = bracket_term * susceptibility
                        , Rc = beta*mult_t
                        , Rc_lwr = beta_lwr * mult_t
                        , Rc_upr = beta_upr * mult_t
           ) 
           |> filter(date >= omicron_start, date <= last_date) 
)


cutoff_time_check <- 232
cutoff_date_check <- burn_in_start + (cutoff_time_check - 1)
base_beta_dates <- burn_in_start + (seq_along(base_beta) - 1)
base_beta_df <- data.frame(date = base_beta_dates, base_beta_value = base_beta)
verify_alignment <- (baseRc
                     |> select(date, beta) 
                     |> inner_join(base_beta_df, by = "date")
                     |> mutate(diff = beta - base_beta_value)
                     |> filter(date < cutoff_date_check)
)

new_constant_value <- 0.3675
cutoff_time <- 232
cutoff_date <- burn_in_start + (cutoff_time - 1)
past_cutoff <- baseRc$date >= cutoff_date

baseRc$beta[past_cutoff]     <- new_constant_value
baseRc$beta_lwr[past_cutoff] <- new_constant_value
baseRc$beta_upr[past_cutoff] <- new_constant_value

baseRc$Rc[past_cutoff]     <- baseRc$beta[past_cutoff]     * baseRc$mult_t[past_cutoff]
baseRc$Rc_lwr[past_cutoff] <- baseRc$beta_lwr[past_cutoff] * baseRc$mult_t[past_cutoff]
baseRc$Rc_upr[past_cutoff] <- baseRc$beta_upr[past_cutoff] * baseRc$mult_t[past_cutoff]


beta_vec_to_df <- function(beta_vec) {
  tibble(
    time = seq_along(beta_vec),
    date = burn_in_start + (time - 1),
    beta_scenario = beta_vec
  )
}

newRc <- function(beta_vec,scenario){
  beta_df <- beta_vec_to_df(beta_vec)
  Rcdat <- (baseRc
            |> select(date,beta,Rc,Rc_lwr,Rc_upr)
            |> inner_join(beta_df |> select(date, beta_scenario), by = "date")
            |> mutate(beta_ratio = beta_scenario/beta
                      , new_beta = beta_ratio* beta
                      , Rc = beta_ratio * Rc
                      , Rc_lwr = beta_ratio * Rc_lwr
                      , Rc_upr = beta_ratio * Rc_upr
                      , scenario = scenario
            )
            |> select(-beta_scenario)
  )
  if (nrow(Rcdat) != nrow(baseRc)) {
    warning(sprintf(
      "newRc('%s'): joined %d of %d baseRc rows -- %d date(s) in baseRc had no matching entry in beta_vec.",
      scenario, nrow(Rcdat), nrow(baseRc), nrow(baseRc) - nrow(Rcdat)
    ))
  }
  Rcdat
}

noALS_Rc <- newRc(beta_noALS,scenario="noALS")
ALS4only_Rc <- newRc(beta_ALS4only,scenario="strict")

comboRc <- (baseRc
            |> select(date,Rc,Rc_lwr,Rc_upr)
            |> mutate(scenario = "base")
            |> bind_rows(noALS_Rc)
            |> bind_rows(ALS4only_Rc)
)

print(comboRc)

scenario_levels <- c("base", "noALS", "strict")
scenario_labels <- c(
  "base"   = "Actual NPIs",
  "noALS"  = "No NPIs",
  "strict" = "Strict NPIs"
)
scenario_colors <- c(
  "base"   = "black",
  "noALS"  = "red",
  "strict" = "#00C5CD"
)

gg <- (ggplot(comboRc, aes(date, Rc, color = scenario, fill = scenario))
       + geom_ribbon(aes(ymin = Rc_lwr, ymax = Rc_upr), alpha = 0.2, colour = NA)
       + geom_line(linewidth = 2)
       + geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1.2)
       + geom_vline(
         xintercept = intervention_removed,
         linetype = "dashed",
         linewidth = 3
       )
       + scale_color_manual(
         name   = NULL,
         breaks = scenario_levels,
         values = scenario_colors,
         labels = scenario_labels
       )
       + scale_fill_manual(
         name   = NULL,
         breaks = scenario_levels,
         values = scenario_colors,
         labels = scenario_labels
       )
       + scale_x_date(
         limits = c(omicron_start, last_date),
         date_breaks = "2 weeks",
         date_labels = "%d %b %Y",
         expand = c(0.01, 0)
       )
       + labs(
         title = "Reproduction number for counterfactual scenarios (timevarying prop)",
         y = "Reproduction number",
         x = NULL,
         tag = "C"
       )
       + theme_clean()
       + theme(
         axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
         axis.text.y = element_text(size = 30),
         axis.title.y = element_text(size = 30, color = "black"),
         legend.position = c(0.85, 0.65),
         legend.justification = c(0.5, 0.5),
         legend.direction = "vertical",
         legend.box = "vertical",
         legend.title = element_text(size = 25),
         legend.text = element_text(size = 25),
         legend.background = element_blank(),
         plot.title = element_text(size = 30, hjust = 0.5, face = "plain"),
         plot.background = element_blank(),
         plot.tag = element_text(size = 35, face = "bold"),
         plot.tag.position = c(0.01, 0.98)
       )
)

png("../figures/reprod_numb.png", width = 5000, height = 2000, res = 300, bg = "white",type = "cairo")
print(gg)
dev.off()

baseRc_with_scenarios <- (baseRc
                          |> inner_join(beta_vec_to_df(beta_noALS)    |> select(date, beta_noALS_value = beta_scenario), by = "date")
                          |> inner_join(beta_vec_to_df(beta_ALS4only) |> select(date, beta_ALS4only_value = beta_scenario), by = "date")
)

gg1 <- (ggplot(baseRc_with_scenarios, aes(date, beta))
        + geom_ribbon(
          aes(ymin = beta_lwr, ymax = beta_upr, fill = "Actual NPIs"),
          alpha = 0.2, colour = NA
        )
        + geom_line(aes(colour = "Actual NPIs"), linewidth = 2)
        + geom_line(
          aes(y = beta_noALS_value, colour = "No NPIs"),
          linewidth = 2
        )
        + geom_line(
          aes(y = beta_ALS4only_value, colour = "Strict NPIs"),
          linewidth = 2
        )
        + geom_vline(
          xintercept = intervention_removed,
          linetype = "dashed",
          linewidth = 3
        )
        + scale_color_manual(
          name = NULL,
          breaks = c("Actual NPIs", "No NPIs", "Strict NPIs"),
          values = c(
            "Actual NPIs" = "black",
            "No NPIs" = "red",
            "Strict NPIs" = "#00C5CD"
          )
        )
        + scale_fill_manual(
          name = NULL,
          breaks = c("Actual NPIs", "No NPIs", "Strict NPIs"),
          values = c(
            "Actual NPIs" = "black",
            "No NPIs" = "red",
            "Strict NPIs" = "#00C5CD"
          ),
          guide = "none"
        )
        + scale_x_date(
          limits = c(omicron_start, last_date),
          date_breaks = "2 weeks",
          date_labels = "%d %b %Y",
          expand = c(0.01, 0)
        )
        + labs(
          title = "Transmission rate for counterfactual scenarios",
          y = "Transmission rate",
          x = NULL,
          tag = "A"
        )
        + theme_clean()
        + theme(
          axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 30),
          axis.title.y = element_text(size = 30, color = "black"),
          legend.position = c(0.85, 0.5),
          legend.justification = c(0.5, 0.5),
          legend.direction = "vertical",
          legend.box = "vertical",
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 25),
          legend.background = element_blank(),
          plot.title = element_text(size = 30, hjust = 0.5, face = "plain"),
          plot.background = element_blank(),
          plot.tag = element_text(size = 35, face = "bold"),
          plot.tag.position = c(0.01, 0.98)
        )
)
png("../figures/transmission_rate.png", width = 5000, height = 2000, res = 300, bg = "white", type = "cairo")
print(gg1)
dev.off()

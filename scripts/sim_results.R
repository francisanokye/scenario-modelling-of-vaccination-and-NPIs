suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggthemes)
  library(shellpipes)
  library(macpan2)
  library(scales)
  library(lubridate)
})

set.seed(2025)
options(macpan2_log_dir = ".")
options(dplyr.summarise.inform = FALSE)
options(scipen = 9999)

loadEnvironments()

all_sim_results <- vector("list", nrow(sim_grid))

for (i in seq_len(nrow(sim_grid))) {
  
  policy_name <- sim_grid$policy_name[i]
  vax_name    <- sim_grid$vax_name[i]
  
  spec_i <- update_vax_spec_beta(vax_scenario_specs[[vax_name]], policy_list[[policy_name]])
  
  all_sim_results[[i]] <- (spec_i
    |> mp_simulator(time_steps = sim_horizon, outputs = model_outputs)
    |> mp_trajectory()
    |> mutate(date = burn_in_start + (time - 1), policy = policy_name, vax_scenario = vax_name)
  )
  
  # diagnostic <- (spec_i
  #   |> mp_simulator(sim_horizon, outputs = diagnostic_outputs)
  #   |> mp_trajectory_sd(conf.int = TRUE, back_transform = TRUE)
  #   |> dplyr::mutate(date = burn_in_start + (time - 1))
  #   |> ggplot()
  #   + geom_rect(data = data.frame(xmin = burn_in_start, xmax = omicron_start),
  #              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
  #              inherit.aes = FALSE, fill = "grey30", alpha = 0.25)
  #   + geom_vline(xintercept = burn_in_start + (231 - 1), color = "black", linetype = 4, linewidth = 1.0)
  #   + geom_vline(xintercept = burn_in_start + (300 - 1), colour = "purple", linetype = 4, linewidth = 1.5)
  #   + aes(date, value)
  #   + geom_line()
  #   + labs(title = paste0("Diagnostic compartment trajectories - ", policy_name, " | ", vax_name))
  #   + scale_x_date(expand = c(0, 0), date_breaks = "month", date_labels = "%b %d %Y")
  #   + facet_wrap(~matrix, scales = "free")
  #   + theme_clean()
  #   + theme(
  #    axis.text.x       = element_text(size = 12,  hjust = 1, angle = 45, margin = margin(t = 8)),
  #    axis.text.y       = element_text(size = 18),
  #    axis.title.x      = element_text(size = 18),
  #    axis.title.y      = element_text(size = 18),
  #    plot.title        = element_text(size = 24, hjust = 0.5, face = "plain"),
  #    strip.text        = element_text(size = 22, face = "plain", hjust = 0.5),
  #    strip.background  = element_blank(),
  #    legend.position   = "top",
  #    legend.direction  = "horizontal",
  #    legend.box        = "horizontal",
  #    legend.spacing.x  = unit(1.2, "cm"),
  #    legend.spacing.y  = unit(0.5, "cm"),
  #    legend.background = element_blank(),
  #    legend.text       = element_text(size = 22, face = "plain"),
  #    panel.spacing     = unit(1.5, "cm"),
  #    plot.background   = element_blank()
  #   )
  # )
  # 
  # diagnostic_tag <- paste0(panel_tag(policy_name), "_", panel_tag(vax_name))
  # 
  # save_png(
  #   file.path(fig_dir, paste0("diagnostic_", diagnostic_tag, ".png")),
  #   diagnostic,
  #   width  = 10000,
  #   height = 5000
  # )
}

results_df <- (bind_rows(all_sim_results) 
  |> mutate(date = as.Date(date)
    , value = as.numeric(value)
    , matrix = as.character(matrix)
    , policy = factor(policy
    , levels = policy_order)
    , vax_scenario = factor(as.character(vax_scenario)
    , levels = vax_order)
  )
)

required_outputs <- c(
  "cum_inc_leq2dose",
  "cum_inc_booster",
  "cum_hosp_total",
  "cum_death_total"
)

saveVars(results_df)
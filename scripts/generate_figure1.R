suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggthemes)
  library(shellpipes)
  library(macpan2)
  library(scales)
  library(grid)
  library(patchwork)
  library(lubridate)
  library(forcats)
})

set.seed(2025)

options(macpan2_log_dir = ".")
options(dplyr.summarise.inform = FALSE)
options(scipen = 9999)

loadEnvironments()

offset       <- 150L
sim_horizon  <- 530L
NL_POP       <- 510550

CA_PRIMARY_ONLY_PROP_START <- 0.750
new_V2_init <- NA_real_

omicron_start  <- as.Date("2021-12-15")
burn_in_start  <- omicron_start - offset

intervention_removed <- as.Date("2022-03-14")
# analysis_end   <- burn_in_start + sim_horizon - 1L

booster_start    <- as.Date("2021-11-08")
early_start_date <- as.Date("2021-10-04")

fig_dir <- "../figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

save_png <- function(path, plot_obj, width, height, res = 300, bg = "white") {
  png(path, width = width, height = height, res = res, bg = bg, type = "cairo")
  print(plot_obj)
  dev.off()
  invisible(path)
}

max_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}

panel_tag <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

calibrator <- rdsRead("calibrate")

padded_beta_noals <- rdsRead("policy_neutral")
padded_beta_base  <- rdsRead("base_beta")
padded_beta_als4  <- rdsRead("beta_als4")

scenario_cols <- c(
  "National Rate"   = "springgreen3",
  "Slower Boosters" = "magenta",
  "Actual Boosters" = "#4D4D4D",
  "Early Boosters"  = "blue",
  "Faster Boosters" = "orange"
)

scenario_dose_cols <- list(
  "Actual Boosters" = c("≤2 Doses" = "grey70",   "Booster" = "#4D4D4D"),
  "Slower Boosters" = c("≤2 Doses" = "#E8B6D2",  "Booster" = "magenta"),
  "Faster Boosters" = c("≤2 Doses" = "#FDBB84",  "Booster" = "orange"),
  "Early Boosters"  = c("≤2 Doses" = "#9ECAE1",  "Booster" = "blue")
)

booster_start_dates <- c(
  "Actual Boosters" = booster_start,
  "Early Boosters"  = early_start_date,
  "Faster Boosters" = booster_start,
  "Slower Boosters" = booster_start,
  "National Rate"   = booster_start
)

############# actual boosters ##############################

spec_actual <- mp_optimized_spec(calibrator, "modified")

############# set slower booster ##############################
slower_booster_changepoints <- as.integer(vax_slower_rate$day)
slower_booster_values       <- as.numeric(vax_slower_rate$daily_rate)

if (!any(slower_booster_changepoints == 0)) {
  slower_booster_changepoints <- c(0L, slower_booster_changepoints)
  slower_booster_values       <- c(0,  slower_booster_values)
}

spec_slower <- mp_optimized_spec(calibrator, "modified")
spec_slower$integers$booster_vac_changepoints <- as.integer(slower_booster_changepoints)
spec_slower$default$booster_vac_values        <- as.numeric(slower_booster_values)

############# set faster booster ##############################
faster_booster_changepoints <- as.integer(vax_faster_rate$day)
faster_booster_values       <- as.numeric(vax_faster_rate$daily_rate)

if (!any(faster_booster_changepoints == 0)) {
  faster_booster_changepoints <- c(0L, faster_booster_changepoints)
  faster_booster_values       <- c(0,  faster_booster_values)
}

spec_faster <- mp_optimized_spec(calibrator, "modified")
spec_faster$integers$booster_vac_changepoints <- as.integer(faster_booster_changepoints)
spec_faster$default$booster_vac_values        <- as.numeric(faster_booster_values)

############# set early booster ##############################
early_booster_changepoints <- as.integer(vax_early_rate$day)
early_booster_values       <- as.numeric(vax_early_rate$daily_rate)

if (!any(early_booster_changepoints == 0)) {
  early_booster_changepoints <- c(0L, early_booster_changepoints)
  early_booster_values       <- c(0,  early_booster_values)
}

spec_early <- mp_optimized_spec(calibrator, "modified")
spec_early$integers$booster_vac_changepoints <- as.integer(early_booster_changepoints)
spec_early$default$booster_vac_values        <- as.numeric(early_booster_values)

############# set National-NL booster ##############################
national_booster_changepoints <- as.integer(can_boosterframe$day)
national_booster_values       <- as.numeric(can_boosterframe$daily_rate)

if (!any(national_booster_changepoints == 0)) {
  national_booster_changepoints <- c(0L, national_booster_changepoints)
  national_booster_values       <- c(0,  national_booster_values)
}

spec_natl <- mp_optimized_spec(calibrator, "modified")
spec_natl$integers$booster_vac_changepoints <- as.integer(national_booster_changepoints)
spec_natl$default$booster_vac_values        <- as.numeric(national_booster_values)

# insert new natioanl initial vaccinated
spec_natl <- (spec_natl |>
  mp_tmb_insert(default = list(V2_init = new_V2_init))
  )

############# update the betas ##############################

update_vax_spec_beta <- function(spec, new_beta) {
  stopifnot(length(new_beta) == spec$default$time_steps)
  
  spec |>
    mp_tmb_delete(phase = "before", at = 26L) |>
    mp_tmb_insert(default = list(timeseries_beta = new_beta))
}

vax_scenario_specs <- list(
  "National Rate"   = spec_natl,
  "Slower Boosters" = spec_slower,
  "Actual Boosters" = spec_actual,
  "Early Boosters"  = spec_early,
  "Faster Boosters" = spec_faster
)

policy_list <- list(
  "Actual ALS"  = padded_beta_base,
  "ALS-4 Only"  = padded_beta_als4,
  "No-ALS Only" = padded_beta_noals
)

policy_order <- c("No-ALS Only", "Actual ALS", "ALS-4 Only")

vax_order <- c(
  "National Rate",
  "Slower Boosters",
  "Actual Boosters",
  "Early Boosters",
  "Faster Boosters"
)

model_outputs <- c(
  "incS", "incv2", "incv3", "inc",
  "inc_leq2dose", "inc_booster",
  "cum_inc_leq2dose", "cum_inc_booster", "cum_inc_total",
  "hosp_leq2dose", "hosp_booster", "hosp_total",
  "death_leq2dose", "death_booster", "death_total",
  "cum_hosp_leq2dose", "cum_hosp_booster", "cum_hosp_total",
  "cum_death_leq2dose", "cum_death_booster", "cum_death_total",
  "cum_booster"
)

diagnostic_outputs <- c(
  "S1", "V2", "V3", "E1", "E2", "E3", "A1", "A2",
  "A3", "I1", "I2", "I3", "R1", "R2", "R3", "inc",
  "inc_leq2dose", "inc_booster", "p1", "p2", "p3",
  "cum_booster", "cum_double_vac"
)

sim_grid <- expand.grid(
  policy_name = names(policy_list),
  vax_name    = names(vax_scenario_specs),
  stringsAsFactors = FALSE
)

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
  
  diagnostic <- (spec_i
                 |> mp_simulator(sim_horizon,outputs = diagnostic_outputs)
                 |> mp_trajectory_sd(conf.int = TRUE, back_transform = TRUE)
                 |> dplyr::mutate(date = burn_in_start + (time - 1))
                 |> ggplot() 
                 + geom_rect(data = data.frame(xmin = burn_in_start, xmax = omicron_start),
                             aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "grey30", alpha = 0.25) 
                 + geom_vline(xintercept = burn_in_start + (231 - 1), color = "black", linetype = 4, linewidth = 1.0) 
                 + geom_vline(xintercept = burn_in_start + (300 - 1), colour = "purple", linetype = 4, linewidth = 1.5) 
                 + aes(date, value) 
                 + geom_line() 
                 + labs(title = paste0("Diagnostic compartment trajectories - ", policy_name, " | ", vax_name)) 
                 + scale_x_date(expand = c(0, 0), date_breaks = "month", date_labels = "%b %d %Y") 
                 + facet_wrap(~matrix, scales = "free") 
                 + theme_clean() 
                 + theme(
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
  
  diagnostic_tag <- paste0(panel_tag(policy_name), "_", panel_tag(vax_name))
  
  save_png(
    file.path(fig_dir, paste0("diagnostic_", diagnostic_tag, ".png")),
    diagnostic,
    width  = 10000,
    height = 5000
  )
}

results_df <- bind_rows(all_sim_results) |>
  mutate(
    date         = as.Date(date),
    value        = as.numeric(value),
    matrix       = as.character(matrix),
    policy       = factor(policy, levels = policy_order),
    vax_scenario = factor(as.character(vax_scenario), levels = vax_order)
  )

saveVars(results_df)

required_outputs <- c(
  "cum_inc_leq2dose",
  "cum_inc_booster",
  "cum_hosp_total",
  "cum_death_total"
)

missing_outputs <- setdiff(required_outputs, unique(results_df$matrix))

if (length(missing_outputs) > 0) {
  stop(
    "Missing model outputs: ",
    paste(missing_outputs, collapse = ", "),
    "\nRe-run upstream after updating flows and params."
  )
}

incidence_df  <- results_df |> filter(date >= omicron_start) #, date <= analysis_end)
cumulative_df <- results_df |> filter(date >= omicron_start) #, date <= analysis_end)

ab_ymax <- incidence_df |>
  filter(
    vax_scenario == "Actual Boosters",
    policy %in% c("No-ALS Only", "ALS-4 Only"),
    matrix %in% c("inc_leq2dose", "inc_booster")
  ) |>
  group_by(policy, date) |>
  summarise(total_inc = sum(value, na.rm = TRUE), .groups = "drop") |>
  summarise(ymax = max_or_na(total_inc)) |>
  pull(ymax)

gh_ymax <- incidence_df |>
  filter(
    policy == "ALS-4 Only",
    vax_scenario %in% c("Early Boosters", "Actual Boosters"),
    matrix %in% c("inc_leq2dose", "inc_booster")
  ) |>
  group_by(vax_scenario, date) |>
  summarise(total_inc = sum(value, na.rm = TRUE), .groups = "drop") |>
  summarise(ymax = max_or_na(total_inc)) |>
  pull(ymax)

de_ymax <- incidence_df |>
  filter(
    policy == "ALS-4 Only",
    vax_scenario %in% c("Slower Boosters", "Faster Boosters"),
    matrix %in% c("inc_leq2dose", "inc_booster")
  ) |>
  group_by(vax_scenario, date) |>
  summarise(total_inc = sum(value, na.rm = TRUE), .groups = "drop") |>
  summarise(ymax = max_or_na(total_inc)) |>
  pull(ymax)


incidence_ymax <- max(ab_ymax, de_ymax, gh_ymax, na.rm = TRUE)

get_panel_cols <- function(scenario_name) {
  cols <- scenario_dose_cols[[as.character(scenario_name)]]
  if (is.null(cols)) stop("No dose colours defined for scenario: ", scenario_name)
  cols
}

make_shaded_incidence_plot <- function(
    df,
    scenario_name,
    policy_name,
    panel_tag,
    panel_title,
    legend_position   = c(0.50, 0.85),
    y_limit           = NULL,
    boundary_linetype = "solid",
    show_legend       = TRUE
) {
  
  panel_cols <- get_panel_cols(scenario_name)
  
  plot_df <- df |>
    filter(
      policy       == policy_name,
      vax_scenario == scenario_name,
      matrix       %in% c("inc_leq2dose", "inc_booster")
    ) |>
    select(date, matrix, value) |>
    pivot_wider(names_from = matrix, values_from = value, values_fill = 0) |>
    mutate(total_inc = inc_leq2dose + inc_booster)
  
  y_upper <- if (is.null(y_limit)) max_or_na(plot_df$total_inc) else y_limit
  
  ggplot(plot_df, aes(x = date)) +
    geom_area(aes(y = total_inc,  fill = "≤2 Doses")) +
    geom_area(aes(y = inc_booster, fill = "Booster")) +
    geom_line(
      aes(y = total_inc),
      colour = "black",
      linewidth = 2,
      linetype = boundary_linetype
    ) +
    geom_line(
      aes(y = inc_booster),
      colour = "black",
      linewidth = 2,
      linetype = boundary_linetype
    ) +
    geom_vline(xintercept = intervention_removed, linetype = "dashed", linewidth = 2) +
    scale_fill_manual(
      values = panel_cols,
      breaks = c("Booster", "≤2 Doses")
    ) +
    scale_x_date(
      # limits       = c(omicron_start, analysis_end),
      date_breaks  = "1 month",
      date_labels  = "%d %b %Y",
      expand       = c(0.01, 0)
    ) +
    scale_y_continuous(
      labels = comma,
      expand = expansion(mult = c(0, 0.04))
    ) +
    coord_cartesian(ylim = c(0, y_upper * 1.06), clip = "off") +
    labs(
      title = panel_title,
      x = NULL,
      y = "Daily incidence",
      fill = NULL,
      tag = panel_tag
    ) +
    theme_clean() +
    theme(
      axis.text.x = element_text(size = 35, hjust = 1, angle = 45),
      axis.text.y = element_text(size = 35),
      axis.title.y = element_text(size = 35),
      legend.position = if (show_legend) legend_position else "none",
      legend.justification = "center",
      legend.direction = "vertical",
      legend.text = element_text(size = 30),
      legend.background = element_blank(),
      plot.title = element_text(size = 45, hjust = 0.5),
      panel.grid.major.x = element_blank(),
      plot.tag = element_text(size = 45, face = "bold"),
      plot.background = element_blank(),
      plot.margin = margin(12, 35, 12, 12)
    )
}

make_cumulative_policy_plot <- function(df, scenario_name, panel_tag) {
  
  panel_cols <- get_panel_cols(scenario_name)
  
  plot_df <- df |>
    filter(
      vax_scenario == scenario_name,
      policy %in% c("No-ALS Only", "ALS-4 Only"),
      matrix %in% c("cum_inc_leq2dose", "cum_inc_booster")
    ) |>
    mutate(
      dose_group = recode(
        matrix,
        "cum_inc_leq2dose" = "≤2 Doses",
        "cum_inc_booster"  = "Booster"
      ),
      policy_label = recode(
        as.character(policy),
        "No-ALS Only" = "No NPIs",
        "ALS-4 Only"  = "Strict NPIs"
      ),
      dose_group = factor(dose_group, levels = c("≤2 Doses", "Booster")),
      policy_label = factor(policy_label, levels = c("No NPIs", "Strict NPIs"))
    )
  
  ggplot(
    plot_df,
    aes(
      x = date,
      y = value,
      colour = dose_group,
      linetype = policy_label,
      group = interaction(dose_group, policy_label)
    )
  ) +
    geom_line(linewidth = 3) +
    scale_colour_manual(
      values = panel_cols,
      breaks = c("≤2 Doses", "Booster"),
      name = NULL
    ) +
    scale_linetype_manual(
      values = c("No NPIs" = "solid", "Strict NPIs" = "twodash"),
      name = NULL
    ) +
    scale_x_date(
      #limits = c(omicron_start, analysis_end),
      date_breaks = "1 month",
      date_labels = "%d %b %Y",
      expand = c(0.01, 0.08)
    ) +
    scale_y_continuous(
      labels = comma,
      expand = expansion(mult = c(0, 0.12))
    ) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = "Cumulative infections", tag = panel_tag) +
    guides(
      colour = guide_legend(
        order = 1,
        override.aes = list(linewidth = 4, linetype = "solid")
      ),
      linetype = guide_legend(
        order = 2,
        override.aes = list(linewidth = 4, colour = "black")
      )
    ) +
    theme_clean() +
    theme(
      axis.text.x = element_text(size = 35, hjust = 1, angle = 45),
      axis.text.y = element_text(size = 35),
      axis.title.y = element_text(size = 35),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.text = element_text(size = 30),
      legend.key.width = unit(2, "cm"),
      legend.background = element_blank(),
      legend.spacing.x = unit(0.20, "cm"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
      panel.grid.major.x = element_blank(),
      plot.tag = element_text(size = 45, face = "bold"),
      plot.background = element_blank(),
      plot.margin = margin(12, 55, 12, 12)
    )
}

make_cumulative_speed_plot <- function(df, panel_tag) {
  
  plot_df <- df |>
    filter(
      policy == "ALS-4 Only",
      vax_scenario %in% c("Slower Boosters", "Faster Boosters"),
      matrix %in% c("cum_inc_leq2dose", "cum_inc_booster")
    ) |>
    mutate(
      dose_group = recode(
        matrix,
        "cum_inc_leq2dose" = "≤2 Doses",
        "cum_inc_booster"  = "Booster"
      ),
      vax_scenario = factor(
        as.character(vax_scenario),
        levels = c("Slower Boosters", "Faster Boosters")
      ),
      dose_group = factor(dose_group, levels = c("≤2 Doses", "Booster"))
    )
  
  ggplot(
    plot_df,
    aes(
      x = date,
      y = value,
      colour = vax_scenario,
      linetype = dose_group,
      group = interaction(vax_scenario, dose_group)
    )
  ) +
    geom_line(linewidth = 3) +
    scale_colour_manual(
      values = scenario_cols[c("Slower Boosters", "Faster Boosters")],
      name = NULL
    ) +
    scale_linetype_manual(
      values = c("≤2 Doses" = "solid", "Booster" = "twodash"),
      name = NULL
    ) +
    scale_x_date(
      #limits = c(omicron_start, analysis_end),
      date_breaks = "1 month",
      date_labels = "%d %b %Y",
      expand = c(0.01, 0.08)
    ) +
    scale_y_continuous(
      labels = comma,
      expand = expansion(mult = c(0, 0.12))
    ) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = "Cumulative infections", tag = panel_tag) +
    guides(
      colour = guide_legend(
        order = 1,
        override.aes = list(linewidth = 4, linetype = "solid")
      ),
      linetype = guide_legend(
        order = 2,
        override.aes = list(linewidth = 4, colour = "black")
      )
    ) +
    theme_clean() +
    theme(
      axis.text.x = element_text(size = 35, hjust = 1, angle = 45),
      axis.text.y = element_text(size = 35),
      axis.title.y = element_text(size = 35),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.text = element_text(size = 30),
      legend.key.width = unit(2, "cm"),
      legend.background = element_blank(),
      legend.spacing.x = unit(0.20, "cm"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
      panel.grid.major.x = element_blank(),
      plot.tag = element_text(size = 45, face = "bold"),
      plot.background = element_blank(),
      plot.margin = margin(12, 55, 12, 12)
    )
}

p_a <- make_shaded_incidence_plot(
  incidence_df, "Actual Boosters", "No-ALS Only",
  "(A)", "No NPIs", c(0.60, 0.85), incidence_ymax, "solid"
)

p_b <- make_shaded_incidence_plot(
  incidence_df, "Actual Boosters", "ALS-4 Only",
  "(B)", "Strict NPIs", c(0.15, 0.85), incidence_ymax, "twodash"
)

p_c <- make_cumulative_policy_plot(cumulative_df, "Actual Boosters", "(C)")

p_d <- make_shaded_incidence_plot(
  incidence_df, "Slower Boosters", "ALS-4 Only",
  "(D)", "Slower Boosters", c(0.15, 0.85), incidence_ymax, "twodash"
)

p_e <- make_shaded_incidence_plot(
  incidence_df, "Faster Boosters", "ALS-4 Only",
  "(E)", "Faster Boosters", c(0.15, 0.85), incidence_ymax, "twodash"
)

p_f <- make_cumulative_speed_plot(cumulative_df, "(F)")

p_g <- make_shaded_incidence_plot(
  incidence_df, "Early Boosters", "ALS-4 Only",
  "(G)", "Early Boosters", c(0.15, 0.85), incidence_ymax, "twodash"
)

p_h <- make_shaded_incidence_plot(
  incidence_df, "Actual Boosters", "ALS-4 Only",
  "(H)", "Actual Boosters", c(0.15, 0.85), incidence_ymax, "twodash"
)

plot_i_df <- incidence_df |>
  filter(
    policy == "ALS-4 Only",
    vax_scenario %in% c("Early Boosters", "Actual Boosters"),
    matrix %in% c("inc_leq2dose", "inc_booster")
  ) |>
  group_by(vax_scenario, matrix) |>
  summarise(
    peak_daily_incidence = max_or_na(value),
    .groups = "drop"
  ) |>
  mutate(
    dose_group = recode(
      matrix,
      "inc_leq2dose" = "≤2 Doses",
      "inc_booster"  = "Booster"
    ),
    scenario_label = factor(
      as.character(vax_scenario),
      levels = c("Early Boosters", "Actual Boosters")
    ),
    dose_group = factor(dose_group, levels = c("≤2 Doses", "Booster"))
  )

p_i <- ggplot(
  plot_i_df,
  aes(x = dose_group, y = peak_daily_incidence, fill = scenario_label)
) +
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.68,
    colour = "#4D4D4D",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = scenario_cols[c("Early Boosters", "Actual Boosters")]
  ) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(x = NULL, y = "Peak daily incidence", fill = NULL, tag = "(I)") +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 35, hjust = 0.5),
    axis.text.y = element_text(size = 35),
    axis.title.y = element_text(size = 35),
    legend.position = c(0.80, 0.90),
    legend.justification = "center",
    legend.direction = "vertical",
    legend.text = element_text(size = 35),
    legend.background = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.tag = element_text(size = 45, face = "bold"),
    plot.background = element_blank(),
    plot.margin = margin(12, 35, 12, 12)
  )

figure1_3x3 <- (
  (p_a | p_b | p_c) /
    (p_d | p_e | p_f) /
    (p_g | p_h | p_i)
) +
  plot_layout(heights = c(1, 1, 1), widths = c(1, 1, 1)) &
  theme(plot.margin = margin(8, 35, 8, 8))

save_png(
  file.path(fig_dir, "figure2.png"),
  figure1_3x3,
  width  = 13500,
  height = 8000
)

severity_fill_cols <- scenario_cols

vax_short_labs <- c(
  "National Rate"   = "National\nRate",
  "Slower Boosters" = "Slower\nBoosters",
  "Actual Boosters" = "Actual\nBoosters",
  "Early Boosters"  = "Early\nBoosters",
  "Faster Boosters" = "Faster\nBoosters"
)

figure2_df <- bind_rows(
  results_df |>
    filter(
      date >= omicron_start,
      #date <= analysis_end,
      matrix == "cum_hosp_total"
    ) |>
    group_by(policy, vax_scenario) |>
    summarise(severe_outcome = max_or_na(value), .groups = "drop") |>
    mutate(outcome = "Hospitalizations"),
  
  results_df |>
    filter(
      date >= omicron_start,
      #date <= analysis_end,
      matrix == "cum_death_total"
    ) |>
    group_by(policy, vax_scenario) |>
    summarise(severe_outcome = max_or_na(value), .groups = "drop") |>
    mutate(outcome = "Deaths")
) |>
  filter(is.finite(severe_outcome)) |>
  mutate(
    policy_label = recode(
      as.character(policy),
      "No-ALS Only" = "No NPIs",
      "Actual ALS"  = "Actual NPIs",
      "ALS-4 Only"  = "Strict NPIs"
    ),
    policy_label = factor(
      policy_label,
      levels = c("No NPIs", "Actual NPIs", "Strict NPIs")
    ),
    outcome = factor(outcome, levels = c("Hospitalizations", "Deaths")),
    vax_scenario = factor(as.character(vax_scenario), levels = vax_order)
  )

x_order_df <- figure2_df |>
  filter(outcome == "Hospitalizations") |>
  group_by(policy_label) |>
  arrange(desc(severe_outcome), .by_group = TRUE) |>
  mutate(x_pos = row_number()) |>
  ungroup() |>
  select(policy_label, vax_scenario, x_pos)

figure2_df <- figure2_df |>
  left_join(x_order_df, by = c("policy_label", "vax_scenario")) |>
  mutate(
    x_lab = vax_short_labs[as.character(vax_scenario)],
    value_label = comma(round(severe_outcome, 0))
  ) |>
  filter(!is.na(x_pos))

figure2 <- ggplot(
  figure2_df,
  aes(x = x_pos, y = severe_outcome, fill = vax_scenario)
) +
  geom_col(colour = "black", width = 0.72, linewidth = 0.25) +
  geom_text(aes(label = value_label), vjust = -0.35, size = 10) +
  facet_grid(outcome ~ policy_label, scales = "free_y") +
  scale_fill_manual(
    values = severity_fill_cols,
    breaks = vax_order,
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = 1:5,
    labels = NULL,
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0.08, 0.18))
  ) +
  labs(
    title = "Severe outcomes across NPI and booster scenarios",
    x = NULL,
    y = "Severe outcomes",
    fill = NULL
  ) +
  theme_clean() +
  theme(
    axis.text.x = element_text(size = 30),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 35),
    axis.title.y = element_text(size = 35, colour = "black"),
    strip.text = element_text(size = 35),
    strip.background = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 25),
    legend.background = element_blank(),
    legend.spacing.x = unit(0.45, "cm"),
    plot.title = element_text(size = 40, hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(2, "cm"),
    plot.margin = margin(t = 10, r = 20, b = 55, l = 20),
    plot.background = element_blank()
  ) +
  geom_text(
    data = figure2_df |> filter(outcome == "Deaths"),
    aes(y = -Inf, label = x_lab),
    vjust = 1.15,
    size = 8,
    lineheight = 0.85,
    inherit.aes = TRUE
  ) +
  coord_cartesian(clip = "off")

save_png(
  file.path(fig_dir, "figure3.png"),
  figure2,
  width  = 10000,
  height = 5000
)



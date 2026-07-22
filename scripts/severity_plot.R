suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggthemes)
  library(shellpipes)
  library(scales)
  library(grid)
  library(patchwork)
  library(lubridate)
})

options(dplyr.summarise.inform = FALSE)
options(scipen = 9999)

loadEnvironments()

scenario_order <- c(
  "Slower Boosters",
  "National Rate",
  "Actual Boosters",
  "Faster Boosters",
  "Early Boosters"
)

policy_order <- c("No NPIs", "Actual NPIs", "Strict NPIs")

scenario_cols <- c(
  "Slower Boosters" = "magenta",
  "National Rate"   = "#33a02c",
  "Actual Boosters" = "#4D4D4D",
  "Faster Boosters" = "#ff7f00",
  "Early Boosters"  = "#1f78b4"
)

vax_short_labs <- c(
  "Slower Boosters" = "Slower\nBoosters",
  "National Rate"   = "National\nRate",
  "Actual Boosters" = "Actual\nBoosters",
  "Faster Boosters" = "Faster\nBoosters",
  "Early Boosters"  = "Early\nBoosters"
)

study_start <- as.Date(omicron_start)
intervention_date <- as.Date(intervention_removed)

combosims <- bind_rows(rdsReadList()) |>
  mutate(
    vax_scenario = recode(
      as.character(vac),
      "vax_actual_boosters" = "Actual Boosters",
      "vax_early_boosters"  = "Early Boosters",
      "vax_faster_rate"     = "Faster Boosters",
      "vax_slower_rate"     = "Slower Boosters",
      "vax_national_rate"   = "National Rate"
    ),
    policy = recode(
      as.character(beta),
      "no_npis"     = "No NPIs",
      "base"        = "Actual NPIs",
      "strict_npis" = "Strict NPIs"
    ),
    date = as.Date(burn_in_start) + (time - 1)
  ) |>
  filter(vax_scenario %in% scenario_order, policy %in% policy_order)

cumulative_infections_df <- combosims |>
  filter(matrix %in% c("cum_inc_total", "cum_inc_booster")) |>
  select(vax_scenario, policy, date, matrix, value) |>
  pivot_wider(names_from = matrix, values_from = value, values_fill = 0) |>
  group_by(vax_scenario, policy) |>
  arrange(date, .by_group = TRUE) |>
  mutate(
    baseline_total = {
      x <- cum_inc_total[date < study_start]
      if (length(x) == 0) 0 else tail(x, 1)
    },
    baseline_booster = {
      x <- cum_inc_booster[date < study_start]
      if (length(x) == 0) 0 else tail(x, 1)
    },
    cum_inc_total = pmax(cum_inc_total - baseline_total, 0),
    cum_inc_booster = pmax(cum_inc_booster - baseline_booster, 0),
    cum_inc_leq2dose = pmax(cum_inc_total - cum_inc_booster, 0)
  ) |>
  ungroup() |>
  filter(date >= study_start) |>
  mutate(
    vax_scenario = factor(vax_scenario, levels = scenario_order),
    policy = factor(policy, levels = policy_order)
  )

cumulative_ymax <- max(cumulative_infections_df$cum_inc_leq2dose, na.rm = TRUE)
if (!is.finite(cumulative_ymax) || cumulative_ymax <= 0) cumulative_ymax <- 1

calculate_study_period_outcome <- function(data, matrix_name, outcome_name) {
  outcome_data <- data |>
    filter(matrix == matrix_name, is.finite(value)) |>
    arrange(policy, vax_scenario, date)
  
  baseline_df <- outcome_data |>
    filter(date < study_start) |>
    group_by(policy, vax_scenario) |>
    slice_max(date, n = 1, with_ties = FALSE) |>
    ungroup() |>
    transmute(policy, vax_scenario, baseline_value = value)
  
  endpoint_df <- outcome_data |>
    filter(date >= study_start) |>
    group_by(policy, vax_scenario) |>
    slice_max(date, n = 1, with_ties = FALSE) |>
    ungroup() |>
    transmute(policy, vax_scenario, endpoint_value = value)
  
  endpoint_df |>
    left_join(baseline_df, by = c("policy", "vax_scenario")) |>
    mutate(
      baseline_value = coalesce(baseline_value, 0),
      severe_outcome = pmax(endpoint_value - baseline_value, 0),
      outcome = outcome_name
    ) |>
    select(policy, vax_scenario, outcome, severe_outcome)
}

severe_df <- bind_rows(
  calculate_study_period_outcome(combosims, "cum_hosp_total", "Hospitalizations"),
  calculate_study_period_outcome(combosims, "cum_death_total", "Deaths")
) |>
  filter(is.finite(severe_outcome)) |>
  mutate(
    policy = factor(policy, levels = policy_order),
    vax_scenario = factor(vax_scenario, levels = scenario_order),
    outcome = factor(outcome, levels = c("Hospitalizations", "Deaths"))
  )

make_cumulative_panel <- function(policy_name, panel_tag, panel_title) {
  panel_df <- cumulative_infections_df |>
    filter(policy == policy_name)
  
  ggplot(panel_df, aes(
    x = date,
    y = cum_inc_leq2dose,
    colour = vax_scenario,
    group = vax_scenario
  )) +
    geom_line(linewidth = 2.2) +
    geom_vline(
      xintercept = intervention_date,
      linetype = "dashed",
      linewidth = 1.3
    ) +
    scale_colour_manual(
      values = scenario_cols,
      breaks = scenario_order,
      drop = FALSE,
      name = NULL
    ) +
    scale_x_date(
      date_breaks = "1 month",
      date_labels = "%d %b %Y",
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_y_continuous(
      labels = comma,
      limits = c(0, cumulative_ymax * 1.04),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = panel_title,
      x = NULL,
      y = if (policy_name == "No NPIs") "Cumulative infections\n(≤2 doses)" else NULL,
      tag = panel_tag
    ) +
    theme_clean() +
    theme(
      axis.text.x = element_text(size = 30, hjust = 1, angle = 45),
      axis.text.y = element_text(size = 35),
      axis.title.x = element_text(size = 30),
      axis.title.y = element_text(size = 40),
      plot.title = element_text(size = 35, hjust = 0.5, face = "plain"),
      plot.subtitle = element_text(size = 25, hjust = 0.5),
      legend.position = "none",
      legend.justification = "center",
      legend.direction = "horizontal",
      legend.text = element_text(size = 35),
      legend.key.width = unit(1.2, "cm"),
      legend.key.height = unit(0.35, "cm"),
      legend.spacing.x = unit(0.12, "cm"),
      legend.box.spacing = unit(0.1, "cm"),
      legend.key.spacing.y = unit(0.5, "cm"),
      legend.background = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.background = element_blank(),
      plot.tag = element_text(size = 45, face = "plain"),
      plot.margin = margin(12, 35, 12, 12)
    )
}

make_severe_panel <- function(policy_name, panel_tag) {
  panel_df <- severe_df |>
    filter(policy == policy_name) |>
    select(vax_scenario, outcome, severe_outcome) |>
    complete(
      vax_scenario = factor(scenario_order, levels = scenario_order),
      outcome = factor(
        c("Hospitalizations", "Deaths"),
        levels = c("Hospitalizations", "Deaths")
      ),
      fill = list(severe_outcome = 0)
    ) |>
    mutate(
      x = as.numeric(vax_scenario),
      x_plot = x + if_else(outcome == "Hospitalizations", -0.19, 0.19)
    )
  
  hosp_max <- max(
    panel_df$severe_outcome[panel_df$outcome == "Hospitalizations"],
    na.rm = TRUE
  )
  death_max <- max(
    panel_df$severe_outcome[panel_df$outcome == "Deaths"],
    na.rm = TRUE
  )
  
  scale_factor <- if (is.finite(death_max) && death_max > 0) hosp_max / death_max else 1
  
  panel_df <- panel_df |>
    mutate(
      plotted_value = if_else(
        outcome == "Deaths",
        severe_outcome * scale_factor,
        severe_outcome
      ),
      value_label = comma(round(severe_outcome))
    )
  
  ggplot(panel_df, aes(x = x_plot, y = plotted_value)) +
    geom_col(
      data = panel_df |> filter(outcome == "Hospitalizations"),
      aes(fill = vax_scenario),
      width = 0.34,
      colour = "black",
      linewidth = 0.45,
      alpha = 1
    ) +
    geom_col(
      data = panel_df |> filter(outcome == "Deaths"),
      aes(fill = vax_scenario),
      width = 0.3,
      colour = "black",
      linewidth = 0.65,
      alpha = 0.35
    ) +
    geom_text(aes(label = value_label), vjust = -0.35, size = 8) +
    scale_fill_manual(values = scenario_cols, breaks = scenario_order, drop = FALSE) +
    scale_x_continuous(
      breaks = seq_along(scenario_order),
      labels = unname(vax_short_labs[scenario_order]),
      expand = expansion(mult = c(0.06, 0.06))
    ) +
    scale_y_continuous(
      name = if (policy_name == "No NPIs") "Hospitalizations" else NULL,
      labels = comma,
      expand = expansion(mult = c(0, 0.18)),
      sec.axis = sec_axis(
        trans = ~ . / scale_factor,
        name = if (policy_name == "Strict NPIs") "Deaths" else NULL,
        labels = comma
      )
    ) +
    labs(x = NULL, tag = panel_tag) +
    guides(fill = "none") +
    coord_cartesian(clip = "off") +
    theme_clean() +
    theme(
      axis.text.x = element_text(size = 30, hjust = 0.5),
      axis.text.y = element_text(size = 35),
      axis.title.x = element_text(size = 30),
      axis.title.y = element_text(size = 40),
      plot.title = element_text(size = 35, hjust = 0.5, face = "plain"),
      plot.subtitle = element_text(size = 25, hjust = 0.5),
      legend.position = "top",
      legend.justification = "center",
      legend.direction = "horizontal",
      legend.text = element_text(size = 35),
      legend.spacing.y = unit(0.8, "cm"),
      legend.key.spacing.y = unit(0.5, "cm"),
      legend.background = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.background = element_blank(),
      plot.tag = element_text(size = 45, face = "plain"),
      plot.margin = margin(12, 35, 12, 12)
    )
}

p_a <- make_cumulative_panel("No NPIs", "(A)", "No NPIs")
p_b <- make_cumulative_panel("Actual NPIs", "(B)", "Actual NPIs")
p_c <- make_cumulative_panel("Strict NPIs", "(C)", "Strict NPIs")
p_d <- make_severe_panel("No NPIs", "(D)")
p_e <- make_severe_panel("Actual NPIs", "(E)")
p_f <- make_severe_panel("Strict NPIs", "(F)")

scenario_legend_grob <- grobTree(
  segmentsGrob(
    x0 = unit(c(0.21, 0.37, 0.53, 0.69, 0.85), "npc"),
    x1 = unit(c(0.245, 0.405, 0.565, 0.725, 0.885), "npc"),
    y0 = unit(rep(0.5, 5), "npc"),
    y1 = unit(rep(0.5, 5), "npc"),
    gp = gpar(
      col = unname(scenario_cols[scenario_order]),
      lwd = 6
    )
  ),
  textGrob(
    label = scenario_order,
    x = unit(c(0.252, 0.412, 0.572, 0.732, 0.892), "npc"),
    y = unit(rep(0.5, 5), "npc"),
    just = "left",
    gp = gpar(fontsize = 35)
  )
)

scenario_legend <- wrap_elements(full = scenario_legend_grob)

outcome_legend_grob <- grobTree(
  rectGrob(
    x = unit(c(0.44, 0.56), "npc"),
    y = unit(c(0.5, 0.5), "npc"),
    width = unit(c(0.01, 0.01), "npc"),
    height = unit(c(0.36, 0.36), "npc"),
    gp = gpar(
      fill = c("grey35", "grey75"),
      col = "black",
      lwd = 1
    )
  ),
  textGrob(
    label = c("Hospitalizations", "Deaths"),
    x = unit(c(0.457, 0.577), "npc"),
    y = unit(c(0.5, 0.5), "npc"),
    just = "left",
    gp = gpar(fontsize = 35)
  )
)

outcome_legend <- wrap_elements(full = outcome_legend_grob)

cumulative_row <- (p_a | p_b | p_c) + plot_layout(widths = c(1, 1, 1))
severe_row <- (p_d | p_e | p_f) + plot_layout(widths = c(1, 1, 1))

figure3 <- scenario_legend /
  cumulative_row /
  outcome_legend /
  severe_row +
  plot_layout(
    widths = 1,
    heights = c(0.045, 1.03, 0.05, 1)
  )

png(
  "../figures/figure3.png",
  width = 13500,
  height = 7200,
  res = 300,
  bg = "white",
  type = "cairo"
)
print(figure3)
dev.off()
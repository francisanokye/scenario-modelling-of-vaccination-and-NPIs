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
  "national rate",
  "national rate (nl init)",
  "actual (can init)",
  "actual"
)

policy_order <- c("No NPIs", "Actual NPIs", "Strict NPIs")

scenario_cols <- c(
  "national rate" = "#33A02C",
  "national rate (nl init)" = "goldenrod",
  "actual (can init)" = "#006400",
  "actual" = "#000000"
)

vax_short_labs <- c(
  "national rate" = "Canadian\nrate",
  "national rate (nl init)" = "Canadian\nrate\n(NL initial)",
  "actual (can init)" = "Actual NL\n(Canadian\ninitial)",
  "actual" = "Actual\nNL"
)

vax_legend_labs <- c(
  "national rate" = "Canadian rate",
  "national rate (nl init)" = "Canadian rate (NL initial)",
  "actual (can init)" = "Actual NL (Canadian initial)",
  "actual" = "Actual NL"
)

study_start <- as.Date(omicron_start)
intervention_date <- as.Date(intervention_removed)
plot_end <- as.Date("2022-12-15")

combosims <- bind_rows(rdsReadList()) |>
  mutate(
    vax_scenario = recode(
      as.character(vac),
      "vax_national_rate" = "national rate",
      "nl_init_national_rate" = "national rate (nl init)",
      "vax_actual_rate_can_init" = "actual (can init)",
      "vax_actual_boosters" = "actual"
    ),
    policy = recode(
      as.character(beta),
      "no_npis" = "No NPIs",
      "base" = "Actual NPIs",
      "strict_npis" = "Strict NPIs"
    ),
    date = as.Date(burn_in_start) + (time - 1)
  ) |>
  filter(
    vax_scenario %in% scenario_order,
    policy %in% policy_order
  )

# ------------------------------------------------------------
# Cumulative infections without booster doses for Figure A-C
# ------------------------------------------------------------

cumulative_infections_df <- combosims |>
  filter(matrix %in% c("cum_inc_total", "cum_inc_booster")) |>
  select(vax_scenario, policy, date, matrix, value) |>
  pivot_wider(
    names_from = matrix,
    values_from = value,
    values_fill = 0
  )

if (!"cum_inc_total" %in% names(cumulative_infections_df)) {
  stop("Matrix 'cum_inc_total' was not found in the simulation output")
}

if (!"cum_inc_booster" %in% names(cumulative_infections_df)) {
  stop("Matrix 'cum_inc_booster' was not found in the simulation output")
}

cumulative_infections_df <- cumulative_infections_df |>
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
  filter(
    date >= study_start,
    date <= plot_end
  ) |>
  mutate(
    vax_scenario = factor(vax_scenario, levels = scenario_order),
    policy = factor(policy, levels = policy_order)
  )

no_booster_ymax <- max(
  cumulative_infections_df$cum_inc_leq2dose,
  na.rm = TRUE
)

if (!is.finite(no_booster_ymax) || no_booster_ymax <= 0) {
  no_booster_ymax <- 1
}

make_cumulative_panel <- function(
    policy_name,
    panel_tag,
    y_label = NULL
) {
  
  panel_df <- cumulative_infections_df |>
    filter(policy == policy_name)
  
  ggplot(
    panel_df,
    aes(
      x = date,
      y = cum_inc_leq2dose,
      colour = vax_scenario,
      group = vax_scenario
    )
  ) +
    geom_line(linewidth = 3.5) +
    geom_vline(
      xintercept = intervention_date,
      linetype = "dashed",
      linewidth = 2
    ) +
    scale_colour_manual(
      values = scenario_cols,
      breaks = scenario_order,
      labels = unname(vax_legend_labs[scenario_order]),
      drop = FALSE,
      name = NULL
    ) +
    scale_x_date(
      date_breaks = "1 month",
      date_labels = "%d %b %Y",
      limits = c(study_start, plot_end),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_y_continuous(
      labels = comma,
      limits = c(0, no_booster_ymax * 1.04),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = policy_name,
      x = NULL,
      y = y_label,
      tag = panel_tag
    ) +
    theme_clean() +
    theme(
      axis.text.x = element_text(size = 50, hjust = 1, angle = 45),
      axis.text.y = element_text(size = 50),
      axis.title.x = element_text(size = 50),
      axis.title.y = element_text(size = 50),
      plot.title = element_text(size = 50, hjust = 0.5, face = "plain"),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      plot.background = element_blank(),
      plot.tag = element_text(size = 50, face = "plain"),
      plot.margin = margin(12, 35, 12, 12)
    )
}

# ------------------------------------------------------------
# Severe outcomes for Figure D-F
# ------------------------------------------------------------

calculate_study_period_outcome <- function(
    data,
    matrix_name,
    outcome_name
) {
  
  outcome_data <- data |>
    filter(
      matrix == matrix_name,
      is.finite(value)
    ) |>
    arrange(policy, vax_scenario, date)
  
  baseline_df <- outcome_data |>
    filter(date < study_start) |>
    group_by(policy, vax_scenario) |>
    slice_max(date, n = 1, with_ties = FALSE) |>
    ungroup() |>
    transmute(
      policy,
      vax_scenario,
      baseline_value = value
    )
  
  endpoint_df <- outcome_data |>
    filter(
      date >= study_start,
      date <= plot_end
    ) |>
    group_by(policy, vax_scenario) |>
    slice_max(date, n = 1, with_ties = FALSE) |>
    ungroup() |>
    transmute(
      policy,
      vax_scenario,
      endpoint_value = value
    )
  
  endpoint_df |>
    left_join(
      baseline_df,
      by = c("policy", "vax_scenario")
    ) |>
    mutate(
      baseline_value = coalesce(baseline_value, 0),
      severe_outcome = pmax(endpoint_value - baseline_value, 0),
      outcome = outcome_name
    ) |>
    select(
      policy,
      vax_scenario,
      outcome,
      severe_outcome
    )
}

severe_df <- bind_rows(
  calculate_study_period_outcome(
    combosims,
    "cum_hosp_total",
    "Hospitalizations"
  ),
  calculate_study_period_outcome(
    combosims,
    "cum_death_total",
    "Deaths"
  )
) |>
  filter(is.finite(severe_outcome)) |>
  mutate(
    policy = factor(policy, levels = policy_order),
    vax_scenario = factor(vax_scenario, levels = scenario_order),
    outcome = factor(
      outcome,
      levels = c("Hospitalizations", "Deaths")
    )
  )

severe_hosp_max <- max(
  severe_df$severe_outcome[
    severe_df$outcome == "Hospitalizations"
  ],
  na.rm = TRUE
)

severe_death_max <- max(
  severe_df$severe_outcome[
    severe_df$outcome == "Deaths"
  ],
  na.rm = TRUE
)

severe_scale_factor <- if (
  is.finite(severe_death_max) &&
  severe_death_max > 0
) {
  severe_hosp_max / severe_death_max
} else {
  1
}

# Add headroom for labels
severe_ymax <- severe_hosp_max * 1.12

make_severe_panel <- function(
    policy_name,
    panel_tag,
    y_label = NULL
) {
  
  panel_df <- severe_df |>
    filter(policy == policy_name) |>
    select(
      vax_scenario,
      outcome,
      severe_outcome
    ) |>
    complete(
      vax_scenario = factor(
        scenario_order,
        levels = scenario_order
      ),
      outcome = factor(
        c("Hospitalizations", "Deaths"),
        levels = c("Hospitalizations", "Deaths")
      ),
      fill = list(severe_outcome = 0)
    ) |>
    mutate(
      x = as.numeric(vax_scenario),
      x_plot = x + if_else(
        outcome == "Hospitalizations",
        -0.21,
        0.21
      ),
      plotted_value = if_else(
        outcome == "Deaths",
        severe_outcome * severe_scale_factor,
        severe_outcome
      ),
      value_label = comma(round(severe_outcome))
    )
  
  hosp_df <- panel_df |>
    filter(outcome == "Hospitalizations")
  
  death_df <- panel_df |>
    filter(outcome == "Deaths")
  
  ggplot(
    panel_df,
    aes(
      x = x_plot,
      y = plotted_value
    )
  ) +
    geom_col(
      data = hosp_df,
      aes(fill = vax_scenario),
      width = 0.34,
      colour = "black",
      linewidth = 0.45,
      alpha = 1
    ) +
    geom_col(
      data = death_df,
      aes(fill = vax_scenario),
      width = 0.30,
      colour = "black",
      linewidth = 0.65,
      alpha = 0.35
    ) +
    
    # Hospitalization labels: move slightly left
    geom_text(
      data = hosp_df,
      aes(
        x = x_plot - 0.035,
        y = plotted_value,
        label = value_label
      ),
      inherit.aes = FALSE,
      hjust = 0.5,
      vjust = -0.45,
      size = 15
    ) +
    
    # Death labels: move slightly right
    geom_text(
      data = death_df,
      aes(
        x = x_plot + 0.035,
        y = plotted_value,
        label = value_label
      ),
      inherit.aes = FALSE,
      hjust = 0.5,
      vjust = -0.45,
      size = 15
    ) +
    
    scale_fill_manual(
      values = scenario_cols,
      breaks = scenario_order,
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = seq_along(scenario_order),
      labels = unname(vax_short_labs[scenario_order]),
      limits = c(0.55, 4.45),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      name = y_label,
      labels = comma,
      limits = c(0, severe_ymax),
      expand = expansion(mult = c(0, 0)),
      sec.axis = sec_axis(
        trans = ~ . / severe_scale_factor,
        name = "Deaths",
        labels = comma
      )
    ) +
    labs(
      x = NULL,
      tag = panel_tag
    ) +
    guides(fill = "none") +
    coord_cartesian(clip = "off") +
    theme_clean() +
    theme(
      axis.text.x = element_text(
        size = 40,
        hjust = 0.5,
        lineheight = 0.9
      ),
      axis.text.y = element_text(size = 50),
      axis.title.x = element_text(size = 50),
      axis.title.y = element_text(size = 50),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      plot.background = element_blank(),
      plot.tag = element_text(size = 50, face = "plain"),
      plot.margin = margin(20, 45, 12, 25)
    )
}

# ------------------------------------------------------------
# Panels
# ------------------------------------------------------------

p_a <- make_cumulative_panel(
  "No NPIs",
  "(A)",
  "Cumulative infections\n(no booster dose)"
)

p_b <- make_cumulative_panel(
  "Actual NPIs",
  "(B)"
)

p_c <- make_cumulative_panel(
  "Strict NPIs",
  "(C)"
)

p_d <- make_severe_panel(
  "No NPIs",
  "(D)",
  "Hospitalizations"
)

p_e <- make_severe_panel(
  "Actual NPIs",
  "(E)"
)

p_f <- make_severe_panel(
  "Strict NPIs",
  "(F)"
)

# ------------------------------------------------------------
# Vaccination scenario legend
# ------------------------------------------------------------

scenario_legend_grob <- grobTree(
  segmentsGrob(
    x0 = unit(c(0.07, 0.29, 0.54, 0.79), "npc"),
    x1 = unit(c(0.105, 0.325, 0.575, 0.825), "npc"),
    y0 = unit(rep(0.5, 4), "npc"),
    y1 = unit(rep(0.5, 4), "npc"),
    gp = gpar(
      col = unname(scenario_cols[scenario_order]),
      lwd = 9
    )
  ),
  textGrob(
    label = unname(vax_legend_labs[scenario_order]),
    x = unit(c(0.112, 0.332, 0.582, 0.832), "npc"),
    y = unit(rep(0.5, 4), "npc"),
    just = "left",
    gp = gpar(fontsize = 50)
  )
)

scenario_legend <- wrap_elements(
  full = scenario_legend_grob
)

# ------------------------------------------------------------
# Severe outcome legend
# ------------------------------------------------------------

outcome_legend_grob <- grobTree(
  rectGrob(
    x = unit(c(0.425, 0.565), "npc"),
    y = unit(c(0.5, 0.5), "npc"),
    width = unit(c(0.014, 0.014), "npc"),
    height = unit(c(0.28, 0.28), "npc"),
    gp = gpar(
      fill = c("grey35", "grey75"),
      col = "black",
      lwd = 1
    )
  ),
  textGrob(
    label = c("Hospitalizations", "Deaths"),
    x = unit(c(0.442, 0.582), "npc"),
    y = unit(c(0.5, 0.5), "npc"),
    just = "left",
    gp = gpar(fontsize = 50)
  )
)

outcome_legend <- wrap_elements(
  full = outcome_legend_grob
)

# ------------------------------------------------------------
# Assemble SI Figure 3
# ------------------------------------------------------------

cumulative_row <- (p_a | p_b | p_c) +
  plot_layout(widths = c(1, 1, 1))

severe_row <- (p_d | p_e | p_f) +
  plot_layout(widths = c(1, 1, 1))

figureS3 <- scenario_legend /
  cumulative_row /
  outcome_legend /
  severe_row +
  plot_layout(
    widths = 1,
    heights = c(
      0.10,
      0.90,
      0.09,
      1.08
    )
  )

png(
  "../figures/SI_Fig3.png",
  width = 17000,
  height = 7600,
  res = 300,
  bg = "white",
  type = "cairo"
)

print(figureS3)
dev.off()
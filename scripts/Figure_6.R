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

scenario_order <- c("slower rate", "national rate", "actual", "faster rate", "earlier")
bar_scenario_order <- c("slower rate", "national rate", "actual", "earlier", "faster rate")
policy_order <- c("No NPIs", "Actual NPIs", "Strict NPIs")

scenario_cols <- c(
  "actual" = "#4D4D4D",
  "earlier" = "#1F78B4",
  "faster rate" = "#FF7F00",
  "slower rate" = "#CC79A7",
  "national rate" = "#33A02C"
)

# Keep the original two-line labels where needed
vax_short_labs <- c(
  "slower rate" = "Slower\nrate",
  "national rate" = "Canadian\nrate",
  "actual" = "Actual\nNL",
  "earlier" = "Earlier",
  "faster rate" = "Faster\nrate"
)

study_start <- as.Date(omicron_start)

combosims <- bind_rows(rdsReadList()) |>
  mutate(
    vax_scenario = recode(
      as.character(vac),
      "vax_actual_boosters" = "actual",
      "vax_early_boosters" = "earlier",
      "vax_faster_rate" = "faster rate",
      "vax_slower_rate" = "slower rate",
      "vax_national_rate" = "national rate"
    ),
    policy = recode(
      as.character(beta),
      "no_npis" = "No NPIs",
      "base" = "Actual NPIs",
      "strict_npis" = "Strict NPIs"
    ),
    date = as.Date(burn_in_start) + (time - 1)
  ) |>
  filter(vax_scenario %in% scenario_order, policy %in% policy_order)

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

severe_hosp_max <- max(
  severe_df$severe_outcome[severe_df$outcome == "Hospitalizations"],
  na.rm = TRUE
)

severe_death_max <- max(
  severe_df$severe_outcome[severe_df$outcome == "Deaths"],
  na.rm = TRUE
)

severe_scale_factor <- if (is.finite(severe_death_max) && severe_death_max > 0) {
  severe_hosp_max / severe_death_max
} else {
  1
}

# Add headroom for labels above the tallest bars
severe_ymax <- severe_hosp_max * 1.12

make_severe_panel <- function(policy_name, panel_tag, panel_title) {
  panel_df <- severe_df |>
    filter(policy == policy_name) |>
    mutate(vax_scenario = factor(as.character(vax_scenario), levels = bar_scenario_order)) |>
    select(vax_scenario, outcome, severe_outcome) |>
    complete(
      vax_scenario = factor(bar_scenario_order, levels = bar_scenario_order),
      outcome = factor(
        c("Hospitalizations", "Deaths"),
        levels = c("Hospitalizations", "Deaths")
      ),
      fill = list(severe_outcome = 0)
    ) |>
    mutate(
      vax_scenario = factor(as.character(vax_scenario), levels = bar_scenario_order),
      x = as.numeric(vax_scenario),
      
      # Slightly increase separation between paired bars
      x_plot = x + if_else(
        outcome == "Hospitalizations",
        -0.20,
        0.20
      ),
      
      plotted_value = if_else(
        outcome == "Deaths",
        severe_outcome * severe_scale_factor,
        severe_outcome
      ),
      
      value_label = comma(round(severe_outcome)),
      
      # Move labels slightly outward from one another
      label_x = x_plot + if_else(
        outcome == "Hospitalizations",
        -0.025,
        0.045
      ),
      
      # Slightly stagger the two labels vertically
      label_y = plotted_value + if_else(
        outcome == "Hospitalizations",
        severe_hosp_max * 0.025,
        severe_hosp_max * 0.050
      )
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
      width = 0.30,
      colour = "black",
      linewidth = 0.65,
      alpha = 0.35
    ) +
    geom_text(
      data = panel_df |> filter(outcome == "Hospitalizations"),
      aes(
        x = label_x,
        y = label_y,
        label = value_label
      ),
      size = 18,
      hjust = 0.6,
      vjust = 0
    ) +
    geom_text(
      data = panel_df |> filter(outcome == "Deaths"),
      aes(
        x = label_x,
        y = label_y,
        label = value_label
      ),
      size = 18,
      hjust = 0.4,
      vjust = 0
    ) +
    scale_fill_manual(
      values = scenario_cols,
      breaks = bar_scenario_order,
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = seq_along(bar_scenario_order),
      labels = unname(vax_short_labs[bar_scenario_order]),
      expand = expansion(mult = c(0.07, 0.07))
    ) +
    scale_y_continuous(
      name = "Hospitalizations",
      labels = comma,
      limits = c(0, severe_ymax),
      expand = expansion(mult = c(0, 0)),
      sec.axis = sec_axis(
        trans = ~ . / severe_scale_factor,
        name = "Deaths",
        labels = comma
      )
    ) +
    labs(title = panel_title,x = NULL,tag = panel_tag) +
    guides(fill = "none") +
    coord_cartesian(clip = "off") +
    theme_clean() +
    theme(
      axis.text.x = element_text(size = 50,hjust = 0.5, lineheight = 0.90, margin = margin(t = 8)),
      axis.text.y = element_text(size = 60),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 60),
      plot.title = element_text(size = 60,hjust = 0.5,face = "plain"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_blank(),
      plot.tag = element_text(size = 60,face = "plain"),
      plot.margin = margin(12, 30, 12, 12)
    )
}

p_a <- make_severe_panel(
  "No NPIs",
  "(A)",
  "No NPIs"
)

p_b <- make_severe_panel(
  "Actual NPIs",
  "(B)",
  "Actual NPIs"
)

p_c <- make_severe_panel(
  "Strict NPIs",
  "(C)",
  "Strict NPIs"
)

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
    label = c(
      "Hospitalizations",
      "Deaths"
    ),
    x = unit(c(0.457, 0.577), "npc"),
    y = unit(c(0.5, 0.5), "npc"),
    just = "left",
    gp = gpar(fontsize = 60)
  )
)

outcome_legend <- wrap_elements(
  full = outcome_legend_grob
)

severe_row <- (p_a | p_b | p_c) +
  plot_layout(
    widths = c(1, 1, 1)
  )

figure3b <- (outcome_legend / severe_row) +
  plot_layout(
    heights = c(0.07, 1)
  )

figure3b <- figure3b &
  theme(
    plot.background = element_blank()
  )

png(
  "../figures/figure6.png",
  width = 19000,
  height = 5500,
  res = 300,
  bg = "white",
  type = "cairo"
)

print(figure3b)

dev.off()
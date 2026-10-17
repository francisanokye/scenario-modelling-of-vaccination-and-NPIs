
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
policy_order <- c("No NPIs", "Actual NPIs", "Strict NPIs")

scenario_cols <- c(
  "actual" = "#4D4D4D",
  "earlier" = "#1F78B4",
  "faster rate" = "#FF7F00",
  "slower rate" = "#CC79A7",
  "national rate" = "#33A02C"
)

scenario_labs <- c(
  "slower rate" = "Slower rate",
  "national rate" = "Canadian rate",
  "actual" = "Actual NL",
  "faster rate" = "Faster rate",
  "earlier" = "Earlier"
)

row_scenarios <- list(
  "Speed" = c("slower rate", "actual", "faster rate"),
  "Timing" = c("earlier", "actual"),
  "Actual" = c("national rate", "actual")
)

study_start <- as.Date(omicron_start)
intervention_date <- as.Date(intervention_removed)

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

cumulative_infections_df <- combosims |>
  filter(matrix %in% c("cum_inc_total", "cum_inc_booster")) |>
  select(vax_scenario, policy, date, matrix, value) |>
  pivot_wider(names_from = matrix, values_from = value, values_fill = 0)

if (!"cum_inc_total" %in% names(cumulative_infections_df)) {
  stop("Matrix 'cum_inc_total' was not found in the simulation output.")
}

if (!"cum_inc_booster" %in% names(cumulative_infections_df)) {
  cumulative_infections_df$cum_inc_booster <- 0
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
  filter(date >= study_start) |>
  mutate(
    vax_scenario = factor(vax_scenario, levels = scenario_order),
    policy = factor(policy, levels = policy_order)
  )

cumulative_ymax <- max(cumulative_infections_df$cum_inc_leq2dose, na.rm = TRUE)
if (!is.finite(cumulative_ymax) || cumulative_ymax <= 0) cumulative_ymax <- 1


make_cumulative_panel <- function(
    scenario_names,
    policy_name,
    panel_tag,
    panel_title = NULL,
    show_x_text = FALSE
) {
  
  panel_df <- cumulative_infections_df |>
    filter(vax_scenario %in% scenario_names, policy == policy_name) |>
    mutate(vax_scenario = factor(as.character(vax_scenario), levels = scenario_names))
  
  if (nrow(panel_df) == 0) {
    stop("No cumulative infection data found for ", paste(scenario_names, collapse = ", "), " under ", policy_name)
  }
  
  panel_cols <- scenario_cols[scenario_names]
  panel_labs <- scenario_labs[scenario_names]
  
  ggplot(panel_df, aes(x = date, y = cum_inc_leq2dose, colour = vax_scenario, group = vax_scenario)) +
    geom_line(linewidth = 3.5) +
    geom_vline(xintercept = intervention_date, linetype = "dashed", linewidth = 1.8) +
    scale_colour_manual(values = panel_cols, breaks = scenario_names, labels = panel_labs, name = NULL, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(linewidth = 2.5))) +
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
    labs(title = panel_title, x = NULL, y = NULL, tag = panel_tag) +
    theme_clean() +
    theme(
      axis.text.x = if (show_x_text) element_text(size = 45, hjust = 1, angle = 45) else element_blank(),
      axis.ticks.x = if (show_x_text) element_line() else element_blank(),
      axis.text.y = element_text(size = 50),
      axis.title = element_blank(),
      plot.title = element_text(size = 50, hjust = 0.5, face = "plain", margin = margin(b = 2)),
      legend.position = "top",
      legend.justification = "center",
      legend.direction = "horizontal",
      legend.text = element_text(size = 45),
      legend.key = element_blank(),
      legend.background = element_blank(),
      legend.key.width = unit(0.9, "cm"),
      legend.key.height = unit(0.25, "cm"),
      legend.spacing.x = unit(0.08, "cm"),
      legend.box.spacing = unit(0, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      panel.grid.major.x = element_blank(),
      plot.background = element_blank(),
      plot.tag = element_text(size = 50, face = "plain"),
      plot.margin = margin(3, 12, 3, 6)
    )
}


make_row_label <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, size = 15, hjust = 0.5, vjust = 0.5) +
    xlim(0, 1) + ylim(0, 1) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 2, 0, 2))
}

make_y_label <- function() {
  ggplot() +
    annotate(
      "text", x = 0.5, y = 0.5,
      label = "Cumulative infections\n(no booster dose)",
      angle = 90, size = 18, hjust = 0.5, vjust = 0.5
    ) +
    xlim(0, 1) + ylim(0, 1) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0))
}


p_a <- make_cumulative_panel(row_scenarios[["Speed"]], "No NPIs", "(A)", "No NPIs")
p_b <- make_cumulative_panel(row_scenarios[["Speed"]], "Actual NPIs", "(B)", "Actual NPIs")
p_c <- make_cumulative_panel(row_scenarios[["Speed"]], "Strict NPIs", "(C)", "Strict NPIs")

p_d <- make_cumulative_panel(row_scenarios[["Timing"]], "No NPIs", "(D)")
p_e <- make_cumulative_panel(row_scenarios[["Timing"]], "Actual NPIs", "(E)")
p_f <- make_cumulative_panel(row_scenarios[["Timing"]], "Strict NPIs", "(F)")

p_g <- make_cumulative_panel(row_scenarios[["Actual"]], "No NPIs", "(G)", show_x_text = TRUE)
p_h <- make_cumulative_panel(row_scenarios[["Actual"]], "Actual NPIs", "(H)", show_x_text = TRUE)
p_i <- make_cumulative_panel(row_scenarios[["Actual"]], "Strict NPIs", "(I)", show_x_text = TRUE)

row_1 <- (make_row_label("Speed") | p_a | p_b | p_c) +
  plot_layout(widths = c(0.22, 1, 1, 1))

row_2 <- (make_row_label("Timing") | p_d | p_e | p_f) +
  plot_layout(widths = c(0.22, 1, 1, 1))

row_3 <- (make_row_label("Actual") | p_g | p_h | p_i) +
  plot_layout(widths = c(0.22, 1, 1, 1))

figure3a_body <- (row_1 / row_2 / row_3) +
  plot_layout(heights = c(1, 1, 1))

figure3a <- (make_y_label() | figure3a_body) +
  plot_layout(widths = c(0.065, 1))

figure3a <- figure3a &
  theme(plot.background = element_blank())

png("../figures/SI_Fig6.png", width = 14500, height = 7200, res = 300, bg = "white", type = "cairo")
print(figure3a)
dev.off()


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

max_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else max(x)
}

scenario_order <- c("slower rate", "national rate", "actual", "faster rate", "earlier")
policy_order <- c("No NPIs", "Actual NPIs", "Strict NPIs")

scenario_cols <- c(
  "slower rate" = "magenta",
  "national rate" = "#33A02C",
  "actual" = "#4D4D4D",
  "faster rate" = "#ff7f00",
  "earlier" = "#1f78b4"
)

scenario_labs <- c(
  "slower rate" = "Slower rate",
  "national rate" = "Canadian rate",
  "actual" = "Actual NL",
  "faster rate" = "Faster rate",
  "earlier" = "Earlier"
)

dose_alpha <- c(
  "no booster dose" = 0.28,
  "booster" = 0.72
)

row_scenarios <- list(
  "Speed" = c("slower rate", "actual", "faster rate"),
  "Timing" = c("earlier", "actual"),
  "Actual" = c("national rate", "actual")
)

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

incidence_df <- combosims |>
  filter(matrix %in% c("inc", "inc_booster"), date >= as.Date(omicron_start)) |>
  select(vax_scenario, policy, date, matrix, value) |>
  pivot_wider(names_from = matrix, values_from = value, values_fill = 0)

if (!"inc" %in% names(incidence_df)) {
  stop("Matrix 'inc' was not found in the simulation output.")
}

if (!"inc_booster" %in% names(incidence_df)) {
  incidence_df$inc_booster <- 0
}

incidence_df <- incidence_df |>
  mutate(
    inc = pmax(inc, 0),
    inc_booster = pmax(inc_booster, 0),
    inc_leq2dose = pmax(inc - inc_booster, 0),
    total_inc = inc_leq2dose + inc_booster,
    vax_scenario = factor(vax_scenario, levels = scenario_order),
    policy = factor(policy, levels = policy_order)
  )

global_ymax <- max_or_na(incidence_df$total_inc)

if (!is.finite(global_ymax) || global_ymax <= 0) {
  global_ymax <- 1
}

make_incidence_panel <- function(
    scenario_names,
    policy_name,
    panel_tag,
    panel_title = NULL,
    show_x_text = FALSE
) {
  
  panel_df <- incidence_df |>
    filter(vax_scenario %in% scenario_names, policy == policy_name) |>
    mutate(
      vax_scenario = factor(
        as.character(vax_scenario),
        levels = scenario_names
      )
    )
  
  if (nrow(panel_df) == 0) {
    stop(
      "No incidence data found for ",
      paste(scenario_names, collapse = ", "),
      " under ",
      policy_name
    )
  }
  
  panel_cols <- scenario_cols[scenario_names]
  panel_labs <- scenario_labs[scenario_names]
  
  ggplot(panel_df, aes(x = date, group = vax_scenario)) +
    geom_ribbon(
      aes(
        ymin = 0,
        ymax = inc_leq2dose,
        fill = vax_scenario,
        alpha = "no booster dose"
      ),
      colour = NA
    ) +
    geom_ribbon(
      aes(
        ymin = inc_leq2dose,
        ymax = total_inc,
        fill = vax_scenario,
        alpha = "booster"
      ),
      colour = NA
    ) +
    geom_line(
      aes(y = inc_leq2dose, colour = vax_scenario),
      linewidth = 1.5
    ) +
    geom_line(
      aes(y = total_inc, colour = vax_scenario),
      linewidth = 1.5
    ) +
    geom_vline(
      xintercept = as.Date(intervention_removed),
      linetype = "dashed",
      linewidth = 1.3
    ) +
    scale_fill_manual(
      values = panel_cols,
      breaks = scenario_names,
      labels = panel_labs,
      name = NULL,
      drop = FALSE
    ) +
    scale_colour_manual(
      values = panel_cols,
      breaks = scenario_names,
      labels = panel_labs,
      name = NULL,
      drop = FALSE
    ) +
    scale_alpha_manual(
      values = dose_alpha,
      guide = "none"
    ) +
    guides(
      fill = "none",
      colour = guide_legend(
        order = 1,
        override.aes = list(
          linewidth = 2.5,
          alpha = 1
        )
      )
    ) +
    scale_x_date(
      date_breaks = "1 month",
      date_labels = "%d %b %Y",
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
      trans = "sqrt",
      labels = comma,
      breaks = c(0, 250, 1000, 2500, 5000, 7500),
      expand = expansion(mult = c(0, 0.04))
    ) +
    coord_cartesian(
      ylim = c(0, global_ymax * 1.06),
      clip = "off"
    ) +
    labs(
      title = panel_title,
      x = NULL,
      y = NULL,
      tag = panel_tag
    ) +
    theme_clean() +
    theme(
      axis.text.x = if (show_x_text) {element_text(size = 45, hjust = 1, angle = 45)
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x_text) element_line() else element_blank(),
      axis.text.y = element_text(size = 48),
      axis.title = element_blank(),
      plot.title = element_text(size = 48, hjust = 0.5,face = "plain", margin = margin(b = 2)),
      strip.text = element_text(size = 50),
      legend.position = c(0.8, 0.5),
      legend.justification = "center",
      legend.direction = "vertical",
      legend.text = element_text(size = 48),
      legend.key = element_blank(),
      legend.background = element_blank(),
      legend.key.width = unit(0.9, "cm"),
      legend.key.height = unit(0.45, "cm"),
      legend.spacing.x = unit(0.08, "cm"),
      legend.spacing.y = unit(0.1, "cm"),
      legend.box.spacing = unit(0, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      panel.grid.major.x = element_blank(),
      plot.background = element_blank(),
      plot.tag = element_text(size = 48, face = "plain"),
      plot.margin = margin(3, 12, 3, 6)
    )
}

make_row_label <- function(label) {
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = label,
      size = 20,
      fontface = "plain",
      hjust = 0.5,
      vjust = 0.5,
      lineheight = 0.95
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 2, 0, 2))
}

make_y_label <- function() {
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = "Daily incidence",
      angle = 90,
      size = 20,
      fontface = "plain",
      hjust = 0.5,
      vjust = 0.5
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

p_a <- make_incidence_panel(row_scenarios[["Speed"]], "No NPIs", "(A)", "No NPIs", FALSE)
p_b <- make_incidence_panel(row_scenarios[["Speed"]], "Actual NPIs", "(B)", "Actual NPIs", FALSE)
p_c <- make_incidence_panel(row_scenarios[["Speed"]], "Strict NPIs", "(C)", "Strict NPIs", FALSE)

p_d <- make_incidence_panel(row_scenarios[["Timing"]], "No NPIs", "(D)", NULL, FALSE)
p_e <- make_incidence_panel(row_scenarios[["Timing"]], "Actual NPIs", "(E)", NULL, FALSE)
p_f <- make_incidence_panel(row_scenarios[["Timing"]], "Strict NPIs", "(F)", NULL, FALSE)

p_g <- make_incidence_panel(row_scenarios[["Actual"]], "No NPIs", "(G)", NULL, TRUE)
p_h <- make_incidence_panel(row_scenarios[["Actual"]], "Actual NPIs", "(H)", NULL, TRUE)
p_i <- make_incidence_panel(row_scenarios[["Actual"]], "Strict NPIs", "(I)", NULL, TRUE)

row_1 <- (make_row_label("Speed") | p_a | p_b | p_c) +
  plot_layout(widths = c(0.22, 1, 1, 1))

row_2 <- (make_row_label("Timing") | p_d | p_e | p_f) +
  plot_layout(widths = c(0.22, 1, 1, 1))

row_3 <- (make_row_label("Actual") | p_g | p_h | p_i) +
  plot_layout(widths = c(0.22, 1, 1, 1))

figure2_body <- (row_1 / row_2 / row_3) +
  plot_layout(heights = c(1, 1, 1))

figure2 <- (make_y_label() | figure2_body) +
  plot_layout(widths = c(0.055, 1))

figure2 <- figure2 &
  theme(plot.background = element_blank())

png(
  "../figures/figure4.png",
  width = 14500,
  height = 7200,
  res = 300,
  bg = "white",
  type = "cairo"
)

print(figure2)

dev.off()
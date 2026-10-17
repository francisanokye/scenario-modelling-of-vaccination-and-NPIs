
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

scenario_order <- c("national rate", "national rate (nl init)", "actual (can init)", "actual")
policy_order <- c("No NPIs", "Actual NPIs", "Strict NPIs")

scenario_dose_cols <- list(
  "national rate" = c("no booster dose" = "#b2df8a", "booster" = "#33A02C"),
  "national rate (nl init)" = c("no booster dose" = "#F2D79B", "booster" = "goldenrod"),
  "actual (can init)" = c("no booster dose" = "#A9CFA9", "booster" = "#006400"),
  "actual" = c("no booster dose" = "grey70", "booster" = "#000000")
)

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
  filter(vax_scenario %in% scenario_order, policy %in% policy_order)

incidence_df <- combosims |>
  filter(matrix %in% c("inc", "inc_booster"), date >= as.Date(omicron_start)) |>
  select(vax_scenario, policy, date, matrix, value) |>
  pivot_wider(names_from = matrix, values_from = value, values_fill = 0)

if (!"inc" %in% names(incidence_df)) incidence_df$inc <- 0
if (!"inc_booster" %in% names(incidence_df)) incidence_df$inc_booster <- 0

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
if (!is.finite(global_ymax) || global_ymax <= 0) global_ymax <- 1

get_panel_cols <- function(scenario_name) {
  cols <- scenario_dose_cols[[scenario_name]]
  if (is.null(cols)) stop("No colours are defined for: ", scenario_name)
  cols
}

make_incidence_panel <- function(scenario_name, policy_name, panel_tag, panel_title = NULL, show_x_text = FALSE, show_y_title = FALSE) {
  panel_df <- incidence_df |>
    filter(vax_scenario == scenario_name, policy == policy_name)
  
  if (nrow(panel_df) == 0) {
    stop("No incidence data found for ", scenario_name, " under ", policy_name)
  }
  
  p <- ggplot(panel_df, aes(x = date)) +
    geom_ribbon(aes(ymin = 0, ymax = inc_leq2dose, fill = "no booster dose")) +
    geom_ribbon(aes(ymin = inc_leq2dose, ymax = total_inc, fill = "booster")) +
    geom_line(aes(y = inc_leq2dose), colour = "black", linewidth = 1.5) +
    geom_line(aes(y = total_inc), colour = "black", linewidth = 1.5) +
    geom_vline(
      xintercept = as.Date(intervention_removed),
      linetype = "dashed",
      linewidth = 1.3
    ) +
    scale_fill_manual(
      values = get_panel_cols(scenario_name),
      breaks = c("booster", "no booster dose"),
      name = NULL
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
      y = if (show_y_title) "Daily incidence" else NULL,
      tag = panel_tag
    ) +
    theme_clean() +
    theme(
      axis.text.y = element_text(size = 40),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 40),
      plot.title = element_text(size = 40, hjust = 0.5, face = "plain"),
      plot.subtitle = element_text(size = 25, hjust = 0.5),
      legend.position = "top",
      legend.justification = "center",
      legend.direction = "horizontal",
      legend.text = element_text(size = 40),
      legend.spacing.y = unit(0.8, "cm"),
      legend.key.spacing.y = unit(0.5, "cm"),
      legend.background = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.background = element_blank(),
      plot.tag = element_text(size = 45, face = "plain"),
      plot.margin = margin(12, 35, 12, 12)
    )
  
  if (show_x_text) {
    p <- p +
      theme(
        axis.text.x = element_text(size = 40, hjust = 1, angle = 45),
        axis.ticks.x = element_line()
      )
  } else {
    p <- p +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }
  
  p
}

make_row_label <- function(label) {
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = label,
      size = 15,
      fontface = "plain",
      hjust = 0.5,
      vjust = 0.5,
      lineheight = 0.95
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 5, 0, 5))
}

p_a <- make_incidence_panel("national rate", "No NPIs", "(A)", "No NPIs", FALSE, TRUE)
p_b <- make_incidence_panel("national rate", "Actual NPIs", "(B)", "Actual NPIs")
p_c <- make_incidence_panel("national rate", "Strict NPIs", "(C)", "Strict NPIs")

p_d <- make_incidence_panel("national rate (nl init)", "No NPIs", "(D)", show_y_title = TRUE)
p_e <- make_incidence_panel("national rate (nl init)", "Actual NPIs", "(E)")
p_f <- make_incidence_panel("national rate (nl init)", "Strict NPIs", "(F)")

p_g <- make_incidence_panel("actual (can init)", "No NPIs", "(G)", show_y_title = TRUE)
p_h <- make_incidence_panel("actual (can init)", "Actual NPIs", "(H)")
p_i <- make_incidence_panel("actual (can init)", "Strict NPIs", "(I)")

p_j <- make_incidence_panel("actual", "No NPIs", "(J)", show_x_text = TRUE, show_y_title = TRUE)
p_k <- make_incidence_panel("actual", "Actual NPIs", "(K)", show_x_text = TRUE)
p_l <- make_incidence_panel("actual", "Strict NPIs", "(L)", show_x_text = TRUE)

row_1 <- (
  make_row_label("Canadian rate") |
    p_a | p_b | p_c
) +
  plot_layout(widths = c(0.42, 1, 1, 1))

row_2 <- (
  make_row_label("Canadian rate\n(NL initial)") |
    p_d | p_e | p_f
) +
  plot_layout(widths = c(0.42, 1, 1, 1))

row_3 <- (
  make_row_label("Actual NL\n(Canadian initial)") |
    p_g | p_h | p_i
) +
  plot_layout(widths = c(0.42, 1, 1, 1))

row_4 <- (
  make_row_label("Actual NL") |
    p_j | p_k | p_l
) +
  plot_layout(widths = c(0.42, 1, 1, 1))

figureS1_body <- row_1 / row_2 / row_3 / row_4 +
  plot_layout(heights = c(1, 1, 1, 1))

figureS1 <- figureS1_body &
  theme(plot.background = element_blank())

png(
  "../figures/SI_Fig2.png",
  width = 14500,
  height = 9500,
  res = 300,
  bg = "white",
  type = "cairo"
)

print(figureS1)
dev.off()
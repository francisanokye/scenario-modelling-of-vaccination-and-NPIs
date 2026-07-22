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

scenario_order <- c(
  "Slower Boosters",
  "National Rate",
  "Actual Boosters",
  "Faster Boosters",
  "Early Boosters"
)

policy_order <- c("No NPIs", "Actual NPIs", "Strict NPIs")

scenario_dose_cols <- list(
  "Actual Boosters" = c("≤2 Doses" = "grey70", "Booster" = "#4D4D4D"),
  "Early Boosters"  = c("≤2 Doses" = "#a6cee3", "Booster" = "#1f78b4"),
  "Faster Boosters" = c("≤2 Doses" = "#fdbf6f", "Booster" = "#ff7f00"),
  "Slower Boosters" = c("≤2 Doses" = "#E8B6D2", "Booster" = "magenta"),
  "National Rate"   = c("≤2 Doses" = "#b2df8a", "Booster" = "#33a02c")
)

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

incidence_df <- combosims |>
  filter(matrix %in% c("inc", "inc_booster"), date >= as.Date(omicron_start)) |>
  select(vax_scenario, policy, date, matrix, value) |>
  pivot_wider(names_from = matrix, values_from = value, values_fill = 0) |>
  mutate(
    inc = pmax(inc, 0),
    inc_booster = pmax(inc_booster, 0),
    inc_leq2dose = pmax(inc - inc_booster, 0),
    total_inc = inc_leq2dose + inc_booster,
    vax_scenario = factor(vax_scenario, levels = scenario_order),
    policy = factor(policy, levels = policy_order)
  )

row_ymax <- incidence_df |>
  group_by(vax_scenario) |>
  summarise(ymax = max_or_na(total_inc), .groups = "drop") |>
  mutate(ymax = if_else(is.finite(ymax) & ymax > 0, ymax, 1))

get_panel_cols <- function(scenario_name) {
  cols <- scenario_dose_cols[[scenario_name]]
  if (is.null(cols)) stop("No colours are defined for: ", scenario_name)
  cols
}

make_incidence_panel <- function(
    scenario_name,
    policy_name,
    panel_tag,
    panel_title = NULL,
    show_x_text = FALSE,
    show_y_title = FALSE
) {
  panel_df <- incidence_df |>
    filter(vax_scenario == scenario_name, policy == policy_name)
  
  if (nrow(panel_df) == 0) {
    stop("No incidence data found for ", scenario_name, " under ", policy_name)
  }
  
  y_upper <- row_ymax |>
    filter(vax_scenario == scenario_name) |>
    pull(ymax)
  
  ggplot(panel_df, aes(x = date)) +
    geom_ribbon(aes(ymin = 0, ymax = inc_leq2dose, fill = "≤2 Doses")) +
    geom_ribbon(aes(ymin = inc_leq2dose, ymax = total_inc, fill = "Booster")) +
    geom_line(aes(y = inc_leq2dose), colour = "black", linewidth = 1.5) +
    geom_line(aes(y = total_inc), colour = "black", linewidth = 1.5) +
    geom_vline(
      xintercept = as.Date(intervention_removed),
      linetype = "dashed",
      linewidth = 1.3
    ) +
    scale_fill_manual(
      values = get_panel_cols(scenario_name),
      breaks = c("Booster", "≤2 Doses"),
      name = NULL
    ) +
    scale_x_date(
      date_breaks = "1 month",
      date_labels = "%d %b %Y",
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.04))) +
    coord_cartesian(ylim = c(0, y_upper * 1.06), clip = "off") +
    labs(
      title = panel_title,
      x = NULL,
      y = if (show_y_title) "Daily incidence" else NULL,
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

make_row_label <- function(label) {
  ggplot() +
    annotate(
      "text",
      x = 1,
      y = 0.5,
      label = label,
      size = 13,
      fontface = "plain",
      hjust = 1
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void() +
    theme(plot.margin = margin(0, 16, 0, 10))
}

p_a <- make_incidence_panel("Slower Boosters", "No NPIs", "(A)", "No NPIs", FALSE, TRUE)
p_b <- make_incidence_panel("Slower Boosters", "Actual NPIs", "(B)", "Actual NPIs")
p_c <- make_incidence_panel("Slower Boosters", "Strict NPIs", "(C)", "Strict NPIs")

p_d <- make_incidence_panel("National Rate", "No NPIs", "(D)", show_y_title = TRUE)
p_e <- make_incidence_panel("National Rate", "Actual NPIs", "(E)")
p_f <- make_incidence_panel("National Rate", "Strict NPIs", "(F)")

p_g <- make_incidence_panel("Actual Boosters", "No NPIs", "(G)", show_y_title = TRUE)
p_h <- make_incidence_panel("Actual Boosters", "Actual NPIs", "(H)")
p_i <- make_incidence_panel("Actual Boosters", "Strict NPIs", "(I)")

p_j <- make_incidence_panel("Faster Boosters", "No NPIs", "(J)", show_y_title = TRUE)
p_k <- make_incidence_panel("Faster Boosters", "Actual NPIs", "(K)")
p_l <- make_incidence_panel("Faster Boosters", "Strict NPIs", "(L)")

p_m <- make_incidence_panel("Early Boosters", "No NPIs", "(M)", show_x_text = TRUE, show_y_title = TRUE)
p_n <- make_incidence_panel("Early Boosters", "Actual NPIs", "(N)", show_x_text = TRUE)
p_o <- make_incidence_panel("Early Boosters", "Strict NPIs", "(O)", show_x_text = TRUE)

row_1 <- (make_row_label("Slower Boosters") | p_a | p_b | p_c) +
  plot_layout(widths = c(0.32, 1, 1, 1))
row_2 <- (make_row_label("National Rate") | p_d | p_e | p_f) +
  plot_layout(widths = c(0.32, 1, 1, 1))
row_3 <- (make_row_label("Actual Boosters") | p_g | p_h | p_i) +
  plot_layout(widths = c(0.32, 1, 1, 1))
row_4 <- (make_row_label("Faster Boosters") | p_j | p_k | p_l) +
  plot_layout(widths = c(0.32, 1, 1, 1))
row_5 <- (make_row_label("Early Boosters") | p_m | p_n | p_o) +
  plot_layout(widths = c(0.32, 1, 1, 1))


figure2_body <- row_1 / row_2 / row_3 / row_4 / row_5 +
  plot_layout(heights = c(1, 1, 1, 1, 1))

figure2 <- figure2_body &
  theme(plot.background = element_blank())

png(
  "../figures/figure2.png",
  width = 15500,
  height = 11500,
  res = 300,
  bg = "white",
  type = "cairo"
)
print(figure2)
dev.off()
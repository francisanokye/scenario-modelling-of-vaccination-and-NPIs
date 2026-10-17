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

# ------------------------------------------------------------
# Scenario definitions
# ------------------------------------------------------------

scenario_order <- c("national rate", "national rate (nl init)", "actual (can init)", "actual")
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
primary_target_date <- as.Date("2021-11-06")
plot_start <- primary_target_date
plot_end <- as.Date(analysis_end)

vaccination_matrices <- c(
  "S1", "E1", "A1", "I1", "R1",
  "V2", "E2", "A2", "I2", "R2",
  "V3", "E3", "A3", "I3", "R3"
)

# ------------------------------------------------------------
# Combine simulations
# ------------------------------------------------------------

combosims <- bind_rows(rdsReadList()) |>
  mutate(
    vax_scenario = recode(as.character(vac),
                          "vax_national_rate" = "national rate",
                          "nl_init_national_rate" = "national rate (nl init)",
                          "vax_actual_rate_can_init" = "actual (can init)",
                          "vax_actual_boosters" = "actual"
    ),
    policy = recode(as.character(beta),
                    "no_npis" = "No NPIs",
                    "base" = "Actual NPIs",
                    "strict_npis" = "Strict NPIs"
    ),
    date = as.Date(burn_in_start) + (time - 1)
  ) |>
  filter(vax_scenario %in% scenario_order, policy %in% policy_order)

# ------------------------------------------------------------
# Vaccination-state trajectories
# ------------------------------------------------------------

vaccination_numbers <- combosims |>
  filter(policy == "Actual NPIs", matrix %in% vaccination_matrices,
         date >= plot_start, date <= plot_end) |>
  mutate(
    vaccination_status = case_when(
      matrix %in% c("S1", "E1", "A1", "I1", "R1") ~ "No primary series",
      matrix %in% c("V2", "E2", "A2", "I2", "R2") ~ "Primary series",
      matrix %in% c("V3", "E3", "A3", "I3", "R3") ~ "Booster dose"
    )
  ) |>
  group_by(vax_scenario, date, vaccination_status) |>
  summarise(individuals = sum(value, na.rm = TRUE), .groups = "drop") |>
  mutate(
    vax_scenario = factor(vax_scenario, levels = scenario_order),
    vaccination_status = factor(vaccination_status,
                                levels = c("No primary series", "Primary series", "Booster dose"))
  )

population_check <- vaccination_numbers |>
  group_by(vax_scenario, date) |>
  summarise(total = sum(individuals), .groups = "drop") |>
  mutate(difference = total - NL_POP)

if (max(abs(population_check$difference), na.rm = TRUE) > 1) {
  warning("Vaccination strata do not sum to NL_POP. Maximum difference = ",
          round(max(abs(population_check$difference), na.rm = TRUE), 2))
}

vaccination_ymax <- max(vaccination_numbers$individuals, na.rm = TRUE)

# ------------------------------------------------------------
# Vaccination panels A-C
# ------------------------------------------------------------

make_vaccination_panel <- function(status_name, panel_tag, panel_title, y_label = NULL) {
  panel_df <- vaccination_numbers |> filter(vaccination_status == status_name)
  
  ggplot(panel_df, aes(x = date, y = individuals, colour = vax_scenario, group = vax_scenario)) +
    geom_line(linewidth = 3.5) +
    scale_colour_manual(values = scenario_cols, breaks = scenario_order,
                        labels = unname(vax_legend_labs[scenario_order]), drop = FALSE, name = NULL) +
    scale_x_date(date_breaks = "1 month", date_labels = "%d %b %Y",
                 expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(labels = comma, limits = c(0, vaccination_ymax * 1.04),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(title = panel_title, x = NULL, y = y_label, tag = panel_tag) +
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
# Severe outcomes
# ------------------------------------------------------------

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

severe_hosp_max <- max(severe_df$severe_outcome[severe_df$outcome == "Hospitalizations"], na.rm = TRUE)
severe_death_max <- max(severe_df$severe_outcome[severe_df$outcome == "Deaths"], na.rm = TRUE)
severe_scale_factor <- if (is.finite(severe_death_max) && severe_death_max > 0) severe_hosp_max / severe_death_max else 1
severe_ymax <- severe_hosp_max * 1.12

# ------------------------------------------------------------
# Severe-outcome panels D-F
# ------------------------------------------------------------

make_severe_panel <- function(policy_name, panel_tag) {
  panel_df <- severe_df |>
    filter(policy == policy_name) |>
    select(vax_scenario, outcome, severe_outcome) |>
    complete(
      vax_scenario = factor(scenario_order, levels = scenario_order),
      outcome = factor(c("Hospitalizations", "Deaths"), levels = c("Hospitalizations", "Deaths")),
      fill = list(severe_outcome = 0)
    ) |>
    mutate(
      x = as.numeric(vax_scenario),
      x_plot = x + if_else(outcome == "Hospitalizations", -0.21, 0.21),
      plotted_value = if_else(outcome == "Deaths", severe_outcome * severe_scale_factor, severe_outcome),
      value_label = comma(round(severe_outcome))
    )
  
  hosp_df <- panel_df |> filter(outcome == "Hospitalizations")
  death_df <- panel_df |> filter(outcome == "Deaths")
  
  ggplot(panel_df, aes(x = x_plot, y = plotted_value)) +
    geom_col(data = hosp_df, aes(fill = vax_scenario),
             width = 0.34, colour = "black", linewidth = 0.45, alpha = 1) +
    geom_col(data = death_df, aes(fill = vax_scenario),
             width = 0.30, colour = "black", linewidth = 0.65, alpha = 0.35) +
    
    # Separate hospitalization and death labels
    geom_text(data = hosp_df,
              aes(x = x_plot - 0.035, y = plotted_value, label = value_label),
              inherit.aes = FALSE, hjust = 0.5, vjust = -0.45, size = 15) +
    geom_text(data = death_df,
              aes(x = x_plot + 0.035, y = plotted_value, label = value_label),
              inherit.aes = FALSE, hjust = 0.5, vjust = -0.45, size = 15) +
    
    scale_fill_manual(values = scenario_cols, breaks = scenario_order, drop = FALSE) +
    scale_x_continuous(
      breaks = seq_along(scenario_order),
      labels = unname(vax_short_labs[scenario_order]),
      limits = c(0.55, 4.45),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      name = "Hospitalizations",
      labels = comma,
      limits = c(0, severe_ymax),
      expand = expansion(mult = c(0, 0)),
      sec.axis = sec_axis(~ . / severe_scale_factor, name = "Deaths", labels = comma)
    ) +
    labs(x = NULL, title = policy_name, tag = panel_tag) +
    guides(fill = "none") +
    coord_cartesian(clip = "off") +
    theme_clean() +
    theme(
      axis.text.x = element_text(size = 44, hjust = 0.5, lineheight = 0.90),
      axis.text.y = element_text(size = 50),
      axis.title.x = element_text(size = 50),
      axis.title.y = element_text(size = 50),
      plot.title = element_text(size = 50, hjust = 0.5, face = "plain"),
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

p_a <- make_vaccination_panel("No primary series", "(A)", "No primary series", "Number of individuals")
p_b <- make_vaccination_panel("Primary series", "(B)", "Primary series")
p_c <- make_vaccination_panel("Booster dose", "(C)", "Booster doses")

p_d <- make_severe_panel("No NPIs", "(D)")
p_e <- make_severe_panel("Actual NPIs", "(E)")
p_f <- make_severe_panel("Strict NPIs", "(F)")

# ------------------------------------------------------------
# Vaccination scenario legend
# ------------------------------------------------------------

scenario_legend_grob <- grobTree(
  segmentsGrob(
    x0 = unit(c(0.07, 0.29, 0.54, 0.79), "npc"),
    x1 = unit(c(0.105, 0.325, 0.575, 0.825), "npc"),
    y0 = unit(rep(0.5, 4), "npc"),
    y1 = unit(rep(0.5, 4), "npc"),
    gp = gpar(col = unname(scenario_cols[scenario_order]), lwd = 9)
  ),
  textGrob(
    label = unname(vax_legend_labs[scenario_order]),
    x = unit(c(0.112, 0.332, 0.582, 0.832), "npc"),
    y = unit(rep(0.5, 4), "npc"),
    just = "left",
    gp = gpar(fontsize = 50)
  )
)

scenario_legend <- wrap_elements(full = scenario_legend_grob)

# ------------------------------------------------------------
# Severe-outcome legend
# ------------------------------------------------------------

outcome_legend_grob <- grobTree(
  rectGrob(
    x = unit(c(0.425, 0.565), "npc"),
    y = unit(c(0.5, 0.5), "npc"),
    width = unit(c(0.014, 0.014), "npc"),
    height = unit(c(0.28, 0.28), "npc"),
    gp = gpar(fill = c("grey35", "grey75"), col = "black", lwd = 1)
  ),
  textGrob(
    label = c("Hospitalizations", "Deaths"),
    x = unit(c(0.442, 0.582), "npc"),
    y = unit(c(0.5, 0.5), "npc"),
    just = "left",
    gp = gpar(fontsize = 50)
  )
)

outcome_legend <- wrap_elements(full = outcome_legend_grob)

# ------------------------------------------------------------
# Assemble figure
# ------------------------------------------------------------

vaccination_row <- (p_a | p_b | p_c) + plot_layout(widths = c(1, 1, 1))
severe_row <- (p_d | p_e | p_f) + plot_layout(widths = c(1, 1, 1))

figure4 <- scenario_legend /
  vaccination_row /
  outcome_legend /
  severe_row +
  plot_layout(widths = 1, heights = c(0.10, 0.90, 0.09, 1.08))

# ------------------------------------------------------------
# Export
# ------------------------------------------------------------

png("../figures/SI_Fig5.png", width = 17500, height = 8000, res = 300, bg = "white", type = "cairo")
print(figure4)
dev.off()
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
  "national rate (nl init)" = "Canadian rate\n(NL initial)",
  "actual (can init)" = "Actual NL\n(Canadian initial)",
  "actual" = "Actual\nNL"
)

study_start <- as.Date(omicron_start)
intervention_date <- as.Date(intervention_removed)

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

cumulative_infections_df <- combosims |>
  filter(matrix %in% c("cum_inc_total", "cum_inc_booster")) |>
  select(vax_scenario, policy, date, matrix, value) |>
  pivot_wider(names_from = matrix, values_from = value, values_fill = 0)

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
  filter(date >= study_start) |>
  mutate(
    vax_scenario = factor(vax_scenario, levels = scenario_order),
    policy = factor(policy, levels = policy_order)
  )

no_booster_ymax <- max(cumulative_infections_df$cum_inc_leq2dose, na.rm = TRUE)
booster_ymax <- max(cumulative_infections_df$cum_inc_booster, na.rm = TRUE)

if (!is.finite(no_booster_ymax) || no_booster_ymax <= 0) no_booster_ymax <- 1
if (!is.finite(booster_ymax) || booster_ymax <= 0) booster_ymax <- 1

make_cumulative_panel <- function(policy_name, outcome_name, panel_tag, panel_title, y_label, y_max) {
  
  panel_df <- cumulative_infections_df |>
    filter(policy == policy_name)
  
  ggplot(panel_df, aes(
    x = date,
    y = .data[[outcome_name]],
    colour = vax_scenario,
    group = vax_scenario
  )) +
    geom_line(linewidth = 3.5) +
    geom_vline(
      xintercept = intervention_date,
      linetype = "dashed",
      linewidth = 2
    ) +
    scale_colour_manual(
      values = scenario_cols,
      breaks = scenario_order,
      labels = unname(vax_short_labs[scenario_order]),
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
      limits = c(0, y_max * 1.04),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = panel_title,
      x = NULL,
      y = y_label,
      tag = panel_tag
    ) +
    theme_clean() +
    theme(
      axis.text.x = element_text(size = 40, hjust = 1, angle = 45),
      axis.text.y = element_text(size = 40),
      axis.title.x = element_text(size = 40),
      axis.title.y = element_text(size = 45),
      plot.title = element_text(size = 45, hjust = 0.5, face = "plain"),
      plot.subtitle = element_text(size = 45, hjust = 0.5),
      legend.position = "none",
      legend.justification = "center",
      legend.direction = "horizontal",
      legend.text = element_text(size = 45),
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

p_a <- make_cumulative_panel(
  "No NPIs",
  "cum_inc_leq2dose",
  "(A)",
  "No NPIs",
  "Cumulative infections\n(no booster dose)",
  no_booster_ymax
)

p_b <- make_cumulative_panel(
  "Actual NPIs",
  "cum_inc_leq2dose",
  "(B)",
  "Actual NPIs",
  NULL,
  no_booster_ymax
)

p_c <- make_cumulative_panel(
  "Strict NPIs",
  "cum_inc_leq2dose",
  "(C)",
  "Strict NPIs",
  NULL,
  no_booster_ymax
)

p_d <- make_cumulative_panel(
  "No NPIs",
  "cum_inc_booster",
  "(D)",
  "No NPIs",
  "Cumulative infections\n(booster doses)",
  booster_ymax
)

p_e <- make_cumulative_panel(
  "Actual NPIs",
  "cum_inc_booster",
  "(E)",
  "Actual NPIs",
  NULL,
  booster_ymax
)

p_f <- make_cumulative_panel(
  "Strict NPIs",
  "cum_inc_booster",
  "(F)",
  "Strict NPIs",
  NULL,
  booster_ymax
)

scenario_legend_grob <- grobTree(
  segmentsGrob(
    x0 = unit(c(0.10, 0.32, 0.55, 0.79), "npc"),
    x1 = unit(c(0.135, 0.355, 0.585, 0.825), "npc"),
    y0 = unit(rep(0.5, 4), "npc"),
    y1 = unit(rep(0.5, 4), "npc"),
    gp = gpar(
      col = unname(scenario_cols[scenario_order]),
      lwd = 10
    )
  ),
  textGrob(
    label = gsub(
      "\n",
      " ",
      unname(vax_short_labs[scenario_order])
    ),
    x = unit(c(0.142, 0.362, 0.592, 0.832), "npc"),
    y = unit(rep(0.5, 4), "npc"),
    just = "left",
    gp = gpar(fontsize = 45)
  )
)

scenario_legend <- wrap_elements(full = scenario_legend_grob)

no_booster_row <- (p_a | p_b | p_c) +
  plot_layout(widths = c(1, 1, 1))

booster_row <- (p_d | p_e | p_f) +
  plot_layout(widths = c(1, 1, 1))

figureS2 <- scenario_legend /
  no_booster_row /
  booster_row +
  plot_layout(
    widths = 1,
    heights = c(0.05, 1, 1)
  )

png("../figures/SI_Fig4.png",width = 17500,height = 9000,res = 300,bg = "white",type = "cairo")
print(figureS2)
dev.off()
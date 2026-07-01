suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggthemes)
  library(shellpipes)
  library(scales)
  library(lubridate)
})

options(dplyr.summarise.inform = FALSE)
options(scipen = 9999)

loadEnvironments()

vax_short_labs <- c(
  "National Rate"   = "National\nRate",
  "Slower Boosters" = "Slower\nBoosters",
  "Actual Boosters" = "Actual\nBoosters",
  "Early Boosters"  = "Early\nBoosters",
  "Faster Boosters" = "Faster\nBoosters"
)

figure2_df <- bind_rows(
  results_df |>
    filter(date >= omicron_start, matrix == "cum_hosp_total") |>
    group_by(policy, vax_scenario) |>
    summarise(severe_outcome = max_or_na(value), .groups = "drop") |>
    mutate(outcome = "Hospitalizations"),
  
  results_df |>
    filter(date >= omicron_start, matrix == "cum_death_total") |>
    group_by(policy, vax_scenario) |>
    summarise(severe_outcome = max_or_na(value), .groups = "drop") |>
    mutate(outcome = "Deaths")
) |>
  filter(is.finite(severe_outcome)) |>
  mutate(
    policy_label = recode(as.character(policy),
                          "No-ALS Only" = "No NPIs",
                          "Actual ALS"  = "Actual NPIs",
                          "ALS-4 Only"  = "Strict NPIs"),
    policy_label = factor(policy_label, levels = c("No NPIs", "Actual NPIs", "Strict NPIs")),
    outcome      = factor(outcome, levels = c("Hospitalizations", "Deaths")),
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
    x_lab       = vax_short_labs[as.character(vax_scenario)],
    value_label = comma(round(severe_outcome, 0))
  ) |>
  filter(!is.na(x_pos))

figure2 <- ggplot(figure2_df, aes(x = x_pos, y = severe_outcome, fill = vax_scenario)) +
  geom_col(colour = "black", width = 0.72, linewidth = 0.25) +
  geom_text(aes(label = value_label), vjust = -0.35, size = 10) +
  facet_grid(outcome ~ policy_label, scales = "free_y") +
  scale_fill_manual(values = scenario_cols, breaks = vax_order, drop = FALSE) +
  scale_x_continuous(breaks = 1:5, labels = NULL, expand = expansion(mult = c(0.04, 0.04))) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0.08, 0.18))) +
  labs(
    title = "Severe outcomes across NPI and booster scenarios",
    x = NULL, y = "Severe outcomes", fill = NULL
  ) +
  theme_clean() +
  theme(
    axis.text.x        = element_text(size = 30),
    axis.ticks.x       = element_blank(),
    axis.text.y        = element_text(size = 35),
    axis.title.y       = element_text(size = 35, colour = "black"),
    strip.text         = element_text(size = 35),
    strip.background   = element_blank(),
    legend.position    = "top",
    legend.text        = element_text(size = 25),
    legend.background  = element_blank(),
    legend.spacing.x   = unit(0.45, "cm"),
    plot.title         = element_text(size = 40, hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.spacing      = unit(2, "cm"),
    plot.margin        = margin(t = 10, r = 20, b = 55, l = 20),
    plot.background    = element_blank()
  ) +
  geom_text(
    data = figure2_df |> filter(outcome == "Deaths"),
    aes(y = -Inf, label = x_lab),
    vjust = 1.15, size = 8, lineheight = 0.85, inherit.aes = TRUE
  ) +
  coord_cartesian(clip = "off")

save_png(file.path(fig_dir, "figure3.png"), figure2, width = 10000, height = 5000)
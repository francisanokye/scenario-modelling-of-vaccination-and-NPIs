suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(scales)
  library(ggthemes)
  library(shellpipes)
  library(grid)
  library(patchwork)
  library(stringr)
})

loadEnvironments()

# ============================================================
# Load scenario objects
# ============================================================

actual_daily   <- scenario_daily_df_actual
actual_weekly  <- scenario_weekly_df_actual
actual_inputs  <- scenario_booster_inputs_df_actual

early_daily    <- scenario_daily_df_early
early_weekly   <- scenario_weekly_df_early
early_inputs   <- scenario_booster_inputs_df_early

faster_daily   <- scenario_daily_df_faster
faster_weekly  <- scenario_weekly_df_faster
faster_inputs  <- scenario_booster_inputs_df_faster

slower_daily   <- scenario_daily_df_slower
slower_weekly  <- scenario_weekly_df_slower
slower_inputs  <- scenario_booster_inputs_df_slower

national_daily  <- scenario_daily_df_national
national_weekly <- scenario_weekly_df_national
national_inputs <- scenario_booster_inputs_df_national

# ============================================================
# Dates
# ============================================================

analysis_start <- as.Date("2021-12-15")
offset <- 150L
init_date <- analysis_start - offset

sim_horizon <- 530L
analysis_end <- init_date + sim_horizon - 1L

actual_booster_start <- as.Date("2021-11-08")

scenario_name_early <- unique(early_weekly$scenario)
stopifnot(length(scenario_name_early) == 1)

early_shift_days <- str_extract(
  scenario_name_early,
  "-?\\d+(?= Days\\))"
) |>
  as.integer()

if (is.na(early_shift_days)) {
  early_shift_days <- -35L
}

early_start_date <- actual_booster_start + days(early_shift_days)

# ============================================================
# Standardize display labels
# ============================================================

rename_scenario <- function(x) {
  x |>
    mutate(
      scenario = case_when(
        scenario == "Early Boosters (-35 Days)" ~ "Early Boosters",
        TRUE ~ scenario
      )
    )
}

actual_daily    <- rename_scenario(actual_daily)
actual_weekly   <- rename_scenario(actual_weekly)
actual_inputs   <- rename_scenario(actual_inputs)

early_daily     <- rename_scenario(early_daily)
early_weekly    <- rename_scenario(early_weekly)
early_inputs    <- rename_scenario(early_inputs)

faster_daily    <- rename_scenario(faster_daily)
faster_weekly   <- rename_scenario(faster_weekly)
faster_inputs   <- rename_scenario(faster_inputs)

slower_daily    <- rename_scenario(slower_daily)
slower_weekly   <- rename_scenario(slower_weekly)
slower_inputs   <- rename_scenario(slower_inputs)

national_daily  <- rename_scenario(national_daily)
national_weekly <- rename_scenario(national_weekly)
national_inputs <- rename_scenario(national_inputs)

# ============================================================
# Checks
# ============================================================

stopifnot(
  exists("NL_POP"),
  exists("eligible_pop"),
  init_date + offset == analysis_start,
  analysis_start + (0 - offset) == init_date,
  early_start_date >= init_date
)

# ============================================================
# Add missing logical columns if absent
# ============================================================

for (obj in c(
  "actual_weekly",
  "early_weekly",
  "faster_weekly",
  "slower_weekly",
  "national_weekly"
)) {
  
  x <- get(obj)
  
  if (!"is_target" %in% names(x)) x$is_target <- FALSE
  if (!"is_interp" %in% names(x)) x$is_interp <- FALSE
  if (!"is_phac" %in% names(x)) x$is_phac <- FALSE
  
  assign(obj, x)
}

# ============================================================
# Combine daily data
# ============================================================

daily_scenarios_df <- bind_rows(
  actual_daily,
  early_daily,
  faster_daily,
  slower_daily,
  national_daily
) |>
  mutate(date = as.Date(date)) |>
  filter(
    date >= init_date,
    date <= analysis_end
  ) |>
  arrange(scenario, date)

# ============================================================
# Combine weekly data
# ============================================================

weekly_scenarios_df <- bind_rows(
  actual_weekly,
  early_weekly,
  faster_weekly,
  slower_weekly,
  national_weekly
) |>
  mutate(week_end = as.Date(week_end)) |>
  filter(
    week_end >= init_date,
    week_end <= analysis_end
  ) |>
  arrange(scenario, week_end) |>
  group_by(scenario, week_end) |>
  summarise(
    weekly_rate = sum(replace_na(weekly_rate, 0), na.rm = TRUE),
    cum_val = max(replace_na(cum_val, 0), na.rm = TRUE),
    is_target = any(replace_na(is_target, FALSE)),
    is_interp = any(replace_na(is_interp, FALSE)),
    is_phac = any(replace_na(is_phac, FALSE)),
    .groups = "drop"
  ) |>
  group_by(scenario) |>
  arrange(week_end, .by_group = TRUE) |>
  mutate(
    cum_val = pmax(cum_val, 0),
    cum_val = pmin(cum_val, eligible_pop),
    cum_val = cummax(cum_val),
    cum_boosters = cum_val,
    cum_prop = clip01(cum_boosters / NL_POP)
  ) |>
  ungroup()

# ============================================================
# Booster inputs
# ============================================================

booster_inputs_df <- bind_rows(
  actual_inputs,
  early_inputs,
  faster_inputs,
  slower_inputs,
  national_inputs
) |>
  mutate(date = as.Date(date)) |>
  filter(
    date >= init_date,
    date <= analysis_end
  ) |>
  arrange(scenario, day)


# ============================================================
# Plot theme
# ============================================================

theme_pub <- theme_clean() +
  theme(
    axis.text.x = element_text(size = 25, hjust = 1, angle = 45),
    axis.text.y = element_text(size = 25),
    axis.title.y = element_text(size = 30, color = "black"),
    legend.position = c(0.88, 0.22),
    legend.justification = "center",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.title = element_blank(),
    legend.text = element_text(size = 25),
    legend.background = element_blank(),
    legend.margin = margin(t = 8, r = 8, b = 8, l = 8),
    legend.key.height = unit(1.2, "cm"),
    legend.spacing.y = unit(1.5, "cm"),
    plot.title = element_text(size = 35, face = "plain", hjust = 0.5),
    plot.background = element_blank(),
    plot.margin = margin(t = 25, r = 35, b = 90, l = 35),
    plot.tag = element_text(size = 35, face = "bold", colour = "black"),
    plot.tag.position = c(0.01, 0.98)
  )

# ============================================================
# Legend order
# ============================================================

vax_order <- c(
  "Actual Boosters",
  "Early Boosters",
  "Faster Rate",
  "Slower Rate",
  "National Rate"
)

# ============================================================
# Colors
# ============================================================

vax_cols <- c(
  "Actual Boosters" = "#4D4D4D",
  "Early Boosters"  = "blue",
  "Faster Rate"     = "orange",
  "Slower Rate"     = "magenta",
  "National Rate"   = "springgreen3"
)

# ============================================================
# Plot data
# ============================================================

plot_data <- weekly_scenarios_df |>
  filter(week_end <= analysis_end) |>
  mutate(
    scenario = factor(scenario, levels = vax_order)
  )

# ============================================================
# HOA target point
# ============================================================

target_date <- as.Date("2022-01-11")
target_cumulative <- 106730

target_point <- tibble(
  week_end = target_date,
  cum_prop = target_cumulative / NL_POP
)

# ============================================================
# Plot limits
# ============================================================

x_min <- early_start_date - days(2)
x_max <- as.Date("2022-06-01")

plot_data <- plot_data |>
  filter(
    week_end >= x_min,
    week_end <= x_max
  )

y_max <- max(
  plot_data$cum_prop,
  target_point$cum_prop,
  na.rm = TRUE
)

y_upper <- ceiling((y_max + 0.03) * 10) / 10

# ============================================================
# Plot
# ============================================================

pp <- ggplot(plot_data) +
  geom_line(
    aes(x = week_end, y = cum_prop, color = scenario),
    linewidth = 1.5,
    na.rm = TRUE
  ) +
  geom_point(
    data = plot_data |> filter(scenario != "Actual Boosters"),
    aes(x = week_end, y = cum_prop, color = scenario),
    size = 7,
    shape = 16,
    na.rm = TRUE
  ) +
  geom_point(
    data = plot_data |>
      filter(scenario == "Actual Boosters", is_phac),
    aes(x = week_end, y = cum_prop),
    size = 7,
    shape = 16,
    color = "black",
    na.rm = TRUE
  ) +
  geom_point(
    data = plot_data |>
      filter(scenario == "Actual Boosters", is_interp),
    aes(x = week_end, y = cum_prop),
    size = 7,
    shape = 17,
    color = "black",
    na.rm = TRUE
  ) +
  geom_point(
    data = target_point,
    aes(x = week_end, y = cum_prop),
    shape = 8,
    size = 12,
    color = "red",
    stroke = 2
  ) +
  geom_text(
    data = target_point,
    aes(x = week_end, y = cum_prop, label = "HOA"),
    color = "red",
    nudge_x = 9,
    size = 7,
    fontface = "bold"
  ) +
  scale_colour_manual(
    values = vax_cols,
    breaks = vax_order,
    drop = FALSE
  ) +
  scale_x_date(
    limits = c(x_min, x_max),
    date_breaks = "1 month",
    date_labels = "%b %Y",
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.1),
    breaks = scales::breaks_pretty(n = 5),
    limits = c(0, y_upper),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    title = "Proportions of booster vaccination coverage in NL",
    y = "Proportion of total population",
    x = NULL,
    color = NULL,
    tag = "B"
  ) +
  theme_pub +
  guides(
    color = guide_legend(ncol = 1)
  )

# ============================================================
# Save figure
# ============================================================

dir.create(
  "../figures",
  showWarnings = FALSE,
  recursive = TRUE
)

png(
  "../figures/cumulative_vax_rate.png",
  width = 5000,
  height = 2800,
  res = 300,
  bg = "white"
)

print(pp)

dev.off()


# ============================================================
# Plot daily booster-dose rates used in the model
# ============================================================

daily_rate_plot_df <- booster_inputs_df |>
  mutate(
    scenario = recode(
      scenario,
      "Early Boosters (-35 Days)" = "Earlier boosters",
      "Early Boosters" = "Earlier boosters",
      "Actual Boosters" = "Actual boosters",
      "Slower Rate" = "Slower boosters",
      "Faster Rate" = "Faster boosters",
      "National Rate" = "National rate"
    )
  ) |>
  filter(
    date >= x_min,
    date <= x_max
  ) |>
  mutate(
    scenario = factor(
      scenario,
      levels = c(
        "Actual boosters",
        "Earlier boosters",
        "Faster boosters",
        "Slower boosters",
        "National rate"
      )
    )
  )

daily_rate_cols <- c(
  "Actual boosters"  = "#4D4D4D",
  "Earlier boosters" = "blue",
  "Faster boosters"  = "orange",
  "Slower boosters"  = "magenta",
  "National rate"    = "springgreen3"
)

daily_rate_plot <- ggplot(
  daily_rate_plot_df,
  aes(
    x = date,
    y = booster_vac_value,
    colour = scenario
  )
) +
  geom_line(
    linewidth = 2,
    na.rm = TRUE
  ) +
  scale_colour_manual(
    values = daily_rate_cols,
    drop = FALSE
  ) +
  scale_x_date(
    limits = c(x_min, x_max),
    date_breaks = "1 month",
    date_labels = "%b %Y",
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    title = "Daily booster-dose rates used in the model",
    y = "Daily booster doses",
    x = NULL,
    colour = NULL
  ) +
  theme_clean() +
  theme(
    axis.text.x = element_text(
      size = 30,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 30),
    axis.title.y = element_text(size = 35),
    legend.position = c(0.83, 0.80),
    legend.justification = "center",
    legend.direction = "vertical",
    legend.text = element_text(size = 28),
    legend.background = element_blank(),
    legend.title = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(
      size = 38,
      hjust = 0.5
    ),
    plot.background = element_blank(),
    plot.margin = margin(15, 35, 15, 15)
  )

png(
  "../figures/daily_booster_rates_used_in_model.png",
  width = 5000,
  height = 2800,
  res = 300,
  bg = "white"
)
print(daily_rate_plot)

dev.off()

# ============================================================
# Save variables
# ============================================================

saveVars(
  weekly_scenarios_df,
  daily_scenarios_df,
  booster_inputs_df
)
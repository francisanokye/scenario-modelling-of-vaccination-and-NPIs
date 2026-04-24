
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
})

loadEnvironments()

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

# ensure flags exist for all scenarios
for (obj in c("actual_weekly", "early_weekly", "faster_weekly", "slower_weekly", "national_weekly")) {
  x <- get(obj)
  if (!"is_target" %in% names(x)) x$is_target <- FALSE
  if (!"is_interp" %in% names(x)) x$is_interp <- FALSE
  if (!"is_phac" %in% names(x)) x$is_phac <- FALSE
  assign(obj, x)
}

daily_scenarios_df <- bind_rows(
  actual_daily,
  early_daily,
  faster_daily,
  slower_daily,
  national_daily
)

weekly_scenarios_df <- bind_rows(
  actual_weekly,
  early_weekly,
  faster_weekly,
  slower_weekly,
  national_weekly
) |>
  arrange(scenario, week_end) |>
  group_by(scenario, week_end) |>
  summarise(
    weekly_rate = sum(replace_na(weekly_rate, 0), na.rm = TRUE),
    cum_prop    = max(replace_na(cum_prop, 0), na.rm = TRUE),
    cum_val     = max(replace_na(cum_val, 0), na.rm = TRUE),
    is_target   = any(replace_na(is_target, FALSE)),
    is_interp   = any(replace_na(is_interp, FALSE)),
    is_phac     = any(replace_na(is_phac, FALSE)),
    .groups = "drop"
  ) |>
  group_by(scenario) |>
  arrange(week_end, .by_group = TRUE) |>
  mutate(
    cum_val      = pmax(cum_val, 0),
    cum_val      = pmin(cum_val, eligible_pop),
    cum_val      = cummax(cum_val),
    cum_boosters = cum_val,
    cum_prop     = clip01(cum_boosters / eligible_pop),
    weekly_rate  = c(cum_boosters[1], diff(cum_boosters)),
    weekly_rate  = pmax(weekly_rate, 0)
  ) |>
  ungroup()

booster_inputs_df <- bind_rows(
  actual_inputs,
  early_inputs,
  faster_inputs,
  slower_inputs,
  national_inputs
) |>
  arrange(scenario, day)

dup_check <- weekly_scenarios_df |>
  count(scenario, week_end) |>
  filter(n > 1)

stopifnot(nrow(dup_check) == 0)
stopifnot(all(weekly_scenarios_df$cum_prop >= 0))

theme_pub <- theme_clean() +
  theme(
    axis.text.x = element_text(size = 25, hjust = 1, angle = 45),
    axis.text.y = element_text(size = 25),
    axis.title.y = element_text(size = 30, color = "black"),
    legend.position = c(0.8, 0.35),
    legend.justification = "center",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.title = element_blank(),
    legend.text = element_text(size = 25),
    legend.background = element_blank(),
    legend.margin = margin(t = 8, r = 8, b = 8, l = 8),
    legend.key.height = unit(1.2, "cm"),   
    legend.spacing.y  = unit(1.5, "cm"), 
    plot.title = element_text(size = 35, face = "plain", hjust = 0.5),
    plot.background = element_blank()
  ) 

vax_cols <- c(
  "Slower Rate"               = "magenta",
  "Early Boosters (-35 Days)" = "blue",
  "Actual Boosters"           = "black",
  "Faster Rate"               = "orange",
  "National Rate"             = "springgreen3"
)

plot_data <- weekly_scenarios_df |>
  filter(week_end <= as.Date("2022-05-22"))

# force the HOA assumed point directly onto the plot
target_date       <- as.Date("2022-01-11")
target_cumulative <- 106730

target_point <- tibble(
  week_end = target_date,
  cum_prop = target_cumulative / eligible_pop
)

pp <- ggplot(plot_data) +
  geom_line(aes(x = week_end, y = cum_prop, color = scenario), linewidth = 1.5) +
  
  # all non-actual scenario points
  geom_point(
    data = plot_data |> filter(scenario != "Actual Boosters"),
    aes(x = week_end, y = cum_prop, color = scenario),
    size = 7,
    shape = 16
  ) +
  
  # actual PHAC points: black circles
  geom_point(
    data = plot_data |> filter(scenario == "Actual Boosters", is_phac),
    aes(x = week_end, y = cum_prop),
    size = 7,
    shape = 16,
    color = "black"
  ) +
  
  # actual interpolated points: black triangles
  geom_point(
    data = plot_data |> filter(scenario == "Actual Boosters", is_interp),
    aes(x = week_end, y = cum_prop),
    size = 7,
    shape = 17,
    color = "black"
  ) +
  
  # HOA assumed point: red star at exact known location
  geom_point(
    data = target_point,
    aes(x = week_end, y = cum_prop),
    shape = 8,
    size = 12,
    color = "red",
    stroke = 2
  ) +
  
  # optional label to make it easy to spot
  geom_text(
    data = target_point,
    aes(x = week_end, y = cum_prop, label = "HOA"),
    color = "red",
    nudge_x = 9,
    nudge_y = 0.0,
    size = 7,
    fontface = "bold"
  ) +
  
  geom_vline(xintercept = booster_start, color = "maroon", linetype = 2, linewidth = 2) +
  scale_colour_manual(values = vax_cols) +
  scale_x_date(
    limits = c(booster_start, as.Date("2022-05-22")),
    date_breaks = "2 weeks",
    date_labels = "%b %Y",
    expand = c(0.05, 0)
  ) +
  labs(
    title = "Proportions of booster vaccination coverage in NL",
    y = "Proportion of NL eligible population",
    x = NULL
  ) +
  theme_pub +
  guides(color = guide_legend(ncol = 1)) +
  labs(tag = "B") +
  theme(
    plot.tag = element_text(size = 35, face = "bold", colour = "black"),
    plot.tag.position = c(0.01, 0.98)
  )

png("../figures/cumulative_vax_rate.png", width = 5000, height = 2800, res = 300, bg = "white")
print(pp)
dev.off()

saveVars(weekly_scenarios_df, daily_scenarios_df, booster_inputs_df)
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

# max that returns NA (not -Inf) when all values are missing
if (!exists("max_or_na")) {
  max_or_na <- function(x, ...) {
    x <- x[!is.na(x)]
    if (length(x) == 0) NA_real_ else max(x, ...)
  }
}

if (!exists("scenario_cols")) {
  scenario_cols <- c(
    "Actual Boosters" = "#4D4D4D",
    "Early Boosters"  = "#1f78b4",
    "Faster Boosters" = "#ff7f00",
    "Slower Boosters" = "magenta",
    "National Rate"   = "#33a02c"
  )
}

# per-scenario dose-group colors: each element names "≤2 Doses" and "Booster"
if (!exists("scenario_dose_cols")) {
  scenario_dose_cols <- list(
    "Actual Boosters" = c("≤2 Doses" = "grey70", "Booster" = "#4D4D4D"),
    "Early Boosters"  = c("≤2 Doses" = "#a6cee3", "Booster" = "#1f78b4"),
    "Faster Boosters" = c("≤2 Doses" = "#fdbf6f", "Booster" = "#ff7f00"),
    "Slower Boosters" = c("≤2 Doses" = "#E8B6D2",  "Booster" = "magenta"),
    "National Rate"   = c("≤2 Doses" = "#b2df8a", "Booster" = "#33a02c")
  )
}
# "Slower Boosters" = c("≤2 Doses" = "#E8B6D2",  "Booster" = "magenta"),
# "Slower Boosters" = c("≤2 Doses" = "#fb9a99", "Booster" = "#e31a1c"),
combosims <- bind_rows(rdsReadList())

combosims <- (combosims 
  |> mutate(vax_scenario = recode(as.character(vac)
    , "vax_actual_boosters" = "Actual Boosters"
    , "vax_early_boosters"  = "Early Boosters"
    , "vax_faster_rate"     = "Faster Boosters"
    , "vax_slower_rate"     = "Slower Boosters"
    , "vax_national_rate"   = "National Rate"
    )
    , policy = recode(as.character(beta)
    , "strict_npis" = "Strict NPIs"
    , "no_npis"   = "No NPIs"
    , "base" = "Actual"
    )
    , date = burn_in_start+ (time - 1)
  )
)

incidence_df  <- (combosims 
  |> filter(matrix %in% c("inc_booster", "inc_leq2dose", "inc")
    , date >= omicron_start
    )
  )

cumulative_df <- (combosims 
  |> filter(matrix %in% c("cum_inc_total", "cum_inc_leq2dose", "cum_inc_booster")
    , date >= omicron_start
    )
  )

# shared y-axis limits
ab_ymax <- (incidence_df
  |> filter(vax_scenario == "Actual Boosters"
    , policy %in% c("No NPIs", "Strict NPIs")
    , matrix %in% c("inc_leq2dose", "inc_booster")
  )
  |> group_by(policy, date)
  |> summarise(total_inc = sum(value, na.rm = TRUE), .groups = "drop")
  |> summarise(ymax = max_or_na(total_inc))
  |> pull(ymax)
)

gh_ymax <- (incidence_df
  |> filter(policy == "Strict NPIs"
    , vax_scenario %in% c("Early Boosters", "Actual Boosters")
    , matrix %in% c("inc_leq2dose", "inc_booster")
  )
  |> group_by(vax_scenario, date)
  |> summarise(total_inc = sum(value, na.rm = TRUE), .groups = "drop")
  |> summarise(ymax = max_or_na(total_inc))
  |> pull(ymax)
)

de_ymax <- (incidence_df
  |> filter(policy == "Strict NPIs"
    , vax_scenario %in% c("Slower Boosters", "Faster Boosters")
    , matrix %in% c("inc_leq2dose", "inc_booster")
  )
  |> group_by(vax_scenario, date)
  |> summarise(total_inc = sum(value, na.rm = TRUE), .groups = "drop")
  |> summarise(ymax = max_or_na(total_inc))
  |> pull(ymax)
)

incidence_ymax <- max(ab_ymax, de_ymax, gh_ymax, na.rm = TRUE)

# plot helpers ----
get_panel_cols <- function(scenario_name) {
  cols <- scenario_dose_cols[[as.character(scenario_name)]]
  if (is.null(cols)) stop("No dose colours defined for scenario: ", scenario_name)
  cols
}

make_shaded_incidence_plot <- function(
    df,
    scenario_name,
    policy_name,
    panel_tag,
    panel_title,
    legend_position   = c(0.50, 0.85),
    y_limit           = NULL,
    boundary_linetype = "solid",
    show_legend       = TRUE
) {
  panel_cols <- get_panel_cols(scenario_name)
  
  plot_df <- df |>
    filter(
      policy       == policy_name,
      vax_scenario == scenario_name,
      matrix       %in% c("inc", "inc_booster")
    ) |>
    select(date, matrix, value) |>
    pivot_wider(names_from = matrix, values_from = value, values_fill = 0)
  
  if (!"inc"         %in% names(plot_df)) plot_df$inc         <- numeric(nrow(plot_df))
  if (!"inc_booster" %in% names(plot_df)) plot_df$inc_booster <- numeric(nrow(plot_df))
  plot_df$total_inc <- plot_df$inc
  
  y_upper <- if (is.null(y_limit)) max_or_na(plot_df$total_inc) else y_limit
  
  ggplot(plot_df, aes(x = date)) +
    geom_area(aes(y = total_inc,   fill = "≤2 Doses")) +
    geom_area(aes(y = inc_booster, fill = "Booster")) +
    geom_line(aes(y = total_inc),   colour = "black", linewidth = 2, linetype = boundary_linetype) +
    geom_line(aes(y = inc_booster), colour = "black", linewidth = 2, linetype = boundary_linetype) +
    geom_vline(xintercept = intervention_removed, linetype = "dashed", linewidth = 2) +
    scale_fill_manual(values = panel_cols, breaks = c("Booster", "≤2 Doses")) +
    scale_x_date(date_breaks = "1 month", date_labels = "%d %b %Y", expand = c(0.01, 0)) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.04))) +
    coord_cartesian(ylim = c(0, y_upper * 1.06), clip = "off") +
    labs(title = panel_title, x = NULL, y = "Daily incidence", fill = NULL, tag = panel_tag) +
    theme_clean() +
    theme(
      axis.text.x          = element_text(size = 35, hjust = 1, angle = 45),
      axis.text.y          = element_text(size = 35),
      axis.title.y         = element_text(size = 35),
      legend.position      = if (show_legend) legend_position else "none",
      legend.justification = "center",
      legend.direction     = "vertical",
      legend.text          = element_text(size = 30),
      legend.background    = element_blank(),
      plot.title           = element_text(size = 45, hjust = 0.5),
      panel.grid.major.x   = element_blank(),
      plot.tag             = element_text(size = 45, face = "bold"),
      plot.background      = element_blank(),
      plot.margin          = margin(12, 35, 12, 12)
    )
}

make_cumulative_policy_plot <- function(df, scenario_name, panel_tag) {
  panel_cols <- get_panel_cols(scenario_name)
  
  plot_df <- df |>
    filter(
      vax_scenario == scenario_name,
      policy %in% c("No NPIs", "Strict NPIs"),
      matrix %in% c("cum_inc_leq2dose", "cum_inc_booster")
    ) |>
    mutate(
      dose_group = recode(matrix,
                          "cum_inc_leq2dose" = "≤2 Doses",
                          "cum_inc_booster"  = "Booster"),
      policy_label = recode(as.character(policy),
                            "No NPIs" = "No NPIs",
                            "Strict NPIs"  = "Strict NPIs"),
      dose_group   = factor(dose_group,   levels = c("≤2 Doses", "Booster")),
      policy_label = factor(policy_label, levels = c("No NPIs", "Strict NPIs"))
    )
  
  ggplot(plot_df, aes(x = date, y = value, colour = dose_group,
                      linetype = policy_label, group = interaction(dose_group, policy_label))) +
    geom_line(linewidth = 3) +
    scale_colour_manual(values = panel_cols, breaks = c("≤2 Doses", "Booster"), name = NULL) +
    scale_linetype_manual(values = c("No NPIs" = "solid", "Strict NPIs" = "twodash"), name = NULL) +
    scale_x_date(date_breaks = "1 month", date_labels = "%d %b %Y", expand = c(0.01, 0.08)) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = "Cumulative infections", tag = panel_tag) +
    guides(
      colour   = guide_legend(order = 1, override.aes = list(linewidth = 4, linetype = "solid")),
      linetype = guide_legend(order = 2, override.aes = list(linewidth = 4, colour = "black"))
    ) +
    theme_clean() +
    theme(
      axis.text.x        = element_text(size = 35, hjust = 1, angle = 45),
      axis.text.y        = element_text(size = 35),
      axis.title.y       = element_text(size = 35),
      legend.position    = "top",
      legend.direction   = "horizontal",
      legend.text        = element_text(size = 30),
      legend.key.width   = unit(2, "cm"),
      legend.background  = element_blank(),
      legend.spacing.x   = unit(0.20, "cm"),
      legend.margin      = margin(t = 0, r = 0, b = 0, l = 0),
      panel.grid.major.x = element_blank(),
      plot.tag           = element_text(size = 45, face = "bold"),
      plot.background    = element_blank(),
      plot.margin        = margin(12, 55, 12, 12)
    )
}

make_cumulative_speed_plot <- function(df, panel_tag) {
  plot_df <- df |>
    filter(
      policy == "Strict NPIs",
      vax_scenario %in% c("Slower Boosters", "Faster Boosters"),
      matrix %in% c("cum_inc_leq2dose", "cum_inc_booster")
    ) |>
    mutate(
      dose_group = recode(matrix,
                          "cum_inc_leq2dose" = "≤2 Doses",
                          "cum_inc_booster"  = "Booster"),
      vax_scenario = factor(as.character(vax_scenario), levels = c("Slower Boosters", "Faster Boosters")),
      dose_group   = factor(dose_group, levels = c("≤2 Doses", "Booster"))
    )
  
  ggplot(plot_df, aes(x = date, y = value, colour = vax_scenario,
                      linetype = dose_group, group = interaction(vax_scenario, dose_group))) +
    geom_line(linewidth = 3) +
    scale_colour_manual(values = scenario_cols[c("Slower Boosters", "Faster Boosters")], name = NULL) +
    scale_linetype_manual(values = c("≤2 Doses" = "solid", "Booster" = "twodash"), name = NULL) +
    scale_x_date(date_breaks = "1 month", date_labels = "%d %b %Y", expand = c(0.01, 0.08)) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = "Cumulative infections", tag = panel_tag) +
    guides(
      colour   = guide_legend(order = 1, override.aes = list(linewidth = 4, linetype = "solid")),
      linetype = guide_legend(order = 2, override.aes = list(linewidth = 4, colour = "black"))
    ) +
    theme_clean() +
    theme(
      axis.text.x        = element_text(size = 35, hjust = 1, angle = 45),
      axis.text.y        = element_text(size = 35),
      axis.title.y       = element_text(size = 35),
      legend.position    = "top",
      legend.direction   = "horizontal",
      legend.text        = element_text(size = 30),
      legend.key.width   = unit(2, "cm"),
      legend.background  = element_blank(),
      legend.spacing.x   = unit(0.20, "cm"),
      legend.margin      = margin(t = 0, r = 0, b = 0, l = 0),
      panel.grid.major.x = element_blank(),
      plot.tag           = element_text(size = 45, face = "bold"),
      plot.background    = element_blank(),
      plot.margin        = margin(12, 55, 12, 12)
    )
}

# ---- build panels ----
p_a <- make_shaded_incidence_plot(incidence_df, "Actual Boosters", "No NPIs",
                                  "(A)", "No NPIs",         c(0.60, 0.85), incidence_ymax, "solid")
p_b <- make_shaded_incidence_plot(incidence_df, "Actual Boosters", "Strict NPIs",
                                  "(B)", "Strict NPIs",     c(0.15, 0.85), incidence_ymax, "twodash")
p_c <- make_cumulative_policy_plot(cumulative_df, "Actual Boosters", "(C)")
p_d <- make_shaded_incidence_plot(incidence_df, "Slower Boosters", "Strict NPIs",
                                  "(D)", "Slower Boosters", c(0.15, 0.85), incidence_ymax, "twodash")
p_e <- make_shaded_incidence_plot(incidence_df, "Faster Boosters", "Strict NPIs",
                                  "(E)", "Faster Boosters", c(0.15, 0.85), incidence_ymax, "twodash")
p_f <- make_cumulative_speed_plot(cumulative_df, "(F)")
p_g <- make_shaded_incidence_plot(incidence_df, "Early Boosters",  "Strict NPIs",
                                  "(G)", "Early Boosters",  c(0.15, 0.85), incidence_ymax, "twodash")
p_h <- make_shaded_incidence_plot(incidence_df, "Actual Boosters", "Strict NPIs",
                                  "(H)", "Actual Boosters", c(0.15, 0.85), incidence_ymax, "twodash")

plot_i_df <- incidence_df |>
  filter(
    policy == "Strict NPIs",
    vax_scenario %in% c("Early Boosters", "Actual Boosters"),
    matrix %in% c("inc_leq2dose", "inc_booster")
  ) |>
  group_by(vax_scenario, matrix) |>
  summarise(peak_daily_incidence = max_or_na(value), .groups = "drop") |>
  mutate(
    dose_group     = recode(matrix, "inc_leq2dose" = "≤2 Doses", "inc_booster" = "Booster"),
    scenario_label = factor(as.character(vax_scenario), levels = c("Early Boosters", "Actual Boosters")),
    dose_group     = factor(dose_group, levels = c("≤2 Doses", "Booster"))
  )

p_i <- ggplot(plot_i_df, aes(x = dose_group, y = peak_daily_incidence, fill = scenario_label)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.68,
           colour = "#4D4D4D", linewidth = 0.25) +
  scale_fill_manual(values = scenario_cols[c("Early Boosters", "Actual Boosters")]) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Peak daily incidence", fill = NULL, tag = "(I)") +
  theme_clean() +
  theme(
    axis.text.x        = element_text(size = 35, hjust = 0.5),
    axis.text.y        = element_text(size = 35),
    axis.title.y       = element_text(size = 35),
    legend.position    = c(0.80, 0.90),
    legend.justification = "center",
    legend.direction   = "vertical",
    legend.text        = element_text(size = 35),
    legend.background  = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.tag           = element_text(size = 45, face = "bold"),
    plot.background    = element_blank(),
    plot.margin        = margin(12, 35, 12, 12)
  )

figure1_3x3 <- (
  (p_a | p_b | p_c) /
    (p_d | p_e | p_f) /
    (p_g | p_h | p_i)
) +
  plot_layout(heights = c(1, 1, 1), widths = c(1, 1, 1)) &
  theme(plot.margin = margin(8, 35, 8, 8))

png("../figures/incidence_plot.png", width = 13500, height = 8000, res = 300, bg = "white", type = "cairo")
print(figure1_3x3)
dev.off()


library(tidyverse)
library(lubridate)
library(shellpipes)
library(ggthemes)

loadEnvironments()

NL_POP <- 510550

hoa_date <- as.Date("2022-01-11")
hoa_count <- round(0.209 * NL_POP)
hoa_proportion <- hoa_count / NL_POP

snapshot_beta <- "base"

time_origin <- if (exists("burn_in_start")) {
  as.Date(burn_in_start)
} else {
  as.Date("2021-12-15")
}

vax_colours <- c(
  "actual" = "#4D4D4D",
  "earlier" = "#1F78B4",
  "faster rate" = "#FF7F00",
  "slower rate" = "#CC79A7",
  "national rate" = "#33A02C",
  "national rate (nl init)" = "goldenrod",
  "actual (can init)" = "#006400"
)

scenario_lookup <- tibble(
  vac_code = c(
    "vax_actual_boosters",
    "vax_early_boosters",
    "vax_faster_rate",
    "vax_slower_rate",
    "vax_national_rate",
    "nl_init_national_rate",
    "vax_actual_rate_can_init"
  ),
  vax_scenario = c(
    "actual",
    "earlier",
    "faster rate",
    "slower rate",
    "national rate",
    "national rate (nl init)",
    "actual (can init)"
  )
)

required_objects <- c(
  "vax_actual_boosters",
  "vax_early_boosters",
  "vax_faster_rate",
  "vax_slower_rate",
  "vax_national_rate",
  "vax_nl_init_national_rate",
  "vax_actual_rate_can_init"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_objects) > 0) {
  stop(
    "Missing vaccination scenario objects: ",
    paste(missing_objects, collapse = ", ")
  )
}

prepare_scenario <- function(data, scenario_name) {
  
  required_columns <- c("date", "daily_rate")
  missing_columns <- setdiff(required_columns, names(data))
  
  if (length(missing_columns) > 0) {
    stop(
      scenario_name,
      " is missing: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  data |>
    transmute(
      date = as.Date(date)
      , daily_rate = replace_na(
        as.numeric(daily_rate)
        , 0
      )
      , scenario = scenario_name
    ) |>
    filter(
      !is.na(date)
      , is.finite(daily_rate)
      , daily_rate >= 0
    ) |>
    group_by(
      scenario
      , date
    ) |>
    summarise(
      daily_rate = sum(
        daily_rate
        , na.rm = TRUE
      )
      , .groups = "drop"
    ) |>
    arrange(
      scenario
      , date
    ) |>
    group_by(scenario) |>
    mutate(
      cumulative_boosters = cumsum(daily_rate)
      , booster_proportion =
        cumulative_boosters / NL_POP
    ) |>
    ungroup()
}

booster_vac_props <- bind_rows(
  prepare_scenario(
    vax_actual_boosters
    , "actual"
  ),
  prepare_scenario(
    vax_early_boosters
    , "earlier"
  ),
  prepare_scenario(
    vax_faster_rate
    , "faster rate"
  ),
  prepare_scenario(
    vax_slower_rate
    , "slower rate"
  ),
  prepare_scenario(
    vax_national_rate
    , "national rate"
  ),
  prepare_scenario(
    vax_nl_init_national_rate
    , "national rate (nl init)"
  ),
  prepare_scenario(
    vax_actual_rate_can_init
    , "actual (can init)"
  )
) |>
  mutate(
    scenario = factor(
      scenario
      , levels = names(vax_colours)
    )
  )

duplicate_dates <- booster_vac_props |>
  count(
    scenario
    , date
  ) |>
  filter(n > 1)

if (nrow(duplicate_dates) > 0) {
  stop(
    "Duplicate dates remain in one or more vaccination scenarios."
  )
}

plot_coverage_summary <- booster_vac_props |>
  group_by(scenario) |>
  summarise(
    final_boosters = max(
      cumulative_boosters
      , na.rm = TRUE
    )
    , `Final booster coverage (%)` =
      100 * max(
        booster_proportion
        , na.rm = TRUE
      )
    , .groups = "drop"
  ) |>
  mutate(
    `Final booster coverage (%)` =
      round(
        `Final booster coverage (%)`
        , 1
      )
  )

required_matrices <- c(
  "S1"
  , "V2"
  , "V3"
  , "cum_booster"
)

trajectories <- bind_rows(rdsReadList()) |>
  filter(
    beta == snapshot_beta
    , matrix %in% required_matrices
  ) |>
  transmute(
    vac_code = as.character(vac)
    , matrix = as.character(matrix)
    , date = time_origin + as.integer(time)
    , value = as.numeric(value)
  ) |>
  filter(
    vac_code %in% scenario_lookup$vac_code
    , !is.na(date)
    , is.finite(value)
  )

missing_matrices <- setdiff(
  required_matrices
  , unique(trajectories$matrix)
)

missing_scenarios <- setdiff(
  scenario_lookup$vac_code
  , unique(trajectories$vac_code)
)

if (length(missing_scenarios) > 0) {
  stop(
    "Missing simulation outputs for: "
    , paste(
      missing_scenarios
      , collapse = ", "
    )
  )
}

booster_summary <- trajectories |>
  filter(matrix == "cum_booster") |>
  group_by(vac_code) |>
  summarise(
    booster_start = {
      
      positive_dates <- date[
        value > 0 &
          is.finite(value)
      ]
      
      if (length(positive_dates) > 0) {
        min(positive_dates)
      } else {
        as.Date(NA)
      }
    }
    , final_boosters = {
      
      valid <- is.finite(value) & !is.na(date)
      
      if (any(valid)) {
        value[valid][which.max(date[valid])]
      } else {
        NA_real_
      }
    }
    , .groups = "drop"
  )

missing_booster_scenarios <- booster_summary |>
  filter(is.na(booster_start)) |>
  pull(vac_code)


valid_booster_summary <- (booster_summary 
  |> filter(!is.na(booster_start)) 
  |> mutate(snapshot_date = booster_start)
)
state_at_start <- (trajectories 
  |> filter(matrix %in% c("S1", "V2", "V3")) 
  |> inner_join(valid_booster_summary 
    |> select(vac_code, snapshot_date)
    , by = c("vac_code", "date" = "snapshot_date")
  ) 
  |> group_by(vac_code, matrix) 
  |> summarise(value = sum(value, na.rm = TRUE)
    , .groups = "drop"
  ) |>
  pivot_wider(names_from = matrix, values_from = value, values_fill = 0)
)

vaccination_coverage <- (valid_booster_summary 
  |> select(vac_code, booster_start, final_boosters) 
  |> left_join(state_at_start, by = "vac_code") 
  |> left_join(scenario_lookup, by = "vac_code") 
  |> mutate(
    total_population = NL_POP
    , `No primary series (%)` =
      100 * S1 / NL_POP
    , `Primary series only (%)` =
      100 * V2 / NL_POP
    , `Boosters at rollout start (%)` =
      100 * V3 / NL_POP
    , `Final booster coverage (%)` =
      100 * final_boosters / NL_POP
  ) 
  |>
  select(
    vax_scenario
    , vac_code
    , booster_start
    , total_population
    , `No primary series (%)`
    , `Primary series only (%)`
    , `Boosters at rollout start (%)`
    , `Final booster coverage (%)`
  )
)

actual_coverage <- (vaccination_coverage 
  |> filter(vac_code == "vax_actual_boosters")
)

national_coverage <- (vaccination_coverage 
  |> filter(vac_code == "vax_national_rate")
)

vaccination_coverage <- (vaccination_coverage |>
  mutate(
    `No primary series (%)` = case_when(
      vac_code == "vax_actual_rate_can_init" ~
        national_coverage$`No primary series (%)`[1]
      , vac_code == "nl_init_national_rate" ~
        actual_coverage$`No primary series (%)`[1]
      , TRUE ~ `No primary series (%)`
    )
    , `Primary series only (%)` = case_when(
      vac_code == "vax_actual_rate_can_init" ~
        national_coverage$`Primary series only (%)`[1]
      , vac_code == "nl_init_national_rate" ~
        actual_coverage$`Primary series only (%)`[1]
      , TRUE ~ `Primary series only (%)`
    )
    , `Boosters at rollout start (%)` = case_when(
      vac_code == "vax_actual_rate_can_init" ~
        national_coverage$`Boosters at rollout start (%)`[1]
      , vac_code == "nl_init_national_rate" ~
        actual_coverage$`Boosters at rollout start (%)`[1]
      , TRUE ~ `Boosters at rollout start (%)`
    )
    , `Final booster coverage (%)` = case_when(
      vac_code == "vax_actual_rate_can_init" ~
        actual_coverage$`Final booster coverage (%)`[1]
      , vac_code == "nl_init_national_rate" ~
        national_coverage$`Final booster coverage (%)`[1]
      , TRUE ~ `Final booster coverage (%)`
    )
  ) |>
  mutate(
    vax_scenario = factor(
      vax_scenario
      , levels = names(vax_colours)
    )
  ) |>
  arrange(vax_scenario) |>
  mutate(
    vax_scenario = as.character(vax_scenario)
    , across(
      c(
        `No primary series (%)`
        , `Primary series only (%)`
        , `Boosters at rollout start (%)`
        , `Final booster coverage (%)`
      )
      , \(x) round(x, 1)
    )
  )
)

hybrid_check <- (vaccination_coverage 
  |> filter(vac_code %in% c(
    "vax_actual_boosters"
    , "vax_national_rate"
    , "vax_actual_rate_can_init"
    , "nl_init_national_rate")
  ) 
  |> select(vax_scenario
    , vac_code
    , `No primary series (%)`
    , `Primary series only (%)`
    , `Boosters at rollout start (%)`
    , `Final booster coverage (%)`
  )
)

# write.csv(vaccination_coverage, "../output/vaccination_coverage_table.csv", row.names = FALSE)

plot_start <- as.Date("2021-10-01")
plot_end <- as.Date("2022-06-01")

booster_vac_plot_data <- (booster_vac_props 
  |> filter(date >= plot_start, date <= plot_end)
)

booster_vac_point_data <- (booster_vac_plot_data 
  |> mutate(week = floor_date(
    date
    , unit = "week"
    , week_start = 1
    )
  ) 
  |> group_by(scenario, week) 
  |> slice_max(order_by = date
    , n = 1
    , with_ties = FALSE
  ) 
  |> ungroup()
)

actual_triangle_points <- (booster_vac_point_data 
  |> filter(scenario == "actual", date <= as.Date("2022-01-31"))
)

actual_circle_points <- (booster_vac_point_data 
  |> filter(scenario == "actual", date > as.Date("2022-01-31"))
)

counterfactual_points <- (booster_vac_point_data 
  |> filter(scenario != "actual")
)

booster_vac_props_plot <- (ggplot(booster_vac_plot_data, aes(x = date, y = booster_proportion, colour = scenario, group = scenario)) 
  + geom_line(linewidth = 1.05, lineend = "round") 
  + geom_point(data = counterfactual_points, size = 5, shape = 16) 
  + geom_point(data = actual_triangle_points, size = 5, shape = 17) 
  + geom_point(data = actual_circle_points, size = 5, shape = 16) 
  + annotate(geom = "point", x = hoa_date, y = hoa_proportion
    , shape = 8, size = 8, stroke = 1.2, colour = "red") 
  + annotate(
    geom = "text"
    , x = hoa_date + days(4)
    , y = hoa_proportion
    , label = "HOA"
    , colour = "red"
    , size = 8
    , hjust = 0
    , vjust = 0.45
  ) +
  scale_colour_manual(
    values = vax_colours
    , breaks = names(vax_colours)
    , drop = FALSE
  ) +
  scale_x_date(
    limits = c(
      plot_start
      , plot_end
    )
    , breaks = seq(
      as.Date("2021-10-01")
      , as.Date("2022-06-01")
      , by = "1 month"
    )
    , date_labels = "%b %Y"
    , expand = expansion(
      mult = c(0.01, 0.01)
    )
  ) +
  scale_y_continuous(
    limits = c(0, 0.72)
    , breaks = seq(
      0
      , 0.6
      , by = 0.2
    )
    , labels = scales::label_number(
      accuracy = 0.1
    )
    , expand = expansion(
      mult = c(0, 0)
    )
  ) +
  labs(
    title = "Booster dose vaccination coverage in NL"
    , x = NULL
    , y = "Proportion of NL population"
    , colour = NULL
    , tag = "B"
  ) +
  guides(
    colour = guide_legend(override.aes = list(linewidth = 1.05, shape = 16, size = 3))
  ) 
  +
  theme_clean() +
  theme(axis.text.x = element_text(size = 22, angle = 45, hjust = 1)
    , axis.text.y = element_text(size = 30)
    , axis.title.y = element_text(size = 29, color = "black")
    , legend.position = c(0.85, 0.25)
    , legend.justification = c(0.5, 0.5)
    , legend.direction = "vertical"
    , legend.box = "vertical"
    , legend.title = element_text(size = 25)
    , legend.text = element_text(size = 25)
    , legend.background = element_blank()
    , plot.title = element_text(size = 30, hjust = 0.5, face = "plain")
    , plot.background = element_blank()
    , plot.tag = element_text(size = 35, face = "bold")
    , plot.tag.position = c(0.01, 0.98)
  )
)

booster_vac_props_plot <- booster_vac_props_plot +
  annotate(
    geom = "text"
    , x = plot_start - days(17)
    , y = 0.72
    , label = "B"
    , size = 6
    , fontface = "bold"
    , hjust = 0
    , vjust = 1
  ) +
  coord_cartesian(
    ylim = c(-0.012, 0.72)
    , clip = "off"
  )

# png(filename = "../figures/booster_vac_props_plot.png", width = 5000, height = 2000, res = 300, bg = "white", type = "cairo")
# print(booster_vac_props_plot)
# dev.off()

saveEnvironment()

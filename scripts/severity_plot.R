
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

if (!exists("scenario_cols")) {
  scenario_cols <- c(
    "Actual Boosters" = "#4D4D4D",
    "Early Boosters"  = "#1f78b4",
    "Faster Boosters" = "#ff7f00",
    "Slower Boosters" = "magenta",
    "National Rate"   = "#33a02c"
  )
}

if (!exists("vax_order")) {
  vax_order <- c(
    "Actual Boosters",
    "Early Boosters",
    "Faster Boosters",
    "Slower Boosters",
    "National Rate"
  )
}

vax_short_labs <- c(
  "National Rate"   = "National\nRate",
  "Slower Boosters" = "Slower\nBoosters",
  "Actual Boosters" = "Actual\nBoosters",
  "Early Boosters"  = "Early\nBoosters",
  "Faster Boosters" = "Faster\nBoosters"
)

study_start <- as.Date(omicron_start)

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
      "strict_npis" = "Strict NPIs",
      "no_npis"     = "No NPIs",
      "base"        = "Actual NPIs"
    ),
    date = as.Date(burn_in_start) + (time - 1)
  )

calculate_study_period_outcome <- function(data, matrix_name, outcome_name) {
  
  outcome_data <- (data 
    |> filter(matrix == matrix_name
     , is.finite(value)) 
    |> arrange(policy, vax_scenario, date)
  )
  
  baseline_df <- (outcome_data 
    |> filter(date < study_start) 
    |> group_by(policy, vax_scenario) 
    |> slice_max(order_by = date
     , n = 1
     , with_ties = FALSE) 
    |> ungroup() 
    |> transmute(policy
     , vax_scenario
     , baseline_value = value
     )
  )
  
  endpoint_df <- (outcome_data 
    |> filter(date >= study_start)
    |> group_by(policy
      , vax_scenario) 
    |> slice_max(order_by = date
      , n = 1
      , with_ties = FALSE) 
    |> ungroup() 
    |> transmute(policy
      , vax_scenario
      , endpoint_value = value
      )
  )
  
  endpoint_df |>
    left_join(
      baseline_df,
      by = c("policy", "vax_scenario")
    ) |>
    mutate(
      baseline_value = coalesce(baseline_value, 0),
      severe_outcome = endpoint_value - baseline_value,
      outcome = outcome_name
    ) |>
    select(
      policy,
      vax_scenario,
      outcome,
      baseline_value,
      endpoint_value,
      severe_outcome
    )
}

hospitalization_df <- calculate_study_period_outcome(
  data = combosims
  , matrix_name = "cum_hosp_total"
  , outcome_name = "Hospitalizations"
  )

death_df <- calculate_study_period_outcome(
  data = combosims,
  matrix_name = "cum_death_total",
  outcome_name = "Deaths"
)

severe_df <- bind_rows(
  hospitalization_df,
  death_df
) |>
  filter(
    is.finite(severe_outcome),
    severe_outcome >= 0
  ) |>
  mutate(
    policy_label = factor(
      as.character(policy),
      levels = c(
        "No NPIs"
        , "Actual NPIs"
        , "Strict NPIs"
      )
    )
    , outcome = factor( outcome, levels = c("Hospitalizations", "Deaths"))
    , vax_scenario = factor(as.character(vax_scenario),levels = vax_order)
  )

x_order_df <- (severe_df 
  |> filter(outcome == "Hospitalizations") 
  |> group_by(policy_label) 
  |> arrange(desc(severe_outcome), .by_group = TRUE) 
  |> mutate(x_pos = row_number()) 
  |> ungroup() 
  |> select(policy_label, vax_scenario, x_pos)
  )

severe_df <- (severe_df 
  |> left_join(x_order_df , by = c("policy_label", "vax_scenario")) 
  |> mutate(x_lab = vax_short_labs[as.character(vax_scenario)], value_label = comma(round(severe_outcome, 0))) 
  |> filter(!is.na(x_pos)) 
  )

severity_plot <- (ggplot(severe_df, aes(x = x_pos, y = severe_outcome, fill = vax_scenario)) 
  + geom_col(colour = "black", width = 0.72, linewidth = 0.25) 
  + geom_text(aes(label = value_label), vjust = -0.35, size = 10) 
  + facet_grid(outcome ~ policy_label, scales = "free_y") 
  + scale_fill_manual(values = scenario_cols, breaks = vax_order, drop = FALSE) 
  + scale_x_continuous(breaks = 1:5, labels = NULL,expand = expansion(mult = c(0.04, 0.04))) 
  + scale_y_continuous(labels = comma, expand = expansion(mult = c(0.08, 0.18))) 
  + labs(
    title = "Severe outcomes across NPI and booster scenarios",
    x = NULL,
    y = "Severe outcomes",
    fill = NULL) 
  + theme_clean() 
  + theme(
      axis.text.x = element_text(size = 30),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 35),
      axis.title.y = element_text(size = 35, colour = "black"),
      strip.text = element_text(size = 35),
      legend.position = "top",
      legend.text = element_text(size = 25),
      legend.background = element_blank(),
      legend.spacing.x = unit(0.45, "cm"),
      plot.title = element_text(size = 40, hjust = 0.5),
      panel.grid.major.x = element_blank(),
      panel.spacing = unit(2, "cm"),
      plot.margin = margin(t = 10,r = 20,b = 55,l = 20),
      plot.background = element_blank()
    ) 
  + geom_text(data = severe_df |>
      filter(outcome == "Deaths"),aes(y = -Inf,label = x_lab),vjust = 1.15,size = 8,lineheight = 0.85,inherit.aes = TRUE) 
  + coord_cartesian(clip = "off")
)

png("../figures/figure3.png", width = 10000, height = 5000, res = 300, bg = "white",type = "cairo")
print(severity_plot)
dev.off()
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(shellpipes)
  library(macpan2)
  library(readr)
  library(purrr)
})

loadEnvironments()

# ============================================================
# Constants
# ============================================================

NL_POP <- 510550

analysis_start <- as.Date("2021-12-15")
offset <- 150L
init_date <- analysis_start - offset
analysis_end <- as.Date("2022-05-22")

actual_booster_start <- as.Date("2021-11-08")
early_booster_start  <- as.Date("2021-10-04")

sim_horizon <- as.integer(analysis_end - init_date + 1L)

new_V2_init <- 277707

# ============================================================
# Load calibrated model
# ============================================================

cal <- rdsRead("calibrate")

spec_base <- mp_optimized_spec(
  cal,
  spec_structure = "modified"
)

stopifnot(exists("booster_inputs_df"))

# ============================================================
# Scenario information
# ============================================================

scenario_info <- tibble(
  scenario_input = c(
    "Actual Boosters",
    "Slower Rate",
    "Faster Rate",
    "Early Boosters",
    "National Rate"
  ),
  Scenario = c(
    "Actual",
    "Slower rate",
    "Faster rate",
    "Earlier",
    "Canadian rate"
  ),
  start_date = as.Date(c(
    "2021-11-08",
    "2021-11-08",
    "2021-11-08",
    "2021-10-04",
    "2021-11-08"
  ))
)

# ============================================================
# Handle older early-scenario names
# ============================================================

available_scenarios <- unique(booster_inputs_df$scenario)

if (!"Early Boosters" %in% available_scenarios) {
  
  early_name <- available_scenarios[
    grepl("^Early Boosters", available_scenarios, ignore.case = TRUE)
  ]
  
  if (length(early_name) == 1) {
    
    scenario_info$scenario_input[
      scenario_info$Scenario == "Earlier"
    ] <- early_name
  }
}

# ============================================================
# Helper: extract booster inputs
# ============================================================

get_booster_inputs <- function(scenario_name) {
  
  booster_df <- booster_inputs_df |>
    
    filter(scenario == scenario_name) |>
    
    transmute(
      changepoints = as.integer(day),
      vac_values   = as.numeric(booster_vac_value)
    ) |>
    
    group_by(changepoints) |>
    
    summarise(
      vac_values = sum(vac_values, na.rm = TRUE),
      .groups = "drop"
    ) |>
    
    arrange(changepoints)
  
  if (nrow(booster_df) == 0) {
    stop(
      paste0(
        "No booster rows found for scenario: ",
        scenario_name
      )
    )
  }
  
  booster_df <- bind_rows(
    tibble(
      changepoints = 0L,
      vac_values = 0
    ),
    booster_df
  ) |>
    
    group_by(changepoints) |>
    
    summarise(
      vac_values = sum(vac_values, na.rm = TRUE),
      .groups = "drop"
    ) |>
    
    arrange(changepoints)
  
  booster_df
}

# ============================================================
# Helper: create scenario-specific specification
# ============================================================

set_scenario_spec <- function(scenario_name) {
  
  booster_df <- get_booster_inputs(scenario_name)
  
  spec_i <- spec_base
  
  spec_i$integers$booster_vac_changepoints <-
    as.integer(booster_df$changepoints)
  
  spec_i$default$booster_vac_values <-
    as.numeric(booster_df$vac_values)
  
  # National-rate scenario
  if (scenario_name == "National Rate") {
    
    spec_i <- spec_i |>
      mp_tmb_insert(
        default = list(
          V2_init = new_V2_init
        )
      )
  }
  
  spec_i
}

# ============================================================
# Helper: extract vaccination states
# ============================================================

extract_state_row <- function(
    scenario_name,
    scenario_label,
    start_date
) {
  
  cat("\nRunning scenario:", scenario_label, "\n")
  
  spec_i <- set_scenario_spec(scenario_name)
  
  traj <- spec_i |>
    
    mp_simulator(
      time_steps = sim_horizon,
      outputs = c(
        "S1",
        "V2",
        "V3",
        "cum_booster"
      )
    ) |>
    
    mp_trajectory() |>
    
    mutate(
      date = init_date + time,
      matrix = as.character(matrix),
      value = as.numeric(value)
    )
  
  # ----------------------------------------------------------
  # Vaccination states at scenario start date
  # ----------------------------------------------------------
  
  state_start <- traj |>
    
    filter(date == start_date) |>
    
    select(matrix, value) |>
    
    pivot_wider(
      names_from = matrix,
      values_from = value
    )
  
  if (nrow(state_start) != 1) {
    
    stop(
      paste0(
        "Could not extract vaccination states for ",
        scenario_label
      )
    )
  }
  
  # ----------------------------------------------------------
  # Final boosters
  # ----------------------------------------------------------
  
  final_boosters <- traj |>
    
    filter(
      date <= analysis_end,
      matrix == "cum_booster"
    ) |>
    
    summarise(
      value = max(value, na.rm = TRUE)
    ) |>
    
    pull(value)
  
  tibble(
    Scenario = scenario_label,
    
    `Start date` = format(
      start_date,
      "%d-%b-%Y"
    ),
    
    `No primary series` =
      100 * state_start$S1 / NL_POP,
    
    `Primary series only` =
      100 * state_start$V2 / NL_POP,
    
    Boosters =
      100 * state_start$V3 / NL_POP,
    
    `Final boosters` =
      100 * final_boosters / NL_POP
  )
}

# ============================================================
# Build vaccination table
# ============================================================

vaccination_table <- pmap_dfr(
  
  list(
    scenario_info$scenario_input,
    scenario_info$Scenario,
    scenario_info$start_date
  ),
  
  extract_state_row
  
) |>
  
  mutate(
    
    across(
      c(
        `No primary series`,
        `Primary series only`,
        Boosters,
        `Final boosters`
      ),
      
      ~ round(.x, 1)
    )
  )

knitr::kable(vaccination_table)
# options(width = 1000)
# print(vaccination_table, width = Inf) 

# ============================================================
# Save outputs safely
# ============================================================

output_dir <- "../output"

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

output_csv <- file.path(
  output_dir,
  "vaccination_scenarios_summary.csv"
)

# remove old file if locked/corrupted
if (file.exists(output_csv)) {
  try(unlink(output_csv), silent = TRUE)
}

write.csv(
  vaccination_table,
  output_csv,
  row.names = FALSE
)

# ============================================================
# Save environment objects
# ============================================================

saveVars(vaccination_table)
# suppressPackageStartupMessages({
#   library(dplyr)
#   library(tidyr)
#   library(lubridate)
#   library(shellpipes)
#   library(macpan2)
#   library(readr)
#   library(purrr)
#   library(knitr)
# })
# 
# loadEnvironments()
# 
# # ============================================================
# # constants
# # ============================================================
# 
# NL_POP <- 510550
# 
# analysis_start <- as.Date("2021-12-15")
# analysis_end   <- as.Date("2022-05-22")
# 
# offset    <- 150L
# init_date <- analysis_start - offset
# 
# actual_booster_start <- as.Date("2021-11-08")
# early_booster_start  <- as.Date("2021-10-04")
# 
# sim_horizon <- as.integer(analysis_end - init_date) + 1L
# 
# new_V2_init <- 277707
# 
# # ============================================================
# # load calibrated model
# # ============================================================
# 
# cal <- rdsRead("calibrate")
# 
# spec_base <- mp_optimized_spec(
#   cal,
#   spec_structure = "modified"
# )
# 
# stopifnot(exists("booster_inputs_df"))
# 
# # ============================================================
# # standardize scenario names
# # ============================================================
# 
# booster_inputs_df <- booster_inputs_df |>
#   mutate(
#     scenario = recode(
#       scenario,
#       "Actual Rate" = "Actual Boosters",
#       "Early Rate" = "Early Boosters",
#       "Earlier Boosters" = "Early Boosters",
#       "Earlier boosters" = "Early Boosters",
#       "Early boosters" = "Early Boosters",
#       "Faster Rate" = "Faster Boosters",
#       "Slower Rate" = "Slower Boosters",
#       "National rate" = "National Rate",
#       .default = scenario
#     ),
#     date = as.Date(date),
#     day = as.integer(day),
#     booster_vac_value = as.numeric(booster_vac_value)
#   )
# 
# print(
#   tibble(
#     available_scenarios = sort(unique(booster_inputs_df$scenario))
#   )
# )
# 
# # ============================================================
# # scenario information
# # ============================================================
# 
# scenario_info <- tibble(
#   scenario_input = c(
#     "Actual Boosters",
#     "Slower Boosters",
#     "Faster Boosters",
#     "Early Boosters",
#     "National Rate"
#   ),
#   Scenario = c(
#     "Actual",
#     "Slower rate",
#     "Faster rate",
#     "Earlier",
#     "Canadian rate"
#   ),
#   start_date = as.Date(c(
#     "2021-11-08",
#     "2021-11-08",
#     "2021-11-08",
#     "2021-10-04",
#     "2021-11-08"
#   ))
# )
# 
# missing_scenarios <- setdiff(
#   scenario_info$scenario_input,
#   unique(booster_inputs_df$scenario)
# )
# 
# if (length(missing_scenarios) > 0) {
#   stop(
#     "Missing booster rows for scenario(s): ",
#     paste(missing_scenarios, collapse = ", "),
#     "\nAvailable scenarios are: ",
#     paste(sort(unique(booster_inputs_df$scenario)), collapse = ", ")
#   )
# }
# 
# # ============================================================
# # helper: extract booster inputs
# # ============================================================
# 
# get_booster_inputs <- function(scenario_name) {
#   
#   booster_df <- booster_inputs_df |>
#     filter(
#       scenario == scenario_name,
#       date <= analysis_end,
#       day >= 0
#     ) |>
#     transmute(
#       changepoints = as.integer(day),
#       vac_values = as.numeric(booster_vac_value)
#     ) |>
#     group_by(changepoints) |>
#     summarise(
#       vac_values = sum(vac_values, na.rm = TRUE),
#       .groups = "drop"
#     ) |>
#     arrange(changepoints)
#   
#   if (nrow(booster_df) == 0) {
#     stop(
#       "No booster rows found for scenario: ",
#       scenario_name,
#       "\nAvailable scenarios are: ",
#       paste(sort(unique(booster_inputs_df$scenario)), collapse = ", ")
#     )
#   }
#   
#   booster_df |>
#     bind_rows(
#       tibble(
#         changepoints = 0L,
#         vac_values = 0
#       )
#     ) |>
#     group_by(changepoints) |>
#     summarise(
#       vac_values = sum(vac_values, na.rm = TRUE),
#       .groups = "drop"
#     ) |>
#     arrange(changepoints)
# }
# 
# # ============================================================
# # helper: create scenario-specific specification
# # ============================================================
# 
# set_scenario_spec <- function(scenario_name) {
#   
#   booster_df <- get_booster_inputs(scenario_name)
#   
#   spec_i <- spec_base
#   
#   spec_i$integers$booster_vac_changepoints <-
#     as.integer(booster_df$changepoints)
#   
#   spec_i$default$booster_vac_values <-
#     as.numeric(booster_df$vac_values)
#   
#   if (scenario_name == "National Rate") {
#     spec_i <- spec_i |>
#       mp_tmb_insert(
#         default = list(
#           V2_init = new_V2_init
#         )
#       )
#   }
#   
#   spec_i
# }
# 
# # ============================================================
# # helper: extract vaccination states
# # ============================================================
# 
# extract_state_row <- function(
#     scenario_name,
#     scenario_label,
#     start_date
# ) {
#   
#   cat("\nRunning scenario:", scenario_label, "\n")
#   
#   spec_i <- set_scenario_spec(scenario_name)
#   
#   traj <- spec_i |>
#     mp_simulator(
#       time_steps = sim_horizon,
#       outputs = c(
#         "S1",
#         "V2",
#         "V3",
#         "cum_booster"
#       )
#     ) |>
#     mp_trajectory() |>
#     mutate(
#       date = init_date + time,
#       matrix = as.character(matrix),
#       value = as.numeric(value)
#     )
#   
#   state_start <- traj |>
#     filter(date == start_date) |>
#     select(matrix, value) |>
#     pivot_wider(
#       names_from = matrix,
#       values_from = value
#     )
#   
#   if (nrow(state_start) != 1) {
#     stop(
#       "Could not extract vaccination states for ",
#       scenario_label,
#       " on ",
#       start_date
#     )
#   }
#   
#   final_boosters <- traj |>
#     filter(
#       date <= analysis_end,
#       matrix == "cum_booster"
#     ) |>
#     summarise(
#       value = max(value, na.rm = TRUE)
#     ) |>
#     pull(value)
#   
#   tibble(
#     Scenario = scenario_label,
#     `Start date` = format(start_date, "%d-%b-%Y"),
#     `No primary series` = 100 * state_start$S1 / NL_POP,
#     `Primary series only` = 100 * state_start$V2 / NL_POP,
#     Boosters = 100 * state_start$V3 / NL_POP,
#     `Final boosters` = 100 * final_boosters / NL_POP
#   )
# }
# 
# # ============================================================
# # build vaccination table
# # ============================================================
# 
# vaccination_table <- pmap_dfr(
#   list(
#     scenario_info$scenario_input,
#     scenario_info$Scenario,
#     scenario_info$start_date
#   ),
#   extract_state_row
# ) |>
#   mutate(
#     across(
#       c(
#         `No primary series`,
#         `Primary series only`,
#         Boosters,
#         `Final boosters`
#       ),
#       ~ round(.x, 1)
#     )
#   )
# 
# print(knitr::kable(vaccination_table))
# 
# # ============================================================
# # save outputs
# # ============================================================
# 
# output_dir <- "../output"
# 
# if (!dir.exists(output_dir)) {
#   dir.create(
#     output_dir,
#     recursive = TRUE,
#     showWarnings = FALSE
#   )
# }
# 
# output_csv <- file.path(
#   output_dir,
#   "vaccination_scenarios_summary.csv"
# )
# 
# if (file.exists(output_csv)) {
#   try(unlink(output_csv), silent = TRUE)
# }
# 
# write.csv(
#   vaccination_table,
#   output_csv,
#   row.names = FALSE
# )
# 
# saveVars(vaccination_table)






suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(shellpipes)
  library(macpan2)
  library(readr)
  library(purrr)
  library(knitr)
})

loadEnvironments()

# ============================================================
# constants
# ============================================================

NL_POP <- 510550

analysis_start <- as.Date("2021-12-15")
analysis_end   <- as.Date("2022-05-22")

offset    <- 150L
init_date <- analysis_start - offset

actual_booster_start <- as.Date("2021-11-08")
early_booster_start  <- as.Date("2021-10-04")

sim_horizon <- as.integer(analysis_end - init_date) + 1L

CA_PRIMARY_ONLY_PROP_START <- 0.750
new_V2_init <- NA_real_

# ============================================================
# load calibrated model
# ============================================================

cal <- rdsRead("calibrate")

spec_base <- mp_optimized_spec(
  cal,
  spec_structure = "modified"
)

stopifnot(exists("booster_inputs_df"))

# ============================================================
# standardize scenario names
# ============================================================

booster_inputs_df <- booster_inputs_df |>
  mutate(
    scenario = recode(
      scenario,
      "Actual Rate" = "Actual Boosters",
      "Early Rate" = "Early Boosters",
      "Earlier Boosters" = "Early Boosters",
      "Earlier boosters" = "Early Boosters",
      "Early boosters" = "Early Boosters",
      "Faster Rate" = "Faster Boosters",
      "Slower Rate" = "Slower Boosters",
      "National rate" = "National Rate",
      .default = scenario
    ),
    date = as.Date(date),
    day = as.integer(day),
    booster_vac_value = as.numeric(booster_vac_value)
  )

print(
  tibble(
    available_scenarios = sort(unique(booster_inputs_df$scenario))
  )
)

# ============================================================
# scenario information
# ============================================================

scenario_info <- tibble(
  scenario_input = c(
    "Actual Boosters",
    "Slower Boosters",
    "Faster Boosters",
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

missing_scenarios <- setdiff(
  scenario_info$scenario_input,
  unique(booster_inputs_df$scenario)
)

if (length(missing_scenarios) > 0) {
  stop(
    "Missing booster rows for scenario(s): ",
    paste(missing_scenarios, collapse = ", "),
    "\nAvailable scenarios are: ",
    paste(sort(unique(booster_inputs_df$scenario)), collapse = ", ")
  )
}

# ============================================================
# helper: extract booster inputs
# ============================================================

get_booster_inputs <- function(scenario_name) {
  
  booster_df <- booster_inputs_df |>
    filter(
      scenario == scenario_name,
      date <= analysis_end,
      day >= 0
    ) |>
    transmute(
      changepoints = as.integer(day),
      vac_values = as.numeric(booster_vac_value)
    ) |>
    group_by(changepoints) |>
    summarise(
      vac_values = sum(vac_values, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(changepoints)
  
  if (nrow(booster_df) == 0) {
    stop(
      "No booster rows found for scenario: ",
      scenario_name,
      "\nAvailable scenarios are: ",
      paste(sort(unique(booster_inputs_df$scenario)), collapse = ", ")
    )
  }
  
  booster_df |>
    bind_rows(
      tibble(
        changepoints = 0L,
        vac_values = 0
      )
    ) |>
    group_by(changepoints) |>
    summarise(
      vac_values = sum(vac_values, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(changepoints)
}

# ============================================================
# helper: create scenario-specific specification
# ============================================================

set_scenario_spec_raw <- function(scenario_name, V2_init_override = NULL) {
  
  booster_df <- get_booster_inputs(scenario_name)
  
  spec_i <- spec_base
  
  spec_i$integers$booster_vac_changepoints <-
    as.integer(booster_df$changepoints)
  
  spec_i$default$booster_vac_values <-
    as.numeric(booster_df$vac_values)
  
  if (!is.null(V2_init_override)) {

  if (scenario_name == "National Rate") {
    spec_i <- spec_i |>
      mp_tmb_insert(
        default = list(
          V2_init = as.numeric(V2_init_override)
        )
      )
  }
  
  spec_i
  }
}

get_primary_only_prop_at_start <- function(V2_init_candidate) {
  
  spec_i <- set_scenario_spec_raw(
    scenario_name = "National Rate",
    V2_init_override = V2_init_candidate
  )
  
  traj <- spec_i |>
    mp_simulator(
      time_steps = sim_horizon,
      outputs = c("V2")
    ) |>
    mp_trajectory() |>
    mutate(
      date = init_date + time,
      matrix = as.character(matrix),
      value = as.numeric(value)
    )
  
  out <- traj |>
    filter(date == actual_booster_start, matrix == "V2") |>
    summarise(value = first(value)) |>
    pull(value)
  
  if (length(out) != 1 || !is.finite(out)) {
    stop("Could not evaluate National Rate V2 at booster start.")
  }
  
  out / NL_POP
}

calibrate_national_V2_init <- function(
    target_prop = CA_PRIMARY_ONLY_PROP_START,
    lower = 0,
    upper = NL_POP) 
  {
  
  objective <- function(x) {
    get_primary_only_prop_at_start(x) - target_prop
  }
  
  f_lower <- objective(lower)
  f_upper <- objective(upper)
  
  if (f_lower * f_upper > 0) {
    stop(
      "Could not bracket National Rate V2_init. ",
      "f(lower)=", round(f_lower, 4),
      "; f(upper)=", round(f_upper, 4)
    )
  }
  
  round(uniroot(objective, lower = lower, upper = upper)$root)
}

new_V2_init <- calibrate_national_V2_init()

set_scenario_spec <- function(scenario_name) {
  
  if (scenario_name == "National Rate") {
    return(
      set_scenario_spec_raw(
        scenario_name = scenario_name,
        V2_init_override = new_V2_init
      )
    )
  }
  
  set_scenario_spec_raw(scenario_name)
}

# ============================================================
# helper: extract vaccination states
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
  
  state_start <- traj |>
    filter(date == start_date) |>
    select(matrix, value) |>
    pivot_wider(
      names_from = matrix,
      values_from = value
    )
  
  if (nrow(state_start) != 1) {
    stop(
      "Could not extract vaccination states for ",
      scenario_label,
      " on ",
      start_date
    )
  }
  
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
    `Start date` = format(start_date, "%d-%b-%Y"),
    `No primary series` = 100 * state_start$S1 / NL_POP,
    `Primary series only` = 100 * state_start$V2 / NL_POP,
    Boosters = 100 * state_start$V3 / NL_POP,
    `Final boosters` = 100 * final_boosters / NL_POP
  )
}

# ============================================================
# build vaccination table
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

print(knitr::kable(vaccination_table))

# ============================================================
# save outputs
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

if (file.exists(output_csv)) {
  try(unlink(output_csv), silent = TRUE)
}

write.csv(
  vaccination_table,
  output_csv,
  row.names = FALSE
)

saveVars(vaccination_table)
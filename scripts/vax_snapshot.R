suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(shellpipes)
  library(lubridate)
})

options(dplyr.summarise.inform = FALSE)
options(scipen = 9999)

loadEnvironments()

## ---- population denominators ----
NL_POP <- 510550
CA_POP <- 36991981

## scenarios reported against the Canadian total population (rate taken from
## national/Canadian data); everything else is against the NL population.
national_vac_codes <- c("vax_national_rate")

# read the pre-computed simulations ----
combosims <- bind_rows(rdsReadList())

needed_matrices <- c("S1", "V2", "V3", "cum_booster")
missing_matrices <- setdiff(needed_matrices, unique(as.character(combosims$matrix)))

time_origin <- if (exists("burn_in_start")) burn_in_start else as.Date("2021-12-15")

traj_all <- combosims |>
  filter(beta == snapshot_beta, matrix %in% needed_matrices) |>
  mutate(matrix = as.character(matrix),
         value  = as.numeric(value),
         date   = time_origin + (time - 1))

vac_to_label <- c(
  vax_actual_boosters = "Actual Boosters",
  vax_early_boosters  = "Early Boosters",
  vax_faster_rate     = "Faster Boosters",
  vax_slower_rate     = "Slower Boosters",
  vax_national_rate   = "National Rate"
)

booster_start_date_for <- function(traj) {
  d <- traj |>
    filter(matrix == "cum_booster", value > 0) |>
    summarise(d = suppressWarnings(min(date))) |>
    pull(d)
  if (length(d) == 0 || is.infinite(as.numeric(d)) || is.na(d)) return(as.Date(NA))
  d
}

# build one row per vaccination scenario 
vax_codes <- sort(unique(traj_all$vac))

rows <- lapply(vax_codes, function(vc) {
  
  traj <- traj_all |> filter(vac == vc)
  if (nrow(traj) == 0) return(NULL)
  
  start_date <- booster_start_date_for(traj)
  if (is.na(start_date)) {
    warning("No booster administration (cum_booster > 0) found for '", vc, "'; skipping.")
    return(NULL)
  }
  
  state_at_start <- traj |>
    filter(date == start_date - 1, matrix %in% c("S1", "V2", "V3")) |>
    select(matrix, value) |>
    pivot_wider(names_from = matrix, values_from = value)
  
  final_boosters <- traj |>
    filter(matrix == "cum_booster", date == max(date)) |>
    pull(value)
  
  #raw compartment counts (model scale)
  no_primary   <- state_at_start$S1
  primary_only <- state_at_start$V2
  boosted      <- state_at_start$V3
  
  # denominator for the "% of total population" column
  denom_total <- if (vc %in% national_vac_codes) CA_POP else NL_POP
  
  tibble(
    vax_scenario       = unname(vac_to_label[vc]),
    vac_code           = vc,
    booster_start_date = start_date,
    total_population   = denom_total,
    
    # % of the scenario's own total population
    `No primary series (%total)`      = 100 * no_primary   / denom_total,
    `Primary series only (%total)`    = 100 * primary_only / denom_total,
    `Boosters (%total)`               = 100 * boosted      / denom_total,
    `Final booster coverage (%total)` = 100 * final_boosters / denom_total,
    
    # % of the NL total population (all scenarios, incl. National Rate)
    `No primary series (%NL)`         = 100 * no_primary   / NL_POP,
    `Primary series only (%NL)`       = 100 * primary_only / NL_POP,
    `Boosters (%NL)`                  = 100 * boosted      / NL_POP,
    `Final booster coverage (%NL)`    = 100 * final_boosters / NL_POP
  )
})

vaccination_coverage <- bind_rows(rows) |>
  mutate(across(where(is.numeric) & !all_of("total_population"), ~ round(.x, 1)))

print(as.data.frame(vaccination_coverage))

write.csv(vaccination_coverage, "../output/vaccination_coverage_table.csv", row.names = FALSE)

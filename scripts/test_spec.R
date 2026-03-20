library(macpan2)
library(shellpipes)
library(tidyverse)

rm(list = ls())
gc()
loadEnvironments()

spec <- rdsRead("timevar_spec.rds")

cat("\n--- names(spec) ---\n")
print(names(spec))

cat("\n--- names(spec$default) ---\n")
print(names(spec$default))

# Test the smallest possible run
test_outputs <- c("beta_thing", "inc_safe", "newR_safe", "serop_total")

sim <- tryCatch(
  {
    spec |>
      mp_rk4() |>
      mp_simulator(
        time = 5,
        outputs = test_outputs
      ) |>
      mp_trajectory()
  },
  error = function(e) {
    message("\nMinimal simulator failed: ", conditionMessage(e))
    stop(e)
  }
)

print(sim)
print(sim |> dplyr::filter(!is.finite(value)))
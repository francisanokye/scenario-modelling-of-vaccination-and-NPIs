
library(macpan2)
library(tidyverse); theme_set(theme_bw())
library(shellpipes)

spec <- rdsRead()
time_steps <- 300L

## piece-wise beta: 0.2 early, stepping up to 0.4 at t = 200
beta_sched <- c(0.2, 0.2, 0.4)
beta_times <- c(1, 2, 200)

nu <- 0.01  # baseline per-capita vaccination rate

# scenario 1: constant vax rate across the whole horizon
vax1_sched <- c(nu, nu, nu)

vax2_sched <- c(2 * nu, 2 * nu, 0)

vax_times  <- c(1, 2, 150)

## insert beta and vax_rate
make_scenario <- function(spec, vax_sched, vax_times) {
  spec |>
    mp_tmb_insert(
      phase       = "during"
      , at          = 1
      , expressions = list(beta ~ time_var(beta_sched, beta_times))
      , default     = list(beta_sched = beta_sched)
      , integers    = list(beta_times = beta_times)
    ) |>
    mp_tmb_insert(
      phase       = "during"
      , at          = 1
      , expressions = list(vax_rate ~ time_var(vax_sched, vax_times))
      , default     = list(vax_sched = vax_sched)
      , integers    = list(vax_times = vax_times)
    )
}

# simulation of each scenario
outputs <- c("inc", "cinc", "beta_thing", "vax_thing", "Vbox")
sim_scenario <- function(scenario_spec, label) {
  (mp_simulator(
    model = scenario_spec
    , time_steps = time_steps
    , outputs = outputs
  )
  |> mp_trajectory()
  |> mutate(scenario = label)
  )
}

# bisection method to check when scenario 2 has to stop
calibrate <- TRUE

total_vax <- function(traj) {
  traj |> filter(matrix == "Vbox", time == max(time)) |> pull(value)
}

target <- total_vax(sim_scenario(make_scenario(spec, vax1_sched, vax_times), "cal"))

if (calibrate) {
  lo <- 2; hi <- time_steps
  for (i in 1:20) {
    mid <- round((lo + hi) / 2)
    got <- total_vax(sim_scenario(
      make_scenario(spec, vax2_sched, c(1, 2, mid)), "cal"))
    if (got < target) lo <- mid else hi <- mid
  }
  vax_end   <- round((lo + hi) / 2)
  vax_times <- c(1, 2, vax_end)
}

spec_s1 <- make_scenario(spec, vax1_sched, c(1, 2, time_steps))
spec_s2 <- make_scenario(spec, vax2_sched, vax_times)

sim <- bind_rows(sim_scenario(spec_s1, "constant vax rate"), sim_scenario(spec_s2, sprintf("double, ends t=%d", vax_end)))

# check that these two numbers should match
check <- (sim
          |> filter(matrix == "Vbox", time == max(time))
          |> select(scenario, total_vaccinated = value)
)
print(as.data.frame(check))

gg <- (ggplot(
  filter(sim,
         matrix %in% c("cinc", "inc", "beta_thing", "vax_thing", "Vbox"))
  , aes(time, value, colour = scenario))
  + geom_line()
  + facet_wrap(~matrix, scales = "free_y")
  + theme(legend.position = "bottom")
)
print(gg)

rdsSave(sim)

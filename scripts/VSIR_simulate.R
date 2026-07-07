library(macpan2)
library(tidyverse); theme_set(theme_bw())
library(shellpipes)

spec <- rdsRead()

time_steps <- 200L

tvec       <- seq_len(time_steps) 

# piece-wise beta
beta_break <- 150
beta_sched <- ifelse(tvec < beta_break, 0.2, 0.4)

nu      <- 0.01                     # baseline per-capita vaccination rate
vax_end <- time_steps %/% 2L        # scenario 2 finishes halfway (t = 100)

# scenario 1: constant vax rate across the whole horizon
vax1_sched <- rep(nu, time_steps)
# scenario 2: double the daily rate
vax2_sched <- ifelse(tvec < vax_end, 2 * nu, 0)

# insert beta and vax_rate 
make_scenario <- function(spec, beta_sched, vax_sched, times) {
  spec |>
    mp_tmb_insert(
      phase       = "during"
      , at          = 1
      , expressions = list(beta ~ time_var(beta_sched, beta_times))
      , default     = list(beta_sched = beta_sched)
      , integers    = list(beta_times = times)
    ) |>
    mp_tmb_insert(
      phase       = "during"
      , at          = 1
      , expressions = list(vax_rate ~ time_var(vax_sched, vax_times))
      , default     = list(vax_sched = vax_sched)
      , integers    = list(vax_times = times)
    )
}
spec_s1 <- make_scenario(spec, beta_sched, vax1_sched, tvec)
spec_s2 <- make_scenario(spec, beta_sched, vax2_sched, tvec)

# simulation of each scenario
outputs <- c("incidence", "vaccination", "S", "V", "I", "R", "beta_thing", "vax_thing")
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
  
  sim <- bind_rows(sim_scenario(spec_s1, "constant vax rate") , sim_scenario(spec_s2, "double ends halfway"))

gg <- (ggplot(
  filter(sim, 
    matrix %in% c("S", "V", "I", "R", "incidence", "beta_thing", "vax_thing"))
    , aes(time, value, colour = scenario))
  + geom_line()
  + facet_wrap(~matrix, scales = "free_y")
)

print(gg)

rdsSave(sim)

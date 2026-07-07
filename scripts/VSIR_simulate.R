library(macpan2)
library(tidyverse); theme_set(theme_bw())
library(shellpipes)

spec <- rdsRead()

time_steps <- 300L

# piece-wise beta
beta_sched <- c(0.2,0.2,0.4)
beta_times <- c(1,2, 200)

nu      <- 0.01                     # baseline per-capita vaccination rate
vax_end <- time_steps %/% 2L        # scenario 2 finishes halfway (t = 100)

# scenario 1: constant vax rate across the whole horizon
vax1_sched <- c(nu ,nu, nu)

vax2_sched <- c(2*nu,2*nu,0)

vax_times <- c(1,2, 150)

# insert beta and vax_rate 
make_scenario <- function(spec, vax_sched) {
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
spec_s1 <- make_scenario(spec, vax1_sched)
spec_s2 <- make_scenario(spec, vax2_sched)

# simulation of each scenario
outputs <- c("inc", "cinc","beta_thing", "vax_thing","Vbox")
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
    matrix %in% c("cinc", "inc", "beta_thing", "vax_thing","Vbox"))
    , aes(time, value, colour = scenario))
  + geom_line()
  + facet_wrap(~matrix, scales = "free_y")
  + theme(legend.position="bottom")
)

print(gg)

rdsSave(sim)

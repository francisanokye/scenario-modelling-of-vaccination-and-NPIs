library(macpan2)
library(shellpipes)

flows = list(
  mp_per_capita_flow("S", "I", "beta * I / N", "incidence")
  , mp_per_capita_flow("S", "V", "vax_rate",     "vaccination")
  , mp_per_capita_flow("I", "R", "gamma",        "recovery")
  , beta_thing ~ beta        
  , vax_thing  ~ vax_rate    
)

default = list(
    beta     = 0.4    
  , gamma    = 0.1    
  , vax_rate = 0.0    
  , N        = 1000   
  , I        = 1      
  , R        = 0      
  , V        = 0      
)

initialize_state = list(S ~ N - I - R - V)

# model specification
spec = mp_tmb_model_spec(
  before = initialize_state
  , during = flows
  , default = default
)

rdsSave(spec)

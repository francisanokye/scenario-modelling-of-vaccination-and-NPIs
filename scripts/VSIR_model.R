library(macpan2)
library(shellpipes)

flows = list(
  mp_per_capita_flow("S", "I", "beta * kappa * (I + VI) / N","inc_S")
  , mp_per_capita_flow("S", "V", "vax_rate","vaccination")
  , mp_per_capita_flow("V", "VI", "beta * kappa2 * (I + VI)/N","inc_V")
  , mp_per_capita_flow("I", "R", "gamma","recovery_S")
  , mp_per_capita_flow("VI", "VR", "gamma","recovery_V")
  , beta_thing ~ beta        
  , vax_thing  ~ vax_rate    
  , inc ~ inc_S + inc_V
  , cinc ~ cinc + inc
  , Vbox ~ V + VI + VR
)

default = list(
    beta     = 0.4    
  , gamma    = 0.1    
  , vax_rate = 0.0    
  , kappa = 1
  , kappa2 = 0.5
  , N        = 10000   
  , I        = 1
  , VI = 0
  , R        = 0      
  , VR        = 0      
  , V        = 0      
	, inc = 0
	, cinc = 0
	, Vbox = 0
)

initialize_state = list(S ~ N - I - R - V - VI - VR)

# model specification
spec = mp_tmb_model_spec(
  before = initialize_state
  , during = flows
  , default = default
)

rdsSave(spec)

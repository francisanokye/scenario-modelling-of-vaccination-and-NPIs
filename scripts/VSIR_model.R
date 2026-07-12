library(macpan2)
library(shellpipes)

flows = list(
  mp_per_capita_flow("S", "I", "beta * (I) / N","inc_S")
  , mp_per_capita_flow("S", "V", "vax_rate","vaccination")
  , mp_per_capita_flow("V", "I", "(1-ve)*beta*(I/N)","inc_V")
  , mp_per_capita_flow("I", "R", "gamma","recovery_S")
  , beta_thing ~ beta        
  , vax_thing  ~ vax_rate    
  , inc ~ inc_S + inc_V
  , cinc ~ cinc + inc
  , cincprop ~ cinc/N
  , cvax ~ cvax + vaccination
)

default = list(
    beta     = 0.4    
  , gamma    = 0.1    
  , vax_rate = 0.0    
  , ve = 0.5
  , N        = 10000   
  , I        = 1
  , VI = 0
  , R        = 0      
  , VR        = 0      
  , V        = 0      
	, inc = 0
	, cinc = 0
	, Vbox = 0
	, cvax = 0
)

initialize_state = list(S ~ N - I - R - V - VI - VR)

# model specification
spec = mp_tmb_model_spec(
  before = initialize_state
  , during = flows
  , default = default
)

rdsSave(spec)

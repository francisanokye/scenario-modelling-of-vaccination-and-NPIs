library(macpan2)
library(dplyr)

time_steps <- 100

## flow diagram specification
flows = list(
    mp_per_capita_flow("S", "I", "beta * I / N", "infection")
  , mp_per_capita_flow("I", "R", "gamma", "recovery")
  , mp_inflow("cumInc","infection" ,"lol")
)

## default values for quantities required to run simulations
default = list(
    beta = 0.2   ## transmission rate
  , gamma = 0.1  ## recovery rate
  , N = 100      ## total population size (constant in this model)
  , I = 1        ## initial number of infectious individuals
  , R = 0        ## initial number of recovered individuals
  , cumInc = 0
)

## compute the initial number of susceptible individuals
initialize_state = list(S ~ N - I - R)

## model specification
spec = mp_tmb_model_spec(
    before = initialize_state
  , during = flows
  , default = default
)

det_sim <- mp_trajectory(
	mp_simulator(model = spec
	, time = time_steps
	, outputs = c("infection","cumInc")
	)
)

cumInc <- (det_sim
	|> filter(matrix == "cumInc")
	|> select(-matrix)
)

calccInc <- (det_sim
	|> filter(matrix == "infection")
	|> select(-matrix)
	|> mutate(value = cumsum(c(value)))
)

all.equal(calccInc,cumInc)


suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(ggplot2)
  library(tidyr)
  library(zoo)
  library(ggthemes)
  library(patchwork)
  library(shellpipes)
  library(tidyverse)
  library(macpan2)
  library(grid)
})

theme_set(theme_bw())

set.seed(2025)
options(macpan2_log_dir = ".")
loadEnvironments()

newRc <- function(beta_vec){

	beta_df <- (tibble(time = seq_along(beta_vec)
			, date = burn_in_start + (time - 1)
			, beta_scenario = beta_vec
			)
		|> select(-time)
	)

	Rcdat <- (baseRc
		|> select(date,beta,Rc,Rc_lwr,Rc_upr)
		|> inner_join(beta_df)
		|> mutate(beta_ratio = beta_scenario/beta
			, new_beta = beta_ratio* beta
			, Rc = beta_ratio * Rc
			, Rc_lwr = beta_ratio * Rc_lwr
			, Rc_upr = beta_ratio * Rc_upr
		)
	)
	return(Rcdat)
}


newRcdat <- (newRc(rdsRead())
	|> mutate(scenario = pipeStar())
)

rdsSave(newRcdat)


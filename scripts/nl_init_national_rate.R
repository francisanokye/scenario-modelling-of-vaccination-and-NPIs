
library(tidyverse);theme_set(theme_bw())
library(zoo)
library(shellpipes)

loadEnvironments()

primary_series_target_date <- as.Date("2021-11-06")
primary_series_target <- round(0.813 * NL_POP)

use_new_V2_init <- FALSE

# ------------------------------------------------------------
# Load PHAC vaccination data
# ------------------------------------------------------------

dat <- csvRead()

Canada <- (dat 
  |> filter(prename == "Canada") 
  |> select(week_end, contains("numtotal")) 
  |> pivot_longer(
   -week_end
   , names_to = "type"
   , values_to = "value"
  )
  |> mutate(
   week_end = as.Date(week_end)
   , value = as.numeric(value)
  )
)

# ------------------------------------------------------------
# Canadian booster trajectory scaled to NL
# ------------------------------------------------------------

Canada_boosters <- (Canada 
  |> filter(type == "numtotal_additional") 
  |> filter(
    !is.na(week_end)
    , !is.na(value)
  )
  |> arrange(
    type
    , week_end
  ) 
  |> mutate(
    value = round(
      value * (
        NL_POP / CA_POP
      )
    )
  ) 
  |> transmute(
    NULL
    , type
    , date = week_end
    , value
  ) 
  |> filter(
    date >= booster_start
    & date <= analysis_end
  )
)

filldat <- data.frame(
  date = c(
    as.Date("2021-11-06")
    , as.Date("2021-12-11")
  )
  , value = c(
    1
    , 27322
  )
  , type = "numtotal_additional"
)

can_boosterdat <- (bind_rows(
  filldat
  , Canada_boosters
)
|> group_by(
  type
  , date
)
|> summarise(
  value = max(
    value
    , na.rm = TRUE
  )
  , .groups = "drop"
)
|> arrange(date)
|> mutate(
  value = cummax(value)
)
)

can_boosterframe <- (tibble(
  date = seq.Date(
    min(can_boosterdat$date)
    , max(can_boosterdat$date)
    , by = "day"
  )
  , type = "numtotal_additional"
)
|> left_join(
  can_boosterdat
  , by = c(
    "date"
    , "type"
  )
)
|> arrange(date)
|> mutate(
  cum_interp = round(
    na.approx(
      value
      , x = as.numeric(date)
      , na.rm = FALSE
    )
  )
)
|> filter(!is.na(cum_interp))
|> mutate(
  cum_interp = cummax(cum_interp)
  , daily_rate = pmax(
    cum_interp -
      lag(
        cum_interp
        , default = first(cum_interp)
      )
    , 0
  )
)
|> filter(
  date >= booster_start
  , date <= analysis_end
  , daily_rate > 0
)
|> select(
  date
  , cum_interp
  , daily_rate
)
)

new_booster <- (can_boosterframe 
  |> mutate(day = as.numeric(date - burn_in_start))
)

nl_init_national_rate = new_booster

print(new_booster)

saveEnvironment()
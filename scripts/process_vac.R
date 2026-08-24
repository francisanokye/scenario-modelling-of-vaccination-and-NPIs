
library(tidyverse);theme_set(theme_bw())
library(zoo)
library(shellpipes)

loadEnvironments()

#####################################################################
## Load PHAC vaccination data for NL
#####################################################################

dat <- csvRead()

nl <- (dat 
  |> filter(grepl("ewfoundland", prename, ignore.case = TRUE)) 
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
# Initial two-dose value at warm-up start
# ------------------------------------------------------------

V2_init <- (nl 
  |> filter(
    type == "numtotal_fully"
    , week_end == as.Date("2021-07-17")
  ) 
  |> pull(value)
)

if (
  length(V2_init) != 1
  || !is.finite(V2_init)
  || V2_init <= 0
  || V2_init >= NL_POP
) {
  stop(
    "Invalid V2_init obtained for 2021-07-17: "
    , paste(V2_init, collapse = ", ")
  )
}

nl_vax <- (nl 
  |> filter(
   type %in% c(
     "numtotal_additional"
     , "numtotal_fully"
   )
  ) 
  |> filter(
   !is.na(week_end)
   , !is.na(value)
  ) 
  |> arrange(
   type
   , week_end
  ) 
  |> select(
   type
   , week_end
   , value
  ) 
)

#####################################################################
## Process booster vaccinated doses
#####################################################################

filldat <- (data.frame(
  date = c(
    as.Date("2021-11-08")
    , as.Date("2022-01-11")
  )
  , value = c(1, round(0.209 * NL_POP))
  , type = "numtotal_additional"
)
)

boosterdat <- (nl_vax
  |> filter(type == "numtotal_additional")
  |> transmute(NULL
  , type
  , date = week_end
  , value
  )
)

boosterdat <- (bind_rows(
  filldat
  , boosterdat
)
|> filter(
  !is.na(date)
  , !is.na(value)
)
|> group_by(
  type
  , date
)
|> summarise(
  value = max(value, na.rm = TRUE)
  , .groups = "drop"
)
|> arrange(date)
|> mutate(
  value = cummax(value)
)
)

boosterframe <- (tibble(
  date = seq.Date(
    min(boosterdat$date)
    , max(boosterdat$date)
    , by = "day"
  )
  , type = "numtotal_additional"
)
|> left_join(
  boosterdat
  , by = c("date", "type")
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
        , default = 0
      )
    , 0
  )
)
|> filter(
  date >= booster_start
  , date <= analysis_end
)
|> select(
  date
  , cum_interp
  , daily_rate
)
)

if (nrow(boosterframe) == 0) {
  stop(
    "No booster vaccination trajectory was produced."
  )
}

if (any(!is.finite(boosterframe$daily_rate))) {
  stop(
    "The booster vaccination trajectory contains non-finite values."
  )
}

if (any(boosterframe$daily_rate < 0)) {
  stop(
    "The booster vaccination trajectory contains negative daily rates."
  )
}

if (any(diff(boosterframe$cum_interp) < 0)) {
  stop(
    "The cumulative booster vaccination trajectory is not monotonic."
  )
}

#####################################################################
## Process double vaccinated doses
#####################################################################

nl_double_target_date <- as.Date("2021-11-06")
nl_double_target <- 423344

doublevaxdat <- (nl_vax 
   |> filter(type == "numtotal_fully") 
   |> transmute(NULL
    , type
    , date = week_end
    , value
   ) 
   |> filter(
     !is.na(date)
     , !is.na(value)
     , date != nl_double_target_date
   )
   |> bind_rows(
     data.frame(
       type = "numtotal_fully"
       , date = nl_double_target_date
       , value = nl_double_target
     )
   )
   |> group_by(
     type
     , date
   )
   |> summarise(
     value = max(value, na.rm = TRUE)
     , .groups = "drop"
   )
   |> arrange(date) 
   |> mutate(
     value = cummax(value)
   )
)

doublevaxframe <- (tibble(
  date = seq.Date(
    min(doublevaxdat$date)
    , max(doublevaxdat$date)
    , by = "day"
  )
  , type = "numtotal_fully"
) 
|> left_join(
  doublevaxdat
  , by = c("date", "type")
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
  date >= burn_in_start
  , date <= analysis_end
)
|> select(
  date
  , cum_interp
  , daily_rate
)
)

nl_double_target_check <- (doublevaxframe
  |> filter(date == nl_double_target_date)
  |> pull(cum_interp)
)

if (
  length(nl_double_target_check) != 1
  || !is.finite(nl_double_target_check)
  || nl_double_target_check != nl_double_target
) {
  stop(
    "NL fully vaccinated total on 2021-11-06 is "
    , paste(nl_double_target_check, collapse = ", ")
    , "; expected "
    , nl_double_target
    , "."
  )
}


saveEnvironment()
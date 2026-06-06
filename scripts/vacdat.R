# library(tidyverse)
# library(shellpipes)
# 
# loadEnvironments()
# 
# # number of days before the first data point to start simulations
# offset0 <- 150
# booster_offset <- 107
# dat <- csvRead()
# day0 <- as.Date("2021-12-17") - lubridate::days(offset0)
# 
# ## filtering just the NL section
# 
# nl <- (dat
#        |> filter(grepl("ewfoundland",prename)) # "ewfoundland"
#        |> select(week_end,contains("numtotal"))
#        |> mutate( numtotal_2nd_additional = ifelse(is.na(numtotal_2nd_additional),0,numtotal_2nd_additional)
#                   , additional_tot = numtotal_1additional + numtotal_2nd_additional
#        )
#        |> pivot_longer(-week_end,names_to="type",values_to="value")
# )
# 
# gg <- (ggplot(nl,aes(week_end,value))
#        + geom_point()
#        + facet_wrap(~type,scale="free")
# )
# 
# 
# ## Grabbing 2 doses complete series
# 
# V2_init <- (nl
#             |> filter(type == "numtotal_fully")
#             |> filter(week_end == as.Date("2021-07-24"))
#             |> pull(value)
# )
# 
# print(V2_init)
# 
# ## Grabbing additional boosters
# 
# nl_additional <- (nl
#                   |> filter(type %in% c("numtotal_additional","numtotal_fully"))
#                   |> mutate(value = ifelse((week_end == as.Date("2021-12-18")) & (type == "numtotal_additional") ,10,value))
#                   |> filter(!is.na(value))
#                   |> arrange(type)
#                   |> group_by(type)
#                   |> mutate(daydiff = diff(c(week_end,0))
#                             , boosterdiff = diff(c(value,0))
#                             , daily_rate = boosterdiff/as.numeric(daydiff)
#                             , days = as.numeric(week_end - day0)
#                             , week_end = week_end
#                   )
#                   |> select(days,daily_rate,week_end)
#                   |> filter(between(days,0,159 + offset0))
#                   #	|> filter(between(days,-20,159))
#                   
# )
# 
# #print(nl_additional,n=Inf)
# saveVars(nl_additional,V2_init, offset0, booster_offset)







library(tidyverse)
library(shellpipes)

loadEnvironments()

# number of days before the first data point to start simulations
offset0 <- 150
booster_offset <- 107
NL_POP <- 510550
CA_POP <- 36991981

day0 <- as.Date("2021-12-15") - lubridate::days(offset0)

dat <- csvRead()

# ------------------------------------------------------------
# NL vaccination data
# ------------------------------------------------------------

nl <- dat |>
  filter(grepl("ewfoundland", prename, ignore.case = TRUE)) |>
  select(week_end, contains("numtotal")) |>
  mutate(
    numtotal_2nd_additional = ifelse(
      is.na(numtotal_2nd_additional),
      0,
      numtotal_2nd_additional
    ),
    numtotal_additional = numtotal_1additional + numtotal_2nd_additional
  ) |>
  pivot_longer(
    -week_end,
    names_to = "type",
    values_to = "value"
  )

# ------------------------------------------------------------
# Canada vaccination data
# ------------------------------------------------------------

can <- dat |>
  filter(grepl("^Canada$", prename, ignore.case = TRUE)) |>
  select(week_end, contains("numtotal")) |>
  mutate(
    numtotal_2nd_additional = ifelse(
      is.na(numtotal_2nd_additional),
      0,
      numtotal_2nd_additional
    ),
    numtotal_additional = numtotal_1additional + numtotal_2nd_additional
  ) |>
  pivot_longer(
    -week_end,
    names_to = "type",
    values_to = "value"
  )

# ------------------------------------------------------------
# Initial two-dose values at warm-up start
# ------------------------------------------------------------

V2_init <- nl |>
  filter(
    type == "numtotal_fully",
    week_end == as.Date("2021-07-24")
  ) |>
  pull(value)

new_V2_init <- can |>
  filter(
    type == "numtotal_fully",
    week_end == as.Date("2021-07-24")
  ) |>
  summarise(
    value = value * (NL_POP / CA_POP)
  ) |>
  pull(value) |>
  round()

# ------------------------------------------------------------
# NL actual vaccination inputs
# ------------------------------------------------------------

nl_vax <- nl |>
  filter(type %in% c("numtotal_additional", "numtotal_fully")) |>
  mutate(
    value = ifelse(
      week_end == as.Date("2021-12-18") &
        type == "numtotal_additional",
      10,
      value
    )
  ) |>
  filter(!is.na(value)) |>
  arrange(type, week_end) |>
  group_by(type) |>
  mutate(
    daydiff = as.numeric(diff(c(week_end, 0))),
    boosterdiff = diff(c(value, 0)),
    daily_rate = boosterdiff / daydiff,
    days = as.numeric(week_end - day0)
  ) |>
  ungroup() |>
  select(type, days, daily_rate, week_end, value) |>
  filter(between(days, 0, 159 + offset0))

# ------------------------------------------------------------
# Canada national vaccination inputs converted to NL equivalent
# ------------------------------------------------------------

can_vax_nl_equiv <- can |>
  filter(type %in% c("numtotal_additional", "numtotal_fully")) |>
  filter(!is.na(value)) |>
  arrange(type, week_end) |>
  group_by(type) |>
  mutate(
    value_nl_equiv = value * (NL_POP / CA_POP),
    daydiff = as.numeric(diff(c(week_end, 0))),
    boosterdiff = diff(c(value_nl_equiv, 0)),
    daily_rate = boosterdiff / daydiff,
    days = as.numeric(week_end - day0)
  ) |>
  ungroup() |>
  select(type, days, daily_rate, week_end, value_nl_equiv) |>
  filter(between(days, 0, 159 + offset0))

# ------------------------------------------------------------
# Quick checks
# ------------------------------------------------------------

print(
  tibble(
    V2_init_NL_actual = V2_init,
    V2_init_Canada_scaled_to_NL = new_V2_init,
    NL_prop = V2_init / NL_POP,
    Canada_scaled_prop = new_V2_init / NL_POP
  )
)

saveVars(
  nl_vax,
  can_vax_nl_equiv,
  V2_init,
  new_V2_init,
  offset0,
  booster_offset
)
## What is going on here? We are first fitting to the broken/piecewise serology data and then make it daily to fit in macpan2?
library(tidyverse)
library(splines)
library(shellpipes)

loadEnvironments()

serodat <- rdsRead()

offset <- 150
omicron_start <- as.Date("2021-12-15")
burn_in_start <- omicron_start - offset
analysis_end  <- "2022-05-22"

newdat <- data.frame(time = (min(serodat$time)-1):max(serodat$time))

dat <- (
  serodat
  |> transmute(NULL
      , date
      , time
      , value
      , cumInc = value*params[["N"]]
      )
  )

mm <- lm(cumInc ~ bs(time,6), data = dat)

newdat$predcumInc <- predict(mm,newdata=newdat)

gg <- (ggplot()
       + geom_point(data=dat,aes(x=time,y=cumInc),size=2)
       + geom_line(data=newdat,aes(x=time,y=predcumInc),color="red")
       )

print(gg)

incdat <- (
  newdat
  |> mutate(predInc = round(diff(c(12800,predcumInc))) 
      , date = as.Date(omicron_start) - offset - 1 + time
      )
  )

print(plot(x=incdat$time, y=incdat$predInc))

# convert to seroprevalence matrix format for modeling
fitdat <- (
  incdat
  |>transmute(NULL
    , date
    , time = as.numeric(date - as.Date(omicron_start)) + 1 + offset
    , matrix = "newR"
    , value  = predInc
    ) 
  |> filter(between(date, as.Date(omicron_start), as.Date(analysis_end)))
)

print(fitdat)
fitdat$value[1] <- fitdat$value[2]

rdsSave(fitdat)


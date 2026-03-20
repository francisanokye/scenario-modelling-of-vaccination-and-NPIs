
library(tidyverse)
library(splines)
library(shellpipes)

loadEnvironments()

serodat <- rdsRead()


start_date <- "2021-12-15"
last_date  <- "2022-07-15"  # "2022-05-26"

newdat <- data.frame(time = (min(serodat$time)-1):max(serodat$time))
#newdat <- data.frame(time = 1:400)

dat <- (serodat
        |> transmute(NULL
                     , date
                     , time
                     , value
                     , cumInc = value*params[["N"]]
        )
)

print(dat)

mm <- lm(cumInc ~ bs(time,6)
         , data = dat
)

print(summary(mm))

newdat$predcumInc <- predict(mm,newdata=newdat)

print(dat)

gg <- (ggplot()
       + geom_point(data=dat,aes(x=time,y=cumInc),size=2)
       + geom_line(data=newdat,aes(x=time,y=predcumInc),color="red")
)

print(gg)

incdat <- (newdat
           |> mutate(predInc = round(diff(c(12800,predcumInc))) ## 12780 hacking the first point
                     , date = as.Date(start_date) - offset0 - 1 + time
           )
)

print(plot(x=incdat$time, y=incdat$predInc))

# convert to seroprevalence matrix format for modeling
fitdat <- (incdat
           |>transmute(NULL
                       , date
                       , time = as.numeric(date - as.Date(start_date)) + 1 + offset0
                       , matrix = "newR"
                       , value  = predInc
           ) 
           |> filter(between(date, as.Date(start_date), as.Date(last_date)))
)

print(fitdat)

rdsSave(fitdat)


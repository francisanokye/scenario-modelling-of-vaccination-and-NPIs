suppressPackageStartupMessages({
  
  library(dplyr)
  library(tidyr)
  library(macpan2)
  library(shellpipes)
  
})


loadEnvironments()



snapshot <- lapply(
  names(vax_specs),
  
  function(name){
    
    
    traj <- vax_specs[[name]] |>
      
      mp_simulator(
        time_steps=sim_horizon,
        outputs=c(
          "S1",
          "V2",
          "V3",
          "cum_booster"
        )
      ) |>
      
      mp_trajectory()
    
    
    
    start <- booster_start_dates[[name]]
    
    state <- traj |>
      
      filter(
        date == start-1
      ) |>
      
      select(
        matrix,
        value
      ) |>
      
      pivot_wider(
        names_from=matrix,
        values_from=value
      )
    
    
    
    tibble(
      
      scenario=name,
      
      no_primary =
        100*state$S1/NL_POP,
      
      
      primary_only =
        100*state$V2/NL_POP,
      
      
      boosted =
        100*state$V3/NL_POP,
      
      
      final_booster =
        100*max(
          traj$value[
            traj$matrix=="cum_booster"
          ],
          na.rm=TRUE
        )/NL_POP
      
    )
    
    
  }
  
) |>
  
  bind_rows()



write.csv(
  snapshot,
  "../output/vaccination_snapshot.csv",
  row.names=FALSE
)


print(snapshot)

quit()
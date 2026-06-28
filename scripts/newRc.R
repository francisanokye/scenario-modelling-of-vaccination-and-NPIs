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

print(newRcdat)

gg <- (ggplot(newRcdat, aes(date, y= Rc))
	+ geom_line()
	+ geom_ribbon(aes(ymin=Rc_lwr, ymax=Rc_upr),alpha=0.2)
)

print(gg)

rdsSave(newRcdat)

quit()
ALS4only_Rc <- newRc(beta_ALS4only,scenario="strict")

comboRc <- (baseRc
            |> select(date,Rc,Rc_lwr,Rc_upr)
            |> mutate(scenario = "base")
            |> bind_rows(noALS_Rc)
            |> bind_rows(ALS4only_Rc)
)

print(comboRc)

scenario_levels <- c("base", "noALS", "strict")
scenario_labels <- c(
  "base"   = "Actual NPIs",
  "noALS"  = "No NPIs",
  "strict" = "Strict NPIs"
)
scenario_colors <- c(
  "base"   = "black",
  "noALS"  = "red",
  "strict" = "#00C5CD"
)

gg <- (ggplot(comboRc, aes(date, Rc, color = scenario, fill = scenario))
       + geom_ribbon(aes(ymin = Rc_lwr, ymax = Rc_upr), alpha = 0.2, colour = NA)
       + geom_line(linewidth = 2)
       + geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1.2)
       + geom_vline(
         xintercept = intervention_removed,
         linetype = "dashed",
         linewidth = 3
       )
       + scale_color_manual(
         name   = NULL,
         breaks = scenario_levels,
         values = scenario_colors,
         labels = scenario_labels
       )
       + scale_fill_manual(
         name   = NULL,
         breaks = scenario_levels,
         values = scenario_colors,
         labels = scenario_labels
       )
       + scale_x_date(
         limits = c(omicron_start, analysis_end),
         date_breaks = "2 weeks",
         date_labels = "%d %b %Y",
         expand = c(0.01, 0)
       )
       + labs(
         title = "Reproduction number for counterfactual scenarios (timevarying prop)",
         y = "Reproduction number",
         x = NULL,
         tag = "C"
       )
       + theme_clean()
       + theme(
         axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
         axis.text.y = element_text(size = 30),
         axis.title.y = element_text(size = 30, color = "black"),
         legend.position = c(0.85, 0.65),
         legend.justification = c(0.5, 0.5),
         legend.direction = "vertical",
         legend.box = "vertical",
         legend.title = element_text(size = 25),
         legend.text = element_text(size = 25),
         legend.background = element_blank(),
         plot.title = element_text(size = 30, hjust = 0.5, face = "plain"),
         plot.background = element_blank(),
         plot.tag = element_text(size = 35, face = "bold"),
         plot.tag.position = c(0.01, 0.98)
       )
)

#png("../figures/reprod_numb.png", width = 5000, height = 2000, res = 300, bg = "white",type = "cairo")
print(gg)
#dev.off()

baseRc_with_scenarios <- (baseRc
                          |> inner_join(beta_vec_to_df(beta_noALS)    |> select(date, beta_noALS_value = beta_scenario), by = "date")
                          |> inner_join(beta_vec_to_df(beta_ALS4only) |> select(date, beta_ALS4only_value = beta_scenario), by = "date")
)

gg1 <- (ggplot(baseRc_with_scenarios, aes(date, beta))
        + geom_ribbon(
          aes(ymin = beta_lwr, ymax = beta_upr, fill = "Actual NPIs"),
          alpha = 0.2, colour = NA
        )
        + geom_line(aes(colour = "Actual NPIs"), linewidth = 2)
        + geom_line(
          aes(y = beta_noALS_value, colour = "No NPIs"),
          linewidth = 2
        )
        + geom_line(
          aes(y = beta_ALS4only_value, colour = "Strict NPIs"),
          linewidth = 2
        )
        + geom_vline(
          xintercept = intervention_removed,
          linetype = "dashed",
          linewidth = 3
        )
        + scale_color_manual(
          name = NULL,
          breaks = c("Actual NPIs", "No NPIs", "Strict NPIs"),
          values = c(
            "Actual NPIs" = "black",
            "No NPIs" = "red",
            "Strict NPIs" = "#00C5CD"
          )
        )
        + scale_fill_manual(
          name = NULL,
          breaks = c("Actual NPIs", "No NPIs", "Strict NPIs"),
          values = c(
            "Actual NPIs" = "black",
            "No NPIs" = "red",
            "Strict NPIs" = "#00C5CD"
          ),
          guide = "none"
        )
        + scale_x_date(
          limits = c(omicron_start, analysis_end),
          date_breaks = "2 weeks",
          date_labels = "%d %b %Y",
          expand = c(0.01, 0)
        )
        + labs(
          title = "Transmission rate for counterfactual scenarios",
          y = "Transmission rate",
          x = NULL,
          tag = "A"
        )
        + theme_clean()
        + theme(
          axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 30),
          axis.title.y = element_text(size = 30, color = "black"),
          legend.position = c(0.85, 0.5),
          legend.justification = c(0.5, 0.5),
          legend.direction = "vertical",
          legend.box = "vertical",
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 25),
          legend.background = element_blank(),
          plot.title = element_text(size = 30, hjust = 0.5, face = "plain"),
          plot.background = element_blank(),
          plot.tag = element_text(size = 35, face = "bold"),
          plot.tag.position = c(0.01, 0.98)
        )
)
#png("../figures/transmission_rate.png", width = 5000, height = 2000, res = 300, bg = "white", type = "cairo")
print(gg1)
#dev.off()

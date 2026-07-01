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

loadEnvironments()


comboRc <- bind_rows(rdsReadList())

print(comboRc)

scenario_levels <- c("base_beta", "policy_neutral", "beta_als4")
scenario_labels <- c(
  "base_beta"   = "Actual NPIs",
  "policy_neutral"  = "No NPIs",
  "beta_als4" = "Strict NPIs"
)
scenario_colors <- c(
  "base_beta"   = "black",
  "policy_neutral"  = "red",
  "beta_als4" = "#00C5CD"
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

print(gg)

png("../figures/reprod_numb.png", width = 5000, height = 2000, res = 300, bg = "white",type = "cairo")
print(gg)
dev.off()

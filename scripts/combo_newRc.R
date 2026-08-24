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

scenario_levels <- c("base", "no_npis", "strict_npis")
scenario_labels <- c("base" = "Actual NPIs", "no_npis" = "No NPIs", "strict_npis" = "Strict NPIs")
scenario_colors <- c("base" = "black", "no_npis"  = "red","strict_npis" = "#00C5CD")

print(comboRc[c("date","beta","Rc","scenario")], n = Inf)
write.csv(comboRc[c("date","beta","Rc","scenario")], "../output/rc.csv", row.names = FALSE)

gg <- (ggplot(comboRc, aes(date, Rc, color = scenario, fill = scenario))
  + geom_ribbon(
   data = filter(comboRc, scenario == "base", date <= intervention_removed),
   aes(ymin = Rc_lwr, ymax = Rc_upr), alpha = 0.2, colour = NA
  )
  + geom_line(linewidth = 2)
  + geom_hline(yintercept = 1, linetype = "dashed", linewidth = 1.2)
  + geom_vline(xintercept = intervention_removed, linetype = "dashed", linewidth = 3)
  + scale_color_manual(name = NULL, breaks = scenario_levels, values = scenario_colors, labels = scenario_labels)
  + scale_fill_manual(name = NULL, breaks = scenario_levels, values = scenario_colors, labels = scenario_labels)
  + scale_x_date(limits = c(omicron_start, analysis_end), date_breaks = "2 weeks", date_labels = "%d %b %Y", expand = c(0.01, 0))
  + labs(
   title = "Reproduction number for counterfactual scenarios",
   y = "Reproduction number", x = NULL, tag = "A"
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

# png("../figures/reprod_numb.png", width = 5000, height = 2000, res = 300, bg = "white",type = "cairo")
# print(gg)
# dev.off()

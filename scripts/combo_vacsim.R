library(tidyverse);theme_set(theme_bw())
library(shellpipes)
library(ggthemes)

keep <- c("boosterbox" ,"cum_inc_total", "cum_hosp_total", "cum_death_total",
          "inc_booster", "inc_leq2dose", "inc", 
          "inc_leq2dose", "inc_booster", "cum_inc_leq2dose", 
          "cum_inc_booster","S1", "V2", "V3", "cum_booster"
          )

combosims <- (bind_rows(rdsReadList())
	|> filter(matrix %in% keep)
)

print(combosims)

gg <- (ggplot(combosims,aes(time,value,color=vac))
	+ geom_line(linewidth =2)
	+ scale_x_continuous(breaks=c(0, 50, 100, 150, 200, 250, 300, 350))
	+ facet_grid(matrix~beta,scale="free")
	+ theme_clean() 
	+ theme(
	    axis.text.x        = element_text(size = 35, hjust = 0.5),
	    axis.text.y        = element_text(size = 35),
	    axis.title.y       = element_text(size = 35),
	    legend.position    = "bottom",
	    legend.justification = "center",
	    legend.direction   = "vertical",
	    strip.text         = element_text(size = 35),
	    strip.background   = element_blank(),
	    legend.text        = element_text(size = 35),
	    legend.background  = element_blank(),
	    panel.grid.major.x = element_blank(),
	    plot.background    = element_blank(),
	    plot.margin        = margin(12, 35, 12, 12)
	  )
)

# png("../figures/combovacsim_plot.png", width = 13500, height = 8000, res = 300, bg = "white", type = "cairo")
# print(gg)
# dev.off()

rdsSave(combosims)

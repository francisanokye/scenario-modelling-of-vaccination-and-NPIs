library(tidyverse);theme_set(theme_bw())
library(shellpipes)

keep <- c("beta_thing", "boosterbox" ,"cum_inc_total")

combosims <- (bind_rows(rdsReadList())
	|> filter(matrix %in% keep)
)

print(combosims)


gg <- (ggplot(combosims,aes(time,value,color=vac))
	+ geom_line()
	+ scale_x_continuous(breaks=c(0,50,100,150,200,250,300,350))
	+ facet_grid(matrix~beta,scale="free")
	+ theme(legend.position="bottom")
)

print(gg)


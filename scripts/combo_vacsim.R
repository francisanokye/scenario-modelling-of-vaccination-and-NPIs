library(tidyverse);theme_set(theme_bw())
library(shellpipes)

combosims <- bind_rows(rdsReadList())

print(combosims
	|> filter(time == 530)
	|> filter(matrix == "boosterbox")
)

gg <- (ggplot(combosims,aes(time,value,color=vac))
	+ geom_line()
	+ scale_x_continuous(breaks=c(0,50,100,150,200,250,300,350))
	+ facet_wrap(~matrix,scale="free",ncol=2)
	+ theme(legend.position="bottom")
)

print(gg)


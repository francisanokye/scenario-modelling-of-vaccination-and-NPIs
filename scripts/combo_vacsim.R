library(tidyverse);theme_set(theme_bw())
library(shellpipes)

combosims <- bind_rows(rdsReadList())

print(head(combosims))

gg <- (ggplot(combosims,aes(time,value,color=vac))
	+ geom_line()
	+ facet_wrap(~matrix,scale="free",nrow=2)
	+ theme(legend_position="bottom")
)

print(gg)


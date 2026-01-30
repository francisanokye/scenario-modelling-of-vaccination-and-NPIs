suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(forcats)
  library(scales)
  library(ggthemes)
  library(shellpipes)
})
options(scipen = 9999)

als_vax <- csvRead("comb_als_vax_bar.csv")
early_booster <- csvRead("early_booster_bar.csv")
natl_vax_rate <- csvRead("can_natl_vax_rate_bar.csv")

merged <- rbind(als_vax, early_booster, natl_vax_rate)

write.csv(merged, "../output/compare_bar.csv", row.names = FALSE)

df2 <- merged %>%
  mutate(
    vax = factor(vax, levels = c("0%","25%","50%","75%","100%")),
    als_scen = factor(als_scen, levels = c("No-ALS Only", "ALS-2 Only", "Implemented ALS", "ALS-4 Only")),
    period = factor(period, levels = c("Up to March 14, 2022", "March 15 - Dec 02, 2022")),
    type = fct_recode(
      factor(type), 
      "Combined ALS & vaccination" = "als & vax",
      "Early boosters"            = "early boosters",
      "National vaccination rate" = "natl rate vax"
    ),
    type = fct_relevel(type, "Combined ALS & vaccination", "Early boosters", "Nationalvaccination rate")
  )

# subset only data for "Up to March 14, 2022"

df2 <- df2 |>
  dplyr::filter( period == "Up to March 14, 2022")
  
# ---- plot ----
p <- ggplot(df2, aes(x = vax, y = total_infections, fill = type)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.78, colour = "grey20", linewidth = 0.25) +
  facet_grid(period ~ als_scen, scales = "free_y") +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(
    values = c(
      "Combined ALS & vaccination" = "#0072B2",
      "Early boosters"             = "#009E73",
      "National vaccination rate"  = "#D55E00"
    )
  )+
  labs(title="Total infections under different alert level policies and vaccination rates",
    x = "Vaccination scale (relative to observed daily rate)",
    y = "Total infections",
    fill = NULL
  ) +
  theme_clean() +
  theme(
    axis.text.x       = element_text(size = 25, angle = 0, hjust = 0.5),
    axis.text.y       = element_text(size = 25),
    axis.title.x      = element_text(size = 25),
    axis.title.y      = element_text(size = 25),
    plot.title        = element_text(size = 30, hjust = 0.5, face = "plain"),
    strip.text        = element_text(size = 25, face = "plain"),
    strip.background = element_rect(fill = "grey90", color = NA),
    legend.position   = "top",
    legend.direction  = "horizontal",
    legend.box.just   = "left",
    legend.background = element_blank(),
    legend.text       = element_text(size = 25),
    legend.margin = margin(t = 2, r = 2, b = 2, l = 2),
    plot.background   = element_blank(),
    plot.margin = margin(20, 20, 20, 40),
    panel.spacing.x   = unit(1.2, "cm"),
    panel.grid.major.x = element_blank()
  )

# png("../figures/compare_bars.png", width = 7500, height = 4000, res = 300, bg = "white", type = "cairo")
# print(p)
# dev.off()

png("../figures/compare_bars1.png", width = 7500, height = 3000, res = 300, bg = "white", type = "cairo")
print(p)
dev.off()
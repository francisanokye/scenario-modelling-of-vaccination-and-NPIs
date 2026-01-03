suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(forcats)
  library(scales)
  library(shellpipes)
})


als_vax <- csvRead("als_vax_bar.csv")
early_booster <- csvRead("early_booster_bar.csv")
natl_vax_rate <- csvRead("natl_vax_rate_bar.csv")

merged <- rbind(als_vax, early_booster, natl_vax_rate)

df2 <- merged %>%
  mutate(
    vax = factor(vax, levels = c("0%","25%","50%","75%","100%")),
    als_scen = factor(als_scen, levels = c("No-ALS Only", "ALS-2 Only", "Implemented ALS", "ALS-4 Only")),
    period = factor(period, levels = c("Up to March 14, 2022", "March 15 - Dec 02, 2022")),
    type = fct_recode(
      factor(type), 
      "Observed ALS & vax" = "als_vax",
      "Early boosters"            = "early boosters",
      "National vaccination rate" = "natl rate vax"
    ),
    type = fct_relevel(type, "Observed ALS & vax", "Early boosters", "Nationalvaccination rate")
  )

# ---- plot ----
p <- ggplot(df2, aes(x = vax, y = total_infections, fill = type)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.78, colour = "grey20", linewidth = 0.25) +
  facet_grid(period ~ als_scen, scales = "free_y") +
  scale_y_continuous(labels = comma) +
  labs(
    x = "Vaccination scale (relative to observed NL daily rate)",
    y = "Total infections",
    fill = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.x       = element_text(size = 12,  hjust = 0.5, margin = margin(t = 8)),
    axis.text.y       = element_text(size = 12),
    axis.title.x      = element_text(size = 15),
    axis.title.y      = element_text(size = 15),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.text = element_text(size = 15),
    legend.margin = margin(t = 2, r = 2, b = 2, l = 2),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 15),
    # strip.background  = element_blank(),
    axis.title = element_text(face = "bold"),
    panel.spacing     = unit(0.50, "cm"),
    plot.background   = element_blank()
  )

png("../figures/compare_bars.png", width = 4500, height = 2500, res = 300, bg = "white", type = "cairo")
print(p)
dev.off()
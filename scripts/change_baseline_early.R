suppressPackageStartupMessages({
  library(tidyverse)
  library(shellpipes)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(forcats)
  library(ggthemes)
  library(scales)
  library(grid)   
})

options(scipen = 9999999)

## ------------------------- load & clean ------------------------------------

baseline_file <- csvRead("comb_vax_als_bars.csv")     
early_file    <- csvRead("early_boosters_bars.csv")   

als_vax_raw <- bind_rows(baseline_file, early_file)

# standardize period labels
als_vax_clean <- als_vax_raw |>
  mutate(
    period = recode(
      period,
      "Up to March 14, 2022"     = "Dec 15–Mar 14",
      "March 15 – Nov 30, 2022"  = "Mar 15–Nov 30",
      "March 15 - Nov 30, 2022"  = "Mar 15–Nov 30"
    )
  )

## --------------------- normal vs early side by side ------------------------

als_vax_wide <- als_vax_clean |>
  select(als_scen, vaccination, period, type, total_infections) |>
  pivot_wider(
    names_from  = type,               
    values_from = total_infections
  )

# difference > 0 means more infections under baseline than early boosters
als_vax_wide <- als_vax_wide |>
  mutate(diff = normal - early)

## ------------------- period-specific and total diffs -----------------------

period_diffs <- als_vax_wide |>
  select(als_scen, vaccination, period, diff)

total_diffs <- period_diffs |>
  group_by(als_scen, vaccination) |>
  summarise(diff = sum(diff, na.rm = TRUE), .groups = "drop") |>
  mutate(period = "Total (Dec 15–Nov 30)")

als_levels <- c("No-ALS Only", "ALS-2 Only", "Implemented ALS", "ALS-4 Only")

plot_df <- bind_rows(period_diffs, total_diffs) |>
  mutate(
    period = factor(
      period,
      levels = c("Dec 15–Mar 14",
                 "Mar 15–Nov 30",
                 "Total (Dec 15–Nov 30)")
    ),
    vaccination = factor(
      vaccination,
      levels = c("0%", "25%", "50%", "75%", "100%")
    ),
    als_scen = factor(als_scen, levels = als_levels)
  )

## ----------------------------- colours -------------------------------------

# okabe–ito style: colourblind-safe and high contrast
period_palette <- c(
  "Dec 15–Mar 14"         = "#0072B2",  
  "Mar 15–Nov 30"         = "#009E73",  
  "Total (Dec 15–Nov 30)" = "#D55E00"   
)

## ----------------------------- plot ----------------------------------------

p <- ggplot(plot_df,aes(x = vaccination, y = diff, fill = period)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
  geom_col(position = position_dodge(width = 0.75),width    = 0.6) +
  facet_wrap(~ als_scen, nrow = 2, scales = "free") +
  scale_fill_manual(values = period_palette, name = "Period") +
  scale_y_continuous(labels = comma) +
  labs(
    x     = "Vaccination rates",
    y     = expression(Delta * " total infections (baseline - early boosters)"),
    title = "Difference in cumulative infections: baseline vs. early boosters") +
  theme_clean() +
  theme(
    axis.text.x      = element_text(size = 18, angle = 0, hjust = 0.5),
    axis.text.y      = element_text(size = 18),
    axis.title.x     = element_text(size = 18),
    axis.title.y     = element_text(size = 18),
    plot.title       = element_text(size = 20, hjust = 0.5, face = "plain"),
    strip.text       = element_text(size = 20, face = "plain"),
    legend.position  = "top",
    legend.direction = "horizontal",
    legend.box       = "horizontal",
    legend.text      = element_text(size = 16),
    legend.title     = element_blank(),
    legend.background= element_blank(),
    plot.background  = element_blank(),
    panel.spacing.x  = unit(1.2, "cm"),
    panel.spacing.y  = unit(0.8, "cm")
  )

print(p)


png("../figures/alsvax_early_total_diffs.png",width  = 4500, height = 2500,res    = 300,bg     = "white",type   = "cairo")
print(p)
dev.off()

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(scales)
  library(ggthemes)
  library(grid)
  library(shellpipes)
})

options(scipen = 9999999)

## ------------------------- load -------------------------------------------
baseline_file <- csvRead("comb_vax_als_bar.csv")       # "normal"
early_file    <- csvRead("early_boosters_bar.csv")   # "early"

baseline_file <- baseline_file %>% mutate(type = "normal")
early_file    <- early_file    %>% mutate(type = "early")

als_vax_raw <- bind_rows(baseline_file, early_file)

print(als_vax_raw)

## ------------------------- clean period labels ----------------------------
# normalize dash variants and whitespace once, then map reliably
als_vax_clean <- als_vax_raw %>%
  mutate(
    period = str_squish(period),
    period = str_replace_all(period, "[\u2012\u2013\u2014\u2212]", "-"), # en/em/minus → "-"
    period = str_replace_all(period, "–", "-"),
    period = str_replace_all(period, "\\s*-\\s*", " - "),
    period = case_when(
      str_detect(period, "^Up to March 14, 2022$") ~ "Dec 15, 2021 - Mar 14, 2022",
      str_detect(period, "^March 15\\s*-\\s*Dec 02, 2022$") ~ "Mar 15, 2022 - Dec 02, 2022",
      str_detect(period, "^March 15\\s*-\\s*Nov 30, 2022$") ~ "Mar 15, 2022 - Nov 30, 2022",
      TRUE ~ period
    )
  )

## ------------------------- wide + diffs -----------------------------------
als_vax_wide <- als_vax_clean %>%
  select(type, als_scen, vax, period, total_infections) %>%
  pivot_wider(
    names_from  = type,
    values_from = total_infections
  ) %>%
  mutate(
    diff = normal - early
  )

period_diffs <- als_vax_wide %>%
  select(als_scen, vax, period, diff)

total_diffs <- period_diffs %>%
  group_by(als_scen, vax) %>%
  summarise(diff = sum(diff, na.rm = TRUE), .groups = "drop") %>%
  mutate(period = "Total (Dec 15, 2021 - Dec 02, 2022)")

als_levels <- c("No-ALS Only", "ALS-2 Only", "Implemented ALS", "ALS-4 Only")

plot_df <- bind_rows(period_diffs, total_diffs) %>%
  mutate(
    vax = factor(vax, levels = c("0%", "25%", "50%", "75%", "100%")),
    als_scen = factor(als_scen, levels = als_levels),
    period = factor(
      period,
      levels = c(
        "Dec 15, 2021 - Mar 14, 2022",
        "Mar 15, 2022 - Dec 02, 2022",
        "Total (Dec 15, 2021 - Dec 02, 2022)"
      )
    )
  )

## ------------------------- plot -------------------------------------------
period_palette <- c(
  "Dec 15, 2021 - Mar 14, 2022"          = "#0072B2",
  "Mar 15, 2022 - Dec 02, 2022"          = "#009E73",
  "Total (Dec 15, 2021 - Dec 02, 2022)"  = "#D55E00"
)

p <- ggplot(plot_df, aes(x = vax, y = diff, fill = period)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
  geom_col(position = position_dodge(width = 0.75), width = 0.6) +
  facet_wrap(~als_scen, nrow = 2, scales = "free") +
  scale_fill_manual(values = period_palette, breaks = names(period_palette), name = NULL) +
  scale_y_continuous(labels = comma) +
  labs(
    x = "Vaccination rate",
    y = "Difference in total infections (baseline − early boosters)",
    title = "Difference in cumulative infections: baseline vs early booster schedule"
  ) +
  theme_clean() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    strip.text  = element_text(size = 14),
    panel.spacing.x = unit(1.2, "cm"),
    panel.spacing.y = unit(0.8, "cm")
  )

print(p)

png("../figures/alsvaxrate_diffs.png", width = 4500, height = 2500, res = 300, bg = "white", type = "cairo")
print(p)
dev.off()

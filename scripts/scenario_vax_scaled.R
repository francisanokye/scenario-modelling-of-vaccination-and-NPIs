# vac_scenarios_locked_with_summaries.R
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(ggthemes)
  library(RColorBrewer); library(shellpipes); library(macpan2)
})

set.seed(2025)
options(macpan2_log_dir = ".")
loadEnvironments()

# ------------------------- config -------------------------
date0    <- as.Date("2021-12-15")
date_end <- as.Date("2022-05-22")
scales   <- c(1.00, 0.75, 0.50, 0.25)
labels   <- c("observed", "75%", "50%", "25%")# sprintf("%d%%", as.integer(round(scales * 100)))  
png_out  <- "../figures/Figure_vac_scen.png"
csv_out  <- "../output/vac_scenario_summary.csv"
dir.create("../figures", recursive = TRUE, showWarnings = FALSE)

# load saved inputs 
calibrator   <- readRDS("calibrator_opt.rds")  
seroprevdata <- readRDS("seroprevdata.rds")              
time_steps   <- max(seroprevdata$time)

# load spec with fitted parameters 
spec_fixed <- mp_optimized_spec(calibrator, "modified")

# keep pristine vaccination series to rescale cleanly
if (is.null(spec_fixed$default$double_vac_values_base))
  spec_fixed$default$double_vac_values_base  <- spec_fixed$default$double_vac_values
if (is.null(spec_fixed$default$booster_vac_values_base))
  spec_fixed$default$booster_vac_values_base <- spec_fixed$default$booster_vac_values

# simulate shared-parameter scenarios ----------------
simulate_scale <- function(scale, label, spec_fixed, time_steps, date0, date_end) {
  cf <- spec_fixed
  cf$default$double_vac_values  <- scale * cf$default$double_vac_values_base
  cf$default$booster_vac_values <- scale * cf$default$booster_vac_values_base
  
  cf |>
    mp_simulator(time_steps, outputs = "inc") |>
    mp_trajectory() |>
    dplyr::filter(matrix == "inc") |>
    dplyr::mutate(
      date   = date0 + (time - offset0),
      series = label
    ) |>
    dplyr::filter(date >= date0, date <= date_end) |>
    dplyr::select(date, inc = value, series)
}

sims_all <- dplyr::bind_rows(
  mapply(
    simulate_scale,
    scale = scales,
    label = labels,
    MoreArgs = list(
      spec_fixed = spec_fixed,
      time_steps = time_steps,
      date0 = date0,
      date_end = date_end
    ),
    SIMPLIFY = FALSE
  )
) |>
  dplyr::mutate(series = factor(series, levels = labels))  # keep 100%,75%,50%,25% order

# summaries: cumulative totals & peaks ----------------
summ_tab <- sims_all %>%
  dplyr::group_by(series) %>%
  dplyr::summarise(
    cumulative_infections = round(sum(inc, na.rm = TRUE)),
    peak_infections       = round(max(inc, na.rm = TRUE)),
    peak_date             = date[which.max(inc)],
    .groups = "drop"
  ) %>%
  dplyr::arrange(factor(series, levels = labels))

obs_row <- dplyr::filter(summ_tab, series == "observed")
summ_vs_obs <- summ_tab %>%
  dplyr::mutate(
    excess_cumulative = cumulative_infections - obs_row$cumulative_infections,
    peak_diff         = peak_infections - obs_row$peak_infections,
    peak_advance_days = as.integer(peak_date - obs_row$peak_date) # negative = earlier than 100%
  )

print(summ_vs_obs)

# save CSV 
#write.csv(summ_tab, file = csv_out, row.names = FALSE)

# distinct colours, mapped to percentages
pal <- setNames(RColorBrewer::brewer.pal(n = max(3, length(labels)), "Set1")[seq_along(labels)], labels)

p <- ggplot(sims_all, aes(x = date)) +
  geom_line(aes(y = inc, color = series), linewidth = 1.2, alpha = 1) +
  #geom_point(aes(y = inc, color = series), size = 1.6, alpha = 0.9) +
  scale_color_manual(name = "vaccination rate", values = pal, breaks = labels) +
  labs(y = "Estimated daily infections", x = NULL, title = "Effect of vaccination rate on estimated daily infections") +
  scale_x_date(expand = c(0, 0), limits = c(date0, date_end),date_breaks = "2 week", date_labels = "%b %d") +
  geom_vline(aes(xintercept = as.Date("2021-12-15")), colour = "gold4", linetype = 4, linewidth = 1.0) +
  geom_vline(aes(xintercept = as.Date("2022-01-03")), colour = "gold4", linetype = 4, linewidth = 1.0) +
  geom_vline(aes(xintercept = as.Date("2022-01-24")), colour = "gold4", linetype = 4, linewidth = 1.0) +
  geom_vline(aes(xintercept = as.Date("2022-02-25")), colour = "gold4", linetype = 4, linewidth = 1.0) +
  geom_vline(aes(xintercept = as.Date("2022-03-17")), colour = "gold4", linetype = 4, linewidth = 1.0) +
  theme_clean() +
  theme(
    axis.text.x   = element_text(size = 18),
    axis.title.x  = element_blank(),
    axis.text.y   = element_text(size = 18),
    axis.title.y  = element_text(size = 18),
    plot.title    = element_text(size = 20, hjust = 0.5, face = "plain"),
    legend.position  = c(0.19, 0.75),  
    legend.direction = "vertical",
    legend.box.just  = "left",
    legend.title     = element_text(size = 18, face = "plain"),
    legend.text      = element_text(size = 18, face = "plain"),
    legend.background= element_blank(),
    plot.background  = element_blank()
  )

print(p)

# png("../figures/Figure_vac_scen.png", width = 5000, height = 2500, res = 300, bg = "white", type = "cairo")
# p
# dev.off()
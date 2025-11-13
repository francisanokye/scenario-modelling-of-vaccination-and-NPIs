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
scales   <- c(0.75, 0.50, 0.25)
labels   <- c("75%", "50%", "25%") 
csv_out  <- "../output/vax_scenario_summary.csv"
dir.create("../figures", recursive = TRUE, showWarnings = FALSE)
offset0 <- 150

# load saved inputs 
calibrator   <- rdsRead("calibrate") 
seroprevdata <- rdsRead("seroprevdata")
time_steps   <- max(seroprevdata$time)

# load spec with fitted parameters 
spec_fixed <- mp_optimized_spec(calibrator, "modified")

# Ensure pristine baselines exist once
if (is.null(spec_fixed$default$double_vac_values_base))
  spec_fixed$default$double_vac_values_base  <- spec_fixed$default$double_vac_values
if (is.null(spec_fixed$default$booster_vac_values_base))
  spec_fixed$default$booster_vac_values_base <- spec_fixed$default$booster_vac_values
if (is.null(spec_fixed$integers$booster_vac_changepoints_base))
  spec_fixed$integers$booster_vac_changepoints_base <- spec_fixed$integers$booster_vac_changepoints

######################################### Baseline  ######################################### 

baseline <- (spec_fixed
                     |> mp_simulator(time_steps, "inc")
                     |> mp_trajectory()
                     |> dplyr::mutate(date = date0 + (time - offset0), series = "baseline") 
                     |> dplyr::filter(date >= date0, date <= date_end) 
                     |> dplyr::select(date, inc = value, series))


######################################### Early Booster Shots -28 days  ######################################### 

# clone calibrated spec for the counterfactual
cf_early28 <- spec_fixed

# shift booster timing
cp  <- as.integer(cf_early28$integers$booster_vac_changepoints_base) - 28L   
val <- as.numeric(cf_early28$default$booster_vac_values_base)

# clamp negatives to 0 (cannot start before t=0)
cp <- pmax(0L, cp)

# clean up changepoints and values and sort in ascending order
order_index <- order(cp)
cp  <- cp[order_index]
val <- val[order_index]

# write back into the scenario spec
cf_early28$integers$booster_vac_changepoints <- cp
cf_early28$default$booster_vac_values        <- val

# --- Simulate and summarize (same pattern as your other scenarios) ---
early28_boosters <- (cf_early28
                |> mp_simulator(time_steps, "inc")
                |> mp_trajectory()
                |> dplyr::mutate(date = date0 + (time - offset0), series = "early boosters") 
                |> dplyr::filter(date >= date0, date <= date_end) 
                |> dplyr::select(date, inc = value, series))

######################################### Scaled Vaccination Rates ######################################### 

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
  dplyr::mutate(series = factor(series, levels = labels))  

######################################### No Booster ######################################### 

### No Booster 
no_booster_spec <- spec_fixed
no_booster_spec$default$booster_vac_values <- 0 * no_booster_spec$default$booster_vac_values_base 
no_booster_inc = (no_booster_spec
                  |> mp_simulator(time_steps, "inc")
                  |> mp_trajectory()
                  |> dplyr::mutate(date   = date0 + (time - offset0), series = "no boosters") 
                  |> dplyr::filter(date >= date0, date <= date_end) 
                  |> dplyr::transmute(date, inc = value, series)
) 

######################################### No Vaccination ######################################### 

no_vax_spec <- spec_fixed
no_vax_spec$default$double_vac_values  <- 0 * no_vax_spec$default$double_vac_values_base
no_vax_spec$default$booster_vac_values <- 0 * no_vax_spec$default$booster_vac_values_base 

no_vax_inc = (no_vax_spec
              |> mp_simulator(time_steps, "inc")
              |> mp_trajectory() 
              |> dplyr::mutate(date   = date0 + (time - offset0), series = "no vaccination") 
              |> dplyr::filter(date >= date0, date <= date_end) 
              |> dplyr::transmute(date, inc = value, series)
              ) 

######################################### No Booster & No Vaccination ######################################### 

df_all <- dplyr::bind_rows(sims_all, baseline, no_booster_inc, no_vax_inc, early28_boosters)
df_all <- df_all |>
  mutate(inc = as.integer(inc))

######################################### Summaries  ######################################### 

labels   <- c("baseline", "75%", "50%", "25%","no vaccination","no boosters", "early boosters") 

# summaries: cumulative totals & peaks ----------------
summ_tab <- df_all %>%
  dplyr::group_by(series) %>%
  dplyr::summarise(
    cumulative_infections = round(sum(inc, na.rm = TRUE)),
    peak_infections       = round(max(inc, na.rm = TRUE)),
    peak_date             = date[which.max(inc)],
    .groups = "drop"
  ) %>%
  dplyr::arrange(factor(series, levels = labels))

obs_row <- dplyr::filter(summ_tab, series == "baseline")

summ_vs_obs <- summ_tab %>%
  dplyr::mutate(
    excess_cumulative = cumulative_infections - obs_row$cumulative_infections,
    peak_diff         = peak_infections - obs_row$peak_infections,
    peak_advance_days = as.integer(peak_date - obs_row$peak_date) 
  )

# save CSV 
write.csv(summ_tab, file = csv_out, row.names = FALSE)

# distinct colours, mapped to percentages
#pal <- setNames(RColorBrewer::brewer.pal(n = max(8, length(labels)), "Set1")[seq_along(labels)], labels)

colors <- c(
  "baseline"        = "red",
  "75%"             = "#009E73",
  "50%"             = "#E69F00",
  "25%"             = "#D55E00",
  "no vaccination"  = "#999999",
  "no boosters"     = "#CC79A7",
  "early boosters"  = "#56B4E9"
)

p <- ggplot(df_all, aes(x = date)) +
  geom_line(aes(y = inc, color = series), linewidth = 1.2, alpha = 1) +
  #geom_point(data = no_vax_booster |> dplyr::filter(series == "baseline"), aes(y = inc, color = "red"), size = 1.5) +
  #scale_color_manual(name = "scenario", values = pal, breaks = labels) +
  scale_color_manual(name = "scenario", values = colors, breaks = names(colors)) + 
  labs(y = "Estimated daily infections", x = NULL, title = "Effect of vaccination (boosters/rate) on estimated daily infections") +
  scale_x_date(expand = c(0, 0), limits = c(date0, date_end),date_breaks = "2 week", date_labels = "%b %d") +
  geom_vline(aes(xintercept = as.Date("2021-12-15")), colour = "gold4", linetype = 4, linewidth = 1.0) +
  geom_vline(aes(xintercept = as.Date("2022-01-03")), colour = "gold4", linetype = 4, linewidth = 1.0) +
  geom_vline(aes(xintercept = as.Date("2022-01-24")), colour = "gold4", linetype = 4, linewidth = 1.0) +
  geom_vline(aes(xintercept = as.Date("2022-02-25")), colour = "gold4", linetype = 4, linewidth = 1.0) +
  geom_vline(aes(xintercept = as.Date("2022-03-17")), colour = "gold4", linetype = 4, linewidth = 1.0) +
  theme_clean() +
  theme(
    axis.text.x   = element_text(size = 16),
    axis.title.x  = element_blank(),
    axis.text.y   = element_text(size = 16),
    axis.title.y  = element_text(size = 16),
    plot.title    = element_text(size = 20, hjust = 0.5, face = "plain"),
    legend.position  = c(0.15, 0.75),  
    legend.direction = "vertical",
    legend.box.just  = "left",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 18, face = "plain"),
    legend.background= element_blank(),
    plot.background  = element_blank()
  )

print(p)

png("../figures/scenarios.png", width = 2800, height = 1500, res = 300, bg = "white", type = "cairo")
p
dev.off()
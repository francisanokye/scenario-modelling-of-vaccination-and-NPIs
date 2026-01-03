
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(ggthemes)
  library(shellpipes); library(macpan2); library(grid)
  library(scales); library(tibble)
})

set.seed(2025)
loadEnvironments()
options(dplyr.summarise.inform = FALSE)

## ------------------------- config ------------------------------------------
date0       <- as.Date("2021-12-17")
offset0     <- 150
sim_horizon <- 500

fig_dir  <- "../figures"
diag_dir <- file.path(fig_dir, "vax_only_diagnostics")

csv_out_summary <- "../output/vax_only_summary.csv"
csv_out_bars    <- "../output/vax_only_bar.csv"

## ------------------------- io helpers --------------------------------------
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

slugify <- function(x) {
  x <- tolower(x)
  x <- gsub("%", "pct", x, fixed = TRUE)
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_|_$", "", x)
}

save_png <- function(path, plot_obj, width, height, res = 200, bg = "white") {
  png(path, width = width, height = height, res = res, bg = bg, type = "cairo")
  print(plot_obj)
  dev.off()
  invisible(path)
}

ensure_dir(fig_dir); ensure_dir(diag_dir)
ensure_dir(dirname(csv_out_summary)); ensure_dir(dirname(csv_out_bars))

## ------------------------- inputs ------------------------------------------
calibrator      <- rdsRead("calibrate")
padded_beta0_hat <- rdsRead("policy_neutral")
spec_fixed      <- mp_optimized_spec(calibrator, "modified")

## ------------------------- model vars --------------------------------------
state_var <- c("S1","E1","V2","I1","A1","R1","E2","V3","I2","A2","R2","E3","I3","A3","R3")

## ---------------------- scaling vaccination --------------------------------
set_vax <- function(spec, scale = 1.0) {
  if (is.null(spec$default$double_vac_values_base))
    spec$default$double_vac_values_base <- spec$default$double_vac_values
  if (is.null(spec$default$booster_vac_values_base))
    spec$default$booster_vac_values_base <- spec$default$booster_vac_values
  
  spec$default$double_vac_values  <- scale * spec$default$double_vac_values_base
  spec$default$booster_vac_values <- scale * spec$default$booster_vac_values_base
  spec
}

## ------------------------- beta updater ------------------------------------
update_beta <- function(new_beta) {
  cp_spec <- spec_fixed
  stopifnot(length(new_beta) == cp_spec$default$time_steps)
  
  cp_spec |>
    mp_tmb_delete(phase = "before", at = 24L) |>
    mp_tmb_insert(default = list(timeseries_beta = new_beta))
}

## ------------------------- one-run trajectory ------------------------------
traj_all <- function(spec) {
  out <- c("inc", "total_inc", state_var)
  
  spec |>
    mp_simulator(sim_horizon, outputs = out) |>
    mp_trajectory() |>
    dplyr::mutate(date = date0 + (time - offset0))
}

extract_inc_total <- function(traj_df) {
  inc_df <- traj_df |>
    dplyr::filter(matrix == "inc") |>
    dplyr::transmute(date = as.Date(date), inc = as.integer(round(value))) |>
    dplyr::filter(date >= date0)
  
  tot_df <- traj_df |>
    dplyr::filter(matrix == "total_inc") |>
    dplyr::transmute(date = as.Date(date), total_inc = as.integer(round(value)))|>
    dplyr::filter(date >= date0)
  
  dplyr::left_join(inc_df, tot_df, by = "date")
}

plot_diagnostics <- function(traj_df, vax_label, als_label) {
  diag_df <- traj_df |>
    dplyr::filter(matrix %in% c(state_var, "total_inc"))
  
  ggplot(diag_df, aes(date, value)) +
    geom_line(linewidth = 0.45) +
    facet_wrap(~matrix, scales = "free", ncol = 4) +
    geom_rect(
      data = data.frame(xmin = date0 - offset0, xmax = date0),
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "grey30", alpha = 0.25
    ) +
    geom_vline(xintercept = date0 + (231 - offset0), color = "black",  linetype = 4, linewidth = 0.9) +
    geom_vline(xintercept = date0 + (300 - offset0), color = "purple", linetype = 4, linewidth = 1.1) +
    scale_x_date(expand = c(0, 0), date_breaks = "month", date_labels = "%b %Y") +
    labs(title = paste0("Diagnostics: ", als_label, " | vax = ", vax_label), x = NULL, y = NULL) +
    theme_clean() +
    theme(
      axis.text.x       = element_text(size = 12,  hjust = 1, angle = 45, margin = margin(t = 8)),
      axis.text.y       = element_text(size = 18),
      axis.title.x      = element_text(size = 18),
      axis.title.y      = element_text(size = 18),
      plot.title        = element_text(size = 24, hjust = 0.5, face = "plain"),
      strip.text        = element_text(size = 22, face = "plain", hjust = 0.5),
      strip.background  = element_blank(),
      legend.position   = "top",
      legend.direction  = "horizontal",
      legend.box        = "horizontal",
      legend.spacing.x  = unit(1.2, "cm"),
      legend.spacing.y  = unit(0.5, "cm"),
      legend.background = element_blank(),
      legend.text       = element_text(size = 22, face = "plain"),
      panel.spacing     = unit(1.5, "cm"),
      plot.background   = element_blank()
    )
}

## ------------------------- scenarios ---------------------------------------
als_levels <- c("No-ALS Only")
vax_scales <- c("0%"=0, "25%"=0.25, "50%"=0.50, "75%"=0.75, "100%"=1.0)

vax_cf_levels <- names(vax_scales)
panel_levels  <- paste0(vax_cf_levels, " Vax & No-ALS Only")

# diagnostics are expensive: save only these during development
diag_keep <- c("0%", "25%", "50%", "75%", "100%")    

beta_list <- list("No-ALS Only" = padded_beta0_hat)

## ------------------------- actual (one run) --------------------------------
actual_traj <- traj_all(spec_fixed)
actual <- extract_inc_total(actual_traj) |>
  dplyr::mutate(series = "Actual")

## ------------------------- counterfactuals ---------------------------------
simulate_vax_set_fast <- function(new_beta, als_label) {
  spec_cf0 <- update_beta(new_beta)
  
  purrr::imap_dfr(vax_scales, function(scale_val, vlab) {
    spec_cf <- set_vax(spec_cf0, scale_val)
    
    tr <- traj_all(spec_cf)
    it <- extract_inc_total(tr)
    
    if (vlab %in% diag_keep) {
      p_diag <- plot_diagnostics(tr, vax_label = vlab, als_label = als_label)
      fn <- file.path(diag_dir, paste0("diagnostic_vax_only_", slugify(vlab), ".png"))
      save_png(fn, p_diag, width = 7500, height = 5000, res = 300)
    }
    
    it |>
      dplyr::mutate(
        vax        = factor(vlab, levels = vax_cf_levels),
        als_scen    = factor(als_label, levels = als_levels),
        series      = "Counterfactual",
        panel_title = factor(paste0(vax, " Vax & ", als_scen), levels = panel_levels)
      ) |>
      dplyr::filter(date >= date0) |>
      dplyr::select(date, inc, total_inc, vax, als_scen, series, panel_title)
      
  })
}

cf_df_full <- purrr::imap_dfr(beta_list, simulate_vax_set_fast) |>
  dplyr::arrange(als_scen, vax, date)

## ------------------------- actual replicated to panels ----------------------
actual_panel_df <- actual |>
  tidyr::crossing(tibble(
    vax      = factor(vax_cf_levels, levels = vax_cf_levels),
    als_scen = factor(als_levels, levels = als_levels)
  )) |>
  dplyr::mutate(
    series      = "Actual",
    panel_title = factor(paste0(vax, " Vax & ", als_scen), levels = panel_levels)
  )

## ------------------------- summary table -----------------------------------
summ_df <- dplyr::bind_rows(cf_df_full, actual_panel_df)

summ_tab <- summ_df |>
  dplyr::mutate(title = paste(vax, als_scen, sep = " - ")) |>
  dplyr::group_by(title, series) |>
  dplyr::summarise(
    cumulative_infections = round(sum(inc, na.rm = TRUE)),
    peak_infections       = round(max(inc, na.rm = TRUE)),
    peak_date             = date[which.max(inc)],
    .groups = "drop"
  )

print(summ_tab)
write.csv(summ_tab, file = csv_out_summary, row.names = FALSE)

## ------------------------- plots -------------------------------------------
# multi-panel daily inc
p <- ggplot() +
  geom_line(data = cf_df_full, aes(date, inc, color = "Counterfactual"), linewidth = 1.05) +
  geom_point(data = actual_panel_df, aes(date, inc, color = "Actual"), size = 0.9) +
  facet_wrap(~ panel_title, nrow = 3, ncol = 2, scales = "fixed") +
  geom_vline(xintercept = as.Date("2022-03-14"), color = "black", linetype = 4, linewidth = 1.0) +
  geom_rect(inherit.aes = FALSE,
            xmin = as.Date("2022-05-22"), xmax = as.Date("2022-11-30"),
            ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.3) +
  geom_vline(xintercept = as.Date("2022-05-22"), color = "purple", linetype = 4, linewidth = 1.2) +
  scale_x_date(expand = c(0, 0), date_breaks = "2 weeks", date_labels = "%b %d") +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(values = c("Actual"="red", "Counterfactual"="#1b9e77"), name = NULL) +
  labs(title = "Estimated daily infections under vaccination-only scenario (No-ALS)",
       x = NULL, y = "Estimated infections") +
  theme_clean() +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE))+
  theme(
    axis.text.x       = element_text(size = 16, angle = 45, hjust = 1),
    axis.text.y       = element_text(size = 18),
    axis.title.y      = element_text(size = 18),
    plot.title        = element_text(size = 24, hjust = 0.5, face = "plain"),
    strip.text        = element_text(size = 16, face = "plain", hjust = 0.5),
    legend.position   = "top",
    legend.direction  = "horizontal",
    legend.box        = "horizontal",
    legend.spacing.x  = unit(1.2, "cm"),
    legend.spacing.y  = unit(0.5, "cm"),
    legend.background = element_blank(),
    legend.text       = element_text(size = 14, face = "plain"),
    panel.spacing     = unit(1.5, "cm"),
    plot.background   = element_blank()
  )
save_png(file.path(fig_dir, "vax_only_scenarios.png"), p, width = 4200, height = 3200, res = 300)

# single-panel overlay daily inc (vax colors + actual)
pal <- c(
  "Actual"="red",
  "0%"="gray49",
  "25%"="purple",
  "50%"="#0072B2",
  "75%"="#009E73",
  "100%"="#E69F00"
)

pp <- ggplot() +
  geom_line(data = cf_df_full, aes(date, inc, color = as.character(vax), group = vax), linewidth = 1.2) +
  geom_point(data = actual, aes(date, inc, color = "Actual"), size = 1.2) +
  geom_vline(xintercept = as.Date("2022-03-14"), color="black", linetype=4, linewidth=1.0) +
  geom_rect(inherit.aes=FALSE,
            xmin=as.Date("2022-05-22"), xmax=max(cf_df_full$date, na.rm=TRUE),
            ymin=-Inf, ymax=Inf, fill="grey80", alpha=0.3) +
  geom_vline(xintercept=as.Date("2022-05-22"), color="purple", linetype=4, linewidth=1.2) +
  scale_x_date(expand=c(0,0), date_breaks="2 weeks", date_labels="%b %d %Y") +
  scale_y_continuous(labels=scales::comma) +
  scale_color_manual(values = pal, breaks = c("Actual", vax_cf_levels), name = NULL) +
  labs(title="Estimated daily infections under vaccination-only scenario (No-ALS)",
       x=NULL, y="Estimated infections") +
  theme_clean() +
  theme_clean() +
  theme(
    axis.text.x       = element_text(size = 18, hjust = 1, angle = 45),
    axis.text.y       = element_text(size = 18),
    axis.title.y      = element_text(size = 18),
    plot.title        = element_text(size = 24, hjust = 0.5, face = "plain"),
    legend.position   = c(0.75, 0.65),
    legend.direction  = "vertical",
    legend.box.just   = "center",
    legend.text       = element_text(size = 18, face = "plain"),
    legend.background = element_blank(),
    plot.background   = element_blank()
  )

save_png(file.path(fig_dir, "vax_only_scen1.png"), pp, width = 4800, height = 1600, res = 220)

# cumulative total_inc (single mapping key so colors never “collapse”)
total_inc_plot_df <- dplyr::bind_rows(cf_df_full, actual_panel_df) |>
  dplyr::mutate(
    col_key = if_else(series == "Actual", "Actual", as.character(vax)),
    col_key = factor(col_key, levels = c("Actual", vax_cf_levels))
  )

ppp <- ggplot(total_inc_plot_df |> dplyr::filter(date >= date0)) +
  geom_line(aes(date, total_inc, color = col_key), linewidth = 1.05) +
  facet_wrap(~ als_scen, nrow = 1, scales = "fixed") +
  geom_vline(xintercept = as.Date("2022-03-14"), color="black",  linetype=4, linewidth=1.0) +
  geom_vline(xintercept = as.Date("2022-05-22"), color="purple", linetype=4, linewidth=1.0) +
  scale_x_date(expand=c(0,0), date_breaks="2 weeks", date_labels="%b %d %Y") +
  scale_y_continuous(labels=scales::comma) +
  scale_color_manual(values = pal, breaks = c("Actual", vax_cf_levels), name = NULL) +
  labs(title="Cumulative infections under different vaccination (No-ALS)",
       x=NULL, y="Cumulative infections") +
  theme_clean() +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE))+
  theme(
    axis.text.x       = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y       = element_text(size = 18),
    axis.title.y      = element_text(size = 18),
    plot.title        = element_text(size = 24, hjust = 0.5, face = "plain"),
    strip.text        = element_text(size = 16, face = "plain", hjust = 0.5),
    legend.position   = "top",
    legend.direction  = "horizontal",
    legend.box        = "horizontal",
    legend.spacing.x  = unit(1.2, "cm"),
    legend.spacing.y  = unit(0.5, "cm"),
    legend.background = element_blank(),
    legend.text       = element_text(size = 18, face = "plain"),
    panel.spacing     = unit(1.5, "cm"),
    plot.background   = element_blank()
  )
save_png(file.path(fig_dir, "vax_only_cum.png"), ppp, width = 5200, height = 2800, res = 300)

## ------------------------- stacked bars (inc totals) ------------------------
cutoff_date <- as.Date("2022-03-14")

cf_bar <- cf_df_full |>
  dplyr::filter(date >= date0) |>
  dplyr::mutate(
    vaccination = as.character(vax),
    period      = if_else(date <= cutoff_date, "Up to March 14, 2022", "March 15 - Nov 30, 2022"),
    series      = "Counterfactual",
    als_scen    = "No-ALS Only"
  ) |>
  dplyr::group_by(als_scen, vaccination, period, series) |>
  dplyr::summarise(total_infections = sum(inc, na.rm = TRUE), .groups = "drop")

actual_bar <- actual |>
  dplyr::filter(date >= date0) |>
  dplyr::mutate(
    vaccination = "Actual",
    period      = if_else(date <= cutoff_date, "Up to March 14, 2022", "March 15 - Nov 30, 2022"),
    series      = "Actual",
    als_scen    = "No-ALS Only"
  ) |>
  dplyr::group_by(als_scen, vaccination, period, series) |>
  dplyr::summarise(total_infections = sum(inc, na.rm = TRUE), .groups = "drop")

bar_df <- dplyr::bind_rows(cf_bar, actual_bar) |>
  dplyr::mutate(
    als_scen    = factor(als_scen, levels = "No-ALS Only"),
    vaccination = factor(vaccination, levels = c(vax_cf_levels, "Actual")),
    period      = factor(period, levels = c("March 15 - Nov 30, 2022", "Up to March 14, 2022")),
    fill_group  = interaction(series, period, drop = TRUE)
  )

write.csv(bar_df, csv_out_bars, row.names = FALSE)

fill_vals <- c(
  "Counterfactual.Up to March 14, 2022"    = "#1b9e77",
  "Counterfactual.March 15 - Nov 30, 2022" = "gold",
  "Actual.Up to March 14, 2022"            = "navy",
  "Actual.March 15 - Nov 30, 2022"         = "#f55f74"
)

fill_labs <- c(
  "Counterfactual.Up to March 14, 2022"    = "Counterfactual (≤ Mar 14, 2022)",
  "Counterfactual.March 15 - Nov 30, 2022" = "Counterfactual (> Mar 14, 2022)",
  "Actual.Up to March 14, 2022"            = "Actual (≤ Mar 14, 2022)",
  "Actual.March 15 - Nov 30, 2022"         = "Actual (> Mar 14, 2022)"
)

p_bar <- ggplot(bar_df) +
  geom_col(aes(vaccination, total_infections, fill = fill_group), width = 0.75) +
  scale_fill_manual(values = fill_vals, breaks = names(fill_labs), labels = fill_labs, name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Cumulative infections under vaccination-only scenario without ALS",
       x = NULL, y = "Total infections (people)") +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_clean() +
  theme(
    axis.text.x       = element_text(size = 18, angle = 0, hjust = 0.5),
    axis.text.y       = element_text(size = 18),
    axis.title.y      = element_text(size = 18),
    plot.title        = element_text(size = 24, hjust = 0.5, face = "plain"),
    strip.text        = element_blank(),
    legend.position   = "top",
    legend.direction  = "horizontal",
    legend.box        = "horizontal",
    legend.spacing.x  = unit(1.2, "cm"),
    legend.spacing.y  = unit(0.5, "cm"),
    legend.background = element_blank(),
    legend.text       = element_text(size = 18 , face = "plain"),
    panel.spacing     = unit(1.5, "cm"),
    plot.background   = element_blank()
  )
save_png(file.path(fig_dir, "vax_only_bar.png"), p_bar, width = 4800, height = 1600, res = 300)

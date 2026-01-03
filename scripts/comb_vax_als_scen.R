
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(ggthemes)
  library(shellpipes); library(macpan2)
  library(scales); library(tibble); library(zoo)
  library(grid)  # for unit()
})

set.seed(2025)
options(macpan2_log_dir = ".")
loadEnvironments()
options(dplyr.summarise.inform = FALSE)
options(scipen = 9999999999)

## --------------------------- config ---------------------------------------
date0       <- as.Date("2021-12-17")
date_end    <- as.Date("2022-05-22")
start_date  <- date0
offset0     <- 150
sim_horizon <- 500

fig_dir  <- "../figures"
diag_dir <- file.path(fig_dir, "als_vax_diagnostics")
csv_out  <- "../output/als_vax_summary.csv"
csv_bar <- "../output/als_vax_bar.csv"

## --------------------------- io helpers -----------------------------------
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

save_png <- function(path, plot_obj, width, height, res = 220, bg = "white") {
  png(path, width = width, height = height, res = res, bg = bg, type = "cairo")
  print(plot_obj)
  dev.off()
  invisible(path)
}

ensure_dir(fig_dir); ensure_dir(diag_dir)
ensure_dir(dirname(csv_out)); ensure_dir(dirname(csv_bar))

## --------------------------- inputs ---------------------------------------
calibrator   <- rdsRead("calibrate")
fitserodata  <- rdsRead("fitsero.rds")
time_steps   <- max(fitserodata$time)

padded_beta0_hat    <- rdsRead("policy_neutral")
padded_base_beta    <- rdsRead("base_beta")
padded_beta_cf_ALS2 <- rdsRead("beta_als2")
padded_beta_cf_ALS4 <- rdsRead("beta_als4")

spec_fixed <- mp_optimized_spec(calibrator, "modified")

## cache pristine vaccination series once
if (is.null(spec_fixed$default$double_vac_values_base))
  spec_fixed$default$double_vac_values_base  <- spec_fixed$default$double_vac_values
if (is.null(spec_fixed$default$booster_vac_values_base))
  spec_fixed$default$booster_vac_values_base <- spec_fixed$default$booster_vac_values
if (is.null(spec_fixed$integers$booster_vac_changepoints_base))
  spec_fixed$integers$booster_vac_changepoints_base <- spec_fixed$integers$booster_vac_changepoints

## ------------------------- variables --------------------------------------
state_var <- c("S1","E1","V2","I1","A1","R1","E2","V3","I2","A2","R2","E3","I3","A3","R3")

vax_scales <- c("0%"=0, "25%"=0.25, "50%"=0.50, "75%"=0.75, "100%"=1.0)
vax_order  <- names(vax_scales)

als_levels <- c("No-ALS Only", "ALS-2 Only", "Implemented ALS", "ALS-4 Only")
als_fac    <- factor(als_levels, levels = als_levels)

beta_list <- list(
  "No-ALS Only"     = padded_beta0_hat,
  "ALS-2 Only"      = padded_beta_cf_ALS2,
  "Implemented ALS" = padded_base_beta,
  "ALS-4 Only"      = padded_beta_cf_ALS4
)

## ---------------------- scaling vaccination --------------------------------
set_vax <- function(spec, scale = 1.0) {
  spec$default$double_vac_values  <- scale * spec$default$double_vac_values_base
  spec$default$booster_vac_values <- scale * spec$default$booster_vac_values_base
  spec
}

## ------------------------- beta updater (do once per ALS) -------------------
update_beta <- function(spec, new_beta) {
  stopifnot(length(new_beta) == spec$default$time_steps)
  spec |>
    mp_tmb_delete(phase = "before", at = 24L) |>
    mp_tmb_delete(phase = "before", at = 24L) |>
    mp_tmb_insert(default = list(timeseries_beta = new_beta))
}

## ------------------------ one-run trajectory (FAST) -------------------------
traj_all <- function(spec) {
  out <- c("inc","total_inc", state_var)
  spec |>
    mp_simulator(sim_horizon, outputs = out) |>
    mp_trajectory() |>
    #dplyr::filter(time >= offset0) |>
    dplyr::mutate(date = date0 + (time - offset0))
}

extract_inc_total <- function(tr) {
  inc_df <- tr |> dplyr::filter(matrix=="inc") |>
    transmute(date = as.Date(date), inc = as.integer(round(value)))
  tot_df <- tr |> dplyr::filter(matrix=="total_inc") |>
    transmute(date = as.Date(date), total_inc = as.integer(round(value)))
  left_join(inc_df, tot_df, by="date")
}

plot_diagnostics <- function(tr, vax_label, als_label, scenario_label) {
  df <- tr |> dplyr::filter(matrix %in% c(state_var,"total_inc"))
  
  ggplot(df, aes(date, value)) +
    geom_line(linewidth = 0.4) +
    facet_wrap(~matrix, scales="free", ncol=4) +
    geom_rect(
      data = data.frame(xmin = date0 - offset0, xmax = date0),
      aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf),
      inherit.aes = FALSE,
      fill="grey30", alpha=0.25
    ) +
    labs(
      title = paste0("Diagnostics: ", scenario_label, " | ", als_label, " | vax=", vax_label),
      x=NULL, y=NULL
    ) +
    geom_vline(
      xintercept = date0 + (231 - offset0),
      color = "black",
      linetype = 4,
      linewidth = 1.0
    ) +
    geom_vline(
      xintercept = date0 + (300 - offset0),
      colour = "purple",
      linetype = 4,
      linewidth = 1.5
    ) +
    aes(date, value) +
    geom_line() +
    labs(
      title = paste0("Diagnostics: ", scenario_label, " | ", als_label, " | vax=", vax_label),
      x=NULL, y=NULL
    ) +
    scale_x_date(expand = c(0, 0), date_breaks = "month", date_labels = "%b %d %Y") +
    facet_wrap(~matrix, scales = "free") +
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

## ---------------------- run all scenarios (FAST) ----------------------------
diag_keep <- c("0%","25%","50%", "75%","100%")   # or character(0) to skip all

run_all <- function(spec_base, beta_list, scenario_label="early_booster_28d") {
  
  purrr::imap_dfr(beta_list, function(beta_vec, als_label) {
    
    spec_beta <- update_beta(spec_base, beta_vec)  # once per ALS beta
    
    purrr::imap_dfr(vax_scales, function(scale_val, vlab) {
      
      spec_cf <- set_vax(spec_beta, scale_val)
      tr <- traj_all(spec_cf)
      
      if (vlab %in% diag_keep) {
        pdiag <- plot_diagnostics(tr, vax_label = vlab, als_label = als_label, scenario_label = scenario_label)
        fn <- file.path(diag_dir, paste0("diagnostic_", slugify(als_label), "_vax_", slugify(vlab), ".png"))
        save_png(fn, pdiag, width=3600, height=2400, res=180)
      }
      
      extract_inc_total(tr) |>
        mutate(
          vax      = factor(vlab, levels = vax_order),
          als_scen = factor(als_label, levels = als_levels),
          series   = "counterfactuals"
        ) |>
        select(date, inc, total_inc, vax, als_scen, series)
    })
  }) |>
    arrange(als_scen, vax, date)
}

df_all <- run_all(spec_fixed, beta_list, scenario_label="early_booster_28d")

## ---------------------- observed (actual) series (ONE RUN) -----------------
# observed does NOT depend on ALS here, so simulate once and replicate
tr_obs <- traj_all(spec_fixed)  # or spec_fixed if you want “non-early” observed
obs_one <- extract_inc_total(tr_obs) |>
  mutate(series = "observed")

actual <- obs_one |>
  crossing(tibble(als_scen = factor(als_levels, levels=als_levels))) |>
  mutate(vax = factor("100%", levels = vax_order)) |>
  select(date, inc, total_inc, vax, als_scen, series) |>
  arrange(als_scen, date)

## ---------------------- plotting frame ------------------------------------
panel_levels <- as.vector(outer(vax_order, als_levels, FUN = \(v,a) paste0(v," Vax & ",a)))

plot_df <- bind_rows(df_all, actual) |>
  mutate(
    vax = factor(vax, levels = vax_order),
    als_scen = factor(als_scen, levels = als_levels),
    panel_title = factor(paste0(vax," Vax & ",als_scen), levels = panel_levels)
  )

inc_plot_df <- plot_df |> select(-total_inc)
total_inc_plot_df <- plot_df |> select(-inc)

## --------------------------- summary table --------------------------------
summ_tab <- inc_plot_df |>
  filter(date <= date_end) |>
  mutate(title = paste(vax, als_scen, sep=" - ")) |>
  group_by(title) |>
  summarise(
    cumulative_infections = round(sum(inc, na.rm=TRUE)),
    peak_infections       = round(max(inc, na.rm=TRUE)),
    peak_date             = date[which.max(inc)],
    .groups="drop"
  )

print(summ_tab)
write.csv(summ_tab, csv_out, row.names = FALSE)

## ---------------------- figure: 20-panel grid ------------------------------
# (still big, but now plotting is the main cost, not sims)
obs_only <- inc_plot_df |> filter(series=="observed") |>
  group_by(date, als_scen) |> summarise(inc_obs = first(inc), .groups="drop")

cf_only <- inc_plot_df |> filter(series=="counterfactuals")

ann_pos <- cf_only |> group_by(panel_title) |>
  summarise(
    xpos = min(date, na.rm=TRUE) + (max(date, na.rm=TRUE)-min(date, na.rm=TRUE))*0.35,
    ypos = max(inc, na.rm=TRUE) * 0.95,
    .groups="drop"
  )

p_grid <- ggplot() +
  geom_line(data = cf_only, aes(date, inc, color="Counterfactual"), linewidth=1.1) +
  geom_point(data = obs_only, aes(date, inc_obs, color="Actual"), size=1.0) +
  facet_wrap(~panel_title, nrow=4, ncol=5, scales="fixed") +
  geom_vline(xintercept = as.Date("2022-03-14"), color="black", linetype = 4, linewidth=1.0) +
  geom_vline(xintercept = as.Date("2022-05-22"), color = "purple", linetype = 4, linewidth = 1.0) +
  scale_x_date(expand=c(0,0), date_breaks="2 weeks", date_labels="%b %d") +
  scale_y_continuous(labels=scales::comma) +
  scale_color_manual(values=c("Actual"="red","Counterfactual"="#2ca02c"), name=NULL) +
  labs(
    title="Combined vaccination and ALS: Counterfactual infections vs actual",
    x=NULL, y="Estimated infections"
  ) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE))+
  theme_clean() +
  theme(
    axis.text.x       = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y       = element_text(size = 18),
    axis.title.y      = element_text(size = 18),
    plot.title        = element_text(size = 24, hjust = 0.5, face = "plain"),
    strip.text        = element_text(size = 22, face = "plain", hjust = 0.5),
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

save_png(file.path(fig_dir, "als_vax_scenario.png"), p_grid, width=8000, height=5600, res=300)

## ---------------------- Figure: 2x2 ALS × vax curves (inc) -----------------

# counterfactual curves
cf_df <- inc_plot_df |>
  dplyr::filter(
    series == "counterfactuals",
    date   >= date0,
    vax    %in% vax_order,
    als_scen %in% als_levels
  )

# observed points: one observed series per ALS panel (same observed replicated)
obs_df <- inc_plot_df |>
  dplyr::filter(series == "observed", date >= date0) |>
  dplyr::arrange(date) |>
  dplyr::group_by(date, als_scen) |>
  dplyr::summarise(inc = first(inc), .groups = "drop") |>
  dplyr::mutate(series = "Actual")

# global top for label placement
global_top <- max(c(cf_df$inc, obs_df$inc), na.rm = TRUE)

label_pos <- cf_df |>
  dplyr::distinct(als_scen) |>
  dplyr::arrange(als_scen) |>
  dplyr::mutate(
    x_lab       = min(cf_df$date, na.rm = TRUE) + 3,
    y_lab       = global_top * 0.99,
    panel_label = LETTERS[seq_len(dplyr::n())]
  )

pp <- ggplot() +
  geom_line(data = cf_df,aes(date, inc, color = vax, group = vax),linewidth = 1.2) +
  geom_point(data = obs_df,aes(date, inc, color = series),size = 1.2) +
  facet_wrap(~ als_scen, nrow = 2, ncol = 2, scales = "fixed") +
  geom_text(data = label_pos,aes(x = x_lab, y = y_lab, label = panel_label),hjust = 0, vjust = 1, size = 8, fontface = "bold", inherit.aes = FALSE) +
  geom_vline(xintercept = as.Date("2022-03-14"), color = "black",  linetype = 4, linewidth = 1.0) +
  geom_vline(xintercept = as.Date("2022-05-22"), color = "purple", linetype = 4, linewidth = 1.0) +
  scale_x_date(expand = c(0, 0), date_breaks = "2 weeks", date_labels = "%b %d %Y") +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(
    values = c(
      "Actual" = "red",
      "0%"     = "gray",
      "25%"    = "#9467bf",
      "50%"    = "#2ca02c",
      "75%"    = "#1f77b4",
      "100%"   = "orange"
    ),
    breaks = c("Actual", vax_order),
    name   = NULL
  ) +
  labs(
    title = "Infections by both ALS and vaccination rate scenario",
    x     = NULL,
    y     = "Estimated infections"
  ) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE))+
  theme_clean() +
  theme(
    axis.text.x       = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y       = element_text(size = 18),
    axis.title.y      = element_text(size = 20),
    plot.title        = element_text(size = 24, hjust = 0.5, face = "plain"),
    strip.text        = element_text(size = 22, face = "plain", hjust = 0.5),
    legend.position   = "top",
    legend.direction  = "horizontal",
    legend.box        = "horizontal",
    legend.spacing.x  = unit(1.2, "cm"),
    legend.spacing.y  = unit(0.5, "cm"),
    legend.background = element_blank(),
    legend.text       = element_text(size = 20),
    panel.spacing     = unit(2.5, "cm"),
    plot.background   = element_blank(),
    strip.background  = element_blank()
  )

save_png(file.path(fig_dir, "als_vax_scen1.png"), pp, width = 5500, height = 3500, res = 300)

## ---------------------- Figure: 2x2 ALS × vax curves (total_inc) -----------

total_cf_df <- total_inc_plot_df |>
  dplyr::filter(
    series == "counterfactuals",
    date   >= date0,
    vax    %in% vax_order,
    als_scen %in% als_levels
  )

total_obs_df <- total_inc_plot_df |>
  dplyr::filter(series == "observed", date >= date0) |>
  dplyr::arrange(date) |>
  dplyr::group_by(date, als_scen) |>
  dplyr::summarise(total_inc = first(total_inc), .groups = "drop") |>
  dplyr::mutate(series = "Actual")

# (reuse) label_pos is fine; y_lab should reflect cumulative scale
global_top_tot <- max(c(total_cf_df$total_inc, total_obs_df$total_inc), na.rm = TRUE)

label_pos_tot <- label_pos |>
  dplyr::mutate(y_lab = global_top_tot * 0.99)

ppp <- ggplot() +
  geom_line(data = total_cf_df,aes(date, total_inc, color = vax, group = vax),linewidth = 1.2) +
  geom_line(data = total_obs_df,aes(date, total_inc, color = series),linewidth = 1.2) +
  facet_wrap(~ als_scen, nrow = 2, ncol = 2, scales = "fixed") +
  geom_text(data = label_pos_tot,aes(x = x_lab, y = y_lab, label = panel_label),hjust = 0, vjust = 1, size = 8, fontface = "bold", inherit.aes = FALSE) +
  geom_vline(xintercept = as.Date("2022-03-14"), color = "black",  linetype = 4, linewidth = 1.0) +
  geom_vline(xintercept = as.Date("2022-05-22"), color = "purple", linetype = 4, linewidth = 1.0) +
  scale_x_date(expand = c(0, 0), date_breaks = "2 weeks", date_labels = "%b %d %Y") +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(
    values = c(
      "Actual" = "red",
      "0%"     = "gray",
      "25%"    = "#9467bf",
      "50%"    = "#2ca02c",
      "75%"    = "#1f77b4",
      "100%"   = "orange"
    ),
    breaks = c("Actual", vax_order),
    name   = NULL
  ) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE))+
  labs(
    title = "Cumulative infections under both ALS scenario and vaccination rate",
    x     = NULL,
    y     = "Cumulative infections"
  ) +
  theme_clean() +
  theme(
    axis.text.x       = element_text(size = 18, angle = 45, hjust = 1),
    axis.text.y       = element_text(size = 18),
    axis.title.y      = element_text(size = 20),
    plot.title        = element_text(size = 24, hjust = 0.5, face = "plain"),
    strip.text        = element_text(size = 22, face = "plain", hjust = 0.5),
    legend.position   = "top",
    legend.direction  = "horizontal",
    legend.box        = "horizontal",
    legend.spacing.x  = unit(1.2, "cm"),
    legend.spacing.y  = unit(0.5, "cm"),
    legend.background = element_blank(),
    legend.text       = element_text(size = 20),
    panel.spacing     = unit(2.5, "cm"),
    plot.background   = element_blank(),
    strip.background  = element_blank()
  )

save_png(file.path(fig_dir, "als_vax_cum.png"), ppp, width = 5500, height = 3500, res = 300)


## ------------------ stacked bar: totals split at Mar 14 --------------------
cut_date <- as.Date("2022-03-14")

sum_df <- total_inc_plot_df |>
  filter(series=="counterfactuals") |>
  group_by(vax, als_scen) |>
  arrange(date) |>
  summarise(
    total_cut  = {
      tmp <- total_inc[date <= cut_date]
      if (length(tmp)==0) 0 else tmp[length(tmp)]
    },
    total_last = total_inc[length(total_inc)],
    .groups="drop"
  ) |>
  mutate(
    total_post = pmax(total_last - total_cut, 0)
  ) |>
  select(vax, als_scen, total_cut, total_post) |>
  pivot_longer(c(total_cut, total_post), names_to="period", values_to="total_infections") |>
  mutate(
    period = recode(period,
                    total_cut  = "Up to March 14, 2022",
                    total_post = "March 15 - Dec 02, 2022"
    ),
    period = factor(period, levels=c("March 15 - Dec 02, 2022", "Up to March 14, 2022")),
    vax = factor(vax, levels=vax_order),
    als_scen = factor(als_scen, levels=als_levels),
    type = "als_vax"
  )

write.csv(sum_df, csv_bar, row.names = FALSE)

p_bar <- ggplot(sum_df, aes(vax, total_infections, fill=period)) +
  geom_col(width=0.75) +
  facet_wrap(~als_scen, nrow=1) +
  scale_y_continuous(labels=scales::comma) +
  coord_cartesian(clip = "off") +
  guides(fill = guide_legend(order = 1, reverse = FALSE)) +
  scale_fill_manual(
    values = c(
      "Up to March 14, 2022"    = "blue",
      "March 15 - Dec 02, 2022" = "goldenrod"
    ),
    breaks = c("Up to March 14, 2022", "March 15 - Dec 02, 2022"),
    name = NULL
  ) +
  labs(
    title="Total infections under both ALS and vaccination scenarios",
    x="Vaccination rate", y="Total infections (people)", fill=NULL
  ) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE))+
  theme_clean() +
  theme(
    axis.text.x       = element_text(size = 20, angle = 0, hjust = 0.5),
    axis.text.y       = element_text(size = 20),
    axis.title.x      = element_text(size = 20),
    axis.title.y      = element_text(size = 20),
    plot.title        = element_text(size = 24, hjust = 0.5, face = "plain"),
    strip.text        = element_text(size = 22, face = "plain"),
    strip.background = element_rect(fill = "grey90", color = NA),
    legend.position   = "top",
    legend.direction  = "horizontal",
    legend.box.just   = "left",
    legend.background = element_blank(),
    legend.text       = element_text(size = 20),
    plot.background   = element_blank(),
    plot.margin = margin(20, 20, 20, 40),
    panel.spacing.x   = unit(1.2, "cm"),
    panel.grid.major.x = element_blank()
  )

save_png(file.path(fig_dir, "als_vax_bar.png"), p_bar, width = 6500, height = 3200, res = 300)

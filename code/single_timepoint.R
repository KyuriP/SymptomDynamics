# =========================================================
# Figure E1 : Single-Time-Point Network Comparison
# =========================================================
# Purpose:
#   (1) Extracts 10-step micro-windows centered at t = 500, 1500, 4000
#       for low / medium / high resilience simulations.
#   (2) Estimates GGMs for each phase and resilience level.
#   (3) Plots a 3×3 grid of networks with consistent scaling.
#   (4) Adds a common legend.
# =========================================================

# --- Load dependencies and data ---------------------------------------------
source("code/libraries.R")
source("code/utils.R")

low  <- readRDS("data/aggregated_low.rds")
med  <- readRDS("data/aggregated.rds")
high <- readRDS("data/aggregated_res.rds")

# Representative time points for each phase
t_before <- 500
t_during <- 1500
t_after  <- 4000

# --- Helper functions -------------------------------------------------------

## Estimate GGM from a subset of data
estimate_single_GGM <- function(df) {
  estimateNetwork(df, default = "EBICglasso")
}

## Extract a 10-step micro-window around a target time
get_phase_data <- function(data, t_target, dt = 0.1, n_steps = 10) {
  half_window <- (n_steps / 2) * dt
  subset <- data %>%
    dplyr::filter(t >= (t_target - half_window) & t < (t_target + half_window)) %>%
    dplyr::select(S_anh:S_sui)
  
  message("t = ", t_target, " | Range: [", t_target - half_window, ", ", 
          t_target + half_window, ") | n = ", nrow(subset))
  
  if (nrow(subset) == 0) stop("No data found near t = ", t_target)
  subset
}

# --- Extract single-time-point data ----------------------------------------
low_before  <- get_phase_data(low,  t_before)
low_during  <- get_phase_data(low,  t_during)
low_after   <- get_phase_data(low,  t_after)

med_before  <- get_phase_data(med,  t_before)
med_during  <- get_phase_data(med,  t_during)
med_after   <- get_phase_data(med,  t_after)

high_before <- get_phase_data(high, t_before)
high_during <- get_phase_data(high, t_during)
high_after  <- get_phase_data(high, t_after)

# --- Estimate single-time-point networks -----------------------------------
low_before_net  <- estimate_single_GGM(low_before)
low_during_net  <- estimate_single_GGM(low_during)
low_after_net   <- estimate_single_GGM(low_after)

med_before_net  <- estimate_single_GGM(med_before)
med_during_net  <- estimate_single_GGM(med_during)
med_after_net   <- estimate_single_GGM(med_after)

high_before_net <- estimate_single_GGM(high_before)
high_during_net <- estimate_single_GGM(high_during)
high_after_net  <- estimate_single_GGM(high_after)

# --- Plot setup -------------------------------------------------------------

nodenames <- c("anh", "sad", "slp", "ene", "app", "glt", "con", "mot", "sui")
net_layout <- readRDS("data/net_layout.rds")


cols <- c("Low resilience"    = "#CC3300",   # red
          "Medium resilience" = "#3366CC",   # blue
          "High resilience"   = "#339966")   # green

# Determine consistent edge scaling
max_value_all <- max(
  abs(low_before_net$graph), abs(low_during_net$graph), abs(low_after_net$graph),
  abs(med_before_net$graph), abs(med_during_net$graph), abs(med_after_net$graph),
  abs(high_before_net$graph), abs(high_during_net$graph), abs(high_after_net$graph)
)

# --- Plot function ----------------------------------------------------------
plot_net <- function(net, color) {
  qgraph(net$graph,
         layout = net_layout,
         labels = nodenames,
         edge.color = color,
         border.color = color,
         color = "white",
         label.cex = 1.3,          # larger labels
         node.width = 3.5,
         edge.width = 1.2,
         maximum = max_value_all,
         fade = TRUE,
         gray = FALSE,
         vsize = 3.5,                # slightly larger node circles
         plot.margin = 0.15,
         asp = 1)
}

# --- Create figure ----------------------------------------------------------
png("figures/single_time_networks.png", width = 600, height = 1600, res = 300)

par(mfrow = c(3, 3),
    mar = c(1, 1, 1, 1),  # internal margins inside each panel
    oma = c(8, 5, 5, 5))  # outer margins around the grid

# ------------------ Plot all 9 networks ------------------

# Row 1 – Low resilience
plot_net(low_before_net,  cols["Low resilience"])
plot_net(low_during_net,  cols["Low resilience"])
plot_net(low_after_net,   cols["Low resilience"])

# Row 2 – Medium resilience
plot_net(med_before_net,  cols["Medium resilience"])
plot_net(med_during_net,  cols["Medium resilience"])
plot_net(med_after_net,   cols["Medium resilience"])

# Row 3 – High resilience
plot_net(high_before_net, cols["High resilience"])
plot_net(high_during_net, cols["High resilience"])
plot_net(high_after_net,  cols["High resilience"])

# --- Common legend ----------------------------------------------------------
par(fig = c(0, 1, 0, 0.12), new = TRUE, mar = c(0, 0, 0, 0))  # ⬅ more vertical space
plot.new()
legend("bottom",
       legend = names(cols),
       fill = cols,
       border = "black",
       horiz = TRUE,
       bty = "n",
       cex = 5,       # ⬅ larger font size
       pt.cex = 5,      # ⬅ larger color boxes
       text.col = "black",
       xpd = TRUE,
       y.intersp = 5)  # ⬅ increase line spacing


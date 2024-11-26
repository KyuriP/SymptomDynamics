## =========================================================
## Main Simulation 
##
## Run the simulation 500 times for each of the resilience scenarios (low, baseline, high).
## Setting the scenario to "baseline" generates:
## - Figure 1: Symptom dynamics during a shock event for individuals in the baseline resilience scenario.
## - Figure 4: Example trajectories and median trajectories for each symptom from individuals in the baseline scenario.
## These figures are referenced in the paper.
## =========================================================

## install packages
source("code/libraries.R")

## source necessary functions
source("code/utils.R")
source("code/euler_stochastic2.R")
source("code/mod_specification.R")

## set the seed
set.seed(123)

## define model specifics: choose the scenario and initial value for symptoms
choice <- "base"
choice <- "high"
choice <- "low"

mod <- mod_spec(scenario = choice, init_val = 0.01)


## define "f"
f <- function(x) x^2

## original params
parms1 <- c(mod$Beta_bistable, mod$delta)
## params given shock (beta: increases)  
parms2 <- c(mod$Beta_sick, mod$delta)

## define dt and the number of time steps:
deltaT <- .1 # dt 
timelength <- 4000 # length of simulation
n_steps <- as.integer(timelength / deltaT) # must be a number greater than 1

## specify the magnitude of noise and specifics of shock 
D_stoeq1 <- 0.01  # before shock
t_shock <- 1000 # time that shock begins
shock_duration <- 500 # shock duration time points

## results
sde_out <- euler_stochastic2(
  Amat = mod$A, 
  deterministic_rate = mod$dif_eq,
  stochastic_rate = mod$sto_eq,
  initial_condition = mod$initial_values,
  parameters1 = parms1,
  parameters2 = parms2, 
  deltaT = deltaT,
  timelength = timelength,
  D1 = D_stoeq1,
  shock = TRUE,
  t_shock = t_shock, 
  duration = shock_duration,
  seed = 123  # set seed
)

# shock period 
shock_period <- data.frame(time = 0:timelength) |>
  mutate(shock = ifelse(time >= t_shock & time <= t_shock + shock_duration, TRUE, FALSE))
  

# labeller
labelllername <- c("S_anh" = "anhedonia", "S_app" = "appetite", "S_con" = "concentration", "S_ene" = "energy", "S_glt" = "guilty", "S_mot" = "motor", "S_sad" = "sad", "S_slp" = "sleep", "S_sui" = "suicidal") 

eachsym <- sde_out |>
  tidyr::pivot_longer(!t, names_to = "symptoms") |>
  ggplot() +
  geom_area(data = shock_period, aes(x = time, y = shock*max(sde_out[,-1])
  ), inherit.aes = FALSE, fill = "darkgray", alpha = 0.2) +
  geom_line(aes(x = t, y = value, color = symptoms), alpha = 0.5, lwd= 0.2) +
  scale_color_discrete(name = "Symptoms", labels = c("anhedonia", "sad", "sleep", "energy", "appetite", "guilty", "concentration", "motor", "suicidal"))+
  #labs(x="time", y = "", title ="Trajectory of each symptom") +
  labs(x="time", y = "") +
  lims(y=c(0,1))+
  theme_classic(base_size = 15) +
  theme(legend.position="none") +
  facet_wrap(~symptoms,
  labeller = labeller(symptoms = labelllername) 
  )

# ggsave("eachsym_v2.pdf", plot = eachsym, width = 30, height =15, units = "cm", dpi = 300)


## aggregate symptom level
n_sims <- 500

# run simulation n_sims times
aggregated <- map(1:n_sims, ~ euler_stochastic2(
  Amat = mod$A, 
  deterministic_rate = mod$dif_eq,
  stochastic_rate = mod$sto_eq,
  initial_condition = mod$initial_values,
  parameters1 = parms1,
  parameters2 = parms2, 
  deltaT = deltaT,
  timelength = timelength,
  D1 = D_stoeq1,
  shock = TRUE,
  t_shock = t_shock, 
  duration = shock_duration) |>
  mutate(totalsymptom = rowSums(pick(S_anh:S_sui)))
) |> list_rbind(names_to = "sim")


# read all datasets
aggregated <- readRDS("data/aggregated.rds")
aggregated_high <- readRDS("data/aggregated_res.rds")
aggregated_low <- readRDS("data/aggregated_low.rds")

# aggregated2 <- readRDS("data/aggregated_base2.rds")
# aggregated_high2 <- readRDS("data/aggregated_high2.rds")
# aggregated_low2 <- readRDS("data/aggregated_low2.rds")


## Each symptom plot
each_summ <- aggregated |> 
  select(-c(totalsymptom, sim)) |>
  tidyr::pivot_longer(!t, names_to = "symptoms") |>
  reframe(inter_quantile(value), .by = c(t, symptoms)) |>
  pivot_wider(names_from = "quant", values_from = "val", names_glue = "q{quant}") 

eachsym_quant <- ggplot(data = each_summ) +
  geom_area(data = shock_period, aes(x = time, y = shock*max(sde_out[,-1])
  ), inherit.aes = FALSE, fill = "darkgray", alpha = 0.2) +
  geom_line(aes(x = t, y = q0.5, color = symptoms), alpha = 0.6, lwd = 0.2) +
  geom_ribbon(aes(x=t,ymin=q0.25,ymax=q0.75, fill = symptoms), alpha=0.3) +
  scale_color_discrete(name = "Symptoms", labels = c("anhedonia", "sad", "sleep", "energy", "appetite", "guilty", "concentration", "motor", "suicidal")) +
  labs(x = "time", y= "") +
  theme_classic(base_size=15) +
  theme(legend.position="none") +
  facet_wrap(~symptoms,
             labeller = labeller(symptoms = labelllername) 
  ) 
  
# ggsave("eachsym_quant.pdf", plot = eachsym_quant, width = 30, height = 15, units = "cm", dpi = 300)


## Total symptom level plot
summarized <- aggregated |>
  reframe(inter_quantile(totalsymptom), .by = t) |>
  pivot_wider(names_from = "quant", values_from = "val", names_glue = "q{quant}")

totalsym <- ggplot(data = summarized) +
  geom_area(data = shock_period, aes(x = time, y = shock*max(sde_out[,-1])*ncol(sde_out[,-1])
  ), inherit.aes = FALSE, fill = "darkgray", alpha = 0.2) +
  geom_line(aes(x = t, y = q0.5), col = "red4", lwd = 0.2) +
  geom_ribbon(aes(x=t,ymin=q0.25,ymax=q0.75),fill = "tomato1", alpha=0.3) +
  labs(x = "time", y= "") +
  geom_hline(yintercept = 5/3, linetype = 2, color = "azure4", lwd = 0.1) +
  geom_hline(yintercept = 10/3, linetype = 2, color = "azure4", lwd = 0.1) +
  theme_classic(base_size = 15)

# ggsave("totalsym_v2.pdf", plot = totalsym, width = 25, height =10, units = "cm", dpi = 300)




## fitting statistical networks
burnin <- 100
# before shock
beforeshock <- aggregated |> 
  filter(t > burnin & t < t_shock) |> 
  select(S_anh:S_sui) |> 
  estimateNetwork(default = "EBICglasso")
# during shock
duringshock <- aggregated |> 
  filter(t >= t_shock & t <= t_shock + shock_duration) |> 
  select(S_anh:S_sui) |>
  estimateNetwork(default = "EBICglasso")
# after shock
aftershock <- aggregated |> 
  # filter(t >= 3000) |> # keep the burn-in the same ! BUT use the old layout with t >= 3000 condition so old figure works --- this change practically doesn't change network a bit.
  filter(t > t_shock + shock_duration + burnin) |>
  select(S_anh:S_sui) |>
  estimateNetwork(default = "EBICglasso")


# hyperparameter *max_value*
max_value <- max(
  max(abs(beforeshock$graph)), 
  max(abs(duringshock$graph)),
  max(abs(aftershock$graph))
)

# hyperparameter *net_layout*
# net_layout <- averageLayout(beforeshock,
#                             duringshock,
#                             aftershock)
net_layout <- readRDS("data/net_layout.rds")


# pdf(file = "triplenetworks2_v2.pdf", width = 16, height = 9) 
# pdf(file = "triplenetworks2_v22.pdf", width = 16, height = 9) #add grey background during shock
par(mfrow = c(1, 3), mar = c(3, 0.5, 5, 0.5), xpd = NA)
plot(beforeshock, layout = net_layout, maximum = max_value, labels = nodenames, node.width=2)
title("(a) Before shock", line = 3, cex.main = 3)

plot(duringshock, layout = net_layout, maximum = max_value,  labels = nodenames,node.width=2)
title("(b) During shock", line = 3, cex.main = 3)
rect(-1.3, -2, 1.3, 2.0, col = rgb(0.5, 0.5, 0.5 ,alpha=0.1), border=FALSE) # add gray shade

plot(aftershock, layout = net_layout, maximum = max_value, labels = nodenames, node.width=2)
title("(c) After shock", line = 3, cex.main = 3)
abline(v = -3.9, lty=3, col = "lightgray")
abline(v = -1.3, lty=3, col = "lightgray")

# dev.off()


## fitting aggregated network
avgnet <- aggregated |> dplyr::select(S_anh:S_sui) 
colnames(avgnet) <- nodenames
avgnet_high <- aggregated_high |> dplyr::select(S_anh:S_sui)
colnames(avgnet_high) <- nodenames
avgnet_low <- aggregated_low |> dplyr::select(S_anh:S_sui)
colnames(avgnet_low) <- nodenames

# total average net  (this is fairly big)
totavgnet <- avgnet |> bind_rows(avgnet_high, avgnet_low) 
totavgnetwork <- estimateNetwork(totavgnet, default = "EBICglasso")

# totavgnetwork <- readRDS("data/overallnetwork.Rds")
totavgnet <- plot(totavgnetwork, labels = nodenames)

## adjust network layout
# aggnetLayout <- totavgnet$layout
# saveRDS(aggnetLayout, file = "aggLayout.Rds")
aggnetLayout <- readRDS("data/aggLayout.Rds")
totavgnet$layout <- aggnetLayout
# saveRDS(totavgnet, file = "overallNetworkvar.rds")

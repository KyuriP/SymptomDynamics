
## install packages
library(modelr)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(qgraph)
library(tidyverse)
library(bootnet)

rm(list=ls())
## source necessary functions
source("code/euler_stochastic2.R")
source("code/mod_specification.R")

set.seed(123)

## define model specifics: choose the scenario and initial value for symptoms
mod <- mod_spec(scenario = "base", init_val = 0.01)
mod_high <- mod_spec(scenario = "high", init_val = 0.01)
mod_low <- mod_spec(scenario = "low", init_val = 0.01)

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
  duration = shock_duration
  #seed = 123  # set seed
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

# run sims n_sims times
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

# saveRDS(aggregated, "aggregated2.rds") # delta = 8.5
saveRDS(aggregated, "aggregated3.rds") # delta = 9


inter_quantile <- function(x, probs = c(0.25, 0.5, 0.75)) {
  tibble(
    val = quantile(x, probs, na.rm = TRUE),
    quant = probs
  )
}

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

summarized <- aggregated |>
  reframe(inter_quantile(totalsymptom), .by = t) |>
  pivot_wider(names_from = "quant", values_from = "val", names_glue = "q{quant}")

totalsym <- ggplot(data = summarized) +
  geom_area(data = shock_period, aes(x = time, y = shock*max(sde_out[,-1])*ncol(sde_out[,-1])
  ), inherit.aes = FALSE, fill = "darkgray", alpha = 0.2) +
  geom_line(aes(x = t, y = q0.5), col = "red4", lwd = 0.2) +
  geom_ribbon(aes(x=t,ymin=q0.25,ymax=q0.75),fill = "tomato1", alpha=0.3) +
  #ggtitle("Average level of the aggregated symptoms") +
  labs(x = "time", y= "") +
  geom_hline(yintercept = 5/3, linetype = 2, color = "azure4", lwd = 0.1) +
  geom_hline(yintercept = 10/3, linetype = 2, color = "azure4", lwd = 0.1) +
  theme_classic(base_size = 15)

# ggsave("totalsym_v2.pdf", plot = totalsym, width = 25, height =10, units = "cm", dpi = 300)


# sde_out |>
#   mutate(totalsymptom = rowSums(pick(S_anh:S_sui))) |>
#   ggplot(aes (x = t, y = totalsymptom)) +
#   geom_line(col = "salmon") +
#   geom_area(data = shock_period, aes(x = time, y = shock*max(sde_out[,-1])*ncol(sde_out[,-1])
#   ), inherit.aes = FALSE, fill = "orange", alpha = 0.2) +
#   labs(y = "", title = "Total symptom level") +
#   theme_classic()



## fitting statistical networks:
burnout <- 100
beforeshock <- aggregated |> filter(t > burnout & t < t_shock) |> select(S_anh:S_sui)
# beforeshock <- aggregated |> filter(t < t_shock) |> select(S_anh:S_sui)
duringshock <- aggregated |> filter(t >= t_shock & t <= t_shock + shock_duration) |> select(S_anh:S_sui)
aftershock <- aggregated |> filter(t >= 3000) |> select(S_anh:S_sui)
# aftershock <- aggregated |> filter(t >= t_shock + shock_duration) |> select(S_anh:S_sui)


beforeshock <- estimateNetwork(beforeshock, default = "EBICglasso")
duringshock <- estimateNetwork(duringshock, default = "EBICglasso")
aftershock <- estimateNetwork(aftershock, default = "EBICglasso")


#' Creating hyperparameter *max_value*
max_value <- max(
  max(abs(beforeshock$graph)), 
  max(abs(duringshock$graph)),
  max(abs(aftershock$graph))
)

#' Creating hyperparameter *net_layout*
net_layout <- averageLayout(beforeshock,
                            duringshock,
                            aftershock)


# Creating node names
nodenames <- c("anhedonia", "sad", "sleep", "energy", "appetite", "guilty", "concentration", "motor", "suicidal")

# pdf(file = "triplenetworks2_v2.pdf", width = 16, height = 9)
par(mfrow = c(1, 3), mar = c(3, 0.5, 5, 0.5), xpd = NA)
plot(beforeshock, layout = net_layout, maximum = max_value, labels = colnames(A), node.width=2)
title("(a) Before shock", line = 3, cex.main = 3)

plot(duringshock, layout = net_layout, maximum = max_value,  labels = colnames(A),node.width=2)
title("(b) During shock", line = 3, cex.main = 3)

plot(aftershock, layout = net_layout, maximum = max_value, labels = colnames(A), node.width=2)
title("(c) After shock", line = 3, cex.main = 3)
abline(v = -3.9, lty=3, col = "lightgray")
abline(v = -1.3, lty=3, col = "lightgray")

# dev.off()



## fitting aggregated network
aggregated <- readRDS("data/aggregated.rds")
aggregated_res <- readRDS("data/aggregated_res.rds")
aggregated_low <- readRDS("data/aggregated_low.rds")

avgnet <- aggregated |> dplyr::select(S_anh:S_sui) 
colnames(avgnet) <- colnames(A)
avgnet_res <- aggregated_res |> dplyr::select(S_anh:S_sui)
colnames(avgnet_res) <- colnames(A)
avgnet_low <- aggregated_low |> dplyr::select(S_anh:S_sui)
colnames(avgnet_low) <- colnames(A)

#total average net 
totavgnet <- avgnet |> bind_rows(avgnet_res, avgnet_res) #|> slice_sample(prop = 0.8)
totavgnetwork <- estimateNetwork(totavgnet, default = "EBICglasso")

grp <- list(`cycle` = 2:6, `no cycle` = c(1,7:9))


totavgnetwork2 <- plot(totavgnetwork,  labels = colnames(A))
# totavgnetwork2 <- plot(totavgnetwork,  labels = colnames(A), groups = grp, color = c("#F5651C", "#58B5BC"), legend=F, border.color = "white",border.width = 2, label.color = "white")

#swap their locations (con & glt)
conloc <- totavgnetwork2$layout[7,]
gltloc <- totavgnetwork2$layout[6,]
totavgnetwork2$layout[6,] <- conloc
totavgnetwork2$layout[7,] <- gltloc


# png(file = "net_comp.png", width=2000, height=1000, bg = 'transparent')
# dev.off()
plot(totavgnetwork2)


# centrality
res_totavg <- centrality_auto(totavgnetwork)
rownames(res_totavg$node.centrality) <- colnames(A)

cent_totavg <- res_totavg$node.centrality$Strength |>as_data_frame() |>  mutate(node = factor(rownames(res_totavg$node.centrality), levels= c(rownames(res_totavg$node.centrality)))
)

cent_totavg$dummyB <- "Strength"
cent_totavg_plot <- cent_totavg |>
  ggplot(aes(x=value, y = node, group=1)) + 
  geom_point(color = "deepskyblue4") + 
  geom_path(color = "deepskyblue4") + 
  theme_bw() +
  labs(x = "", y = "", color = "") +
  theme(text=element_text(size=15))+
  facet_grid(. ~ dummyB)

ggsave("cent_sim.png", plot = cent_totavg_plot, width = 3.5, height =5.5, dpi = 300)

library(patchwork)

figpatches <- cent_totavg_plot + ~plot(totavgnetwork,  labels = colnames(A), node.width=2, main = "population network") 

pdf(file = "totalnetwork.pdf", width = 12, height = 8)
old_par <- par(mar = c(0, 2, 0, 0), bg = NA)
wrap_elements(panel = ~plot(totavgnetwork2,  labels = colnames(A), edge.color = "deepskyblue4"), clip = FALSE) + cent_totavg_plot + plot_layout(widths = c(1.7, 1))
par(old_par)
dev.off()

# per resilience level
avgnet <- estimateNetwork(avgnet, default = "EBICglasso", tuning = 0.1)
avgnet_res <- estimateNetwork(avgnet_res, default = "EBICglasso", tuning = 0.1)
avgnet_low <- estimateNetwork(avgnet_low, default = "EBICglasso", tuning = 0.1)

max_value_avg <- max(
  max(abs(avgnet$graph)),
  max(abs(avgnet_res$graph)),
  max(abs(avgnet_low$graph))
)

#' Creating hyperparameter *net_layout*
net_layout_avg <- averageLayout(avgnet,
                                avgnet_res,
                                avgnet_low)


plot(avgnet_low, layout = net_layout_avg, maximum = max_value_avg,  labels = colnames(A), node.width=2)
plot(avgnet, layout = net_layout_avg, maximum = max_value_avg,  labels = colnames(A), node.width=2)
plot(avgnet_res, layout = net_layout_avg, maximum = max_value_avg,  labels = colnames(A), node.width=2)


## centrality
reslow <- centrality_auto(avgnet_low)
resmid <- centrality_auto(avgnet)
reshigh <- centrality_auto(avgnet_res)

centrality_sim <- tibble(reslow2 = reslow$node.centrality$Strength, resmid2 = resmid$node.centrality$Strength, reshigh2 = reshigh$node.centrality$Strength) |> 
  mutate(node = factor(rownames(reslow$node.centrality), levels= c(rownames(reslow$node.centrality)))
  )

centrality_sim %<>%  pivot_longer(!node, names_to = "resilience", values_to = "value")

cent_sim <- centrality_sim |>
  ggplot(aes(x=value, y = node, group=1)) + 
  geom_point() + 
  geom_path() + 
  theme_bw() +
  labs(x = "", y = "") +
  theme(text=element_text(size=20))+
  facet_wrap(~resilience) 




reslow2 <- centrality_auto(duringshock_low)
resmid2 <- centrality_auto(duringshock)
reshigh2 <- centrality_auto(duringshock_res)

centrality_sim2 <- tibble(reslow3 = reslow2$node.centrality$Strength, resmid3 = resmid2$node.centrality$Strength, reshigh3 = reshigh2$node.centrality$Strength) |> 
  mutate(node = factor(rownames(reslow$node.centrality), levels= c(rownames(reslow$node.centrality)))
  )

centrality_sim2 %<>%  pivot_longer(!node, names_to = "resilience", values_to = "value")

cent_sim2 <- centrality_sim2 |>
  mutate(resilience = factor(resilience, levels = c("reslow3", "resmid3", "reshigh3"))) |>
  ggplot(aes(x=value, y = node, group=1)) + 
  geom_point(color = "deepskyblue4") + 
  geom_path(color = "deepskyblue4") + 
  theme_bw() +
  labs(x = "", y = "") +
  theme(text=element_text(size=20))+
  facet_wrap(~resilience, labeller = as_labeller(c("reslow3" = "low resilience", "resmid3" = "mid resilience", "reshigh3" = "high resilience"))) 

# centrality_sim2 |>
#   mutate(resilience = factor(resilience, levels = c("reslow3", "resmid3", "reshigh3"))) |>
#   ggplot(aes(x=value, y = node, group=resilience, color = resilience)) + 
#   geom_point() + 
#   geom_path() + 
#   theme_bw() +
#   labs(x = "", y = "") +
#   theme(text=element_text(size=20),
#         legend.position = "bottom") + 
#   coord_flip() 

ggsave("centrality_sim.pdf", plot = cent_sim2, width = 25, height = 16, units = "cm", dpi = 300)


# beforeGGM <- qgraph(cor_auto(beforeshock),
#                     graph = 'glasso',
#                     layout = 'spring',
#                     theme = 'colorblind',
#                     sampleSize = nrow(beforeshock),
#                     title = "Before shock")
# # fix layout
# L <- beforeGGM$layout
# duringGGM <- qgraph(cor_auto(duringshock),
#                    graph = 'glasso',
#                    theme = 'colorblind',
#                    sampleSize = nrow(duringshock),
#                    title = "During shock",
#                    layout = L)
# 
# afterGGM <- qgraph(cor_auto(aftershock),
#                    graph = 'glasso',
#                    theme = 'colorblind',
#                    sampleSize = nrow(aftershock),
#                    title = "After shock",
#                    layout = L)





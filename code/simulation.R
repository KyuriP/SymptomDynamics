
# install packages
library(modelr)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(qgraph)
#library(psychonetrics)
library(bootnet)

rm(list=ls())


# source the function
source("code/euler_stochastic2.R")

# change suicidal alpha to be the same
A <- matrix(c( .30, 0, 0, 0, 0, 0, 0, 0, 0,
               .33, .30, .14, .15, 0, .13, 0, 0, .15,
               .13, .14, .30, .22, .23, 0, 0, 0, 0,
               .21, .15, .22, .30, 0, 0, .12, 0, 0,
               0, 0, 0, .17, .30, 0, 0, 0, 0,
               0, .13, 0, 0, .15, .30, .2, .15, .22,
               0, 0, 0, 0, 0, 0, .30, .17, 0,
               0, 0, 0, 0, 0, 0, 0, .30, 0,
               0, 0, 0, 0, 0, 0, 0, .3, 0.30), 9, 9, byrow = T)


rownames(A) <- colnames(A) <- c("anh", "sad", "slp", "ene", "app", "glt", "con", "mot", "sui")

# define differential equations:
#equation <- "dS[i] ~ S[i]*(1-S[i])*(Beta[i] + A_i[i]*S[i] + A_[j][i]*S_[j]*(1+ delta[i]*f(S[i])))"
# dif_eq <- 1:9 |> map(function(x){
#   lhs <- paste0("dS", paste0("_",colnames(A)[x]))
#   rhs <- stringr::str_replace_all(
#     "Sk*(1-Sk)*(Betak + A_ik*Sk + A_jk*(1+ deltak*f(Sk)))", "k", paste0("_",colnames(A)[x]))
#   form <- as.formula(paste(lhs,"~", rhs))
# })
# 
# 
# dif_eq <- c(dS_anh ~ S_anh * (1 -  S_anh) * (Beta_anh + A[1,1] *  S_anh + 
#                                                          A[-1,1] %*% c(S_sad, S_slp, S_ene, S_app, S_glt, S_con, S_mot, S_sui)*(1 + delta_anh * f(S_anh))),
#             dS_sad ~ S_sad * (1 -  S_sad) * (Beta_sad + A[2,2] *  S_sad + 
#                                                          A[-2,2] %*% c(S_anh, S_slp, S_ene, S_app, S_glt, S_con, S_mot, S_sui)*(1 + delta_sad * f(S_sad))),
#             dS_slp ~ S_slp * (1 -  S_slp) * (Beta_slp + A[3,3] *  S_slp + 
#                                                          A[-3,3] %*% c(S_anh, S_app, S_ene, S_app, S_glt, S_con, S_mot, S_sui)*(1 + delta_slp * f(S_slp))),
#             dS_ene ~ S_ene * (1 -  S_ene) * (Beta_ene + A[4,4] *  S_ene + 
#                                                          A[-4,4] %*% c(S_anh, S_app, S_slp, S_app, S_glt, S_con, S_mot, S_sui)*(1 + delta_ene * f(S_ene))),
#             dS_app ~ S_app * (1 -  S_app) * (Beta_app + A[5,5] *  S_app + 
#                                                          A[-5,5] %*% c(S_anh, S_sad, S_slp, S_ene, S_glt, S_con, S_mot, S_sui)*(1 + delta_app * f(S_app))),
#             dS_glt ~ S_glt * (1 - S_glt) * (Beta_glt + A[6,6] * S_glt + 
#                                                          A[-6,6] %*% c(S_anh, S_sad, S_slp, S_ene, S_app, S_con, S_mot, S_sui)*(1 + delta_glt * f(S_glt))),
#             dS_con ~ S_con * (1 - S_con) * (Beta_con + A[7,7] *  S_con + 
#                                                          A[-7,7] %*% c(S_anh, S_sad, S_slp, S_ene, S_app, S_glt, S_mot, S_sui)*(1 + delta_con * f(S_con))),
#             dS_mot ~ S_mot * (1 -  S_mot) * (Beta_mot + A[8,8] *  S_mot + 
#                                                          A[-8,8] %*% c(S_anh, S_sad, S_slp, S_ene, S_app, S_glt, S_con, S_sui)*(1 + delta_mot * f(S_mot))),
#             dS_sui ~ S_sui * (1 -  S_sui) * (Beta_sui + A[9,9] *  S_sui + 
#                                                          A[-9,9] %*% c(S_anh, S_sad, S_slp, S_ene, S_app, S_glt, S_con, S_mot)*(1 + delta_sui * f(S_sui)))
#             
# )

## define "f"
f <- function(x) x^2

dif_eq <- 1:9 |> map(function(x){
  lhs <- paste0("dS", paste0("_",colnames(A)[x]))
  sumAj <- paste0("S", paste0("_",colnames(A)[-x]), collapse = ",") 
  rhs <- stringr::str_replace_all(
    "Sk * (1-Sk) * (Betak + A[q,q] * Sk + (1+ deltak * f(Sk))) * A[-q,q] %*% ", c(k = paste0("_",colnames(A)[x]), q = x)) |> paste0(paste0("c(", sumAj, ")"))
  form <- as.formula(paste(lhs,"~", rhs))
})


# give noise:
sto_eq <-  1:9 |> map(function(x){
  lhs <- paste0("dS", paste0("_",colnames(A)[x]))
  rhs <- 1  # change as you need
  form <- as.formula(paste(lhs,"~", rhs))
})

# sto_eq <-  c(dSA ~ .1,
#              dSB ~ .1,
#              dSC ~ .1)

# define the parameters (as a named vector):
A_i <- diag(A) |> set_names(paste0("A_i_", colnames(A))) 
A_j <- colSums(A) |> set_names(paste0("A_j_", colnames(A))) 
delta <- rep(9.5, 9) |> set_names(paste0("delta_", colnames(A))) 

## Healthy: sum(A_ji)*(1 + delta(S_i^4)) + beta_i) < -A_ii
# beta_i < -A_ii - sum(A_ji)*(1 + delta(S_i^2))


mid <- ((- A_i - A_j*(1+delta)) - A_j)/2
dist <- (- A_j - (- A_i - A_j*(1+delta))) / 100
mid + dist*45.5

dit <- (- A_j - (- A_i - A_j*(1+delta))) 

# Beta_healthy <- c(-6.68, -6.84, -6.88, -7.74, -4.56, -5.70, -4.60, -9.72, -10.90) |> set_names(paste0("Beta_", colnames(A))) 
Beta_healthy <- c(-1.94, -1.44, -1.32, -1.68, -1.36, -0.86, -1.24, -1.84, -0.74) |> set_names(paste0("Beta_", colnames(A))) 

## bistable: ...< beta_i < - A_j
- A_j*(1+2)
# Beta_bistable <- c(-5.58, -5.34, -5.78, -5.64, -5.46, -4.60, -4.50, -7.62, -6.81) |> set_names(paste0("Beta_", colnames(A)))
# Beta_bistable <- c(-1.07, -0.82, -0.76, -0.94, -0.78, -0.53, -0.72, -1.02, -0.62) |> set_names(paste0("Beta_", colnames (A)))
# Beta_bistable <- c(-1.24, -0.83, -0.73, -1.03, -0.77, -0.35, -0.67, -1.16, -1.01) |> set_names(paste0("Beta_", colnames(A)))
Beta_bistable <- c(-1.355, -0.995, -0.905, -1.165, -0.935, -0.575, -0.845, -1.285, -2.00) |> set_names(paste0("Beta_", colnames(A))) 
# Beta_bistable <- c(-1.246,-0.921, -0.843, -1.077, -0.869, -0.544, -0.791, -1.181, -1.00) |> set_names(paste0("Beta_", colnames(A))) 
rownames(A) <- colnames(A) <- c("anh", "sad", "slp", "ene", "app", "glt", "con", "mot", "sui")

- A_j + dist*20
Beta_bistable*(1 - 0.325)
## sick: beta_i > - A_j
# Beta_sick <- c(1.68, 2.04, 1.88, -2.74, -2.56, 1.95, -1.60, -3.72, -5.90) |> set_names(paste0("Beta_", colnames(A))
# Beta_sick <- c(-0.72, -0.47, -0.41, -0.59, -0.53, -0.18, -0.37, -0.67, -0.41) |> set_names(paste0("Beta_", colnames(A)))
# Beta_sick <- c(-0.90,   -0.65,   -0.59,   -0.77,   -0.61,   -0.36,   -0.55,   -0.85 , -1.25) |> set_names(paste0("Beta_", colnames(A)))
Beta_sick <- c(-0.91, -0.66, -0.60, -0.78, -0.62, -0.37, -0.56, -0.86, -1.30) |> set_names(paste0("Beta_", colnames(A))) #- A_j + dist*5

# Beta <- Beta_sick |> set_names(paste0("Beta", 1:9))

rescaling_factor <- 3
## original params
parms1 <- c(A, Beta_bistable, delta) 
## given shock (beta: stressor increases)  
parms2 <- c(A, Beta_sick, delta) 

# define the initial condition (as a named vector):
init <- c(S_anh = .01, 
          S_sad = .01, 
          S_slp = .01, 
          S_ene =.01, 
          S_app =.01, 
          S_glt =.01, 
          S_con =.01, 
          S_mot =.01, 
          S_sui =.01)

# define deltaT and the number of time steps:
deltaT <- .1 # timestep length
timelength <- 4000
n_steps <- as.integer(timelength / deltaT) # must be a number greater than 1

# specify the standard deviation of the stochastic noise
D_stoeq1 <- 0.01  # before shock
# D_stoeq2 <- 0.0015 / 20 # after shock
t_shock <- 1000
shock_duration <- 500

#set.seed(678) # set the seed (system's behavior varies much by random noise)
# simulate :
# > shock given at t = 500
# > shock --> change b and noise S.D.
# > shock duration = 50
# > after shock --> revert back to original parameters but noise S.D. remains larger (debating... reasonable?)

sde_out <- euler_stochastic2(
  Amat = A,
  deterministic_rate = dif_eq,
  stochastic_rate = sto_eq,
  initial_condition = init,
  parameters1 = parms1,
  deltaT = deltaT,
  timelength = timelength,
  D1 = D_stoeq1,
  shock = TRUE,
  parameters2 = parms2, 
  t_shock = t_shock, 
  duration = shock_duration,
  seed = 123
)

shock_period <- data.frame(time = 0:timelength) |>
  mutate(shock = ifelse(time >= t_shock & time <= t_shock + shock_duration, TRUE, FALSE))
  
  
sde_out |>
  tidyr::pivot_longer(!t, names_to = "symptoms") |>
  ggplot() +
  geom_line(aes(x = t, y = value, color = symptoms), alpha = 0.5) +
  geom_area(data = shock_period, aes(x = time, y = shock*max(sde_out[,-1])
  ), inherit.aes = FALSE, fill = "orange", alpha = 0.2) +
  scale_color_discrete(name = "Symptoms", labels = c("anhedonia", "sad", "sleep", "energy", "appetite", "guilty", "concentration", "motor", "suicidal"))+
  labs(y = "", title ="Each symptom level") +
  theme_classic() +
  theme(legend.position="bottom") 

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
  deterministic_rate = dif_eq,
  stochastic_rate = sto_eq,
  initial_condition = init,
  parameters1 = parms1,
  deltaT = deltaT,
  timelength = timelength,
  D1 = D_stoeq1,
  D2 = D_stoeq1,
  shock = TRUE,
  parameters2 = parms2, 
  t_shock = t_shock, 
  duration = shock_duration) |>
  mutate(totalsymptom = rowSums(pick(S_anh:S_sui)))
) |> list_rbind(names_to = "sim")

# saveRDS(aggregated, "aggregated.rds")

quantile_df <- function(x, probs = c(0.025, 0.5, 0.975)) {
  tibble(
    val = quantile(x, probs, na.rm = TRUE),
    quant = probs
  )
}

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

# beforeshock <- aggregated |> filter( t >=100 & t <= 101) |> select(S_anh:S_sui)
# duringshock <- aggregated |> filter(t >= 490 & t <= 491) |> select(S_anh:S_sui) 
# aftershock <- aggregated |> filter(t  >= 1000 & t <= 1001) |> select(S_anh:S_sui)

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
# net_layout[9,] <- c(0.3, -0.3)
# net_layout[7,1] <- c(-0.5)
# net_layout[5,] <- c(-0.5, -0.2)
# net_layout[4,1] <- 0.17
# net_layout[1,2] <- -0.7



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



## High resilience person

## adjacency matrix
# Ares <- matrix(c(0.3, 1.07, 0, 0, 0, 0, 0, 0, 0,
#                   1.07, 0.3, 1.5, 0, 0, 1.5, 0, 0, 1.13,
#                   0, 0, 0.3, 1.35, 1.38, 0, 0, 0, 0,
#                   1.27, 0, 0, 0.3, 0, 0, 0, 0, 0,
#                   0, 0, 0, 1.52, 0.3, 0, 0, 0, 0,
#                   0, 1.35, 1.44, 0, 0, 0.3, 1.3, 1.19, 1.07,
#                   0, 0, 0, 0, 0, 0, 0.3, 1.64, 0,
#                   0, 0, 0, 0, 0, 0, 0, 0.3, 0,
#                   0, 0, 0, 0, 0, 0, 0, 1.03, 0.1), 9, 9, byrow = T)

Ahigh <- matrix(c( .10, 0, 0, 0, 0, 0, 0, 0, 0,
                   .33, .10, .14, .15, 0, .13, 0, 0, .15,
                   .13, .14, .10, .22, .23, 0, 0, 0, 0,
                   .21, .15, .22, .10, 0, 0, .12, 0, 0,
                   0, 0, 0, .17, .10, 0, 0, 0, 0,
                   0, .13, 0, 0, .15, .10, .2, .15, .22,
                   0, 0, 0, 0, 0, 0, .10, .17, 0,
                   0, 0, 0, 0, 0, 0, 0, .10, 0,
                   0, 0, 0, 0, 0, 0, 0, .3, 0.02), 9, 9, byrow = T)

Ahigh <- matrix(c( .10, 0, 0, 0, 0, 0, 0, 0, 0,
                   .33, .10, .14, .15, 0, .13, 0, 0, .15,
                   .13, .14, .10, .22, .23, 0, 0, 0, 0,
                   .21, .15, .22, .10, 0, 0, .12, 0, 0,
                   0, 0, 0, .17, .10, 0, 0, 0, 0,
                   0, .13, 0, 0, .15, .10, .2, .15, .22,
                   0, 0, 0, 0, 0, 0, .10, .17, 0,
                   0, 0, 0, 0, 0, 0, 0, .10, 0,
                   0, 0, 0, 0, 0, 0, 0, .3, 0.10), 9, 9, byrow = T)

rownames(Ahigh) <- colnames(Ahigh) <- c("anh", "sad", "slp", "ene", "app", "glt", "con", "mot", "sui")


# define the parameters (as a named vector):
Ahigh_i <- diag(Ahigh) |> set_names(paste0("A_i_", colnames(Ahigh))) # for now, set it all to 1
Ahigh_j <- colSums(Ahigh) |> set_names(paste0("A_j_", colnames(Ahigh)))
delta_high <- rep(6, 9) |> set_names(paste0("delta_", colnames(Ahigh)))

## Healthy: sum(A_ji)*(1 + delta(S_i^4)) + beta_i) < -A_ii
# beta_i < -A_ii - sum(A_ji)*(1 + delta(S_i^4))
(- A_i - A_j*2) # as the max of delta(S_i^4) == delta*1

# Beta_healthy <- c(-6.68, -6.84, -6.88, -7.74, -4.56, -5.70, -4.60, -9.72, -10.90) |> set_names(paste0("Beta_", colnames(A))) 
# Beta_bistable <- c(-5.58, -5.34, -5.78, -5.64, -5.46, -4.60, -4.50, -7.62, -6.81) |> set_names(paste0("Beta_", colnames(A)))
# # Beta_sick <- c(5.18, 5.34, 5.38, 6.24, 3.06, 8.50, -3.10, -8.22, -3.90) |> set_names(paste0("Beta", 1:9))
# Beta_sick <- c(1.68, 2.04, 1.88, -2.74, -2.56, 1.95, -1.60, -3.72, -5.90) |> set_names(paste0("Beta_", colnames(A)))
# Beta <- Beta_sick |> set_names(paste0("Beta", 1:9))

## original params
parms1H <- c(Ahigh, Beta_bistable, delta_high) # / rescaling_factor
## given shock (beta: stressor increases)  
parms2H <- c(Ahigh, Beta_sick, delta_high)  #/ rescaling_factor

# define the initial condition (as a named vector):
init <- c(S_anh = .01, 
          S_sad = .01, 
          S_slp = .01, 
          S_ene =.01, 
          S_app =.01, 
          S_glt =.01, 
          S_con =.01, 
          S_mot =.01, 
          S_sui =.01)

# define deltaT and the number of time steps:
# deltaT <- .1 # timestep length
# timelength <- 2000
# n_steps <- as.integer(timelength / deltaT) # must be a number greater than 1

# specify the standard deviation of the stochastic noise
# D_stoeq1 <- 0.01 # before shock
# D_stoeq2 <- 0.02 # after shock
# t_shock <- 300
# shock_duration <- 200

#set.seed(678) # set the seed (system's behavior varies much by random noise)
# simulate :
# > shock given at t = 500
# > shock --> change b and noise S.D.
# > shock duration = 50
# > after shock --> revert back to original parameters but noise S.D. remains larger (debating... reasonable?)

sde_out_res <- euler_stochastic2(
  deterministic_rate = dif_eq,
  stochastic_rate = sto_eq,
  initial_condition = init,
  parameters1 = parms1H,
  deltaT = deltaT,
  timelength = timelength,
  D1 = D_stoeq1,
  D2 = D_stoeq1,
  shock = TRUE,
  parameters2 = parms2H, 
  t_shock = t_shock, 
  duration = shock_duration,
  seed = 123
)

shock_period <- data.frame(time = 0:timelength) |>
  mutate(shock = ifelse(time >= t_shock & time <= t_shock + shock_duration, TRUE, FALSE))

sde_out_res |>
  tidyr::pivot_longer(!t, names_to = "symptoms") |>
  ggplot() +
  geom_line(aes(x = t, y = value, color = symptoms), alpha = 0.5) +
  geom_area(data = shock_period, aes(x = time, y = shock*max(sde_out[,-1])
  ), inherit.aes = FALSE, fill = "orange", alpha = 0.2) +
  scale_color_discrete(name = "Symptoms", labels = c("anhedonia", "sad", "sleep", "energy", "appetite", "guilty", "concentration", "motor", "suicidal"))+
  labs(y = "", title ="Each symptom level") +
  theme_classic() +
  theme(legend.position="bottom") 


# labeller
labelllername <- c("S_anh" = "anhedonia", "S_app" = "appetite", "S_con" = "concentration", "S_ene" = "energy", "S_glt" = "guilty", "S_mot" = "motor", "S_sad" = "sad", "S_slp" = "sleep", "S_sui" = "suicidal") 

eachsym_res <- sde_out_res |>
  tidyr::pivot_longer(!t, names_to = "symptoms") |>
  ggplot() +
  geom_area(data = shock_period, aes(x = time, y = shock*max(sde_out[,-1])
  ), inherit.aes = FALSE, fill = "darkgray", alpha = 0.2) +
  geom_line(aes(x = t, y = value, color = symptoms), alpha = 0.5, lwd= 0.2) +
  scale_color_discrete(name = "Symptoms", labels = c("anhedonia", "sad", "sleep", "energy", "appetite", "guilty", "concentration", "motor", "suicidal"))+
  #labs(x="time", y = "", title ="Trajectory of each symptom") +
  labs(x="time", y = "") +
  theme_classic(base_size = 15) +
  theme(legend.position="none") +
  facet_wrap(~symptoms, labeller = labeller(symptoms = labelllername) 
  )

# ggsave("eachsym_resv2.pdf", plot = eachsym_res, width = 30, height =15, units = "cm", dpi = 300)


## aggregate symptom level
# n_sims <- 300

# run sims n_sims times
aggregated_res <- map(1:n_sims, ~ euler_stochastic2(
  deterministic_rate = dif_eq,
  stochastic_rate = sto_eq,
  initial_condition = init,
  parameters1 = parms1H,
  deltaT = deltaT,
  timelength = timelength,
  D1 = D_stoeq1,
  D2 = D_stoeq1,
  shock = TRUE,
  parameters2 = parms2H, 
  t_shock = t_shock, 
  duration = shock_duration) |>
    mutate(totalsymptom = rowSums(pick(S_anh:S_sui)))
) |> list_rbind(names_to = "sim")

# saveRDS(aggregated_res, "aggregated_res.rds")
# aggregated_res <- readRDS("data/aggregated_res.rds")



quantile_df <- function(x, probs = c(0.025, 0.5, 0.975)) {
  tibble(
    val = quantile(x, probs, na.rm = TRUE),
    quant = probs
  )
}




each_summ_res <- aggregated_res |>
  select(-c(totalsymptom, sim)) |>
  tidyr::pivot_longer(!t, names_to = "symptoms") |>
  reframe(inter_quantile(value), .by = c(t, symptoms)) |>
  pivot_wider(names_from = "quant", values_from = "val", names_glue = "q{quant}") 

eachsym_quant_res <- ggplot(data = each_summ_res) +
  geom_area(data = shock_period, aes(x = time, y = shock*max(sde_out[,-1])
  ), inherit.aes = FALSE, fill = "darkgray", alpha = 0.2) +
  geom_line(aes(x = t, y = q0.5, color = symptoms), alpha = 0.5, lwd = 0.2) +
  geom_ribbon(aes(x=t,ymin=q0.25,ymax=q0.75, fill = symptoms), alpha=0.3) +
  scale_color_discrete(name = "Symptoms", labels = c("anhedonia", "sad", "sleep", "energy", "appetite", "guilty", "concentration", "motor", "suicidal")) +
  labs(x = "time", y= "") +
  theme_classic(base_size=15) +
  theme(legend.position="none") +
  facet_wrap(~symptoms,
             labeller = labeller(symptoms = labelllername) 
  )

# ggsave("eachsym_quant_res.pdf", plot = eachsym_quant_res, width = 30, height = 15, units = "cm", dpi = 300)


summarized_res <- aggregated_res |>
  reframe(inter_quantile(totalsymptom), .by = t) |>
  pivot_wider(names_from = "quant", values_from = "val", names_glue = "q{quant}")


totalsym_res <- ggplot(data = summarized_res) +
  geom_area(data = shock_period, aes(x = time, y = shock*max(sde_out[,-1])*ncol(sde_out[,-1])
  ), inherit.aes = FALSE, fill = "darkgray", alpha = 0.2) +
  geom_line(aes(x = t, y = q0.5), col = "red4", lwd = 0.2) +
  geom_ribbon(aes(x=t,ymin=q0.25,ymax=q0.75),fill = "tomato1", alpha=0.3) +
  #ggtitle("Average level of the aggregated symptoms") +
  labs(x = "time", y= "") +
  geom_hline(yintercept = 5/3, linetype =2, color = "azure4", lwd = 0.1) +
  geom_hline(yintercept = 10/3, linetype = 2, color = "azure4", lwd = 0.1) +
  theme_classic(base_size = 15)

# ggsave("totalsym_resv2.pdf", plot = totalsym_res, width = 25, height =10, units = "cm", dpi = 300)



## fitting statistical networks:
burnout <- 100
beforeshock_res <- aggregated_res |> filter(t > burnout & t < t_shock) |> select(S_anh:S_sui)
# beforeshock_res <- aggregated_res |> filter(t < t_shock) |> select(S_anh:S_sui)
duringshock_res <- aggregated_res |> filter(t >= t_shock & t <= t_shock + shock_duration) |> select(S_anh:S_sui) 
aftershock_res <- aggregated_res |> filter(t >= 3000) |> select(S_anh:S_sui)
# aftershock_res <- aggregated_res |> filter(t >=t_shock + shock_duration) |> select(S_anh:S_sui)

# beforeshock_res <- aggregated_res |> filter(t > burnout & t < 250) |> select(S_anh:S_sui)
# duringshock_res <- aggregated_res |> filter(t >= 400 & t <= 500) |> select(S_anh:S_sui) 
# aftershock_res <- aggregated_res |> filter(t >=1500) |> select(S_anh:S_sui)


beforeshock_res <- estimateNetwork(beforeshock_res, default = "EBICglasso")
duringshock_res <- estimateNetwork(duringshock_res, default = "EBICglasso")
aftershock_res <- estimateNetwork(aftershock_res, default = "EBICglasso")

#' #' Creating hyperparameter *max_value*
#' max_value_res <- max(
#'   max(abs(beforeshock_res$graph)), 
#'   max(abs(duringshock_res$graph)),
#'   max(abs(aftershock_res$graph))
#' )

#' #' Creating hyperparameter *net_layout*
#' net_layout_res <- averageLayout(beforeshock_res,
#'                             duringshock_res,
#'                             aftershock_res)
#' #net_layout[9, ] <- c(0.6, -0.3)

# Creating node names
nodenames <- c("anhedonia", "sad", "sleep", "energy", "appetite", "guilty", "concentration", "motor", "suicidal")

# pdf(file = "triplenetworks_res2_v2.pdf", width = 16, height = 9)
par(mfrow = c(1, 3), mar = c(3, 0.5, 5, 0.5), xpd = NA)
plot(beforeshock_res, layout = net_layout, maximum = max_value, labels = colnames(A), node.width=2)
title("(a) Before shock", line = 3, cex.main = 3)

plot(duringshock_res, layout = net_layout, maximum = max_value,  labels = colnames(A),node.width=2)
title("(b) During shock", line = 3, cex.main = 3)

plot(aftershock_res, layout = net_layout, maximum = max_value, labels = colnames(A), node.width=2)
title("(c) After shock", line = 3, cex.main = 3)

abline(v = -3.9, lty=3, col = "lightgray")
abline(v = -1.3, lty=3, col = "lightgray")

# dev.off()




#########################
## LOW resilience person
#########################

## adjacency matrix
# Alow <- matrix(c(1.3, 1.07, 0, 0, 0, 0, 0, 0, 0,
#                  1.07, 1.3, 1.5, 0, 0, 1.5, 0, 0, 1.13,
#                  0, 0, 1.3, 1.35, 1.38, 0, 0, 0, 0,
#                  1.27, 0, 0, 1.3, 0, 0, 0, 0, 0,
#                  0, 0, 0, 1.52, 1.3, 0, 0, 0, 0,
#                  0, 1.35, 1.44, 0, 0, 1.3, 1.3, 1.19, 1.07,
#                  0, 0, 0, 0, 0, 0, 1.3, 1.64, 0,
#                  0, 0, 0, 0, 0, 0, 0, 1.3, 0,
#                  0, 0, 0, 0, 0, 0, 0, 1.03, 1), 9, 9, byrow = T)

Alow <- matrix(c( .50, 0, 0, 0, 0, 0, 0, 0, 0,
                  .33, .50, .14, .15, 0, .13, 0, 0, .15,
                  .13, .14, .50, .22, .23, 0, 0, 0, 0,
                  .21, .15, .22, .50, 0, 0, .12, 0, 0,
                  0, 0, 0, .17, .50, 0, 0, 0, 0,
                  0, .13, 0, 0, .15, .50, .2, .15, .22,
                  0, 0, 0, 0, 0, 0, .50, .17, 0,
                  0, 0, 0, 0, 0, 0, 0, .50, 0,
                  0, 0, 0, 0, 0, 0, 0, .30, 0.15), 9, 9, byrow = T)

Alow <- matrix(c( .50, 0, 0, 0, 0, 0, 0, 0, 0,
                  .33, .50, .14, .15, 0, .13, 0, 0, .15,
                  .13, .14, .50, .22, .23, 0, 0, 0, 0,
                  .21, .15, .22, .50, 0, 0, .12, 0, 0,
                  0, 0, 0, .17, .50, 0, 0, 0, 0,
                  0, .13, 0, 0, .15, .50, .2, .15, .22,
                  0, 0, 0, 0, 0, 0, .50, .17, 0,
                  0, 0, 0, 0, 0, 0, 0, .50, 0,
                  0, 0, 0, 0, 0, 0, 0, .30, 0.50), 9, 9, byrow = T)

rownames(Alow) <- colnames(Alow) <- c("anh", "sad", "slp", "ene", "app", "glt", "con", "mot", "sui")


# define the parameters (as a named vector):
Alow_i <- diag(Alow) |> set_names(paste0("A_i_", colnames(Alow)))
Alow_j <- colSums(Alow) |> set_names(paste0("A_j_", colnames(Alow)))
delta_low <- rep(11, 9) |> set_names(paste0("delta_", colnames(Alow)))

## Healthy: sum(A_ji)*(1 + delta(S_i^4)) + beta_i) < -A_ii
# beta_i < -A_ii - sum(A_ji)*(1 + delta(S_i^4))
(- A_i - A_j*2) # as the max of delta(S_i^4) == delta*1
# Beta_healthy <- c(-6.68, -6.84, -6.88, -7.74, -4.56, -5.70, -4.60, -9.72, -10.90) |> set_names(paste0("Beta_", colnames(A))) 
# Beta_bistable <- c(-5.58, -5.34, -5.78, -5.64, -5.46, -4.60, -4.50, -7.62, -6.81) |> set_names(paste0("Beta_", colnames(A)))
# # Beta_sick <- c(5.18, 5.34, 5.38, 6.24, 3.06, 8.50, -3.10, -8.22, -3.90) |> set_names(paste0("Beta", 1:9))
# Beta_sick <- c(1.68, 2.04, 1.88, -2.74, -2.56, 1.95, -1.60, -3.72, -5.90) |> set_names(paste0("Beta_", colnames(A)))
# Beta <- Beta_sick |> set_names(paste0("Beta", 1:9))

## original params
parms1L <- c(Alow, Beta_bistable, delta_low) #/ rescaling_factor
## given shock (beta: stressor increases)  
parms2L <- c(Alow, Beta_sick, delta_low) #/ rescaling_factor

# define the initial condition (as a named vector):
init <- c(S_anh = .01, 
          S_sad = .01, 
          S_slp = .01, 
          S_ene =.01, 
          S_app =.01, 
          S_glt =.01, 
          S_con =.01, 
          S_mot =.01, 
          S_sui =.01)

# define deltaT and the number of time steps:
# deltaT <- .1 # timestep length
# timelength <- 2500
# n_steps <- as.integer(timelength / deltaT) # must be a number greater than 1

# specify the standard deviation of the stochastic noise
# D_stoeq1 <- 0.01 # before shock
# D_stoeq2 <- 0.02 # after shock
# t_shock <- 300
# shock_duration <- 200

#set.seed(678) # set the seed (system's behavior varies much by random noise)
# simulate :
# > shock given at t = 500
# > shock --> change b and noise S.D.
# > shock duration = 50
# > after shock --> revert back to original parameters but noise S.D. remains larger (debating... reasonable?)

sde_out_low <- euler_stochastic2(
  deterministic_rate = dif_eq,
  stochastic_rate = sto_eq,
  initial_condition = init,
  parameters1 = parms1L,
  deltaT = deltaT,
  timelength = timelength,
  D1 = D_stoeq1,
  D2 = D_stoeq1,
  shock = TRUE,
  parameters2 = parms2L, 
  t_shock = t_shock, 
  duration = shock_duration,
  seed = 123
)

shock_period <- data.frame(time = 0:timelength) |>
  mutate(shock = ifelse(time >= t_shock & time <= t_shock + shock_duration, TRUE, FALSE))

sde_out_low |>
  tidyr::pivot_longer(!t, names_to = "symptoms") |>
  ggplot() +
  geom_line(aes(x = t, y = value, color = symptoms), alpha = 0.5) +
  geom_area(data = shock_period, aes(x = time, y = shock*max(sde_out[,-1])
  ), inherit.aes = FALSE, fill = "orange", alpha = 0.2) +
  scale_color_discrete(name = "Symptoms", labels = c("anhedonia", "sad", "sleep", "energy", "appetite", "guilty", "concentration", "motor", "suicidal"))+
  labs(y = "", title ="Each symptom level") +
  theme_classic() +
  theme(legend.position="bottom") 


# labeller
labelllername <- c("S_anh" = "anhedonia", "S_app" = "appetite", "S_con" = "concentration", "S_ene" = "energy", "S_glt" = "guilty", "S_mot" = "motor", "S_sad" = "sad", "S_slp" = "sleep", "S_sui" = "suicidal") 

eachsym_low <- sde_out_low |>
  tidyr::pivot_longer(!t, names_to = "symptoms") |>
  ggplot() +
  geom_area(data = shock_period, aes(x = time, y = shock*max(sde_out[,-1])
  ), inherit.aes = FALSE, fill = "darkgray", alpha = 0.2) +
  geom_line(aes(x = t, y = value, color = symptoms), alpha = 0.5, lwd= 0.2) +
  scale_color_discrete(name = "Symptoms", labels = c("anhedonia", "sad", "sleep", "energy", "appetite", "guilty", "concentration", "motor", "suicidal"))+
  #labs(x="time", y = "", title ="Trajectory of each symptom") +
  labs(x="time", y = "") +
  theme_classic(base_size = 15) +
  theme(legend.position="none") +
  facet_wrap(~symptoms,
             labeller = labeller(symptoms = labelllername) 
  )

# ggsave("eachsym_lowv2.pdf", plot = eachsym_low, width = 30, height =15, units = "cm", dpi = 300)


## aggregate symptom level
# n_sims <- 300

# run sims n_sims times
aggregated_low <- map(1:n_sims, ~ euler_stochastic2(
  deterministic_rate = dif_eq,
  stochastic_rate = sto_eq,
  initial_condition = init,
  parameters1 = parms1L,
  deltaT = deltaT,
  timelength = timelength,
  D1 = D_stoeq1,
  D2 = D_stoeq1,
  shock = TRUE,
  parameters2 = parms2L, 
  t_shock = t_shock, 
  duration = shock_duration) |>
    mutate(totalsymptom = rowSums(pick(S_anh:S_sui)))
) |> list_rbind(names_to = "sim")

# saveRDS(aggregated_low, "aggregated_low.rds")

quantile_df <- function(x, probs = c(0.025, 0.5, 0.975)) {
  tibble(
    val = quantile(x, probs, na.rm = TRUE),
    quant = probs
  )
}


each_summ_low <- aggregated_low |>
  select(-c(totalsymptom, sim)) |>
  tidyr::pivot_longer(!t, names_to = "symptoms") |>
  reframe(inter_quantile(value), .by = c(t, symptoms)) |>
  pivot_wider(names_from = "quant", values_from = "val", names_glue = "q{quant}") 

eachsym_quant_low <- ggplot(data = each_summ_low) +
  geom_area(data = shock_period, aes(x = time, y = shock*max(sde_out[,-1])
  ), inherit.aes = FALSE, fill = "darkgray", alpha = 0.1) +
  geom_line(aes(x = t, y = q0.5, color = symptoms), alpha = 0.5, lwd = 0.2) +
  geom_ribbon(aes(x=t,ymin=q0.25,ymax=q0.75, fill = symptoms), alpha=0.3) +
  scale_color_discrete(name = "Symptoms", labels = c("anhedonia", "sad", "sleep", "energy", "appetite", "guilty", "concentration", "motor", "suicidal")) +
  labs(x = "time", y= "") +
  theme_classic(base_size=15) +
  theme(legend.position="none") +
  facet_wrap(~symptoms,
             labeller = labeller(symptoms = labelllername) 
  )

# ggsave("eachsym_quant_low.pdf", plot = eachsym_quant_low, width = 30, height = 15, units = "cm", dpi = 300)


summarized_low <- aggregated_low |>
  reframe(inter_quantile(totalsymptom), .by = t) |>
  pivot_wider(names_from = "quant", values_from = "val", names_glue = "q{quant}")


totalsym_low <- ggplot(data = summarized_low) +
  geom_area(data = shock_period, aes(x = time, y = shock*max(sde_out[,-1])*ncol(sde_out[,-1])
  ), inherit.aes = FALSE, fill = "darkgray", alpha = 0.2) +
  geom_line(aes(x = t, y = q0.5), col = "red4", lwd = 0.2) +
  geom_ribbon(aes(x=t,ymin=q0.25,ymax=q0.75),fill = "tomato1", alpha=0.3) +
  #ggtitle("Average level of the aggregated symptoms") +
  labs(x = "time", y= "") +
  geom_hline(yintercept = 5/3, linetype =2, color = "azure4", lwd = 0.1) +
  geom_hline(yintercept = 10/3, linetype = 2, color = "azure4", lwd = 0.1) +
  theme_classic(base_size = 15)

# ggsave("totalsym_lowv2.pdf", plot = totalsym_low, width = 25, height =10, units = "cm", dpi = 300)



## fitting statistical networks:
burnout <- 100
beforeshock_low <- aggregated_low |> filter(t > burnout & t < t_shock) |> select(S_anh:S_sui)
# beforeshock_low <- aggregated_low |> filter(t < t_shock) |> select(S_anh:S_sui)
duringshock_low <- aggregated_low |> filter(t >= t_shock & t <= t_shock + shock_duration) |> select(S_anh:S_sui) 
# aftershock_low <- aggregated_low |> filter(t >= t_shock + shock_duration) |> select(S_anh:S_sui)
aftershock_low <- aggregated_low |> filter(t >= 3000) |> select(S_anh:S_sui)
# aftershock_low <- aggregated_low |> filter(t >=t_shock + shock_duration + burnout) |> select(S_anh:S_sui)

# 
# beforeshock_low <- aggregated_low |> filter(t > burnout & t < 250) |> select(S_anh:S_sui)
# duringshock_low <- aggregated_low |> filter(t >= 400 & t <= 500) |> select(S_anh:S_sui) 
# aftershock_low <- aggregated_low |> filter(t >=1000) |> select(S_anh:S_sui)


beforeshock_low <- estimateNetwork(beforeshock_low, default = "EBICglasso")
duringshock_low <- estimateNetwork(duringshock_low, default = "EBICglasso")
aftershock_low <- estimateNetwork(aftershock_low, default = "EBICglasso")



#' #' Creating hyperparameter *max_value*
#' max_value_res <- max(
#'   max(abs(beforeshock_res$graph)), 
#'   max(abs(duringshock_res$graph)),
#'   max(abs(aftershock_res$graph))
#' )

#' #' Creating hyperparameter *net_layout*
#' net_layout_res <- averageLayout(beforeshock_res,
#'                             duringshock_res,
#'                             aftershock_res)
#' #net_layout[9, ] <- c(0.6, -0.3)

# Creating node names
nodenames <- c("anhedonia", "sad", "sleep", "energy", "appetite", "guilty", "concentration", "motor", "suicidal")

# pdf(file = "triplenetworks_low2_v2.pdf", width = 16, height = 9)
par(mfrow = c(1, 3), mar = c(3, 0.5, 5, 0.5), xpd = NA)
plot(beforeshock_low, layout = net_layout, maximum = max_value, labels = colnames(A), node.width=2)
title("(a) Before shock", line = 3, cex.main = 3)

plot(duringshock_low, layout = net_layout, maximum = max_value,  labels = colnames(A),node.width=2)
title("(b) During shock", line = 3, cex.main = 3)

plot(aftershock_low, layout = net_layout, maximum = max_value, labels = colnames(A), node.width=2)
title("(c) After shock", line = 3, cex.main = 3)
abline(v = -3.9, lty=3, col = "lightgray")
abline(v = -1.3, lty=3, col = "lightgray")

# dev.off()
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





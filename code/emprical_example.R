
# install packages
library(modelr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(qgraph)
library(psychonetrics)

source("code/euler_stochastic2.R")

f <- function(x) x^2

## adjacency matrix
A <- matrix(c( 1, 1.07, 0, 0, 0, 0, 0, 0, 0,
               1.07, 1, 1.5, 0, 0, 1.5, 0, 0, 1.13,
               0, 0, 1, 1.35, 1.38, 0, 0, 0, 0,
               1.27, 0, 0, 1, 0, 0, 0, 0, 0,
               0, 0, 0, 1.52, 1, 0, 0, 0, 0,
               0, 1.35, 1.44, 0, 0, 1, 1.3, 1.19, 1.07,
               0, 0, 0, 0, 0, 0, 1, 1.64, 0,
               0, 0, 0, 0, 0, 0, 0, 1, 0,
               0, 0, 0, 0, 0, 0, 0, 1.03, 0.5), 9, 9, byrow = T)
rownames(A) <- colnames(A) <- c("anh", "sad", "slp", "ene", "app", "glt", "con", "mot", "sui")


# define differential equations:
#equation <- "dS[i] ~ S[i]*(1-S[i])*(Beta[i] + A_i[i]*S[i] + A_j[i]*(1+ delta[i]*f(S[i])))"
dif_eq <- 1:9 |> map(function(x){
  lhs <- paste0("dS", paste0("_",colnames(A)[x]))
  rhs <- stringr::str_replace_all(
    "Sk*(1-Sk)*(Betak + A_ik*Sk + A_jk*(1+ deltak*f(Sk)))", "k", paste0("_",colnames(A)[x]))
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
A_i <- diag(A) |> set_names(paste0("A_i_", colnames(A))) # for now, set it all to 1
A_j <- colSums(A) |> set_names(paste0("A_j_", colnames(A)))
delta <- rep(1.3, 9) |> set_names(paste0("delta_", colnames(A)))
## Healthy: sum(A_ji)*(1 + delta(S_i^4)) + beta_i) < -A_ii
# beta_i < -A_ii - sum(A_ji)*(1 + delta(S_i^4))
(- A_i - A_j*2) # as the max of delta(S_i^4) == delta*1
Beta_healthy <- c(-6.68, -6.84, -6.88, -7.74, -4.56, -5.70, -4.60, -9.72, -10.90) |> set_names(paste0("Beta_", colnames(A))) 
Beta_bistable <- c(-5.58, -5.34, -5.78, -5.64, -5.46, -4.60, -4.50, -7.62, -6.81) |> set_names(paste0("Beta_", colnames(A)))
# Beta_sick <- c(5.18, 5.34, 5.38, 6.24, 3.06, 8.50, -3.10, -8.22, -3.90) |> set_names(paste0("Beta", 1:9))
Beta_sick <- c(1.68, 2.04, 1.88, -2.74, -2.56, 1.95, -1.60, -3.72, -5.90) |> set_names(paste0("Beta_", colnames(A)))
# Beta <- Beta_sick |> set_names(paste0("Beta", 1:9))

## original params
parms1 <- c(A_i, A_j, Beta_bistable, delta)
## given shock (beta: stressor increases)  
parms2 <- c(A_i, A_j, Beta_sick, delta)

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
timelength <- 1000
n_steps <- as.integer(timelength / deltaT) # must be a number greater than 1

# specify the standard deviation of the stochastic noise
D_stoeq1 <- 0.01 # before shock
D_stoeq2 <- 0.02 # after shock
t_shock <- 300
shock_duration <- 200

#set.seed(678) # set the seed (system's behavior varies much by random noise)
# simulate :
# > shock given at t = 500
# > shock --> change b and noise S.D.
# > shock duration = 50
# > after shock --> revert back to original parameters but noise S.D. remains larger (debating... reasonable?)

sde_out <- euler_stochastic2(
  deterministic_rate = dif_eq,
  stochastic_rate = sto_eq,
  initial_condition = init,
  parameters1 = parms1,
  deltaT = deltaT,
  timelength = timelength,
  D1 = D_stoeq1,
  D2 = D_stoeq2,
  shock = TRUE,
  parameters2 = parms2, 
  t_shock = t_shock, 
  duration = shock_duration,
  seed = 100
)

shock_period <- data.frame(time = 0:timelength) |>
  mutate(shock = ifelse(time >= t_shock & time <= t_shock + shock_duration, TRUE, FALSE))

sde_out |>
  tidyr::pivot_longer(!t, names_to = "symptoms") |>
  ggplot() +
  geom_line(aes(x = t, y = value, color = symptoms), alpha = 0.5) +
  geom_area(data = shock_period, aes(x = time, y = shock*max(sde_out[,-1])
  ), inherit.aes = FALSE, fill = "orange", alpha = 0.2) +
  scale_color_discrete(name = "Symptoms", labels = c("anhedonia", "appetite", "concentration", "energy", "guilty", "motor", "sad", "sleep", "suicidal"))+
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
  scale_color_discrete(name = "Symptoms", labels = c("anhedonia", "appetite", "concentration", "energy", "guilty", "motor", "sad", "sleep", "suicidal"))+
  #labs(x="time", y = "", title ="Trajectory of each symptom") +
  labs(x="time", y = "") +
  theme_classic() +
  theme(legend.position="none") +
  facet_wrap(~symptoms,
  labeller = labeller(symptoms = labelllername) 
  )

ggsave("eachsym.pdf", plot = eachsym, width = 20, height =10, units = "cm", dpi = 300)


## aggregate symptom level
n_sims <- 100

# run sims n_sims times
aggregated <- map(1:n_sims, ~ euler_stochastic2(
  deterministic_rate = dif_eq,
  stochastic_rate = sto_eq,
  initial_condition = init,
  parameters1 = parms1,
  deltaT = deltaT,
  timelength = timelength,
  D1 = D_stoeq1,
  D2 = D_stoeq2,
  shock = TRUE,
  parameters2 = parms2, 
  t_shock = t_shock, 
  duration = shock_duration) |>
  mutate(totalsymptom = rowSums(pick(S_anh:S_sui)))
) |> list_rbind(names_to = "sim")


quantile_df <- function(x, probs = c(0.025, 0.5, 0.975)) {
  tibble(
    val = quantile(x, probs, na.rm = TRUE),
    quant = probs
  )
}

summarized <- aggregated |>
  reframe(quantile_df(totalsymptom), .by = t) |>
  pivot_wider(names_from = "quant", values_from = "val", names_glue = "q{quant}")


totalsym <- ggplot(data = summarized) +
  geom_area(data = shock_period, aes(x = time, y = shock*max(sde_out[,-1])*ncol(sde_out[,-1])
  ), inherit.aes = FALSE, fill = "darkgray", alpha = 0.2) +
  geom_line(aes(x = t, y = q0.5), col = "red4", lwd = 0.2) +
  geom_ribbon(aes(x=t,ymin=q0.025,ymax=q0.975),fill = "tomato1", alpha=0.3) +
  #ggtitle("Average level of the aggregated symptoms") +
  labs(x = "time", y= "") +
  geom_hline(yintercept = 5/3, linetype = 2, color = "azure4", lwd = 0.1) +
  geom_hline(yintercept = 10/3, linetype = 2, color = "azure4", lwd = 0.1) +
  theme_classic()

ggsave("totalsym.pdf", plot = totalsym, width = 20, height =10, units = "cm", dpi = 300)


# sde_out |>
#   mutate(totalsymptom = rowSums(pick(S_anh:S_sui))) |>
#   ggplot(aes (x = t, y = totalsymptom)) +
#   geom_line(col = "salmon") +
#   geom_area(data = shock_period, aes(x = time, y = shock*max(sde_out[,-1])*ncol(sde_out[,-1])
#   ), inherit.aes = FALSE, fill = "orange", alpha = 0.2) +
#   labs(y = "", title = "Total symptom level") +
#   theme_classic()



## fitting statistical networks:
burnout <- 30
beforeshock <- aggregated |> filter(t > burnout & t < t_shock) |> select(S_anh:S_sui)
duringshock <- aggregated |> filter(t >= t_shock & t <= t_shock + shock_duration) |> select(S_anh:S_sui) 
aftershock <- aggregated |> filter(t >=t_shock + shock_duration) |> select(S_anh:S_sui)

library(bootnet)
layout(t(1:3))
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
net_layout[9,2] <- -0.3

# Creating node names
nodenames <- c("anhedonia", "appetite", "concentration", "energy", "guilty", "motor", "sad", "sleep", "suicidal")

plot(beforeshock, layout = net_layout, maximum = max_value, font = 2, title = "(a) Before shock", labels = colnames(A))
plot(duringshock, layout = net_layout, maximum = max_value, font = 2, title = "(b) During shock", labels = colnames(A))
plot(aftershock, layout = net_layout, maximum = max_value, font = 2, title = "(c) After shock", labels = colnames(A))


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





layout(1)
## fake networks for Fig1
fakeA <- matrix(c(0,0,0,0,0,
                  0,0,0,0.1,0,
                  0,0,0,0,0,
                  0,0.3,0,0,0,
                  0.3,0,0,0,0), 5, 5, byrow = T)

fakeA[lower.tri(fakeA)] = t(fakeA)[lower.tri(fakeA)]
rownames(fakeA) <- colnames(fakeA) <- nl

fakeB <-matrix(c(0,1,1,1,0,
                 1,0,1,0.1,1,
                 1,0,0,0,1,
                 1,1,1,0,1,
                 0,1,0,1,0), 5, 5, byrow = T)
rownames(fakeB) <- colnames(fakeB) <- nl
fakeB[lower.tri(fakeB)] = t(fakeB)[lower.tri(fakeB)]

nl <- c("S1", "S2", "S3", "S4", "S5")
qgraph(fakeA, theme = 'colorblind', maximum = 0.3)
qgraph(fakeB, theme = 'colorblind')

png(file="healthynetwork.png", units="in", width=5, height=5, res=300, bg = 'transparent')
qgraph(fakeA, theme = 'colorblind', maximum = 0.3)
dev.off()


png(file="sicknetwork.png", units="in", width=5, height=5, res=300, bg = 'transparent')
qgraph(fakeB, theme = 'colorblind')
dev.off()


## Directed network example
Names <- c("anhedonia", "sadness", "sleep", "energy", "appetite", "guilty", "concentration", "motor", "suicidal")

grp <- list(`cyclic` = 1:6, `no cycle` = 7:9)


manual_layout <- matrix(c( -0.3463598, -1.0000000,
             -0.4144493, -0.2077197,
             0.3317425 ,-0.2414516,
             0.4811538 ,-0.9022438,
             1.0000000 ,-0.4836201,
             -0.2013931,  0.2878382,
             0.6342513 , 0.3000000,
             0.1718597,  0.9427039,
             -1.000000,  0.2813530), 9, 2, byrow=T)


pdf(file = "toymodel.pdf", width=16, height=12, bg = 'transparent')

qgraph(A, theme = 'colorblind',legend = TRUE, groups = grp, color = c("#F5651C", "#58B5BC"), nodeNames = Names,      border.color = "white",border.width = 2, edge.color = "darkgray", edge.width = 0.8, curve = 0.8, curveAll = T, label.color = "white", legend.cex = 0.9, layout=manual_layout, asize= 5)

dev.off()

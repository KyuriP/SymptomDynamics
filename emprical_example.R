
# install packages
library(modelr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(qgraph)
library(psychonetrics)

source("euler_stochastic2.R")

f <- function(x) x^2

## adjacency matrix
A <- matrix(c( 1, 1.07, 0, 0, 0, 0, 0, 0, 0,
               1.07, 1, 1.5, 0, 0, 1.35, 0, 0, 1.13,
               0, 0, 1, 1.35, 1.28, 0, 0, 0, 0,
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
Beta_healthy <- c(-7.68, -7.84, -6.88, -8.74, -5.56, -5.70, -5.60, -9.72, -10.90) |> set_names(paste0("Beta_", colnames(A))) 
Beta_bistable <- c(-5.68, -5.84, -5.88, -5.74, -5.56, -4.70, -4.60, -7.72, -8.90) |> set_names(paste0("Beta_", colnames(A)))
# Beta_sick <- c(5.18, 5.34, 5.38, 6.24, 3.06, 8.50, -3.10, -8.22, -3.90) |> set_names(paste0("Beta", 1:9))
Beta_sick <- c(-2.68, 2.84, -1.88, -2.74, -2.56, 1.70, -1.60, -3.72, -5.90) |> set_names(paste0("Beta_", colnames(A)))
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
n_steps <- as.integer(time / deltaT) # must be a number greater than 1

# specify the standard deviation of the stochastic noise
D_stoeq1 <- 0.01 # before shock
D_stoeq2 <- 0.02 # after shock
t_shock <- 300
shock_duration <- 100

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
  duration = shock_duration
)

shock_period <- data.frame(time = 0:timelength) |>
  mutate(shock = ifelse(time >= t_shock & time <= t_shock + shock_duration, TRUE, FALSE))

sde_out |>
  tidyr::pivot_longer(!t, names_to = "symptoms") |>
  ggplot() +
  geom_line(aes(x = t, y = value, color = symptoms), alpha = 0.5) +
  geom_area(data = shock_period, aes(x = time, y = shock*max(sde_out[,-1])
  ), inherit.aes = FALSE, fill = "orange", alpha = 0.2) +
  labs(y = "", title ="Each symptom level") +
  theme_classic() +
  theme(legend.position="bottom") 

sde_out |>
  mutate(totalsymptom = rowSums(pick(S_anh:S_sui))) |>
  ggplot(aes (x = t, y = totalsymptom)) +
  geom_line(col = "salmon") +
  geom_area(data = shock_period, aes(x = time, y = shock*max(sde_out[,-1])*ncol(sde_out[,-1])
  ), inherit.aes = FALSE, fill = "orange", alpha = 0.2) +
  labs(y = "", title = "Total symptom level") +
  theme_classic()


## fitting statistical networks:
burnout <- 30
beforeshock <- sde_out |> filter(t > burnout & t < t_shock) |> select(!t)
aftershock <- sde_out |> filter(t >=t_shock) |> select(!t)

layout(t(1:2))
beforeGGM <- qgraph(cor_auto(beforeshock),
                    graph = 'glasso',
                    layout = 'spring',
                    theme = 'colorblind',
                    sampleSize = nrow(beforeshock),
                    title = "beforeGGM")
# fix layout
L <- beforeGGM$layout

afterGGM <- qgraph(cor_auto(aftershock),
                   graph = 'glasso',
                   theme = 'colorblind',
                   sampleSize = nrow(aftershock),
                   title = "afterGGM",
                   layout = L)



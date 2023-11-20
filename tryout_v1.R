# install package
library(modelr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)


f <- function(x) {
  if(x >= 0.5 | is.nan(x)) 1 else 0
}

# define diff equation:
det_eq <- c(
  dSA ~ SA*(1-SA)*(bA - aA*SA + SC*(1+ dA*f(SA))),
  dSB ~ SB*(1-SB)*(bB - aB*SB + SA*(1+ dB*f(SB))),
  dSC ~ SC*(1-SC)*(bC - aC*SC + SB*(1+ dC*f(SC)))
)

sto_eq <-  c(dSA ~ .1,
             dSB ~ .1,
             dSC ~ .1)

## HEALTHY STATE ##
# define the parameters (as a named vector):
parms1 <- c(bA = -.6, bB = -.5, bC = -.6, 
           aA = .3, aB = .2, aC = .4, 
           dA = .2, dB = .5, dC = .3)

parms2 <- c(bA = .6, bB = .5, bC = .7, 
            aA = .1, aB = .01, aC = .2, 
            dA = .2, dB = .5, dC = .3)
  
# define the initial condition (as a named vector):
init <- c(SA = 0.3, SB = 0.5, SC= 0.8)

# define deltaT and the number of time steps:
deltaT <- .1 # timestep length
n_steps <- 2300 # must be a number greater than 1

# Identify the standard deviation of the stochastic noise
D_stoeq1 <- .3
D_stoeq2 <- .6

set.seed(45678)
# Do one simulation of this differential equation
sde_out <- euler_stochastic2(
  deterministic_rate = det_eq,
  stochastic_rate = sto_eq,
  initial_condition = init,
  parameters1 = parms1,
  deltaT = deltaT,
  n_steps = n_steps,
  D1 = D_stoeq1,
  D2 = D_stoeq2,
  shock = TRUE,
  parameters2 = parms2, 
  t_shock = 300, 
  duration = 40
)

p1 <- sde_out |>
  mutate(totalsymptom = SA + SB + SC) |>
  tidyr::pivot_longer(!c(t, totalsymptom), names_to = "symptoms") |>
  ggplot() +
  geom_line(aes(x = t, y = value, color = symptoms), alpha = 0.5) +
  # geom_line(aes(x = t, y = totalsymptom)) +
  labs(y = "", title ="Each symptom level") +
  theme_classic() +
  theme(legend.position="bottom")
  

p2 <- sde_out |>
  mutate(totalsymptom = SA + SB + SC) |>
  ggplot(aes (x = t, y = totalsymptom)) +
  geom_line(col = "salmon") +
  labs(y = "", title = "Total symptom level") +
  theme_classic()

# patchwork
p1+p2




## WITHOUT NOISE
# # compute the solution via Runge-Kutta method:
# out_solution <- rk4(system_eq = deq,
#                     parameters = parms,
#                     initial_condition = init,
#                     deltaT = deltaT,
#                     n_steps = n_steps
# )
# 
# # make a plot of the solution:
# out_solution |>
#   tidyr::pivot_longer(!t, names_to = "symptoms") |>
#   ggplot(aes(x = t, y = value, color = symptoms)) +
#   geom_line() 


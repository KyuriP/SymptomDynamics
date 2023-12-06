# install packages
library(modelr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(qgraph)
library(psychonetrics)


f <- function(x) {
  if(x >= 0.5 | is.nan(x)) 1 else 0  
}  # linear function maybe or power of S_i.



# define diff equation:
det_eq <- c(
  dSA ~ SA*(1-SA)*(bA + aA*SA + SC*(1+ dA*f(SA))),
  dSB ~ SB*(1-SB)*(bB + aB*SB + SA*(1+ dB*f(SB))),
  dSC ~ SC*(1-SC)*(bC + aC*SC + SB*(1+ dC*f(SC)))
)

det_eq <- c(
  dSA ~ SA*(1-SA)*(bA + aA*SA + SC*(1 + dA*(SA)^5)),
  dSB ~ SB*(1-SB)*(bB + aB*SB + SA*(1 + dB*(SB)^5)),
  dSC ~ SC*(1-SC)*(bC + aC*SC + SB*(1 + dC*(SC)^5))
)

sto_eq <-  c(dSA ~ .1,
             dSB ~ .1,
             dSC ~ .1)

# define the parameters (as a named vector):
## original params ##
parms1 <- c(bA = -.6, bB = -.3, bC = -.7, 
           aA = -.3, aB = -.2, aC = -.4, 
           dA = .3, dB = .3, dC = .3)

## given shock ##
parms2 <- c(bA = .6, bB = .3, bC = .7, 
            # aA = 0, aB = 0, aC = , 0
            aA = -.3, aB = -.2, aC = -.4, 
            dA = .3, dB = .3, dC = .3)
  
# define the initial condition (as a named vector):
init <- c(SA = 0.3, SB = 0.5, SC= 0.8)

# define deltaT and the number of time steps:
deltaT <- .1 # timestep length
n_steps <- 1500 # must be a number greater than 1

# specify the standard deviation of the stochastic noise
D_stoeq1 <- 0.1 # before shock
D_stoeq2 <- 0.5 # after shock

# set.seed(45678) # set the seed (system behavior varies much by random noise)
# simulate :
# > shock given at t = 300
# > shock --> change b and noise S.D.
# > shock duration = 40 
# > after shock --> revert back to original parameters but noise S.D. remains larger (debating... reasonable?)

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
  t_shock = 500, 
  duration = 50
)


sde_out |>
  mutate(totalsymptom = SA + SB + SC) |>
  tidyr::pivot_longer(!c(t, totalsymptom), names_to = "symptoms") |>
  ggplot() +
  geom_line(aes(x = t, y = value, color = symptoms), alpha = 0.5) +
  # geom_line(aes(x = t, y = totalsymptom)) +
  labs(y = "", title ="Each symptom level") +
  theme_classic() +
  theme(legend.position="bottom")



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
p1 + p2


## fitting network
beforeshock <- sde_out |> filter(t > 25 & t < 49.9) |> select(!t)
aftershock <- sde_out |> filter(t >=49.9) |> select(!t)

layout(matrix(1:4,2,2))
beforeGGM <- qgraph(cor_auto(beforeshock),
       graph = 'glasso',
       layout = 'spring',
       theme = 'colorblind',
       sampleSize = nrow(beforeshock),
       title = "beforeGGM")

afterGGM <- qgraph(cor_auto(aftershock),
       graph = 'glasso',
       layout = 'spring',
       theme = 'colorblind',
       sampleSize = nrow(aftershock),
       title = "afterGGM")

before <- gvar(beforeshock, estimator="FIML") %>% runmodel
beforeGVAR <- getmatrix(before, "PDC") |> 
  qgraph(theme = "colorblind", directed=TRUE, diag=TRUE,
       title = "beforeGVAR")

after <- gvar(aftershock, estimator="FIML") %>% runmodel
afterGVAR <- getmatrix(after, "PDC") |>
  qgraph(theme = "colorblind", directed=TRUE, diag=TRUE,
       title = "afterGVAR")



centralityPlot(beforeGVAR,scale="raw0")
centralityPlot(afterGVAR,scale="raw0")

## Questions:
# - how to show the hysteresis 
# - how to show some symptom comes down while some doesn't (differentiate behavior)
# - the parameter value range conditions -- in order to keep bistabiliy .. ?





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

x <- seq(0, 1, 0.01)
y = x^(3)
y1 = x
plot(y1)
plot(y)
lines(y1)
y1 > y

## ======================
## simulate with A_hel 
## ======================

## install packages
source("code/libraries.R")

# rm(list=ls())

## source necessary functions
source("code/euler_stochastic2.R")

# change suicidal alpha to be the same
A_hel <- matrix(c( .10, 0, 0, 0, .12, 0, 0, 0, 0,
                   .39, .10, .10, 0, 0, .13, 0, .12, .14,
                   0, 0, .10, .3, .12, 0, 0, 0, 0,
                   .23, 0, 0, .10, 0, 0, .12, 0, 0,
                   0, 0, 0, .20, .10, 0, 0, 0, 0,
                   0, .13, 0, 0, .12, .10, .14, .13, .21,
                   0, 0, 0, 0, 0, 0, .10, 0.25, 0,
                   0, 0, 0, 0, 0, 0, 0, .10, 0,
                   0, 0, 0, 0, 0, 0, 0, .12, .10), 9, 9, byrow = T)
rownames(A_hel) <- colnames(A_hel) <- c("anh", "sad", "slp", "ene", "app", "glt", "con", "mot", "sui")
diag(A_hel) <- 0.03
# define differential equations:

dif_eq <- 1:9 |> map(function(x){
  lhs <- paste0("dS", paste0("_",colnames(A_hel)[x]))
  sumAj <- paste0("S", paste0("_",colnames(A_hel)[-x]), collapse = ",") 
  rhs <- stringr::str_replace_all(
    "Sk * (1-Sk) * (Betak + A_hel[q,q] * Sk + (1+ deltak * f(Sk))) * A_hel[-q,q] %*% ", c(k = paste0("_",colnames(A)[x]), q = x)) |> paste0(paste0("c(", sumAj, ")"))
  form <- as.formula(paste(lhs,"~", rhs))
})


# define "f"
f <- function(x) x^2

# give noise:
sto_eq <-  1:9 |> map(function(x){
  lhs <- paste0("dS", paste0("_",colnames(A_hel)[x]))
  rhs <- 1  # change as you need
  form <- as.formula(paste(lhs,"~", rhs))
})

# sto_eq <-  c(dSA ~ .1,
#              dSB ~ .1,
#              dSC ~ .1)

# define the parameters (as a named vector):
A_i <- diag(A_hel) |> set_names(paste0("A_i_", colnames(A_hel))) 
A_j <- colSums(A_hel) |> set_names(paste0("A_j_", colnames(A_hel))) 
delta <- rep(1, 9) |> set_names(paste0("delta_", colnames(A_hel))) 

dist <- (- A_j - (- A_i - A_j*(1+delta))) / 100
- A_j - dist*110
## bistable: ...< beta_i < - A_j # -omega * (1 + 0.41)
- A_j*(1+1.8)

Beta_bistable <- c(-1.6831, -0.4924, -0.4195, -1.3915, -1.0513, -0.4924, -0.8083, -1.6831, -1.0270 ) |> set_names(paste0("Beta_", colnames(A_hel))) 

## sick: beta_i > - A_j
# Beta_sick <- c(1.68, 2.04, 1.88, -2.74, -2.56, 1.95, -1.60, -3.72, -5.90) |> set_names(paste0("Beta_", colnames(A))
# Beta_sick <- c(-0.72, -0.47, -0.41, -0.59, -0.53, -0.18, -0.37, -0.67, -0.41) |> set_names(paste0("Beta_", colnames(A)))
# Beta_sick <- c(-0.90,   -0.65,   -0.59,   -0.77,   -0.61,   -0.36,   -0.55,   -0.85 , -1.25) |> set_names(paste0("Beta_", colnames(A)))
Beta_sick <- c(-0.46800, -0.21075, -0.19500, -0.40500, -0.33150, -0.21075, -0.27900, -0.46800, -0.32625 ) |> set_names(paste0("Beta_", colnames(A_hel))) #- A_j + dist*5

# Beta <- Beta_sick |> set_names(paste0("Beta", 1:9))

## original params
parms1 <- c(Beta_bistable, delta) 
## given shock (beta: stressor increases)  
parms2 <- c(Beta_sick, delta) 

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
timelength <- 3000
n_steps <- as.integer(timelength / deltaT) # must be a number greater than 1

# specify the standard deviation of the stochastic noise
D_stoeq1 <- 0.005  # before shock
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
  Amat = A_hel,
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
  duration = shock_duration
  #seed = 123
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
  scale_color_discrete(name = "Symptoms", labels = c("anhedonia", "sad", "sleep", "energy", "appetite", "guilty", "concentration", "motor", "suicidal")) +
  #labs(x="time", y = "", title ="Trajectory of each symptom") +
  labs(x="time", y = "") +
  lims(y=c(0,1)) +
  theme_classic(base_size = 15) +
  theme(legend.position="none") +
  facet_wrap(~symptoms,
             labeller = labeller(symptoms = labelllername) 
  )

# ggsave("eachsym_v2.pdf", plot = eachsym, width = 30, height =15, units = "cm", dpi = 300)


## aggregate symptom level
n_sims <- 100

# run sims n_sims times
aggregated <- map(1:n_sims, ~ euler_stochastic2(
  Amat = A_hel,
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
  duration = shock_duration) |>
    mutate(totalsymptom = rowSums(pick(S_anh:S_sui)))
) |> list_rbind(names_to = "sim")

# saveRDS(aggregated, "simdat_Ahel.rds")


# .libPaths(c('/gpfs/home1/kpark/rpackages',.libPaths()))

library(modelr)
library(dplyr)
library(purrr)
library(tidyr)
#library(ggplot2)
library(qgraph)
library(bootnet, lib='/gpfs/home1/kpark/rpackages')library(deSolve)
library(magrittr)
library(stringr)

## source necessary functions
euler_stochastic2 <- function (deterministic_rate, stochastic_rate, initial_condition, 
                               Amat, parameters1 = NULL, t_start = 0, deltaT = 1, timelength = 1, 
                               D1 = 1, shock = F, parameters2 = NULL, t_shock = NULL, duration = NULL, seed = NULL) 
{
  set.seed(seed)
  curr_vec <- c(initial_condition, t = t_start)
  vec_names <- names(curr_vec)
  n_vars <- length(vec_names)
  time_eq <- c(dt ~ 1)
  new_rate_eq <- c(deterministic_rate, time_eq) %>% formula.tools::rhs()
  time_eq_stoc <- c(dt ~ 0)
  new_stochastic_rate <- c(stochastic_rate, time_eq_stoc) %>% 
    formula.tools::rhs()
  n_steps = as.integer(timelength / deltaT)
  out_list <- vector("list", length = n_steps)
  out_list[[1]] <- curr_vec
  for (i in 2:n_steps) {
    if(shock){
      if(i < t_shock / deltaT | i > (t_shock + duration)/deltaT){
        in_list <- c(parameters1, curr_vec) %>% as.list() %>% append(list(A = Amat))
      } else {
        in_list <- c(parameters2, curr_vec) %>% as.list() %>% append(list(A = Amat))
      }
    } else {
      in_list <- c(parameters1, curr_vec) %>% as.list()
    }
    curr_rate <- sapply(new_rate_eq, FUN = eval, envir = in_list) %>% 
      purrr::set_names(nm = vec_names)
    curr_stoch_rate <- sapply(new_stochastic_rate, FUN = eval, 
                              envir = in_list) %>% purrr::set_names(nm = vec_names)
    
    v3 <- c(curr_vec, curr_rate * deltaT, curr_stoch_rate *
              sqrt(2 * D1 * deltaT) * rnorm(n_vars))
    
    curr_vec <- tapply(v3, names(v3), sum) 
    # correction for violating boundaries[0,1]
    curr_vec[-length(curr_vec)] %<>% ifelse(. > 1, abs(. - 2), .) # this would not be over 2 with optimal noise magnitutde
    curr_vec[-length(curr_vec)] %<>% ifelse(. < 0, abs(.), .) # this would not be lower than -2 with optimal noise magnitutde
    out_list[[i]] <- curr_vec
  }
  out_results <- out_list %>% dplyr::bind_rows() %>% dplyr::relocate(t)
  return(out_results)
}



# compute density before shock
# sum(weightMatrix_lowresil) - sum(weightMatrix_highresil)

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

## define "f"
f <- function(x) x^2

# 1st critical pt: -colSums(A)*(1+delta)-a_ii; given baseline scenario, -colSums(A)*(1+9)-0.3 
# 2nd critical pt: -colSums(A)
criticalpt1 <- -colSums(A)*(1+9) - 0.3
criticalpt2 <- -colSums(A)
betas <- 1:9 |> purrr::map(function(x) {seq(criticalpt1[x], criticalpt2[x], length.out=10)}) |>
  purrr::set_names(names(criticalpt1)) |> as.data.frame() |> purrr::set_names(paste0("Beta_", colnames(A))) 
sigmas <- seq(0.005, 0.02, length.out = 10) # 0.1 = sqrt(2*0.005) - 0.2 = sqrt(2*0.02)



# Function computes the density of matrix
compute_density <- function(choice = "base", Beta = NULL, Sigma, n_sim = 50){
  
  # A matrix
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
  
  ## define "f"
  f <- function(x) x^2
  
  dif_eq <- 1:9 |> purrr::map(function(x){
    lhs <- paste0("dS", paste0("_",colnames(A)[x]))
    sumAj <- paste0("S", paste0("_",colnames(A)[-x]), collapse = ",") 
    rhs <- stringr::str_replace_all(
      "Sk * (1-Sk) * (Betak + (A[q,q] * Sk) + ((1+ deltak * f(Sk)) * A[-q,q] %*% ", c(k = paste0("_",colnames(A)[x]), q = x)) |> paste0(paste0("c(", sumAj, ")))"))
    form <- as.formula(paste(lhs,"~", rhs))
  })
  
  # define stochastic part:
  sto_eq <-  1:9 |> purrr::map(function(x){
    lhs <- paste0("dS", paste0("_",colnames(A)[x]))
    rhs <- 1  # change as you need
    form <- as.formula(paste(lhs,"~", rhs))
  })
  
  # beta_sick is fixed
  Beta_sick <- c(-0.91, -0.66, -0.60, -0.78, -0.62, -0.37, -0.56, -0.86, -1.30) |> purrr::set_names(paste0("Beta_", colnames(A))) 
  
  ## delta 
  delta <- rep(9, 9) |> purrr::set_names(paste0("delta_", colnames(A)))
  
  if(choice == "high") {
    delta <- rep(8, 9) |> purrr::set_names(paste0("delta_", colnames(A))) 
    diag(A) <- 0.25
  } else if (choice == "low"){
    delta <- rep(10, 9) |> purrr::set_names(paste0("delta_", colnames(A))) 
    diag(A) <- 0.35
  }
  
  ## original params
  if (is.null(Beta)){
    Beta_bistable <- c(-1.373, -1.019, -0.934, -1.189, -0.962, -0.608, -0.877, -1.302, -2.30) |> purrr::set_names(paste0("Beta_", colnames(A))) 
    parms1 <- c(Beta_bistable, delta)
  } else {
    parms1 <- c(purrr::as_vector(Beta), delta)}
  ## params given shock (beta: increases)  
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
  
  ## define dt and the number of time steps:
  deltaT <- .1 # dt 
  timelength <- 4000 # length of simulation
  n_steps <- as.integer(timelength / deltaT) # must be a number greater than 1
  
  ## specify the magnitude of noise and specifics of shock 
  D_stoeq1 <- Sigma  # before shock
  t_shock <- 1000 # time that shock begins
  shock_duration <- 500 # shock duration time points
  n_sims = n_sim
  aggregated <- purrr::map(1:n_sims, ~ euler_stochastic2(
    deterministic_rate = dif_eq,
    stochastic_rate = sto_eq,
    initial_condition = init,
    parameters1 = parms1,
    parameters2 = parms2, 
    Amat = A,
    deltaT = deltaT,
    timelength = timelength,
    D1 = D_stoeq1,
    shock = TRUE,
    t_shock = t_shock, 
    duration = shock_duration)
  ) |> purrr::list_rbind(names_to = "sim")
  
  # before shock
  beforeshock <- aggregated |> 
    filter(t < t_shock) |> 
    select(S_anh:S_sui) |> 
    bootnet::estimateNetwork(default = "EBICglasso") |> 
    suppressMessages() |> suppressWarnings() |>
    {\(x) x$graph}() |> sum()
  # during shock
  duringshock <- aggregated |> 
    filter(t >= t_shock & t <= t_shock + shock_duration) |> 
    select(S_anh:S_sui) |>
    bootnet::estimateNetwork(default = "EBICglasso")|> 
    suppressMessages() |> suppressWarnings() |>
    {\(x) x$graph}() |> sum()
  # after shock
  aftershock <- aggregated |> 
    filter(t >=t_shock + shock_duration) |> 
    select(S_anh:S_sui) |>
    bootnet::estimateNetwork(default = "EBICglasso")|> 
    suppressMessages() |> suppressWarnings() |>
    {\(x) x$graph}() |> sum()
  
  return(list(before = beforeshock
              , during = duringshock, after=aftershock))
}


## results over 30nsims
res_highs <- purrr::map(1:nrow(betas), function(x) purrr::map_dfr(1:length(sigmas), function(y) compute_density(choice = "high", Beta = betas[x,], Sigma = sigmas[y]),.id = "id") |> dplyr::mutate(sigma = sigmas[as.numeric(id)]))
res_lows <- purrr::map(1:nrow(betas), function(x) purrr::map_dfr(1:length(sigmas), function(y) compute_density(choice = "low", Beta = betas[x,], Sigma = sigmas[y]),.id = "id") |> dplyr::mutate(sigma = sigmas[as.numeric(id)]))

saveRDS(res_highs, "res_high100.rds")
saveRDS(res_lows, "res_low100.rds")
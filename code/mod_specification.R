## =========================================================
## Model specification
##
## specify parameter values given the chosen scenario
## define differential equations
## 
## return: list of relevant parameters and equations
## =========================================================

mod_spec <- function(scenario = "base", init_val = 0, ...){
  
  ## Adjacency matrices
  # A <- matrix(c( .30, 0, 0, 0, 0, 0, 0, 0, 0,
  #                .33, .30, .14, .15, 0, .13, 0, 0, .15,
  #                .13, .14, .30, .22, .23, 0, 0, 0, 0,
  #                .21, .15, .22, .30, 0, 0, .12, 0, 0,
  #                0, 0, 0, .17, .30, 0, 0, 0, 0,
  #                0, .13, 0, 0, .15, .30, .2, .15, .22,
  #                0, 0, 0, 0, 0, 0, .30, .17, 0,
  #                0, 0, 0, 0, 0, 0, 0, .30, 0,
  #                0, 0, 0, 0, 0, 0, 0, .3, 0.07), 9, 9, byrow = T)
  
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

    ## beta
  # Beta_bistable <- c(-1.455, -1.080, -0.990, -1.260, -1.020, -0.645, -0.930, -1.380, -2.50) |> set_names(paste0("Beta_", colnames(A))) 
  # Beta_bistable <- c(-1.74, -1.24, -1.12, -1.48, -1.16, -0.66, -1.04, -1.64, -2.50) |> set_names(paste0("Beta_", colnames(A)))
  Beta_bistable <- c(-1.358, -1.008, -0.924, -1.176, -0.952, -0.602, -0.868, -1.288, -2.50) |> set_names(paste0("Beta_", colnames(A))) # -omega * (1 + 0.4)

  Beta_sick <- c(-0.91, -0.66, -0.60, -0.78, -0.62, -0.37, -0.56, -0.86, -1.50) |> set_names(paste0("Beta_", colnames(A))) #- A_j + dist*5
  
  ## delta 
  delta <- rep(9, 9) |> set_names(paste0("delta_", colnames(A)))
  
  ## high resilience case
  if (scenario == "high") {
    diag(A) <- 0.1
    delta <- rep(7, 9) |> set_names(paste0("delta_", colnames(A)))
  } 
  
  ## low resilience case
  if (scenario == "low") {
    diag(A) <- 0.5
    delta <- rep(10, 9) |> set_names(paste0("delta_", colnames(A)))
  } 
  
  
  # define deterministic part:
  # "dS[i] ~ S[i]*(1-S[i])*(Beta[i] + A_i[i]*S[i] + SUM(A_[j][i]*S_[j])*(1+ delta[i]*f(S[i])))"
  
  dif_eq <- 1:9 |> map(function(x){
    lhs <- paste0("dS", paste0("_",colnames(A)[x]))
    sumAj <- paste0("S", paste0("_",colnames(A)[-x]), collapse = ",") 
    rhs <- stringr::str_replace_all(
      "Sk * (1-Sk) * (Betak + (A[q,q] * Sk) + ((1+ deltak * f(Sk)) * A[-q,q] %*% ", c(k = paste0("_",colnames(A)[x]), q = x)) |> paste0(paste0("c(", sumAj, ")))"))
    form <- as.formula(paste(lhs,"~", rhs))
  })

  
  # define stochastic part:
  sto_eq <-  1:9 |> map(function(x){
    lhs <- paste0("dS", paste0("_",colnames(A)[x]))
    rhs <- 1  # change as you need
    form <- as.formula(paste(lhs,"~", rhs))
  })
  
  init <- c(S_anh = init_val, 
            S_sad = init_val, 
            S_slp = init_val, 
            S_ene = init_val, 
            S_app = init_val, 
            S_glt = init_val, 
            S_con = init_val, 
            S_mot = init_val, 
            S_sui = init_val)
  
  return(model = list(delta = delta, Beta_bistable = Beta_bistable, Beta_sick = Beta_sick, A = A, initial_values = init, dif_eq = dif_eq, sto_eq = sto_eq))
}

res <- mod_spec(scenario = "low")


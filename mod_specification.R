## ===============================
## Model specification
## specify parameter vals depending on the scenario
## define differential equations
## 
## ===============================

mod_spec <- function(scenario = "base"){
  
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
  Beta_bistable <- c(-1.355, -0.995, -0.905, -1.165, -0.935, -0.575, -0.845, -1.285, -2.00) |> set_names(paste0("Beta_", colnames(A))) 
  Beta_sick <- c(-0.91, -0.66, -0.60, -0.78, -0.62, -0.37, -0.56, -0.86, -1.30) |> set_names(paste0("Beta_", colnames(A))) #- A_j + dist*5
  
  ## delta 
  delta <- rep(9.5, 9) |> set_names(paste0("delta_", colnames(A)))
  
  ## high resislience case
  if (scenario == "high") {
    diag(A) <- 0.1
    delta <- rep(8, 9) |> set_names(paste0("delta_", colnames(A)))
  } 
  
  ## low resilience case
  if (scenario == "low") {
    diag(A) <- 0.5
    delta <- rep(11, 9) |> set_names(paste0("delta_", colnames(A)))
  } 
  
  ## define "f"
  f <- function(x) x^2
  
  # define differential equations:
  #equation <- "dS[i] ~ S[i]*(1-S[i])*(Beta[i] + A_i[i]*S[i] + A_[j][i]*S_[j]*(1+ delta[i]*f(S[i])))"
  
  dif_eq <- 1:9 |> map(function(x){
    lhs <- paste0("dS", paste0("_",colnames(A)[x]))
    sumAj <- paste0("dS", paste0("_",colnames(A)[-x]), collapse = ",") 
    rhs <- stringr::str_replace_all(
      "Sk * (1-Sk) * (Betak + A[q,q] * Sk + (1+ deltak * f(Sk))) * A[-q,q] %*% ", c(k = paste0("_",colnames(A)[x]), q = x)) |> paste0(paste0("c(", sumAj, ")"))
    form <- as.formula(paste(lhs,"~", rhs))
  })
  
  # dif_eq <- c(dS_anh ~ S_anh * (1 -  S_anh) * (Beta_anh + A[1,1] *  S_anh + 
  #                                                A[-1,1] %*% c(S_sad, S_slp, S_ene, S_app, S_glt, S_con, S_mot, S_sui)*(1 + delta_anh * f(S_anh))),
  #             dS_sad ~ S_sad * (1 -  S_sad) * (Beta_sad + A[2,2] *  S_sad + 
  #                                                A[-2,2] %*% c(S_anh, S_slp, S_ene, S_app, S_glt, S_con, S_mot, S_sui)*(1 + delta_sad * f(S_sad))),
  #             dS_slp ~ S_slp * (1 -  S_slp) * (Beta_slp + A[3,3] *  S_slp + 
  #                                                A[-3,3] %*% c(S_anh, S_app, S_ene, S_app, S_glt, S_con, S_mot, S_sui)*(1 + delta_slp * f(S_slp))),
  #             dS_ene ~ S_ene * (1 -  S_ene) * (Beta_ene + A[4,4] *  S_ene + 
  #                                                A[-4,4] %*% c(S_anh, S_app, S_slp, S_app, S_glt, S_con, S_mot, S_sui)*(1 + delta_ene * f(S_ene))),
  #             dS_app ~ S_app * (1 -  S_app) * (Beta_app + A[5,5] *  S_app + 
  #                                                A[-5,5] %*% c(S_anh, S_sad, S_slp, S_ene, S_glt, S_con, S_mot, S_sui)*(1 + delta_app * f(S_app))),
  #             dS_glt ~ S_glt * (1 - S_glt) * (Beta_glt + A[6,6] * S_glt + 
  #                                               A[-6,6] %*% c(S_anh, S_sad, S_slp, S_ene, S_app, S_con, S_mot, S_sui)*(1 + delta_glt * f(S_glt))),
  #             dS_con ~ S_con * (1 - S_con) * (Beta_con + A[7,7] *  S_con + 
  #                                               A[-7,7] %*% c(S_anh, S_sad, S_slp, S_ene, S_app, S_glt, S_mot, S_sui)*(1 + delta_con * f(S_con))),
  #             dS_mot ~ S_mot * (1 -  S_mot) * (Beta_mot + A[8,8] *  S_mot + 
  #                                                A[-8,8] %*% c(S_anh, S_sad, S_slp, S_ene, S_app, S_glt, S_con, S_sui)*(1 + delta_mot * f(S_mot))),
  #             dS_sui ~ S_sui * (1 -  S_sui) * (Beta_sui + A[9,9] *  S_sui + 
  #                                                A[-9,9] %*% c(S_anh, S_sad, S_slp, S_ene, S_app, S_glt, S_con, S_mot)*(1 + delta_sui * f(S_sui)))
  #             
  # )
  
  # give noise:
  sto_eq <-  1:9 |> map(function(x){
    lhs <- paste0("dS", paste0("_",colnames(A)[x]))
    rhs <- 1  # change as you need
    form <- as.formula(paste(lhs,"~", rhs))
  })
  
  # sto_eq <-  c(dSA ~ .1,
  #              dSB ~ .1,
  #              dSC ~ .1)
  
  return(model = list(A = A, dif_eq = dif_eq))
}

res <- mod_spec(scenario = "high")


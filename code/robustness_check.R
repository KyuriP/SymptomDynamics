
## install packages
source("code/libraries.R")


# compute density before shock
# sum(weightMatrix_lowresil) - sum(weightMatrix_highresil)

# 1st critical pt: -colSums(A)*(1+delta)-a_ii; given baseline scenario, -colSums(A)*(1+9)-0.3 
# 2nd critical pt: -colSums(A)
criticalpt1 <- -colSums(A)*(1+9) - 0.3
criticalpt2 <- -colSums(A)
betas <- 1:9 |> map(function(x) {seq(criticalpt1[x], criticalpt2[x], length.out=10)}) |>
  set_names(names(criticalpt1)) |> as.data.frame() |> set_names(paste0("Beta_", colnames(A))) 

sigmas <- seq(0.0707, 0.1414, length.out = 10) # 0.1/sqrt(2) - 0.2/sqrt(2)


# Function computes the density of matrix
compute_density <- function(choice = "base", Beta = NULL, Sigma, n_sim = 10){
  
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
  
  # beta_sick is fixed
  Beta_sick <- c(-0.91, -0.66, -0.60, -0.78, -0.62, -0.37, -0.56, -0.86, -1.30) |> set_names(paste0("Beta_", colnames(A))) 
  
  ## delta 
  delta <- rep(9, 9) |> set_names(paste0("delta_", colnames(A)))
  
  if(choice == "high") {
    delta <- rep(8, 9) |> set_names(paste0("delta_", colnames(A))) 
    diag(A) <- 0.25
  } else if (choice == "low"){
    delta <- rep(10, 9) |> set_names(paste0("delta_", colnames(A))) 
    diag(A) <- 0.35
  }
  
  ## original params
  if (is.null(Beta)){
    Beta_bistable <- c(-1.373, -1.019, -0.934, -1.189, -0.962, -0.608, -0.877, -1.302, -2.30) |> set_names(paste0("Beta_", colnames(A))) 
    parms1 <- c(Beta_bistable, delta)
  } else {
    parms1 <- c(as_vector(Beta), delta)}
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
  aggregated <- map(1:n_sims, ~ euler_stochastic2(
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
  ) |> list_rbind(names_to = "sim")
  
  # before shock
  beforeshock <- aggregated |> 
    filter(t < t_shock) |> 
    select(S_anh:S_sui) |> 
    estimateNetwork(default = "EBICglasso") |> 
    suppressMessages() |> suppressWarnings() |>
    {\(x) x$graph}() |> sum()
  # during shock
  duringshock <- aggregated |> 
    filter(t >= t_shock & t <= t_shock + shock_duration) |> 
    select(S_anh:S_sui) |>
    estimateNetwork(default = "EBICglasso")|> 
    suppressMessages() |> suppressWarnings() |>
    {\(x) x$graph}() |> sum()
  # after shock
  aftershock <- aggregated |> 
    filter(t >=t_shock + shock_duration) |> 
    select(S_anh:S_sui) |>
    estimateNetwork(default = "EBICglasso")|> 
    suppressMessages() |> suppressWarnings() |>
    {\(x) x$graph}() |> sum()
  
  return(list(before = beforeshock
              , during = duringshock, after=aftershock))
}


## results

res_highs <- map(1:nrow(betas), function(x) map_dfr(1:length(sigmas), function(y) compute_density(choice = "high", Beta = betas[x,], Sigma = sigmas[y]),.id = "id") |> mutate(sigma = sigmas[as.numeric(id)]))
res_lows <- map(1:nrow(betas), function(x) map_dfr(1:length(sigmas), function(y) compute_density(choice = "low", Beta = betas[x,], Sigma = sigmas[y]),.id = "id") |> mutate(sigma = sigmas[as.numeric(id)]))


## over 10 beta and 10 sigma
high_resil <- readRDS("data/res_highs.rds")
low_resil <- readRDS("data/res_lows.rds")

high_resil %<>%  map(~as.data.frame(.x) |>
                    tidyr::pivot_longer(everything())|>
                                          mutate(sigma = rep(sigmas, each=3))
                  )
low_resil %<>%  map(~as.data.frame(.x) |>
                    tidyr::pivot_longer(everything())|>
                    mutate(sigma = rep(sigmas, each=3))
)

diff_result <- map2_dfc(low_resil, high_resil, \(x,y) x$value - y$value) |> 
  rename_with(~paste0("beta", 1:nrow(betas))) |>
  mutate(sigma = rep(sigmas, each=3),
         phase = rep(c("before", "during", "after"), 10)) 

## plot diff results

# given sigma value
p_sigma <- diff_result |>
  tidyr::pivot_longer(!c(sigma, phase), names_to = "beta", values_to = "value") |>
  ggplot(aes(x = as.factor(round(sigma,3)), y = value, color = factor(phase, levels=c("before", "during", "after")), fill = factor(phase, levels=c("before", "during", "after")))) +
  geom_boxplot(alpha =0.3,   outlier.size = 1) +
  labs(color = "", fill = "", x = "sigma", title = "network density difference", subtitle = "high vs. low resilience", y = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme_pubr()

# given beta value
p_beta <- diff_result |>
  tidyr::pivot_longer(!c(sigma, phase), names_to = "beta", values_to = "value") |>
  ggplot(aes(x = factor(beta, levels =c(paste0("beta", 1:10))), y = value, color = factor(phase, levels=c("before", "during", "after")), fill = factor(phase, levels=c("before", "during", "after")))) +
  geom_boxplot(alpha =0.3,   outlier.size = 1) +
  labs(color = "", fill = "", x = "beta", title = "network density difference", subtitle = "high vs. low resilience", y = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme_pubr()

## during shock -- heatmap
interactionplot <- diff_result |>
  tidyr::pivot_longer(!c(sigma, phase), names_to = "beta", values_to = "value") |>
  filter(phase == "during") |>
  ggplot(aes(x = factor(beta, levels =c(paste0("beta", 1:10))), y = sigma, fill = value)) +
  # geom_tile() +
  geom_raster()+
  labs(x = "beta", y = "sigma", fill=expression(Delta*density),title = "network density difference", subtitle = "high vs. low resilience (during the shock)" ) +
  theme_classic() +
  scale_fill_continuous(high = "#132B43", low = "#56B1F7")

## =========================================================
## Assessing Model Flexibility
## 
## This script evaluates model flexibility by applying the 
## lowest and highest beta values within the bistability range 
## for each resilience scenario.
## =========================================================
## install packages
source("code/libraries.R")
## source necessary functions
source("code/utils.R")
source("code/euler_stochastic2.R")
source("code/mod_specification.R")


## set the seed
set.seed(123)

## define model specifics: choose the scenario and initial value for symptoms
choice <- "base"
choice <- "high"
choice <- "low"

mod <- mod_spec(scenario = choice, init_val = 0.01)


## define "f"
f <- function(x) x^2

## different beta sets
criticalpt1 <- -colSums(A)*(1+10) - 0.35 ## low resilient 1st critical pt
criticalpt2 <- -colSums(A)
#suicidal manipulate
criticalpt2[9] <- -1.3 - 0.06


Beta_bistable1 <- criticalpt1 |> set_names(paste0("Beta_", nodenames)) 
rownames(A) <- colnames(A) <- nodenames
Beta_bistable2 <- criticalpt2 |> set_names(paste0("Beta_", nodenames)) 
rownames(A) <- colnames(A) <- nodenames

## original params
parms1 <- c(Beta_bistable1, mod$delta)
parms1 <- c(Beta_bistable2, mod$delta)

## params given shock (beta: increases)  
parms2 <- c(mod$Beta_sick, mod$delta)

## define dt and the number of time steps:
deltaT <- .1 # dt 
timelength <- 4000 # length of simulation
n_steps <- as.integer(timelength / deltaT) # must be a number greater than 1

## specify the magnitude of noise and specifics of shock 
D_stoeq1 <- 0.01  # before shock
t_shock <- 1000 # time that shock begins
shock_duration <- 500 # shock duration time points

# shock period 
shock_period <- data.frame(time = 0:timelength) |>
  mutate(shock = ifelse(time >= t_shock & time <= t_shock + shock_duration, TRUE, FALSE))

## aggregate symptom level
n_sims <- 500

# run simulation n_sims times
aggregated_low_beta2_sim500 <- map(1:n_sims, ~ euler_stochastic2(
  Amat = mod$A, 
  deterministic_rate = mod$dif_eq,
  stochastic_rate = mod$sto_eq,
  initial_condition = mod$initial_values,
  parameters1 = parms1,
  parameters2 = parms2, 
  deltaT = deltaT,
  timelength = timelength,
  D1 = D_stoeq1,
  shock = TRUE,
  t_shock = t_shock, 
  duration = shock_duration) |>
    mutate(totalsymptom = rowSums(pick(S_anh:S_sui)))
) |> list_rbind(names_to = "sim")

## save resutls
# saveRDS(aggregated_base_beta1_sim500, "data/aggregated500_base_beta1.rds") # beta_bistable1, scenario = base, nsim = 500
# saveRDS(aggregated_high_beta1_sim500, "data/aggregated500_high_beta1.rds") # beta_bistable1, scenario = high, nsim = 500
# saveRDS(aggregated_low_beta1_sim500, "data/aggregated500_low_beta1.rds") # beta_bistable1, scenario = low, nsim = 500
# saveRDS(aggregated_base_beta2_sim500, "data/aggregated500_base_beta2.rds") # beta_bistable2, scenario = base, nsim = 500
# saveRDS(aggregated_high_beta2_sim500, "data/aggregated500_high_beta2.rds") # beta_bistable2, scenario = high, nsim = 500
# saveRDS(aggregated_low_beta2_sim500, "data/aggregated500_low_beta2.rds") # beta_bistable2, scenario = low, nsim = 500

## load results
# define scenarios and beta values
scenarios <- c("base", "high", "low")
beta_values <- c("beta1", "beta2")
# initialize a list to store avgnet data
avgnet_data <- list()

# iterate over scenarios and beta values
for (beta in beta_values) {
  for (scenario in scenarios) {
    # construct the file path
    file_path <- paste0("data/aggregated500_", scenario, "_", beta, ".rds")
    
    # load the data, process it, and store only the selected avgnet data
    avgnet_data[[paste0(scenario, "_", beta)]] <- readRDS(file_path) |>
      dplyr::select(S_anh:S_sui) |>
      data.table::data.table()
    
    # update column names
    colnames(avgnet_data[[paste0(scenario, "_", beta)]]) <- nodenames
  }
}

# combine average networks for beta1 and beta2 scenarios
totavgnetb1 <- data.table::rbindlist(
  avgnet_data[grepl("_beta1$", names(avgnet_data))]  # filter datasets with beta1
)
totavgnetb2 <- data.table::rbindlist(
  avgnet_data[grepl("_beta2$", names(avgnet_data))]  # filter datasets with beta2
)
# saveRDS(totavgnetb1, "data/alldat_beta1.rds") # beta_bistable1
# saveRDS(totavgnetb2, "data/alldat_beta2.rds") # beta_bistable2

# then fit the network for all the scenarios aggregated
totavgnetworkb1_sim500 <- estimateNetwork(totavgnetb1, default = "EBICglasso")
totavgnetworkb2_sim500 <- estimateNetwork(totavgnetb2, default = "EBICglasso")
# saveRDS(totavgnetworkb1_sim500, "data/estnetwork_beta1.rds") # beta_bistable1
# saveRDS(totavgnetworkb2_sim500, "data/estnetwork_beta2.rds") # beta_bistable1
totavgnetworkb1_sim500 <- readRDS("data/estnetwork_beta1.rds")
totavgnetworkb2_sim500 <- readRDS("data/estnetwork_beta2.rds")

# compute centrality
cent_b1 <- centrality_auto(totavgnetworkb1_sim500)
rownames(cent_b1$node.centrality) <- nodenames

cent_b2 <- centrality_auto(totavgnetworkb2_sim500)
rownames(cent_b2$node.centrality) <- nodenames

cent_helius <- centrality_auto(totavgnetwork_helius)
rownames(cent_helius$node.centrality) <- nodenames

str_b1 <- cent_b1$node.centrality$Strength |>
  as_data_frame() |>  
  mutate(node = factor(rownames(cent_b1$node.centrality), 
                       levels= c(rownames(cent_b1$node.centrality)))
  )|> mutate(value = scale2(value))

str_b2 <- cent_b2$node.centrality$Strength |>
  as_data_frame() |>  
  mutate(node = factor(rownames(cent_b2$node.centrality), 
                       levels= c(rownames(cent_b2$node.centrality)))
  )|> mutate(value = scale2(value))

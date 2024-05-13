
## install packages
source("code/libraries.R")
## source necessary functions
source("code/euler_stochastic2.R")
source("code/mod_specification.R")

## standardizing function (Sacha qgraph)
scale2 <- function(x) {
  if (all(is.na(x))) return(NA)
  if (sd(x,na.rm=TRUE)!=0){
    return((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))
  } else {
    return(rep(0, length(x)))
  }
}

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


Beta_bistable1 <- criticalpt1 |> set_names(paste0("Beta_", colnames(A))) 
rownames(A) <- colnames(A) <- c("anh", "sad", "slp", "ene", "app", "glt", "con", "mot", "sui")
Beta_bistable2 <- criticalpt2 |> set_names(paste0("Beta_", colnames(A))) 
rownames(A) <- colnames(A) <- c("anh", "sad", "slp", "ene", "app", "glt", "con", "mot", "sui")


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



## save datasets
saveRDS(aggregated_base_beta1_sim500, "aggregated500_base_beta1.rds") # beta_bistable1, scenario = base, nsim = 500
saveRDS(aggregated_high_beta1_sim500, "aggregated500_high_beta1.rds") # beta_bistable1, scenario = high, nsim = 500
saveRDS(aggregated_low_beta1_sim500, "aggregated500_low_beta1.rds") # beta_bistable1, scenario = low, nsim = 500
saveRDS(aggregated_base_beta2_sim500, "aggregated500_base_beta2.rds") # beta_bistable2, scenario = base, nsim = 500
saveRDS(aggregated_high_beta2_sim500, "aggregated500_high_beta2.rds") # beta_bistable2, scenario = high, nsim = 500
saveRDS(aggregated_low_beta2_sim500, "aggregated500_low_beta2.rds") # beta_bistable2, scenario = low, nsim = 500

avgnetb1 <- aggregated_base_beta1_sim500 |> dplyr::select(S_anh:S_sui) 
colnames(aggregated_base_beta1_sim500) <- colnames(A)
avgnet_highb1 <- aggregated_high_beta1_sim500 |> dplyr::select(S_anh:S_sui)
colnames(aggregated_high_beta1_sim500) <- colnames(A)
avgnet_lowb1 <- aggregated_low_beta1_sim500 |> dplyr::select(S_anh:S_sui)
colnames(aggregated_low_beta1_sim500) <- colnames(A)
avgnetb2 <- aggregated_base_beta2_sim500 |> dplyr::select(S_anh:S_sui) 
colnames(aggregated_base_beta2_sim500) <- colnames(A)
avgnet_highb2 <- aggregated_high_beta2_sim500 |> dplyr::select(S_anh:S_sui)
colnames(aggregated_high_beta2_sim500) <- colnames(A)
avgnet_lowb2 <- aggregated_low_beta2_sim500 |> dplyr::select(S_anh:S_sui)
colnames(aggregated_low_beta2_sim500) <- colnames(A)

# saveRDS(aggregated, "aggregated_base_beta1.rds") # beta_bistable1, scenario = base, nsim = 300
# saveRDS(aggregated_high_beta1, "aggregated_high_beta1.rds") # beta_bistable1, scenario = high, nsim = 300
# saveRDS(aggregated_low_beta1, "aggregated_low_beta1.rds") # beta_bistable1, scenario = low, nsim = 300
# saveRDS(aggregated_base_beta2, "aggregated_base_beta2.rds") # beta_bistable2, scenario = base, nsim = 300
# saveRDS(aggregated_high_beta2, "aggregated_high_beta2.rds") # beta_bistable2, scenario = high, nsim = 300
# saveRDS(aggregated_low_beta2, "aggregated_low_beta2.rds") # beta_bistable2, scenario = low, nsim = 300

# then fit the network for all the scenarios aggregated

# avgnetb1 <- aggregated |> dplyr::select(S_anh:S_sui) 
# colnames(avgnetb1) <- colnames(A)
# avgnet_highb1 <- aggregated_high_beta1 |> dplyr::select(S_anh:S_sui)
# colnames(avgnet_highb1) <- colnames(A)
# avgnet_lowb1 <- aggregated_low_beta1 |> dplyr::select(S_anh:S_sui)
# colnames(avgnet_lowb1) <- colnames(A)
# avgnetb2 <- aggregated_base_beta2 |> dplyr::select(S_anh:S_sui) 
# colnames(avgnetb2) <- colnames(A)
# avgnet_highb2 <- aggregated_high_beta2 |> dplyr::select(S_anh:S_sui)
# colnames(avgnet_highb2) <- colnames(A)
# avgnet_lowb2 <- aggregated_low_beta2 |> dplyr::select(S_anh:S_sui)
# colnames(avgnet_lowb2) <- colnames(A)

# totavgnet2 <- bind_rows(avgnet2, avgnet_high2, avgnet_low2)
# total average net 
totavgnetb1 <- avgnetb1 |> bind_rows(avgnet_highb1, avgnet_lowb1)
totavgnetb2 <- avgnetb2 |> bind_rows(avgnet_highb2, avgnet_lowb2)

totavgnetworkb1 <- estimateNetwork(totavgnetb1, default = "EBICglasso")
totavgnetworkb2 <- estimateNetwork(totavgnetb2, default = "EBICglasso")
# saveRDS(totavgnetworkb1, "estnetwork_beta1.rds") # beta_bistable1
# saveRDS(totavgnetworkb2, "estnetwork_beta2.rds") # beta_bistable1
# totavgnetworkb1 <- readRDS("data/estnetwork_beta1.rds")
# totavgnetworkb2 <- readRDS("data/estnetwork_beta2.rds")

totavgnetworkb1_sim500 <- estimateNetwork(totavgnetb1, default = "EBICglasso")
totavgnetworkb2_sim500 <- estimateNetwork(totavgnetb2, default = "EBICglasso")
# saveRDS(totavgnetworkb1_sim500, "estnetwork_beta1.rds") # beta_bistable1
# saveRDS(totavgnetworkb1_sim500, "estnetwork_beta2.rds") # beta_bistable1
totavgnetworkb1_sim500 <- readRDS("data/estnetwork_beta1.rds")
totavgnetworkb2_sim500 <- readRDS("data/estnetwork_beta2.rds")


plot(totavgnetworkb1_sim500)
plot(totavgnetworkb2_sim500)

# compute centrality
cent_b1 <- centrality_auto(totavgnetworkb1_sim500)
rownames(cent_b1$node.centrality) <- colnames(A)

cent_b2 <- centrality_auto(totavgnetworkb2_sim500)
rownames(cent_b2$node.centrality) <- colnames(A)

cent_helius <- centrality_auto(totavgnetwork_helius)
rownames(cent_helius$node.centrality) <- colnames(A)

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

str_helius <- cent_helius$node.centrality$Strength |>
  as_data_frame() |>  
  mutate(node = factor(rownames(cent_helius$node.centrality), 
                       levels= c(rownames(cent_helius$node.centrality)))
  )|> mutate(value = scale2(value))

str_heliusH1 <- cent_heliusH1$node.centrality$Strength |>
  as_data_frame() |>  
  mutate(node = factor(rownames(cent_heliusH1$node.centrality), 
                       levels= c(rownames(cent_heliusH1$node.centrality)))
  )|> mutate(value = scale2(value))

str_poly_helius <- cent_poly_helius$node.centrality$Strength |>
  as_data_frame() |>  
  mutate(node = factor(rownames(cent_poly_helius$node.centrality), 
                       levels= c(rownames(cent_poly_helius$node.centrality)))
  )  |> mutate(value = scale2(value))


str_norm_helius <- cent_norm_helius$node.centrality$Strength |>
  as_data_frame() |>  
  mutate(node = factor(rownames(cent_norm_helius$node.centrality), 
                       levels= c(rownames(cent_norm_helius$node.centrality)))
  )

str_sick_helius <- cent_sickhelius$node.centrality$Strength |>
  as_data_frame() |>  
  mutate(node = factor(rownames(cent_sickhelius$node.centrality), 
                       levels= c(rownames(cent_sickhelius$node.centrality)))
  ) 

str_helius$dummy <- "Strength"


strength_cent <- bind_rows(str_b1, str_b2, .id="id") |> 
  pivot_wider(names_from = id, names_prefix = "beta") |>
  rowwise() |>
  mutate(avg = mean(c(beta1, beta2))) |>
  full_join(str_helius) |>
  rename("helius" = value) |>
  pivot_longer(c(avg, helius), names_to = "grp", values_to = "avg")
strength_cent$dummy <- "Standardized Strength"


cent_helius_compare <- strength_cent |>
  ggplot(aes(x = avg, y = node, group = grp, color = grp)) + 
  geom_point(alpha = 0.7) +
  geom_path() +
  scale_color_manual(values = alpha(c("deepskyblue4", "indianred3"),  0.7), labels = c("Simulated", "HELIUS")) +
  geom_errorbarh(aes(xmin = beta1, xmax = beta2), color = "deepskyblue4", alpha = 0.3, height = 0.5) +
  theme_bw() +
  labs(x = "", y = "", color = "") +
  theme(text=element_text(size=23),
        legend.position = "bottom") +
  facet_grid(. ~ dummy)

ggsave("cent_helius_compare.png", plot = cent_helius_compare, width = 4.5, height =5.5, dpi = 300)

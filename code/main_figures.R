## =========================================================
## Main Figures
##
## This script generates:
## - Figure 2: Median aggregate symptom levels and GGMs over time 
##             for individuals with high and low resilience.
## - Figure 3: Comparison of overall network structures derived from aggregated simulated data, 
##             the reference mechanistic network, and the empirical network constructed from HELIUS data.
## =========================================================

## install packages
source("code/libraries.R")
## source necessary functions
source("code/utils.R")

## ---------------------------------------
## Figure 2: Create combined result figure
## ---------------------------------------
# load data
low <- readRDS("data/aggregated_low.rds") |> data.table::data.table()
med <- readRDS("data/aggregated.rds") |> data.table::data.table()
high <- readRDS("data/aggregated_res.rds") |> data.table::data.table()

# shock 
t_shock <- 1000 # time that shock begins
shock_duration <- 500 # shock duration time points
timelength <- 4000 # length of simulation

# shock period 
shock_period <- data.frame(time = 0:timelength) |>
  mutate(shock = ifelse(time >= t_shock & time <= t_shock + shock_duration, TRUE, FALSE))

# combine data (low & high)
comb <- data.table::rbindlist(list(low,high), idcol="ID") |>
  reframe(inter_quantile(totalsymptom), .by = c(t, ID)) |>
  mutate(ID = ifelse(ID == 1, "low", "high")) |>
  pivot_wider(names_from = "quant", values_from = "val", names_glue = "q{quant}") 


## combined dataset figure for total symptoms
p_comb <- ggplot(data = comb) +
  geom_line(aes(x = t, y = q0.5, col = ID), lwd = 0.2) +
  geom_ribbon(aes(x = t, ymin = q0.25, ymax = q0.75, fill = ID), alpha=0.2) +
  geom_area(data = shock_period, aes(x = time, y = shock*max(low$totalsymptom)), 
            inherit.aes = FALSE, fill = "darkgray", alpha = 0.2) +
  scale_color_manual("resilience", values = c("#006699", "#CC3300"), 
                     aesthetics = c("color", "fill")) +
  labs(x = "time", y= "aggregated symptom level") +
  geom_hline(yintercept = 5/3, linetype = 2, color = "azure4", lwd = 0.1) +
  geom_hline(yintercept = 10/3, linetype = 2, color = "azure4", lwd = 0.1) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")
# ggsave("totalsym_comb.png", plot = p_comb, width = 20, height =10, units = "cm", dpi = 300)


## fitting statistical networks
burnin <- 100
aggregated <- low #high
# before shock
beforeshock <- aggregated |> 
  filter(t > burnin & t < t_shock) |> 
  select(S_anh:S_sui) |> 
  estimateNetwork(default = "EBICglasso")
# during shock
duringshock <- aggregated |> 
  filter(t >= t_shock & t <= t_shock + shock_duration) |> 
  select(S_anh:S_sui) |>
  estimateNetwork(default = "EBICglasso")
# after shock
aftershock <- aggregated |> 
  # filter(t >= 3000) |> 
  filter(t > t_shock + shock_duration + burnin) |> 
  select(S_anh:S_sui) |>
  estimateNetwork(default = "EBICglasso")


# hyperparameter *max_value*
max_value <- max(
  max(abs(beforeshock$graph)), 
  max(abs(duringshock$graph)),
  max(abs(aftershock$graph))
)

# hyperparameter *net_layout*
# net_layout <- averageLayout(beforeshock,
#                             duringshock,
#                             aftershock)
# save layout used in the paper
# saveRDS(net_layout, "data/net_layout.rds")
net_layout <- readRDS("data/net_layout.rds")

# create node names
nodenames <- c("anhedonia", "sad", "sleep", "energy", "appetite", "guilty", "concentration", "motor", "suicidal")

# pdf(file = "triplenetworks2_v2.pdf", width = 16, height = 9)

# png(file = "triplenet_low.png", width = 1800, height = 900)
par(mfrow = c(1, 3), mar = c(3, 9, 6, 0.5), xpd = NA, oma = c(0, 2, 0.5, 0))

plot(beforeshock, layout = net_layout, maximum = max_value, labels = nodenames, node.width=2, border.color = "#CC3300", edge.color = "#CC3300")
#title("(a) Before shock", line = 3, cex.main = 5)
mtext("Low resilience", side=2, padj = 1, cex=3.5, col="black", outer= T, font =2)

plot(duringshock, layout = net_layout, maximum = max_value,  labels = nodenames,node.width=2, border.color = "#CC3300", edge.color = "#CC3300")
#title("(b) During shock", line = 3, cex.main = 5)
rect(-1.3, -2, 1.3, 2.0, col = rgb(0.5, 0.5, 0.5 ,alpha=0.1), border=FALSE) 

plot(aftershock, layout = net_layout, maximum = max_value, labels = nodenames, node.width=2, border.color = "#CC3300", edge.color = "#CC3300")
#title("(c) After shock", line = 3, cex.main = 5)
abline(v = -3.9, lty=3, col = "darkgray")
abline(v = -1.3, lty=3, col = "darkgray")

# dev.off()




## --------------------------------------------------
## Figure 3: Create ref & sim networks comparison fig
## --------------------------------------------------

# set layout for all three networks
manual_layout <- matrix(c( -1.00000000,-0.3697451,
                           -0.25362943,-0.4206165,
                           -0.86442347, 0.4941126,
                           -0.08080896, 0.3245454,
                           0.49177041, 0.9000000,
                           0.43942364,-0.1656119,
                           0.99799343, 0.2837986,
                           1.00000000,-0.7135021,
                           0.41570786,-1.0000000), 9, 2, byrow=T)

# weigthed adjacency matrix for reference network
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
diag(A) <- 0.1 # for plotting purpose

# simulated network (saved result from `main_simulation.R`)
totavgnetwork <- readRDS("data/overallnetwork.Rds")
totavgnet <- plot(totavgnetwork, labels = nodenames)

# data network (saved result from `data_network.R`)
helius_H1net <- readRDS("data/helius_H1net.rds")


## plot networks
# png(file = "network_comparison_v2.png", width = 2300, height = 1000)
# ref network
par(mfrow = c(1, 3), mar = c(1, 9, 4, 0.5), xpd = NA, oma = c(0, 2, 3, 0))

ref_net <- qgraph(A, theme = 'colorblind', color = c("lightsteelblue1", rep("bisque",5), rep("lightsteelblue1",3)), border.color = "white",border.width = 2, edge.width = 0.8, curve = 0.3, curveAll = T, label.color = "black",  asize= 8,  diag = T, layout = manual_layout, node.width=2, edge.color = "deepskyblue4")
title("Mechanistic network", line = 3, cex.main = 5.5, adj = 0.5)
# p_ref <- recordPlot()  

# sim network
sim_net <- qgraph(totavgnet, labels = nodenames, layout = manual_layout, node.width=2, edge.color = "indianred3")
title("Simulated network", line = 3, cex.main = 5.5, adj = 0.5)
# p_networks <- recordPlot()  

# helius network
hel_net <- plot(helius_H1net, labels = colnames(A), layout = manual_layout, minimum = 0.03, node.width=2, edge.color ="darkseagreen4")
title("Empirical network", line = 3, cex.main = 5.5, adj = 0.5)
# p_hel <- recordPlot()  

# dev.off()


## compute centrality
cent_refnet <- centrality_auto(ref_net)
rownames(cent_refnet$node.centrality) <- colnames(A)

cent_simnet <- centrality_auto(totavgnetwork)
rownames(cent_simnet$node.centrality) <- colnames(A)

cent_helius <- centrality_auto(helius_H1net)
rownames(cent_helius$node.centrality) <- colnames(A)

str_refnet <- cent_refnet$node.centrality$InStrength |> # in-strength
  as_data_frame() |>  
  mutate(node = factor(rownames(cent_refnet$node.centrality), 
                       levels= c(rownames(cent_refnet$node.centrality)))
         # standardize
  ) |> mutate(value = scale2(value))

str_simnet <- cent_simnet$node.centrality$Strength |>
  as_data_frame() |>  
  mutate(node = factor(rownames(cent_simnet$node.centrality), 
                       levels= c(rownames(cent_simnet$node.centrality)))
         # standardize
  ) |> mutate(value = scale2(value))

str_heliusnet <- cent_helius$node.centrality$Strength |>
  as_data_frame() |>  
  mutate(node = factor(rownames(cent_helius$node.centrality), 
                       levels= c(rownames(cent_helius$node.centrality)))
         # standardize
  )|> mutate(value = scale2(value))


## b1, b2 strength centrality from "mod_flexibility.R"
totavgnetworkb1_sim500 <- readRDS("data/estnetwork_beta1.rds")
totavgnetworkb2_sim500 <- readRDS("data/estnetwork_beta2.rds")
cent_b1 <- centrality_auto(totavgnetworkb1_sim500)
rownames(cent_b1$node.centrality) <- nodenames
cent_b2 <- centrality_auto(totavgnetworkb2_sim500)
rownames(cent_b2$node.centrality) <- nodenames

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


combined_cent <- bind_rows(str_b1, str_b2, .id="id") |> 
  pivot_wider(names_from = id, names_prefix = "beta") |>
  rowwise() |>
  mutate(sim = mean(c(beta1, beta2))) |>
  full_join(str_heliusnet) |>
  rename("helius" = value) |>
  full_join(str_refnet) |>
  rename("ref" = value) |>
  pivot_longer(c(sim, helius, ref), names_to = "grp", values_to = "val")
combined_cent$dummyA <- "Standardized (In-)Strength"
combined_cent$dummyB <- "Standardized Strength"

## all three centralities in one 
cent_all <- combined_cent |>
  ggplot(aes(x = val, y = node, group = grp, color = grp)) + 
  geom_point(alpha = 0.7) +
  geom_path() +
  scale_color_manual(values = alpha(c("darkseagreen4", "deepskyblue4", "indianred3"),  0.7), labels = c("HELIUS", "Reference", "Simulated")) +
  geom_errorbarh(aes(xmin = beta1, xmax = beta2), color = "indianred3", alpha = 0.3, height = 0.5) +
  theme_bw() +
  labs(x = "", y = "", color = "") +
  theme(text=element_text(size=15),
        legend.position = "bottom") +
  facet_grid(. ~ dummy)
# grab legend from all cent plot
legend <- cowplot::get_legend(cent_all)

## centrality plot for toy vs sim
cent1 <- combined_cent |> filter(grp == 'ref' | grp == 'sim') |>
  ggplot(aes(x=node, y = val, group=grp, color = grp)) + 
  geom_point() + 
  geom_path() + 
  geom_errorbar(aes(ymin = beta1, ymax = beta2), color = "indianred3", alpha = 0.2, width = 0.4) +
  scale_color_manual(values = alpha(c("deepskyblue4", "indianred3"),  0.7), labels = c("Mechanistic", "Simulated")) +
  ylim(-2.5,2) +
  theme_bw() +
  labs(x = "", y = "", color = "") +
  theme(text=element_text(size=15),
        legend.position = "none") +
  facet_wrap(. ~ dummyA, strip.position = "left")
  
## centrality plot for toy vs sim
cent2 <- combined_cent |> filter(grp == 'helius' | grp == 'sim') |>
  ggplot(aes(x=node, y = val, group=grp, color = grp)) + 
  geom_point() + 
  geom_path() + 
  geom_errorbar(aes(ymin = beta1, ymax = beta2), color = "indianred3", alpha = 0.2, width = 0.4) +
  scale_color_manual(values = alpha(c("darkseagreen4", "indianred3"),  0.7), labels = c("HELIUS", "Simulated")) +
  theme_bw() +
  labs(x = "", y = "", color = "") +
  ylim(-2.5,2) +
  theme(text=element_text(size=15),
        legend.position = "none") +
  facet_wrap(. ~ dummyB, strip.position = "left")

# put them together
centrality_comp <- cowplot::plot_grid(cent1, cent2, legend, ncol=1, align='v', rel_heights = c(4,4,1))
# ggsave("centrality_comp.png", plot = centrality_comp, width = 20, height = 17, units = "cm", dpi = 300)


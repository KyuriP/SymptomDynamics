
source("code/libraries.R")
library(data.table)



# load data
low <- readRDS("data/aggregated_low.rds") |> data.table::data.table()
med <- readRDS("data/aggregated.rds") |> data.table::data.table()
high <- readRDS("data/aggregated_res.rds") |> data.table::data.table()

## IQR function
inter_quantile <- function(x, probs = c(0.25, 0.5, 0.75)) {
  tibble(
    val = quantile(x, probs, na.rm = TRUE),
    quant = probs
  )
}

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


ggsave("totalsym_comb.png", plot = p_comb, width = 20, height =10, units = "cm", dpi = 300)

## fitting statistical networks
burnin <- 100
aggregated <- low
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
  filter(t >= 3000) |> 
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
# net_layout <- readRDS("data/net_layout.rds")

# create node names
nodenames <- c("anhedonia", "sad", "sleep", "energy", "appetite", "guilty", "concentration", "motor", "suicidal")

# pdf(file = "triplenetworks2_v2.pdf", width = 16, height = 9)

png(file = "triplenet_low.png", width = 1800, height = 900)
par(mfrow = c(1, 3), mar = c(3, 9, 6, 0.5), xpd = NA, oma = c(0, 2, 0.5, 0))

plot(beforeshock, layout = net_layout, maximum = max_value, labels = colnames(A), node.width=2, border.color = "#CC3300", edge.color = "#CC3300")
#title("(a) Before shock", line = 3, cex.main = 5)
mtext("High resilience", side=2, padj = 1, cex=3.5, col="black", outer= T, font =2)

plot(duringshock, layout = net_layout, maximum = max_value,  labels = colnames(A),node.width=2, border.color = "#CC3300", edge.color = "#CC3300")
#title("(b) During shock", line = 3, cex.main = 5)
rect(-1.3, -2, 1.3, 2.0, col = rgb(0.5, 0.5, 0.5 ,alpha=0.1), border=FALSE) 

plot(aftershock, layout = net_layout, maximum = max_value, labels = colnames(A), node.width=2, border.color = "#CC3300", edge.color = "#CC3300")
#title("(c) After shock", line = 3, cex.main = 5)
abline(v = -3.9, lty=3, col = "darkgray")
abline(v = -1.3, lty=3, col = "darkgray")

dev.off()

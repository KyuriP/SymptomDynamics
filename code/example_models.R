## =========================================================
## Example network models
##
## creating Figure 1 (example bistable landscape with networks),
## Figure 2 (example directed network model) and
## Figure 3 (centrality of the example network model)
## =========================================================
## install packages
source("code/libraries.R")


## Creating Figure 1: example bistable landscape with networks

# bistable landscape for Fig1
x = seq(-5, 5, 0.01)
y = 1/5 * x^4 - 4*x^2 + 10

plot(x, y, type = "l")


# example networks for Fig1
fakeA <- matrix(c(0,0.1,0,0,0.1,
                  0.1,0,0,0.1,0,
                  0,0,0,0,0,
                  0,0.1,0,0,0,
                  0.1,0,0,0.1,0), 5, 5, byrow = T)

fakeA[lower.tri(fakeA)] = t(fakeA)[lower.tri(fakeA)]
rownames(fakeA) <- colnames(fakeA) <- nl

fakeB <-matrix(c(0,1,1,1,0,
                 1,0,1,0.1,1,
                 1,0,0,0,1,
                 1,1,1,0,1,
                 0,1,0,1,0), 5, 5, byrow = T)
rownames(fakeB) <- colnames(fakeB) <- nl
fakeB[lower.tri(fakeB)] = t(fakeB)[lower.tri(fakeB)]

nl <- c("S1", "S2", "S3", "S4", "S5")
qgraph(fakeA, theme = 'colorblind', maximum = 0.3)
qgraph(fakeB, theme = 'colorblind')

png(file="healthynetwork.png", units="in", width=5, height=5, res=300, bg = 'transparent')
qgraph(fakeA, theme = 'colorblind', maximum = 0.3)
dev.off()


png(file="sicknetwork.png", units="in", width=5, height=5, res=300, bg = 'transparent')
qgraph(fakeB, theme = 'colorblind')
dev.off()


## Creating Figure 2: example directed network 
Names <- c("anhedonia", "sadness", "sleep", "energy", "appetite", "guilty", "concentration", "motor", "suicidal")

grp <- list(`cycle` = 2:6, `no cycle` = c(1,7:9))


manual_layout <- matrix(c( -1.00000000,-0.3697451,
                           -0.25362943,-0.4206165,
                           -0.86442347, 0.4941126,
                           -0.08080896, 0.3245454,
                           0.49177041, 0.9000000,
                           0.43942364,-0.1656119,
                           0.99799343, 0.2837986,
                           1.00000000,-0.7135021,
                           0.41570786,-1.0000000), 9, 2, byrow=T)


pdf(file = "toymodel2.pdf", width=16, height=12, bg = 'transparent')

ho <- qgraph(A, theme = 'colorblind', legend = TRUE, groups = grp, color = c("#F5651C", "#58B5BC"), nodeNames = Names, border.color = "white",border.width = 2, edge.color = "darkgray", edge.width = 0.8, curve = 0.3, curveAll = T, label.color = "white", legend.cex = 1.2, asize= 4, layout = manual_layout)

dev.off()


# diag(A) <- 0.1
png(file = "toymodelnet.png", width=1000, height=1000, bg = 'transparent')
qgraph(A, theme = 'colorblind', legend = F, groups = grp, color = c("#F5651C", "#58B5BC"), border.color = "white",border.width = 2, edge.color = "darkgray", edge.width = 0.8, curve = 0.3, curveAll = T, label.color = "white", legend.cex = 1.2, asize= 4, layout = manual_layout, diag = T)
dev.off()


## specifying weight matrix based on empirical networks
# lee et al.(2023) partial correlations
net1 <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0.32, 0, 0, 0, 0, 0, 0, 0, 0,
                 0.24, 0.07, 0, 0, 0, 0, 0, 0, 0,
                 0.17, 0.13, 0.24, 0, 0, 0, 0, 0, 0,
                 0.09, 0.04, 0.13, 0.24, 0, 0, 0, 0, 0,
                 0.08, 0.16, 0.00, 0.00, 0.14, 0, 0, 0, 0,
                 0.00, 0.04, 0.13, 0.06, 0.10, 0.17, 0, 0, 0,
                 0.06, 0.00, 0.10, 0.03, 0.01, 0.23, 0.19, 0, 0,	
                 0.00, 0.29, 0.01, 0.00, 0.00, 0.33, 0.12, 0.18, 0), 9, 9, byrow = T)

net1[upper.tri(net1)] = t(net1)[upper.tri(net1)]

# peruvian Ramos-Vera (2021) regularized partial correlations (GGM)
net2 <- matrix(c(0.000, 0.318, 0.095, 0.170, 0.101, 0.082,	0.088, 0.000,	0.000,
                 0.318, 0.000, 0.209, 0.193, 0.081, 0.042,	0.081, 0.071,	0.082,
                 0.095, 0.209, 0.000, 0.152, 0.152, 0.059,	0.100, 0.067,	0.000,
                 0.170, 0.193, 0.152, 0.000, 0.122, 0.142,	0.104, 0.038,	0.046,
                 0.101, 0.081, 0.152, 0.122, 0.000, 0.141,	0.062, 0.051,	0.033,
                 0.082, 0.042, 0.059, 0.142, 0.141, 0.000,	0.261, 0.000,	0.141,
                 0.088, 0.081, 0.100, 0.104, 0.062, 0.261,	0.000, 0.084,	0.112,
                 0.000, 0.071, 0.067, 0.038, 0.051, 0.000,	0.084, 0.000,	0.396,
                 0.000,	0.082, 0.000,	0.046, 0.033, 0.141,  0.112, 0.396,	0.000), 9, 9, byrow = T)

# adolescents Jing et al.(2023) # regularized partial correlations
net3 <- matrix(c(0    , 0.353, 0.056, 0.295, 0.045, 0 ,0.141, 0.067, 0.014,
                 0.353, 0, 0.155, 0.120, 0.067, 0.180, 0, 0.080, 0.077,
                 0.056, 0.155, 0, 0.281, 0.114, 0.036, 0, 0.014, 0.087,
                 0.295, 0.120, 0.281, 0, 0.156 ,0.036, 0.188, 0.070, -0.066,
                 0.045, 0.067, 0.114, 0.156, 0, 0.168, 0.071, 0.094, 0.110,
                 0    , 0.180, 0.036, 0.036, 0.168, 0, 0.166, 0.209, 0.175,
                 0.141, 0 ,0 ,0.188, 0.071, 0.166, 0 ,0.250, 0.001,
                 0.067, 0.080, 0.014, 0.070, 0.094, 0.209, 0.250, 0 ,0.323,
                 0.014, 0.077, 0.087, -0.066, 0.110, 0.175, 0.001, 0.323, 0), 9, 9, byrow =T)

# Wuhan Zhang et al.(2023) # Ising model weighted adjacency matrix
net4 <- matrix(c(0.00,	1.72,	0.16,	1.99,	0.68,	0.13,	0.78,	0.16,	0.00,
                 1.72,	0.00,	0.40,	1.03,	0.55,	1.36,	0.42,	0.75,	0.86,
                 0.16,	0.40,	0.00,	1.58,	1.02,	0.22,	0.46,	0.31,	0.00,
                 1.99,	1.03,	1.58,	0.00,	0.87,	0.94,	0.16,	0.87,	0.00,
                 0.68,	0.55,	1.02,	0.87,	0.00,	0.58,	0.64,	0.34,	0.58,
                 0.13,	1.36,	0.22,	0.94,	0.58,	0.00,	0.96,	0.89,	1.51,
                 0.78,	0.42,	0.46,	0.16,	0.64,	0.96,	0.00,	1.91,	0.00,
                 0.16,	0.75,	0.31,	0.87,	0.34,	0.89,	1.91,	0.00,	1.44,
                 0.00,	0.86,	0.00,	0.00,	0.58,	1.51,	0.00,	1.44,	0.00), 9, 9, byrow =T)

# Hongkong Cheung et al.(2021) Ising model weighted adjacency matrix
net5 <- matrix(c(0.00,	2.14,	0.32,	1.27,	0.64,	0.50,	0.66,	0.31,	0.30,
                 2.14,	0.00,	0.50,	0.78,	0.52,	1.35,	0.29,	0.59,	1.13,
                 0.32,	0.50,	0.00,	1.35,	1.28,	0.44,	0.44,	0.41,	0.43,
                 1.27,	0.78,	1.35,	0.00,	1.52,	0.74,	0.51,	0.46,	0.00,
                 0.64,	0.52,	1.28,	1.52,	0.00,	0.73,	0.61,	0.80,	0.28,
                 0.50,	1.35,	0.44,	0.74,	0.73,	0.00,	1.30,	1.19,	1.70,
                 0.66,	0.29,	0.44,	0.51,	0.61,	1.30,	0.00,	1.64,	0.67,
                 0.31,	0.59,	0.41,	0.46,	0.80,	1.19,	1.64,	0.00,	1.03,
                 0.30,	1.13,	0.43,	0.00,	0.28,	1.70,	0.67,	1.03,	0.00), 9, 9, byrow=T)



# get the average weight matrices
ggmlist <- list(net1, net2, net3)
isinglist <- list(net4, net5)

averageGGM <- Reduce("+", ggmlist) / length(ggmlist)
averageGGM <- ifelse(averageGGM < 0.1, 0, averageGGM)
averageIsing <- Reduce("+", isinglist) / length(isinglist)

round(averageGGM, 2)
round(averageIsing, 2)

colnames(averageGGM) <- rownames(averageGGM) <- Names
colnames(averageIsing) <- rownames(averageIsing) <- Names

# example directed network 
Names <- c("anhedonia", "sadness", "sleep", "energy", "appetite", "guilt", "concentration", "motor", "suicidal")
layout(1)
qgraph(averageGGM, theme = 'colorblind',legend = TRUE,  nodeNames = Names, layout="spring")

qgraph(averageIsing, layout="spring")


## Creating Figure 3: centrality of example network model

A <- matrix(c( .30, 0, 0, 0, 0, 0, 0, 0, 0,
               .33, .30, .14, .15, 0, .13, 0, 0, .15,
               .13, .14, .30, .22, .23, 0, 0, 0, 0,
               .21, .15, .22, .30, 0, 0, .12, 0, 0,
               0, 0, 0, .17, .30, 0, 0, 0, 0,
               0, .13, 0, 0, .15, .30, .2, .15, .22,
               0, 0, 0, 0, 0, 0, .30, .17, 0,
               0, 0, 0, 0, 0, 0, 0, .30, 0,
               0, 0, 0, 0, 0, 0, 0, .3, 0.30), 9, 9, byrow = T)

cent <- centrality_auto(A)
stn <- (cent$node.centrality$InStrength + cent$node.centrality$OutStrength)/2 

centralityA <- tibble(InStrength = cent$node.centrality$InStrength, OutStrength = cent$node.centrality$OutStrength, Strength = stn) |> 
  mutate(node = factor(rownames(cent$node.centrality), levels= c(rownames(cent$node.centrality)))
  )

centralityA %<>%  pivot_longer(!node, names_to = "centrality", values_to = "value")

cent_A <- centralityA |>
  ggplot(aes(x=value, y = node, group=1, color = centrality)) + 
  geom_point() + 
  geom_path() + 
  theme_bw() +
  labs(x = "", y = "") +
  theme(text=element_text(size=20))+
  facet_wrap(~centrality) 

ggsave("centrality_A.pdf", plot = cent_A, width = 25, height =20, units = "cm", dpi = 300)


centralityA$dummyA <- "In/Out Strength"
centralityA$dummyB <- "Average Strength"


plotA <- centralityA |> filter(centrality == "InStrength" | centrality == "OutStrength") |>
  ggplot(aes(x=value, y = node, group=centrality, color = centrality)) + 
  geom_point() + 
  geom_path() + 
  theme_bw() +
  labs(x = "", y = "", color = "") +
  scale_color_brewer(palette = "Set2") + 
  theme(text=element_text(size=20),
        legend.position = "bottom",
        legend.text = element_text(size=15)) +
  facet_grid(. ~ dummyA)
ggsave("cent_toy.png", plot = plotA, width = 4.5, height = 8, dpi = 300)


plotA <- centralityA |> filter(centrality == "InStrength" ) |>
  ggplot(aes(x=value, y = node, group=centrality, color = centrality)) + 
  geom_point() + 
  geom_path() + 
  theme_bw() +
  labs(x = "", y = "", color = "") +
  lims(x = c(0,0.9)) +
  scale_color_brewer(palette = "Set2") + 
  theme(text=element_text(size=20),
        legend.position = "bottom",
        legend.text = element_text(size=15)) +
  facet_grid(. ~ dummyA)
ggsave("cent_toy_instrength.png", plot = plotA, width = 4.5, height = 8, dpi = 300)



plotB <- centralityA |> filter(centrality == "Strength") |>
  ggplot(aes(x=value, y = node, group=1)) + 
  geom_point(color = "deepskyblue4") + 
  geom_path(color = "deepskyblue4") + 
  theme_bw() +
  labs(x = "", y = "", color = "") +
  theme(text=element_text(size=20),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  facet_grid(. ~ dummyB)

plotcomb <- ggpubr::ggarrange(plotA, plotB, widths = c(1.15,1), common.legend = TRUE, legend = 'bottom')

ggsave("centrality_comb.pdf", plot = plotcomb, width = 25, height =20, units = "cm", dpi = 300)

## a simpler toy model based on helius

round(ifelse(helius_aggnet$graph < 0.05, 0, helius_aggnet$graph),2)
round(ifelse(helius_aggnet$graph>0.05, helius_aggnet$graph, 0),2)

A_hel <- matrix(c( .30, 0, 0, 0, .12, 0, 0, 0, 0,
                   .39, .30, .10, 0, 0, .13, 0, .12, .14,
                   0, 0, .30, .3, .12, 0, 0, 0, 0,
                   .23, 0, 0, .30, 0, 0, .12, 0, 0,
                   0, 0, 0, .20, .30, 0, 0, 0, 0,
                   0, .13, 0, 0, .12, .30, .14, .13, .21,
                   0, 0, 0, 0, 0, 0, .30, 0.25, 0,
                   0, 0, 0, 0, 0, 0, 0, .30, 0,
                   0, 0, 0, 0, 0, 0, 0, .12, .30), 9, 9, byrow = T)
rownames(A_hel) <- colnames(A_hel) <- c("anh", "sad", "slp", "ene", "app", "glt", "con", "mot", "sui")

qgraph(A_hel, theme = 'colorblind')



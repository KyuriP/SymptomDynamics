## =========================================================
## HELIUS data
##
## 
## =========================================================
## install packages
source("code/libraries.R")

# import HELIUS data
helius2 <- read_sav("data/HELIUS_itemscores.sav")
# helius2$H1_WlbvRecent8
# helius2$H1_WlbvRecent9
# helius2$H1_WlbvRecent_89


# clean up depression item scores
dep_scores <- helius2 |> 
  select(contains("WlbvRecent") & !ends_with("Ingevuld"), ID) |>
  pivot_longer(-ID, names_to = "wave1", values_to = "value") |>
  mutate(wave = str_extract(wave1, "[^_]+"), item_num = str_extract(wave1,  "(\\d+)(?!.*\\d)")) |>
  drop_na() |>
  select(-wave1) %>% 
  filter(!(item_num == "8" | item_num == "9")) |>
  mutate(item_num = case_when(item_num == "1" ~ "anh",
            item_num == "2"~ "sad",
            item_num == "3"~ "slp",
            item_num == "4"~ "ene",
            item_num == "5"~ "app",
            item_num == "6"~ "glt",
            item_num == "7"~ "con",
            item_num == "89" ~ "mot",
            item_num == "10" ~ "sui")
            )

dep_list <- dep_scores %>% 
  group_by(wave) %>% 
  {setNames(group_split(.), group_keys(.)[[1]])}  %>% 
  map(~.x |>pivot_wider(id_cols = ID, names_from = item_num, values_from = value))

# aggregate all waves
helius_aggnet <- dep_list |> list_rbind() |> drop_na() |> select(-ID) |>
  relocate(sui, .after = mot) |>
  estimateNetwork(default = "EBICglasso") 

png(file = "helius_totavgnetwork.png",bg = 'transparent', res=300)
plot(helius_aggnet, layout = aggnetLayout)  
dev.off()

# separate each wave 
net_list <- dep_list |>
  map(~.x |>select(-ID) |> relocate(sui, .after = mot) |>
      estimateNetwork(.x, default = "EBICglasso") |>
        plot(.x,  layout = aggnetLayout)
      )


plot(net_list$H1)
plot(net_list$H2)

res_netH1_helius <- centrality_auto(net_list$H1)
rownames(res_netH1_helius$node.centrality) <- colnames(A)

cent_netH1_helius <- res_netH1_helius$node.centrality$Strength |>as_data_frame() |>  mutate(node = factor(rownames(res_netH1_helius$node.centrality), levels= c(rownames(res_netH1_helius$node.centrality)))
)

cent_netH1_helius$dummyB <- "Strength"
cent_netH1_helius |>
  ggplot(aes(x=value, y = node, group=1)) + 
  geom_point(color = "deepskyblue4") + 
  geom_path(color = "deepskyblue4") + 
  theme_bw() +
  labs(x = "", y = "", color = "") +
  theme(text=element_text(size=15))+
  facet_grid(. ~ dummyB)


# all waves together
average_net_allwaves <- dep_list |>
  map(~.x |>select(-ID) |> relocate(sui, .after = mot)
  ) |> list_rbind() |>
  drop_na()

totavgnetwork_helius <- estimateNetwork(average_net_allwaves, default = "EBICglasso")
totavgnetwork_helius2 <- plot(totavgnetwork_helius)
totavgnetwork_helius2$layout <- Layout_simnet
plot(totavgnetwork_helius2)

## bootstrapped samples
# bootstrap sample of helius
helius_boots <- bootnet(totavgnetwork_helius, nBoots=1000, type = 'nonparametric')
plot(helius_boots)
helius_boots_cent <- helius_boots$boots |>
  map(~.$graph |> centrality_auto() |> extract(1) |> 
        as.data.frame() |> select(node.centrality.Strength) |> 
        rename(Strength = node.centrality.Strength) |> t() |> as.data.frame()
      ) |> list_rbind()
# saveRDS(helius_boots, file = "helius_bootstrap.Rds")
# helius_boots <- readRDS("data/helius_bootstrap.Rds")

# bootstrap sample of toy model
sim_boots <- bootnet(totavgnetwork, nBoots=1000, type = 'nonparametric') # using 30%aggdat
sim_boots_cent <- sim_boots$boots |>
  map(~.$graph |> centrality_auto() |> extract(1) |> 
        as.data.frame() |> select(node.centrality.Strength) |> 
        rename(Strength = node.centrality.Strength) |> t() |> as.data.frame()
  ) |> list_rbind()

sim_boots_cent <- readRDS("data/sim_bootstrap1000cent_30%agg.Rds")
sim_boots <- readRDS("data/sim_bootstrap1000_30%agg.Rds")

# sim_boots <- readRDS("data/sim_bootstrap100.Rds")
# sim_boots_cent <- readRDS("data/sim_bootstrap100_cent.Rds")
cent_dat <- bind_rows(helius_boots_cent, sim_boots_cent, .id= "id")
  
cent_stats <- cent_dat |>
  group_by(id) |>
  summarise(across(
    .cols = everything(),
    .fns = list(Mean = mean, SD = sd),
    .names = "{fn}_{col}"
  )) |> pivot_longer(!id, names_to = c(".value", "symptom"), names_sep = "_") |>
  mutate(N = case_when(id == 1 ~ 1000, id ==2 ~ 1000),
         SE = SD / sqrt(N),
         low.ci = Mean - 2*SD, 
         upper.ci = Mean + 2*SD) 
         # low.ci = Mean - qnorm(0.975)*SE, # consider qt() # interval so small
         # upper.ci = Mean + qnorm(0.975)*SE) # consider qt()


## centrality plot over bootstrapped samples
boot_centrailtyPlot <- ggplot(data = cent_stats, 
       aes(x = Mean, y = factor(symptom, levels= c(rownames(res_netH1_helius$node.centrality))), xmin = low.ci, xmax = upper.ci, 
           color = id)) +
  geom_pointrange(alpha=0.5) +
  # geom_crossbar(width=0.1, alpha = 0.1)+
  geom_path(aes(group = id), alpha=0.5)+
  scale_color_discrete(labels = c("1" = "Helius", "2" = "Simulated")) +
  labs(title = "Strength centrality of bootstrapped samples",
       x = "",
       y = "",
       color = "",
       caption = "Error bars show the 2*SD") +
  ggpubr::theme_classic2()

ggsave("boot_centrality_SD.png", boot_centrailtyPlot, width = 5, height = 6)




## Edge weight stability plot over bootstrapped samples
edge_stats <- bind_rows(helius_boots$bootTable, sim_boots$bootTable, .id ="grp") |>
  select(id, value, grp) |>
  group_by(id, grp) |>
  summarise(across(
    .cols = everything(),
    .fns = list(Mean = mean, SD = sd),
    .names = "{fn}")
  ) |>
  filter(stringr::str_detect(id, "--")) |>
  mutate(N = case_when(grp == 1 ~ 1000, grp ==2 ~ 100),
         SE = SD / sqrt(N),
         low.ci = Mean - SD*2, # consider qt() # interval so small
         upper.ci = Mean + SD*2) # consider qt()
         # low.ci = Mean - qnorm(0.975)*SE*100, # consider qt() # interval so small
         # upper.ci = Mean + qnorm(0.975)*SE*100)

# lst_mats <- sim_boots$boots |> map(~.$graph)
# mat_means <- apply(simplify2array(lst_mats), 1:2, mean) |> as.data.frame()
# mat_sds <- apply(simplify2array(lst_mats), 1:2, sd)

edge_stats |> mutate(Mean_grp1 = ifelse(grp == 1, Mean, 0))


boot_edgePlot <- edge_stats |> mutate(Mean_grp1 = ifelse(grp == 1, Mean, 0)) |>
ggplot(aes(x = forcats::fct_reorder(id, desc(Mean_grp1)), y = Mean, ymin = low.ci, ymax = upper.ci, 
                                  color = grp)) +
  # geom_pointrange(alpha=0.5, fatten = 2) +
  geom_errorbar(alpha=0.5, width = 0.4) +
  # geom_crossbar(width=0.2, alpha = 0.5)+
  scale_color_discrete(labels = c("1" = "Helius", "2" = "Simulated")) +
  labs(title = "Edge weight of bootstrapped samples",
       x = "",
       y = "",
       color = "",
       caption = "Error bars show the 2*SD") +
  ggpubr::theme_classic2()+
  theme(axis.text.x = element_text(angle = 80, vjust = 1, hjust=1),
        legend.position = "bottom") 

ggsave("boot_edgePlot_SD.png", boot_edgePlot, width = 12, height = 6)





## plot strength centrality of helius vs. toy model
# weighted adj mat
mat <- round(totavgnetwork_helius$graph, 2)
ifelse(mat < 0.1, 0, mat)
round(totavgnetwork$graph, 2)

res_totavg_helius <- centrality_auto(totavgnetwork_helius2)
rownames(res_totavg_helius$node.centrality) <- colnames(A)


cent_totavg_helius <- res_totavg_helius$node.centrality$Strength |>as_data_frame() |>  mutate(node = factor(rownames(res_totavg_helius$node.centrality), levels= c(rownames(res_totavg_helius$node.centrality)))
)

cent_totavg_helius$dummyB <- "Strength"
cent_totavg_helius_plot <- cent_totavg_helius |>
  ggplot(aes(x=value, y = node, group=1)) + 
  geom_point(color = "deepskyblue4") + 
  geom_path(color = "deepskyblue4") + 
  theme_bw() +
  labs(x = "", y = "", color = "") +
  theme(text=element_text(size=15))+
  facet_grid(. ~ dummyB)

# combine two centralities
cent_totavg$grp <- "toy"
cent_totavg_helius$grp <- "helius"
both_cents <- cent_totavg |> full_join(cent_totavg_helius)
# normalization
maxval <- max(both_cents$value)
minval <- min(both_cents$value)
both_cents$val <- both_cents$value - minval / (maxval - minval)

cent_combined <- both_cents |>
  ggplot(aes(x=value, y = node, group=grp, color = grp)) + 
  geom_point() + 
  geom_path() + 
  theme_bw() +
  labs(x = "", y = "", color = "") +
  theme(text=element_text(size=15))+
  facet_grid(. ~ dummyB)


### ==========================
### Only dep_sum score >= 10
### ==========================

dep_scores_sick <- helius2 |> 
  dplyr::select(contains(c("WlbvRecent","PHQ9_deprsymp")) & !ends_with("Ingevuld"), ID) |>
  pivot_longer(cols = -c(contains("PHQ9_deprsymp"), ID), names_to = "wave1", values_to = "value") |>
  mutate(wave = str_extract(wave1, "[^_]+"), item_num = str_extract(wave1,  "(\\d+)(?!.*\\d)")) |>
  drop_na() |>
  dplyr::select(-wave1) %>% 
  filter(!(item_num == "8" | item_num == "9")) |>
  mutate(item_num = case_when(item_num == "1" ~ "anhedonia",
                              item_num == "2"~ "sadness",
                              item_num == "3"~ "sleep",
                              item_num == "4"~ "energy",
                              item_num == "5"~ "appetite",
                              item_num == "6"~ "guilt",
                              item_num == "7"~ "concentration",
                              item_num == "89" ~ "motor",
                              item_num == "10" ~ "suicidal")
  )


dep_list_sick <- dep_scores_sick %>% 
  group_by(wave) %>% 
  {setNames(group_split(.), group_keys(.)[[1]])}  %>% 
  # select only depresymptom == 1
  imap(~.x |> filter(if_any(starts_with(.y), ~ . == 1))  |> 
         pivot_wider(id_cols = ID, names_from = item_num, values_from = value))


net_list_sick <- dep_list_sick |>
  map(~.x |> dplyr::select(-ID) |>
        estimateNetwork(.x, default = "EBICglasso") |>
        plot(.x)
  )

plot(net_list_sick$H1)
plot(net_list_sick$H2)

# ordinal correction (polychoric corr)
h1sicknet <- qgraph(cor_auto(dep_list_sick$H1[,-1]), graph="glasso", sampleSize=nrow(dep_list_sick$H1), tuning=0.5,layout="spring", details=TRUE, theme = "colorblind")

h2sicknet <- qgraph(cor_auto(dep_list_sick$H2[,-1]), graph="glasso", sampleSize=nrow(dep_list_sick$H2), tuning=0.5, details=TRUE, theme = "colorblind", layout=h1sicknet$layout)

cov1sicknet <- qgraph(cor_auto(dep_list_sick$Cov1[,-1]), graph="glasso", sampleSize=nrow(dep_list_sick$Cov1), tuning=0.5, details=TRUE, theme = "colorblind", layout=h1sicknet$layout)

cov2sicknet <- qgraph(cor_auto(dep_list_sick$Cov2[,-1]), graph="glasso", sampleSize=nrow(dep_list_sick$Cov2), tuning=0.5, details=TRUE, theme = "colorblind", layout=h1sicknet$layout)

cov3sicknet <- qgraph(cor_auto(dep_list_sick$Cov3[,-1]), graph="glasso", sampleSize=nrow(dep_list_sick$Cov3), tuning=0.5, details=TRUE, theme = "colorblind", layout=h1sicknet$layout)






## =============================================================================
## plotAG function (slightly modified from the pcalg package)
## =============================================================================
#' Plot PAG (partial ancestral graph) for FCI and CCI algorithms
#'
#' @param amat adjacency matrix of the resulting estimated graph
#'
#' @details
#' "0": no edge; "1": circle; "2": arrow; "3": tail
#'
#' @return a PAG graph of graphNEL class
plotAG <- function(amat)
{
  # create a graph object
  g <- as(amat,"graphNEL")
  # extract node info
  nn <- nodes(g)
  # extract number of nodes
  p <- numNodes(g)
  # extract number of edges
  n.edges <- numEdges(g)
  ah.list <- at.list <- vector("list", n.edges)
  l.names <- character(n.edges)
  # assign a shape for each edge type
  amat[amat == 1] <- "odot"
  amat[amat == 2] <- "vee"
  amat[amat == 3] <- "none"
  iE <- 0
  for (i in 1:(p-1)) {
    x <- nn[i]
    for (j in (i+1):p) {
      y <- nn[j]
      if (amat[x,y] != 0) {
        iE <- iE + 1
        ah.list[[iE]] <- amat[x,y]
        at.list[[iE]] <- amat[y,x]
        l.names[[iE]] <- paste0(x,"~",y)
      }
    }
  }
  names(ah.list) <- names(at.list) <- l.names
  
  edgeRenderInfo(g) <- list(arrowhead = ah.list, arrowtail = at.list)
  # global features
  graph.par(list(nodes=list(cex = 1)))
  # plot the PAG
  Rgraphviz::renderGraph(Rgraphviz::layoutGraph(g))
}

# CCI algorithm
library(devtools)
install_github("ericstrobl/CCI")
library(CCI)
library(Rgraphviz)


sampledat <- dep_list_sick$H1 |>
  select(-ID) |>
  drop_na()

suffStat=list(); suffStat$C = cor(sampledat); suffStat$n = nrow(sampledat); # get all of the parameters needed by Fisher's z test

G=cci(suffStat,gaussCItest,alpha=0.01,p=ncol(sampledat)) # run CCI

G$maag #print the recovered partially oriented MAAG
plotAG(G$maag)

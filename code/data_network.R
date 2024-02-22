library(haven)
library(purrr)
library(dplyr)
library(stringr)

helius2 <- read_sav("data/HELIUS_itemscores.sav")
helius2$H2_WlbvRecent10

dep_scores <- helius2 |> 
  select(contains("WlbvRecent") & !ends_with("Ingevuld"), ID) |>
  pivot_longer(-ID, names_to = "wave1", values_to = "value") |>
  mutate(wave = str_extract(wave1, "[^_]+"), item_num = str_extract(wave1,  "(\\d+)(?!.*\\d)")) |>
  drop_na() |>
  select(-wave1) %>% 
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

dep_list <- dep_scores %>% 
  group_by(wave) %>% 
  {setNames(group_split(.), group_keys(.)[[1]])}  %>% 
  map(~.x |>pivot_wider(id_cols = ID, names_from = item_num, values_from = value))
  

net_list <- dep_list |>
  map(~.x |>select(-ID) |>
      estimateNetwork(.x, default = "EBICglasso") |>
        plot(.x)
      )

plot(net_list$Cov1)

# H1net <- estimateNetwork(H1, default = "EBICglasso")
# plot(H1net)


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

library(rootSolve)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(deSolve)
library(ggpubr)

## bifurcation plot per symptom


# anhedonia
b <- -8; a <- 0.3; W <- 0.97; d <- 9
# sleep
b <- -8; a <- 0.3; W <- 0.66; d <- 9
# guilty
b <- -8; a <- 0.3; W <- 0.43; d <- 9
# suicide
b <- -8; a <- 0.3; W <- 0.44; d <- 9

# S rate
rate <- function(S, b = 1) S*(1 - S)*(b + a*S + W*(1 + d*(S^2)))

# Stability of a root ~ sign of eigenvalue of Jacobian

stability<-function(b) {
  Eq <- uniroot.all(rate, c(0,1), b = b)
  eig <-vector()
  for(i in 1:length(Eq))
    eig[i]<-sign(gradient(rate,Eq[i],b = b))
  return(list(Eq=Eq,Eigen=eig)) }

# bifurcation diagram 
bseq <- seq(-13, 2, by=0.01)

# plot bifurcation
anh_lst <- bseq |>
  map_df(~ stability(.x), .id = "id") |>
  mutate(beta = as.factor(bseq[as.numeric(id)]),
         Eigen = ifelse(Eigen == 1, "unstable", "stable")) 

anh_stable <- anh_lst |> filter(Eigen =="stable")
anh_unstable <- anh_lst |> filter(Eigen =="unstable")

slp_lst <- bseq |>
  map_df(~ stability(.x), .id = "id") |>
  mutate(beta = as.factor(bseq[as.numeric(id)]),
         Eigen = ifelse(Eigen == 1, "unstable", "stable")) 

slp_stable <- slp_lst |> filter(Eigen =="stable")
slp_unstable <- slp_lst |> filter(Eigen =="unstable")

glt_lst <- bseq |>
  map_df(~ stability(.x), .id = "id") |>
  mutate(beta = as.factor(bseq[as.numeric(id)]),
         Eigen = ifelse(Eigen == 1, "unstable", "stable")) 

glt_stable <- glt_lst |> filter(Eigen =="stable")
glt_unstable <- glt_lst |> filter(Eigen =="unstable")

sui_lst <- bseq |>
  map_df(~ stability(.x), .id = "id") |>
  mutate(beta = as.factor(bseq[as.numeric(id)]),
         Eigen = ifelse(Eigen == 1, "unstable", "stable")) 

sui_stable <- sui_lst |> filter(Eigen =="stable")
sui_unstable <- sui_lst |> filter(Eigen =="unstable")

# find x-intercept
# lst |> filter(beta == -0.6)


# bif_anh <- anh_lst |> ggplot(aes(x = beta, y = Eq, color = Eigen)) +
#   geom_point(size = 0.1) +
#   theme_classic() +
#   scale_color_manual(values=c("royalblue", "white"), labels = c("stable", "unstable")) +
#   geom_line(data= anh_unstable, aes(x = beta, y = Eq), group=1, linetype=2,linewidth=1, color = "royalblue") +
#   annotate('rect', xmin = 301, xmax = 1204, ymin =0, ymax= 1, alpha=0.1, fill="gold1")+
#   geom_vline(xintercept = 1106, linetype="dotted", color = "orange", linewidth=1.5) +
#   geom_vline(xintercept = 1210, linetype="dotted", color = "coral3", linewidth=1.5) +
#   scale_y_continuous(breaks=c(0,1))+
#   labs(title = "(a) anhedonia", y = "", x="")+
#   guides(color="none") +
#   theme(axis.text.x=element_blank(), 
#         axis.ticks.x=element_blank(),
#         plot.title=element_text(size = 20), 
#         axis.text.y=element_text(size = 20))
# #plot.margin = margin(0.5, 0.1, 0.5, 0.1, "cm"))
# 
# bif_slp <- slp_lst |> ggplot(aes(x = beta, y = Eq, color = Eigen)) +
#   geom_point(size = 0.1) +
#   theme_classic() +
#   scale_color_manual(values=c("royalblue", "white"), labels = c("stable", "unstable")) +
#   geom_line(data= slp_unstable, aes(x = beta, y = Eq), group=1, linetype=2,linewidth=1, color = "royalblue") +
#   annotate('rect', xmin = 611, xmax = 1235, ymin =0, ymax= 1, alpha=0.1, fill="gold1")+
#   geom_vline(xintercept = 1166, linetype="dotted", color = "orange", linewidth=1.5) +
#   geom_vline(xintercept = 1241, linetype="dotted", color = "coral3", linewidth=1.5) +
#   scale_y_continuous(breaks=c(0,1))+
#   labs(title = "(c) sleep", y = "", x= "")+
#   guides(color="none") +
#   theme(axis.text.x=element_blank(), 
#         axis.ticks.x=element_blank(),
#         plot.title=element_text(size = 20), 
#         axis.text.y=element_text(size = 20))
#         #plot.margin = margin(0.5, 0.1, 0.5, 0.1, "cm"))


# bif_glt <- glt_lst |> ggplot(aes(x = beta, y = Eq, color = Eigen)) +
#   geom_point(size = 0.1) +
#   theme_classic() +
#   scale_color_manual(values=c("royalblue", "white"), labels = c("stable", "unstable")) +
#   geom_line(data= glt_unstable, aes(x = beta, y = Eq), group=1, linetype=2,linewidth=1, color = "royalblue") +
#   annotate('rect', xmin = 841, xmax = 1258, ymin =0, ymax= 1, alpha=0.1, fill="gold1")+
#   geom_vline(xintercept = 1215, linetype="dotted", color = "orange", linewidth=1.5) +
#   geom_vline(xintercept = 1265, linetype="dotted", color = "coral3", linewidth=1.5) +
#   scale_y_continuous(breaks=c(0,1))+
#   labs(title = "(b) guilty", y = "", x="")+
#   guides(color="none") +
#   theme(axis.text.x=element_blank(), 
#         axis.ticks.x=element_blank(),
#         plot.title=element_text(size = 20), 
#         axis.text.y=element_text(size = 20))
#         #plot.margin = margin(0.5, 0.1, 0.5, 0.1, "cm"))

bif_anh <- anh_lst |> ggplot(aes(x = beta, y = Eq)) +
  geom_point(aes(color = Eigen), size = 0.1) +
  geom_line(data= anh_unstable, aes(x = beta, y = Eq), group=1, linetype=2,linewidth=1, color = "royalblue") +
  annotate('rect', xmin = 301, xmax = 1204, ymin =0, ymax= 1, alpha=0.1, fill="gold1")+
  geom_vline(aes(xintercept = 1106, color = "beta"),linetype="dotted", linewidth=1.5, key_glyph = "path") +
  geom_vline(aes(xintercept = 1210, color = "beta+tau"),linetype="dotted",linewidth=1.5, key_glyph = "path") +
  scale_color_manual(values=c("beta" = "orange", "beta+tau"="coral3", "stable"="royalblue", "unstable"="white"), breaks = c("beta", "beta+tau"), labels = c(expression(beta["bistable"]), expression(beta["shock"]))) +
  guides(color = guide_legend(override.aes = list(linewidth = 2)))+
  scale_y_continuous(breaks=c(0,1))+
  theme_classic() +
  labs(title = "(a) anhedonia", y = "S", color = "Bifurcation", x= expression(beta))+
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(hjust = 1.1, vjust = 3.5, size = 20), 
        axis.title.y=element_text(hjust = 0, vjust = 1.05, size = 20, angle = 360),
        plot.margin = margin(0.5, 1.5, 0, 0.1, "cm"),
        legend.text.align = 0,
        legend.key.width = unit(4, "line"),
        legend.spacing.x = unit(1.0, 'cm'),
        plot.title=element_text(size = 20, hjust = 0.5), 
        axis.text.y=element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        legend.text = element_text(size = 22),
        legend.position = "bottom")


bif_slp <- slp_lst |> ggplot(aes(x = beta, y = Eq)) +
  geom_point(aes(color = Eigen), size = 0.1) +
  geom_line(data= slp_unstable, aes(x = beta, y = Eq), group=1, linetype=2,linewidth=1, color = "royalblue") +
  annotate('rect', xmin = 611, xmax = 1235, ymin =0, ymax= 1, alpha=0.1, fill="gold1")+
  geom_vline(aes(xintercept = 1166, color = "beta"),linetype="dotted", linewidth=1.5, key_glyph = "path") +
  geom_vline(aes(xintercept = 1241, color = "beta+tau"),linetype="dotted",linewidth=1.5, key_glyph = "path") +
  scale_color_manual(values=c("beta" = "orange", "beta+tau"="coral3", "stable"="royalblue", "unstable"="white"), breaks = c("beta", "beta+tau"), labels = c(expression(beta["bistable"]), expression(beta["shock"]))) +
  guides(color = guide_legend(override.aes = list(linewidth = 2)))+
  scale_y_continuous(breaks=c(0,1))+
  theme_classic() +
  labs(title = "(c) sleep", y = "S", color = "Bifurcation", x= expression(beta))+
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(hjust = 1.1, vjust = 3.5, size = 20), 
        axis.title.y=element_text(hjust = 0, vjust = 1.05, size = 20, angle = 360),
        plot.margin = margin(0.5, 1.5, 0, 0.1, "cm"),
        legend.text.align = 0,
        legend.key.width = unit(4, "line"),
        legend.spacing.x = unit(1.0, 'cm'),
        plot.title=element_text(size = 20, hjust = 0.5), 
        axis.text.y=element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        legend.text = element_text(size = 22),
        legend.position = "bottom")


bif_glt <- glt_lst |> ggplot(aes(x = beta, y = Eq)) +
  geom_point(aes(color = Eigen), size = 0.1) +
  geom_line(data= glt_unstable, aes(x = beta, y = Eq), group=1, linetype=2,linewidth=1, color = "royalblue") +
  annotate('rect', xmin = 841, xmax = 1258, ymin =0, ymax= 1, alpha=0.1, fill="gold1")+
  geom_vline(aes(xintercept = 1215, color = "beta"),linetype="dotted", linewidth=1.5, key_glyph = "path") +
  geom_vline(aes(xintercept = 1265, color = "beta+tau"),linetype="dotted",linewidth=1.5, key_glyph = "path") +
  scale_color_manual(values=c("beta" = "orange", "beta+tau"="coral3", "stable"="royalblue", "unstable"="white"), breaks = c("beta", "beta+tau"), labels = c(expression(beta["bistable"]), expression(beta["shock"]))) +
  guides(color = guide_legend(override.aes = list(linewidth = 2)))+
  scale_y_continuous(breaks=c(0,1))+
  theme_classic() +
  labs(title = "(b) guilty", y = "S", color = "Bifurcation", x= expression(beta))+
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(hjust = 1.1, vjust = 3.5, size = 20), 
        axis.title.y=element_text(hjust = 0, vjust = 1.05, size = 20, angle = 360),
        plot.margin = margin(0.5, 1.5, 0, 0.1, "cm"),
        legend.text.align = 0,
        legend.key.width = unit(4, "line"),
        legend.spacing.x = unit(1.0, 'cm'),
        plot.title=element_text(size = 20, hjust = 0.5), 
        axis.text.y=element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        legend.text = element_text(size = 22),
        legend.position = "bottom")

  
bif_sui <- sui_lst |> ggplot(aes(x = beta, y = Eq)) +
  geom_point(aes(color = Eigen), size = 0.1) +
  #scale_color_manual(values=c("royalblue", "white"), labels = c("stable", "unstable")) +
  geom_line(data= sui_unstable, aes(x = beta, y = Eq), group=1, linetype=2,linewidth=1, color = "royalblue") +
  annotate('rect', xmin = 841, xmax = 1257, ymin =0, ymax= 1, alpha=0.1, fill="gold1")+
  geom_vline(aes(xintercept = 1080, color = "beta"),linetype="dotted", linewidth=1.5, key_glyph = "path") +
  geom_vline(aes(xintercept = 1130, color = "beta+tau"),linetype="dotted",linewidth=1.5, key_glyph = "path") +
  scale_color_manual(values=c("beta" = "orange", "beta+tau"="coral3", "stable"="royalblue", "unstable"="white"), breaks = c("beta", "beta+tau"), labels = c(expression(beta["bistable"]), expression(beta["shock"]))) +
  guides(color = guide_legend(override.aes = list(linewidth = 2)))+
  scale_y_continuous(breaks=c(0,1))+
  theme_classic() +
  labs(title = "(d) suicidal", y = "S", color = "Bifurcation", x= expression(beta))+
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(hjust = 1.1, vjust = 3.5, size = 20), 
        axis.title.y=element_text(hjust = 0, vjust = 1.05, size = 20, angle = 360),
        plot.margin = margin(0.5, 1.5, 0, 0.1, "cm"),
        legend.text.align = 0,
        legend.key.width = unit(4, "line"),
        legend.spacing.x = unit(1.0, 'cm'),
        plot.title=element_text(size = 20, hjust = 0.5), 
        axis.text.y=element_text(size = 20),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 17),
        legend.position = "bottom")


x <- seq(-2, 2, by = 0.01)
anh_bi = x^2*(582*x^3 - 615*x^2- 440*x + 455) / 300
anh_sh = x^2*(194*x^3 - 205*x^2- 40*x -15) / 100

slp_bi = x^2*(1608*x^3 - 1560*x^2- 1270*x + 1005) / 1200
slp_sh = x^2*(804*x^3 - 780*x^2- 230*x -105) / 600

glt_bi = x^2*(1032*x^3 - 840*x^2- 890*x + 435) / 1200
glt_sh = x^2*(86*x^3 - 70*x^2- 40*x -15) / 100

sui_bi = x^2*(804*x^3 - 780*x^2- 1630*x + 1990) / 600
sui_sh = x^2*(200*x^3 - 175*x^2- 310*x + 315) / 200

df_landscape <- data.frame(x, anh_bi, anh_sh, slp_bi, slp_sh, glt_bi, glt_sh, sui_bi, sui_sh) |>
  pivot_longer(!x, values_to = "val", names_to = "symptom")


# anh landscape
pointdat_anh <- data.frame(x = 1, y = 0.0075)

anh_ls <- df_landscape |>filter(grepl("anh", symptom)) |>
  ggplot() +
  geom_line(aes(x = x, y = val, col = symptom), linewidth = 1) +
  scale_color_manual(values= c("orange","coral3"), labels = c("bistable", "shock"))+
  geom_point(pointdat_anh, mapping = aes(x = x, y = y),  fill = "gray", size = 5, alpha = 0.2) +
  theme_classic()+
  ylim(-0.66,1) + 
  labs(color = "Landscape", y ="V(s)", x = "S") +
  scale_x_continuous(breaks=c(0,1), limits =c(-0.3, 1.5))+
  theme(axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.x=element_text(hjust = 1.1, vjust = 9, size = 20), 
        axis.title.y=element_text(hjust = 1, vjust = 1, size = 20, angle = 360),
        axis.text.x=element_text(size = 20),
        plot.margin = margin(0, 1.5, 0.5, 1.2, "cm"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18, face ="bold"),
        legend.key.width = unit(4, "line"),
        legend.position='bottom')


# slp landscape
pointdat_slp <- data.frame(x = 1, y = -0.122)

slp_ls <- df_landscape |>filter(grepl("slp", symptom)) |>
  ggplot() +
  geom_line(aes(x = x, y = val, col = symptom), linewidth = 1) +
  scale_color_manual(values= c("orange","coral3"), labels = c("bistable", "shock"))+
  geom_point(pointdat_slp, mapping = aes(x = x, y = y),  fill = "gray", size = 5, alpha = 0.2) +
  theme_classic()+
  ylim(-0.52,1) + 
  labs(color = "Landscape", y ="V(s)", x = "S") +
  scale_x_continuous(breaks=c(0,1), limits =c(-0.3, 1.5))+
  theme(axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.x=element_text(hjust = 1.1, vjust = 9, size = 20), 
        axis.title.y=element_text(hjust = 1, vjust = 1, size = 20, angle = 360),
        axis.text.x=element_text(size = 20),
        plot.margin = margin(0, 1.5, 0.5, 1.2, "cm"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18, face ="bold"),
        legend.key.width = unit(4, "line"),
        legend.position='bottom')

# glt landscape
pointdat_glt <- data.frame(x = 1.0001, y = -0.1645)

glt_ls <- df_landscape |>filter(grepl("glt", symptom)) |>
  ggplot() +
  geom_line(aes(x = x, y = val, col = symptom), linewidth = 1) +
  scale_color_manual(values= c("orange","coral3"), labels = c("bistable", "shock"))+
  geom_point(pointdat_glt, mapping = aes(x = x, y = y),  fill = "gray", size = 5, alpha = 0.2) +
  theme_classic()+
  ylim(-0.4,1) + 
  labs(color = "Landscape", y ="V(s)", x = "S") +
  scale_x_continuous(breaks=c(0,1), limits =c(-0.3, 1.5))+
  theme(axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.x=element_text(hjust = 1.1, vjust = 9, size = 20), 
        axis.title.y=element_text(hjust = 1, vjust = 1, size = 20, angle = 360),
        axis.text.x=element_text(size = 20),
        plot.margin = margin(0, 1.5, 0.5, 1.2, "cm"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18, face ="bold"),
        legend.key.width = unit(4, "line"),
        legend.position='bottom')

# sui landscape
pointdat_sui <- data.frame(x = 0, y = 0.0405)

sui_ls <- df_landscape |>filter(grepl("sui", symptom)) |>
  ggplot() +
  geom_line(aes(x = x, y = val, col = symptom), linewidth = 1) +
  scale_color_manual(values= c("orange","coral3"), labels = c("pre-shock", "shock"))+
  theme_classic()+
  geom_point(pointdat_sui, mapping = aes(x = x, y = y),  fill = "gray", size = 5, alpha = 0.2) +
  ylim(0,1) + 
  labs(color = "Landscape", y ="V(s)", x = "S") +
  scale_x_continuous(breaks=c(0,1), limits =c(-0.3, 1.5))+
  theme(axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.x=element_text(hjust = 1.1, vjust = 9, size = 20), 
        axis.title.y=element_text(hjust = 1, vjust = 1, size = 20, angle = 360),
        axis.text.x=element_text(size = 20),
        plot.margin = margin(0, 1.5, 0.5, 1.2, "cm"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.key.width = unit(4, "line"),
        legend.position='bottom')

legend1 <- ggpubr::get_legend(bif_sui)
legend2 <- ggpubr::get_legend(sui_ls)
legends <- ggarrange(legend2, NULL, legend1, NULL, nrow = 4)

combplot <- ggpubr::ggarrange(bif_anh, bif_glt, anh_ls, glt_ls, bif_slp,bif_sui,slp_ls,sui_ls,
                              #labels = c("A", "B", "C", "D"),
                              ncol = 2, nrow = 4, common.legend = T, legend ="none",
                              heights = c(1, 1.2, 1, 1.2, 1, 1.2, 1, 1),align = "v")

combplot_leg <- ggarrange(combplot, legends, ncol=1, heights = c(10, 1), align = "hv")

ggsave("symptoms_bif1.pdf", plot = combplot_leg, units="in", width=12, height=14, bg = 'transparent')



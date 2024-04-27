## =========================================================
## Bifurcation diagrams
##
## creating Figure 4 and Figure 5 
## =========================================================
## install packages
source("code/libraries.R")


## Create Figure 5 : Bifurcation w.r.t. omega W

# setting 1: 1 1
b <- -8; a <- 1; W <- 1; d <- 1
# setting 2: 1 0
b <- -8; a <- 1; W <- 1; d <- 0
# setting 3: 3 0
b <- -8; a <- 3; W <- 1; d <- 0
# setting 4: 1 3
b <- -8; a <- 1; W <- 1; d <- 3

# S rate
rate <- function(S, W = 3) S*(1 - S)*(b + a*S + W*(1 + d*(S^2)))  

# stability of a root ~ sign of eigenvalue of Jacobian
stability<-function(W) { 
  Eq <- uniroot.all(rate, c(0,1), W = W) 
  eig <-vector() 
  for(i in 1:length(Eq)) 
    eig[i]<-sign(gradient(rate,Eq[i], W = W)) 
  return(list(Eq=Eq,Eigen=eig)) } 


# omega range
wseq <-  seq(0, 10, by=0.01) 

# get the list of stability for each condition
lst1 <- wseq |> 
  map_df(~ stability(.x), .id = "id") |>
  mutate(omega = as.factor(wseq[as.numeric(id)]),
         Eigen = ifelse(Eigen == 1, "unstable", "stable"))

stable1 <- lst |> filter(Eigen =="stable")
unstable1 <- lst |> filter(Eigen =="unstable")

lst2 <- wseq |> 
  map_df(~ stability(.x), .id = "id") |>
  mutate(omega = as.factor(wseq[as.numeric(id)]),
         Eigen = ifelse(Eigen == 1, "unstable", "stable"))

stable2 <- lst |> filter(Eigen =="stable")
unstable2 <- lst |> filter(Eigen =="unstable")

lst3 <- wseq |> 
  map_df(~ stability(.x), .id = "id") |>
  mutate(omega = as.factor(wseq[as.numeric(id)]),
         Eigen = ifelse(Eigen == 1, "unstable", "stable"))

stable3 <- lst |> filter(Eigen =="stable")
unstable3 <- lst |> filter(Eigen =="unstable")

lst4 <- wseq |> 
  map_df(~ stability(.x), .id = "id") |>
  mutate(omega = as.factor(wseq[as.numeric(id)]),
         Eigen = ifelse(Eigen == 1, "unstable", "stable"))

stable4 <- lst |> filter(Eigen =="stable")
unstable4 <- lst |> filter(Eigen =="unstable")

big_lst <- bind_rows(lst2, lst3, lst4, .id = "sets") |>
  mutate(omega = as.numeric(omega))

big_lst_stable <- big_lst |> filter(Eigen == "stable")
big_lst_unstable <- big_lst |> filter(Eigen == "unstable")

# plot bifurcation diagram w.r.t. omega W
bif_W <- big_lst_stable |> ggplot(aes(x = omega, y = Eq)) + 
  geom_point(aes(color = sets)) +
  theme_classic() +
  scale_color_manual(name = "Condition", values= alpha(c("royalblue", "darkgreen", "darkorange"), .5), labels= c(expression(alpha[ii]~" low, "~delta[i]~ " low", alpha[ii]~" high, "~delta[i]~ " low", alpha[ii]~" low, "~delta[i]~ " high"))) +
  geom_line(data= big_lst_unstable, aes(x = omega, y = Eq, group = sets, color = sets), linetype = 2, linewidth=2) +
  geom_point(aes(x=800, y=0), shape = 8, color = "tomato3", size = 8, stroke=2)+
  geom_point(aes(x=700, y = 1), shape = 8, color = "tomato3", size = 8,stroke=2) +
  geom_point(aes(x=5*100, y=1), shape = 8, color = "tomato3", size = 8, stroke=2)+
  geom_point(aes(x=7/4*100, y = 1), shape = 8, color = "tomato3", size = 8,stroke=2)+
  labs(x = expression(omega[i]), y = expression("S"[i]^"*"))+
  scale_y_continuous(breaks=c(0,1))+
  scale_x_continuous(breaks = c(7/4*100, 800), labels =c(expression(frac(-(beta[i] + alpha[ii]), 1 +delta[i])), expression(-beta[i]))) + 
  guides(color = guide_legend(override.aes = list(size = 5, linewidth = 5))) +
  theme(axis.text.x=element_text(size = 35), 
        axis.text.y=element_text(size = 35),
        legend.title = element_text(size=30),
        legend.text = element_text(size = 30),
        #legend.position = 'right',
        legend.text.align = 0,
        axis.title.x = element_text(hjust = 1.05, vjust = 13, size= 40),
        axis.title.y = element_text(hjust = -1, vjust = 1.11, angle = 0, size = 40),
        plot.margin = margin(4, 4, 4, 5, "cm"),
        axis.ticks.length=unit(1, "cm"))

# ggsave("bifurcation_W.pdf", plot = bif_W, units="in", width=25, height=12, bg = 'transparent')


## Create Figure 4: Bifurcation diagram w.r.t. beta 
b <- 1; a <- 1; W <- 1; d <- 3

rate <- function(S, b = 1) S*(1 - S)*(b + a*S + W*(1 + d*(S^2)))

stability<-function(b) {
  Eq <- uniroot.all(rate, c(0,1), b = b)
  eig <-vector()
  for(i in 1:length(Eq))
    eig[i]<-sign(gradient(rate,Eq[i],b = b))
  return(list(Eq=Eq,Eigen=eig)) }

bseq <- seq(-8, 2, by=0.01)

lst <- bseq |>
  map_df(~ stability(.x), .id = "id") |>
  mutate(beta = as.factor(bseq[as.numeric(id)]),
         Eigen = ifelse(Eigen == 1, "unstable", "stable"))
stable <- lst |> filter(Eigen =="stable")
unstable <- lst |> filter(Eigen =="unstable")


bif <- lst |> ggplot(aes(x = beta, y = Eq, color = Eigen)) +
  geom_point(size = 0.3) +
  theme_classic() +
  scale_color_manual(values=c("royalblue", "white"), labels = c("stable", "unstable")) +
  geom_line(data= unstable, aes(x = beta, y = Eq), group=1, linetype=2,linewidth=2, color = "royalblue") +
  geom_point(aes(x=300, y=1), shape = 8, color = "red", size = 8, stroke=2)+
  geom_point(aes(x=700, y = 0), shape = 8, color = "red", size = 8,stroke=2)+
  # annotate('rect', xmin = 205, xmax = 700, ymin =0, ymax= 1, fill = rgb(255,240,220, alpha = 50, maxColorValue=255))+
  annotate('rect', xmin = 300, xmax = 700, ymin =0, ymax= 1, alpha=0.1, fill="orange")+
  labs(x = expression(beta[i]), y = expression("S"[i]^"*"))+
  scale_y_continuous(breaks=c(0,1))+
  scale_x_discrete(breaks = c("-5", "-1"), labels = c(expression(-omega[i] *(1+delta[i]) - alpha[ii]), expression(-omega[i]))) + 
  #scale_x_discrete(breaks = c(205, 700), labels = c(expression(-omega *(delta+1) - alpha), expression(-omega))) + 
  theme(axis.text.x=element_text(size = 40), 
        axis.text.y=element_text(size = 40),
        #legend.title = element_blank(),
        legend.position = 'none',
        axis.title.x = element_text(hjust=1.03, vjust = 9, size= 45),
        axis.title.y = element_text(hjust= -1, vjust = 1.05, angle = 0, size = 45),
        plot.margin = margin(4, 4, 4, 5, "cm"),
        axis.ticks.length=unit(1, "cm")
        # plot.background = element_rect(
        #   fill = "white",
        #   colour = "white")
        ) 
  #guides(color = guide_legend(override.aes = list(shape = NA, linetype  = c(1,2), size = 5, color = "royalblue")))

png(file="bifurcation.png", units="in", width=20, height=10, res=300, bg = 'transparent')
bif
dev.off()



## phase plots
model<-function(t,state,parms)
{ with(as.list(c(state,parms)), 
       { 
         dS <- S*(1 - S)*(b + a*S + W*(1 + d*(S^2)))   
         return(list(dS)) 
       })
} 

# converging to 1
perms <- c(a = 1, b = -1.6, W = 0.5, d = 0.4)
# bistability
perms <- c(a = 1, b = -1, W = 0.5, d = 0.4)
# converging to 0
perms <- c(a = 1, b = -0.54, W = 0.5, d = 0.4)


# converging to 1
perms <- c(a = 1, b = -0.53, W = 0.5, d = 0.4)
# bistability
perms <- c(a = 1, b = -1, W = 0.5, d = 0.4)
# converging to 0
perms <- c(a = 1, b = -1.9, W = 0.8, d = 0.2)

times <-seq(0, 100,by= 0.01)
# converging to 1
possibleS <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3,0.4, 0.6, 0.8, 1)
# bistability
possibleS <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
# converging to 0
possibleS <- c(0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1)

out <- list()
for(i in 1:length(possibleS)){
  state <- c(S=possibleS[i])
  out[[i]] <-ode(y= state, times = times, func= model, parms = perms)[,2]
}

# result1 <- do.call(cbind, out) |> as.data.frame() |> rename_with(~paste0("S_", possibleS)) |> mutate(time = times) |> pivot_longer(!time, names_to = "S", values_to = "value") 

result1 <- do.call(cbind, out) |> as.data.frame() |> 
  rename_with(~paste0(possibleS)) |> 
  mutate(time = times) |> pivot_longer(!time, names_to = "S", values_to = "value") |>
  mutate(S = as.numeric(S))


par(bg =  rgb(255,240,220, alpha = 25, maxColorValue=255))
plot(x = rep(0,length(possibleS)), possibleS, col = cc, xlim = c(0, 100), pch = 19)
for(i in 1:length(possibleS)){
  lines(result1[result1$S == possibleS[i],]$time, result1[result1$S == possibleS[i],]$value, type="l", col = cc[i])
}


cc <- scales::seq_gradient_pal("lightskyblue1", "navy")(seq(0,1,length.out=11))
points <- data.frame(possibleS) |> mutate(time = 0, S = possibleS) 

p1 <- result1 |>
  ggplot() +
  geom_line(aes(x = time, y = value, color = S)) + 
  scale_color_manual(values=cc) +
  geom_point(data = points, aes(x= time, y = S), color = cc, size=2) +
  geom_point(aes(x=1, y = -0.1), color = NA) +
  geom_point(aes(x=1, y = 1.1), color = NA) +
  theme_classic() +
  labs(x = "t", y ="") +
  theme(legend.position = "none",
        axis.text.x=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_text(size = 15),
        panel.background = element_rect(fill = alpha("orange", .1))
        )

png(file="convergebi.png", units="in", width=5, height=5, res=300, bg = 'transparent')
p1
dev.off()




## landscape plots

# converging to 0
perms <- c(a = 1, b = -2, W = 0.5, d = 0.5)
# bistability
perms <- c(a = 1, b = -1, W = 0.5, d = 0.5)
# converging to 1
perms <- c(a = 1, b = -0.2, W = 0.5, d = 0.5)

times <-seq(0, 500, by= 0.01)
possibleS <- seq(0, 1, 0.1)
out <- list()
for(i in 1:length(possibleS)){
  state <- c(S=possibleS[i])
  out[[i]] <-ode(y= state, times = times, func= model, parms = perms)[,2]
}

result2 <- do.call(cbind, out) |> as.data.frame() |> rename_with(~paste0("S_", possibleS)) |> mutate(time = times) |> pivot_longer(!time, names_to = "S", values_to = "value") |> filter(time > 10) 





y <- c(rep(5, 650), rep(4.5, 5), rep(4,  330), rep(3.5,  340), rep(3,  350), rep(2.5,  360), rep(2,  370), rep(1.5,  380), rep(1,  390), rep(0.5, 400))

#y <- c(rep(-0.5, 300), rep(1.5,  200), rep(2, 100), rep(2.5,  140), rep(3,  140), rep(3.5,  150), rep(4,  150), rep(4.5,  150), rep(5,  300))

x <- seq(0, 5, length.out= length(y))
dat <- data.frame(x, y)

ball1 <- ggplot(dat) +
  geom_density(aes(x = y),fill = "gray", color = "black", alpha = .4, bw =0.3) +
  expand_limits(y = 0) +
  theme_void() +
  #geom_point(aes(x=0.57, y=0.0155),color = "lightblue", alpha = .5, size = 5)  + 
  coord_cartesian(xlim =c(1.1,4.7))

png(file="ball1.png", units="in", width=4, height=4, res=300, bg = 'transparent')
ball1
dev.off()


y <- c(rep(-0.5, 400), rep(1.5,  90), rep(2, 95), rep(2.5,  160), rep(3,  170), rep(3.5,  180), rep(4,  185), rep(4.5,  190), rep(5,  191), rep(5.5, 192), rep(6, 193), rep(6.5,194), rep(7, 195), rep(7.5, 196),rep(8, 197),rep(8.5, 198),rep(9, 199),rep(9.5, 200),rep(10, 200))

x <- seq(0, 10, length.out= length(y))
dat <- data.frame(x, y)


ball0 <- ggplot(dat) +
  geom_density(aes(x = y),fill = "gray", color = "black", alpha = .4, bw = 0.4) +
  expand_limits(y = 0) +
  theme_void() +
  coord_cartesian(xlim =c(-0.2, 8.8)) 

# ball0 <- ggplot(dat) +
#   geom_density(aes(x = y),fill = "gray", color = "black", alpha = .4, bw = 0.4) +
#   expand_limits(y = 0) +
#   theme_void() +
#   #geom_point(aes(x=1.105, y=0.0795),color = "lightblue", alpha = .5, size = 5) +
#   #coord_cartesian(xlim =c(0.5,4))

png(file="ball0.png", units="in", width=6, height=4, res=300, bg = 'transparent')
ball0
dev.off()




y <- c(rep(0, 500),rep(0.5, 30), rep(1, 100), rep(1.5, 150), rep(2, 250), rep(2.5, 350), rep(3, 250), rep(3.5, 150), rep(4, 100), rep(4.5,  30), rep(5, 500))
hist(y)
x <- seq(0, 5, length.out= length(y))
dat <- data.frame(x, y)

biball <- ggplot(dat) +
  geom_density(aes(x = y),fill = "gray", color = "black", alpha = .4) +
  expand_limits(y = 0) +
  theme_void() +
  #geom_point(aes(x=1.105, y=0.0795),color = "lightblue", alpha = .5, size = 5) +
  coord_cartesian(xlim =c(0.3,4.7), ylim = c(0.05, 0.3)) 
  # theme(
  #   panel.background = element_rect(fill = alpha("orange", .1))
  # )


png(file="ball_bi.png", units="in", width=6, height=4, res=300, bg = 'transparent')
biball
dev.off()



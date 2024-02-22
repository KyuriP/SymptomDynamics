library(rootSolve)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)

b <- 1
a <- 1
W <- 1
d <- 4


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
bseq <- seq(-8, 1, by=0.01) 

# plot bifurcation
lst <- bseq |> 
  map_df(~ stability(.x), .id = "id") |>
  mutate(beta = as.factor(bseq[as.numeric(id)]),
         Eigen = ifelse(Eigen == 1, "unstable", "stable"))

stable <- lst |> filter(Eigen =="stable")
unstable <- lst |> filter(Eigen =="unstable")


bif <- lst |> ggplot(aes(x = beta, y = Eq, color = Eigen)) +
  geom_point(size = 0.1) +
  theme_classic() +
  scale_color_manual(values=c("royalblue", "white"), labels = c("stable", "unstable")) +
  geom_line(data= unstable, aes(x = beta, y = Eq), group=1, linetype=2,linewidth=0.6, color = "royalblue") +
  geom_point(aes(x=205, y=1), shape = 8, color = "red", size = 6)+
  geom_point(aes(x=700, y = 0), shape = 8, color = "red", size = 6)+
  annotate('rect', xmin = 205, xmax = 700, ymin =0, ymax= 1, alpha=0.1, fill="orange")+
  labs(x = expression(beta[i]), y = expression("S*"[i]))+
  scale_y_continuous(breaks=c(0,1))+
  scale_x_discrete(breaks = c("-6", "-1"), labels = c(expression(-omega[i] *(delta[i]+1) - alpha[ii]), expression(-omega[i]))) + 
  #scale_x_discrete(breaks = c(205, 700), labels = c(expression(-omega *(delta+1) - alpha), expression(-omega))) + 
  theme(axis.text.x=element_text(size = 25), 
        axis.text.y=element_text(size = 25),
        #legend.title = element_blank(),
        legend.position = 'none',
        axis.title.x = element_text(hjust=1.03, vjust = 9, size= 30),
        axis.title.y = element_text(hjust= -1, vjust = 1.05, angle = 0, size = 30),
        plot.margin = margin(4, 4, 4, 5, "cm")
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

# converging to 0
perms <- c(a = 1, b = -1.6, W = 0.5, d = 0.4)
# bistability
perms <- c(a = 1, b = -1, W = 0.5, d = 0.4)
# converging to 1
perms <- c(a = 1, b = -0.55, W = 0.5, d = 0.4)

times <-seq(0, 100,by= 0.01)
possibleS <- seq(0, 1, 0.1)
out <- list()
for(i in 1:length(possibleS)){
  state <- c(S=possibleS[i])
  out[[i]] <-ode(y= state, times = times, func= model, parms = perms)[,2]
}

result1 <- do.call(cbind, out) |> as.data.frame() |> rename_with(~paste0("S_", possibleS)) |> mutate(time = times) |> pivot_longer(!time, names_to = "S", values_to = "value") 

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




# landscape plots

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




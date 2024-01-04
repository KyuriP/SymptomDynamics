
# bistable landscape

x = seq(-5, 5, 0.01)
y = 1/5 * x^4 - 4*x^2 + 10

plot(x, y, type = "l")

rm(list=ls())
# Bifurcation plot

model<-function(t,state,parms)
{ with(as.list(c(state,parms)), 
       { 
         dS <- S*(1 - S)*(b + a*S + W*(1 + d*(S^4)))   
         return(list(dS)) 
       })
} 

beta <- seq(-3, 3, by = 0.01)
states <- seq(0, 1, by = 0.01)

conver <- list()
for (i in 1:length(beta)){
  perms <- c(a = 1, b = beta[i], W = 1, d = 0.5) # pisanamedvectorofparameters 
  times <-seq(0, 100,by= 0.01)
  convergence <- list()
  for(j in 1:length(state)){
  state <- c(S = states[j]) # sisthestate
  out <- ode(y= state, times = times, func= model, parms = perms)
  convergence[[j]] <- unique(round(tail(out[,2], 100),1))
  }
  conver[[i]] <- convergence[j]
}

data.frame(beta, convergence) |>
  ggplot(aes(x = beta, y = convergence)) +
  geom_point()

unique(round(tail(out[,2], 100),1))
class(out)


Bmax <- 5
plot(-1, -1, xlim = c(-3, 5), ylim = c(-2, 2), xlab = "beta", ylab = "N")
alpha <- 1
beta <- seq(-2, Bmax, by = 0.01)
W <- 1
delta <- 1

times <-seq(0, 100,by= 0.01)




library(deSolve)
library(phaaseR)
library(deBif)

model<-function(t,state,parms)
  { with(as.list(c(state,parms)), 
         { 
         dS <- S*(1 - S)*(b + a*S + W*(1 + d*(S^4)))   
         return(list(dS)) 
         })
} 

perms <- c(a = -1, b = -1, W = 0, d = 0.5)#pisanamedvectorofparameters 
state <- c(S = 0.3) #sisthestate

times <-seq(0, 100,by= 0.01)

out<-ode(y= state, times = times, func= model, parms = perms)
head(out)
plot(out)

phaseplane(model, state, perms)
bifurcation(model, state, perms)

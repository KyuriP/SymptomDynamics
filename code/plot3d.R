## =========================================================
## Creating 3d-plots 
##
## make Figure 4 (bifurcation diagram with landscapes) in 3D
## =========================================================
## install packages
source("code/libraries.R")

# open 3d workspace
open3d()
x <-  seq(0, 1, by = 0.0001)
# bistability
x <-  seq(-3, 1, by = 0.0001)
y <-  c(0,1)
# congerging to 1
z <- -dgamma(x, shape = 2, scale = 0.1)
# congerging to 0
z <- -dgamma(rev(x), shape = 2, scale = 0.1)
# bistability
z <- 0.5*x^4 + 2*x^3 + 1*x^2- 2*x + 0.5
plot3d(x, y, z, col = "lightgray", alpha = c(0,1), 
       expand = 1, box = F, type = "s", size = 1, xlab='',ylab='',zlab='', labels = F)
# congerging to 1
spheres3d(x = 0.1, y = 1, z = -3.5, radius = 0.1, color = "#333377")
# congerging to 0
spheres3d(x = 0.9, y = 1, z = -3.5, radius = 0.1, color = "#333377")
# bistability
spheres3d(x = -2.4, y = 1, z = 0.08, radius = 0.1, color = "#333377")
spheres3d(x = 0.4, y = 1, z = 0.08, radius = 0.1, color = "#333377")

axes3d(xlen = 0, ylen= 0, zlen = 0)
title3d('', '', '', '', '')
show2d({
  # converging 0
  # par(mar = c(0,0,2,0))
  # converging 1
  # par(mar = c(2.2,1,2.2,0))
  # bistability
  # par(mar = c(3.9,1,3.9,0), bg =  rgb(255,240,220, alpha = 50, maxColorValue=255))
  par(mar = c(3.9,1,3.9,0))
  plot(x = rep(0,length(possibleS)), possibleS, col = cc, xlim = c(0, 100), pch = 19, cex = 1.5, xaxt='n', yaxt='n', ann=FALSE, bty="n")
  for(i in 1:length(possibleS)){
    lines(result1[result1$S == possibleS[i],]$time, result1[result1$S == possibleS[i],]$value, type="l", col = cc[i], lwd = 2)
  }
}, expand = 1 , texmipmap = F, rotate = 1)
aspect3d(1.5, 2, 1)
#view3d(theta = 100, phi = 30)
par3d(windowRect = c(20, 30, 800, 800))
#par3d(pp)

rgl.snapshot("bistab1.png", fmt="png", top=TRUE)


## Save RGL parameters to a list object
pp_01 <- par3d()$userMatrix
# bistability plot
pp_bi <- par3d(no.readonly = TRUE)

## Save the list to a text file
dput(pp_01, file="3dView_01.R", control = "all")
#dput(pp_bi, file="3dView_bi.R", control = "all")





remotes::install_github("tylermorganwall/rayshader")
library(rayshader)

plot_gg(bif, multicore=TRUE, width=30, height=13, zoom = 0.6, windowsize = c(3000, 3000), shadow = F,   background = "white", raytrace =F, flat_transparent_bg	= T)
render_snapshot("bifurcation3d(2).png", clear = TRUE)

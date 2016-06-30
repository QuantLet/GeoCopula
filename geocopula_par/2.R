# clear variables and close windows
rm(list = ls(all = TRUE))
graphics.off()

# install and load packages
libraries = c("foreign", "gstat", "spacetime", "fossil", "tseries", "mvtnorm", 
              "xts")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)
})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

# load dataset
load("opti.RData")
op = order(opti$value)
opti_10 = opti[op, ][1, ]  #keep the best one, Chapter 5.2
opti_10 = opti_10[4:6]  #keep only parameters needed

# plot
layout(c(1, 2, 3, 4))
plot_v = function(data, vv, par) {
  k         = list()
  k$par     = par
  k$par[4]  = 1  #sigma=1
  data0     = vv
  gammaMod  = k$par[4] - sapply(1:nrow(data0), function(i) 
  Kernel(data0$avgDist[i], data0$timelag[i], nu, k$par[1], k$par[2], k$par[3], k$par[4]))
  vv1       = vv
  vv1$gamma = gammaMod
  plot(vv1, wireframe = T, xlab = list("distance (km)", rot = 30), ylab = list("time lag (days)", rot = -35), 
       scales = list(arrows = F, z = list(distance = 5)), 
       zlim = c(0, 1.2))
}

plot_v(d4, vvd4t, as.matrix(opti_10))

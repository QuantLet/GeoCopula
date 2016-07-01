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

# load dataset of five countries
load("fivecountries.RData")
resultsres = commonans0915mac$vt  #uk, germany, france, italy, sweden

transformC = function(data, method = "uniform") {
  f1    = function(x) {
  fn    = ecdf(x)
  y     = fn(x)
  return(y)
  }
  nr    = nrow(data)
  nc    = ncol(data)
  data1 = sapply(1:nc, function(i) f1(data[, i]))
  colnames(data1) = colnames(data)
  rownames(data1) = rownames(data)
  data1
}


# the basic setting for optimization
nu     = 0.5
reltol = 0.01
trace  = TRUE
maxit  = 50000
sigma  = 1
cutoff = 2000 
tlag   = 20 

namelist = c("UK", "Germany", "France", "Italy", "Sweden")
location = read.csv("location_country.csv")
rownames(location) = location[, 1]

location = location[namelist, ]
location = location[, -1]

## transform to copulae
d3       = t(resultsres)
d4       = qnorm(transformC(d3) * nrow(d3)/(nrow(d3) + 1))  #formula 3.12 in paper
location = SpatialPoints(location, proj4string = CRS("+proj=longlat +datum=WGS84"))
start    = c(2009, 11)
end      = c(2014, 12)
time     = seq(ISOdate(start[1], start[2], 28), ISOdate(end[1], end[2], 28), 
               by = "1 months")
data     = d4
da       = as.data.frame(as.vector(t(data)))

colnames(da) = "inf"
inflation    = STFDF(location, time, da)
summary(sort(spDists(inflation@sp)))

vvd4 = variogram(inf ~ 1, inflation, width = 400, cutoff = 2000, tlags = 0:35)
vvd4[1, 2:3] = 0

na_indi = which(is.na(vvd4), arr.ind = TRUE)
vvd4t = vvd4

for (i in 1:nrow(na_indi)) {
  vvd4t[na_indi[i, 1], na_indi[i, 2]] = 0.5 * 
    (vvd4[na_indi[i, 1] - 1, na_indi[i, 2]] + vvd4[na_indi[i, 1] + 1, na_indi[i, 2]])
}

start      = c(2009, 11)
end        = c(2014, 12)
time = seq(ISOdate(start[1], start[2], 28), ISOdate(end[1], end[2], 28), 
           by = "1 months")
resultsres = commonans0915mac$vt  #uk, germany, france, italy, sweden
d3         = t(resultsres)
d4         = qnorm(transformC(d3) * nrow(d3)/(nrow(d3) + 1))  #formula 3.12 in paper
data       = d4
da         = as.data.frame(as.vector(t(data)))
colnames(da) = "inf"

namelist = c("UK", "Germany", "France", "Italy", "Sweden")
location = read.csv("location_country.csv")
rownames(location) = location[, 1]

location = location[namelist, ]
location = location[, -1]
#inflation    = STFDF(location, time, da)
dd1      = (spDists(inflation@sp))
dd       = dd1
rm(location, time, da, inflation, dd1)

# WLS as, initial parameter, Chapter 5.2 set multiple initial value for
# WLS estimation
a    = c(1e-10, 0.5, 1)
b    = c(1e-20, 1e-05, 0.1, 0.3)
beta = c(1e-10, 0.01, 0.5)

# choice from multiple optimal points of WLS estimates
fit.count = function(vv, a, b, beta) {
  op   = expand.grid(1:length(a), 1:length(b), 1:length(beta))
  fitv = function(vv, a, b, beta) {
  k    = fit.vv(vv, a, b, beta, nu)
  return(c(k$par, k$value))
  }
  k    = sapply(1:nrow(op), function(i) fitv(vv, a = a[op[i, 1]], 
                                             b = b[op[i, 2]], beta = beta[op[i, 3]]))
  xabc = function(x, abc) {
    abc[x]
  }
  op[, 1] = xabc(op[, 1], a)
  op[, 2] = xabc(op[, 2], b)
  op[, 3] = xabc(op[, 3], beta)
  o = cbind(op, t(k))
  colnames(o) = c("a_i", "b_i", "beta_i", "a", "b", "beta", "value")
  o
}

# WLS estimate
fit.vv = function(vv, a = 0.01, b = 0.005, beta = 0.1, nu) {
  sigma = 1
  par0 = c(a, b, beta)
  fit.model = function(data0, par0, nu, sigma) {
    a = par0[1]
    b = par0[2]
    beta = par0[3]
    weightFun = data0$np/data0$gamma^2
    weightFun[1] = 0
    gammaMod = sigma - sapply(1:nrow(data0), function(i) 
      Kernel(data0$avgDist[i], data0$timelag[i], nu, a, b, beta, sigma))
    sum = sum(weightFun * (data0$gamma - gammaMod)^2)
    return(sum)
  }
  vv[1, 2:3] = 0
  k = optim(par = par0, fit.model, data0 = vv, nu = nu, sigma = sigma, 
            method = "L-BFGS-B", lower = c(1e-10, 1e-10, 0.1))
  k
}

Kernel = function(h, u, nu, a, b, beta, sigma) {
  # same as in GeoCopula
  y       = rep(1, length = length(h))
  idx1    = which(h == 0)
  idx2    = which(h != 0)
  c1      = a^2 * u^2 + 1
  c2      = a^2 * u^2 + beta
  y[idx1] = sigma * beta/c1^nu/c2
  y[idx2] = sigma * 2 * beta/c1^nu/c2/gamma(nu) * (b/2 * h[idx2] * (c1/c2)^0.5)^nu * 
    besselK(b * h[idx2] * (c1/c2)^0.5, nu)
  return(y)
}

opti = fit.count(vvd4t, a, b, beta)
save(opti, file = "opti.RData")

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
load("commonans0915mac.RData")
resultsres = commonans0915mac$vt  #uk, germany, france, italy, sweden
nrow(resultsres)
ncol(resultsres)


source("helpfun.r")

# the basic setting for optimization
nu = 0.5
reltol = 0.01
trace = TRUE
maxit = 50000
sigma = 1
cutoff = 2000  # almost include all possible distance
tlag = 20  # tlag is observed from WLS variogram

namelist = c("UK", "Germany", "France", "Italy", "Sweden")
location = read.csv("location_country.csv")
rownames(location) = location[, 1]

location = location[namelist, ]
location = location[, -1]

## transform to copulae
d3 = t(resultsres)
d4 = qnorm(transformC(d3) * nrow(d3)/(nrow(d3) + 1))  #formula 3.12 in paper
location = SpatialPoints(location, proj4string = CRS("+proj=longlat +datum=WGS84"))
start = c(2009, 11)
end = c(2014, 12)
time = seq(ISOdate(start[1], start[2], 28), ISOdate(end[1], end[2], 28), 
           by = "1 months")
data = d4
da = as.data.frame(as.vector(t(data)))
colnames(da) = "inf"
inflation = STFDF(location, time, da)
summary(sort(spDists(inflation@sp)))
par("mar")
par(mar = c(1, 1, 1, 1))
plot(sort(spDists(inflation@sp)))
vvd4 = variogram(inf ~ 1, inflation, width = 400, cutoff = 2000, tlags = 0:35)


vvd4[1, 2:3] = 0

na_indi = which(is.na(vvd4), arr.ind = TRUE)
vvd4t = vvd4

for (i in 1:nrow(na_indi)) {
  
  vvd4t[na_indi[i, 1], na_indi[i, 2]] = 0.5 * (vvd4[na_indi[i, 1] - 1, 
                                                    na_indi[i, 2]] + vvd4[na_indi[i, 1] + 1, na_indi[i, 2]])
  
}

time = seq(ISOdate(start[1], start[2], 28), ISOdate(end[1], end[2], 28), 
           by = "1 months")

data = d4
da = as.data.frame(as.vector(t(data)))
colnames(da) = "inf"
inflation = STFDF(location, time, da)
dd1 = (spDists(inflation@sp))

dd = dd1
rm(location, time, da, inflation, dd1)


# WLS as, initial parameter, Chapter5.2 set multiple initial value for
# WLS estimation
a = c(1e-10, 0.5, 1)
b = c(1e-20, 1e-05, 0.1, 0.3)
beta = c(1e-10, 0.01, 0.5)
system.time(opti <- fit.count(vvd4t, a, b, beta))
op = order(opti$value)
opti_10 = opti[op, ][1, ]  #keep the best one, Chapter 5.2
opti_10 = opti_10[4:6]  #keep only parameters needed

layout(c(1, 2, 3, 4))

plot_v(d4, vvd4t, as.matrix(opti_10))
plot(vvd4t, wireframe = T, xlab = list("distance (km)", rot = 30), ylab = list("time lag (days)", 
                                                                               rot = -35), zlab = "gamma", scales = list(arrows = F, z = list(distance = 5)), 
     zlim = c(0, 1.2))

nu = 0.5


# use c_p function generate a series sets of estimated parameters

namelist = c("UK", "Germany", "France", "Italy", "Sweden")

location = read.csv("location_country.csv")

rownames(location) = location[, 1]

location = location[namelist, ]
location = location[, -1]

tp = c(59:62)  # fit the data from Aug.2002 to May 2008

system.time(pre_y <- sapply(1:length(tp), function(i) c_p(d4, tp[i], start)))
pre_y = t(pre_y)
colnames(pre_y) = c("a", "b", "beta", "value")

system.time(pre_wls <- sapply(1:length(tp), function(i) c_w(d4, tp[i], 
                                                            start)))
colnames(pre_wls) = c("a", "b", "beta")

aa = load("bb.dat")

pre_data = pre_data
## pre_data=pre_y # in case of WCL pre_data=pre_wls # in case of WLS
mse_all = matrix(0, 8, 12)
rownames(mse_all) = colnames(d4)

tp = 62

np = 59
yy_g = matrix(0, np, 5)

for (k in 0:(np - 1)) {
  yy_g[np - k, ] = pre_term1(d4[1:(tp - k), ], tlag = 3, par = pre_data[2, 
                                                                        ][1:3])
  
}




yy = rbind(d4[1:tlag, ], as.matrix(yy_g))
y4 = yy
ncountry = 5

layout(matrix(c(1, 2, 3, 4, 5, 6), 2, 3, byrow = TRUE))
for (i in 1:ncountry) {
  plot(y4[, i], type = "l", xlab = "month", ylab = "inflation rate", 
       col = "blue")
  lines(d3[, i], type = "l", xlab = "month", ylab = "inflation rate", 
        col = "red")
  
  
  
}

# save(y4, d3, file = 'result_final.dat')

res = load("result_final.dat")

par(mar = c(2, 2, 2, 2), cex = 2)
layout(matrix(c(1, 2, 3, 4, 5, 6), 2, 3, byrow = TRUE))
for (i in 1:ncountry) {
  
  plot(time, y4[, i], type = "l", xlab = "", ylab = "", col = "blue", 
       lwd = 2, lty = 1, ylim = c(-0.8, 1.2))
  lines(time, d3[, i], type = "l", xlab = "", ylab = "", col = "red", 
        lwd = 2, lty = 2, ylim = c(-0.8, 1.2))
  
}

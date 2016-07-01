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
  f1       = function(x) {
  fn       = ecdf(x)
  y        = fn(x)
  return(y)
  }
  nr       = nrow(data)
  nc       = ncol(data)
  data1    = sapply(1:nc, function(i) f1(data[, i]))
  colnames(data1) = colnames(data)
  rownames(data1) = rownames(data)
  data1
}

# the basic setting for optimization
nu       = 0.5
reltol   = 0.01
trace    = TRUE
maxit    = 50000
sigma    = 1
cutoff   = 2000 
tlag     = 20 
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
vvd4         = variogram(inf ~ 1, inflation, width = 400, cutoff = 2000, tlags = 0:35)
vvd4[1, 2:3] = 0

na_indi  = which(is.na(vvd4), arr.ind = TRUE)
vvd4t    = vvd4
# save(vvd4, file="vvd4.RData")

# produce plot
plot(vvd4t, wireframe = T, xlab = list("distance (km)", rot = 30), 
     ylab   = list("time lag (days)", rot = -35), zlab = "gamma", 
     scales = list(arrows = F, z = list(distance = 5)), zlim = c(0, 1.2))

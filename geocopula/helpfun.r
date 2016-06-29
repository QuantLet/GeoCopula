transformC = function(data, method = "uniform") {
  f1 = function(x) {
    fn = ecdf(x)
    y = fn(x)
    return(y)
  }
  
  nr = nrow(data)
  nc = ncol(data)
  data1 = sapply(1:nc, function(i) f1(data[, i]))
  colnames(data1) = colnames(data)
  rownames(data1) = rownames(data)
  data1
}


# recover the transformed series back.
recov = function(y, d3, mis) {
  fn = ecdf(d3[, mis])
  
  # uniroot(function(x) fn(x),c(-1,1))
  y1 = pnorm(y)
  # k_1=which(y1>=1) y1[k_1]=1
  f = function(x) fn(x)
  y2 = y
  for (i in 1:length(y1)) {
    y2[i] = uniroot(function(x) f(x) - y1[i], c(-3, max(d3[, mis])))$root
  }
  y2
}

Kernel = function(h, u, nu, a, b, beta, sigma) {
  # same as in GeoCopula
  y = rep(1, length = length(h))
  idx1 = which(h == 0)
  idx2 = which(h != 0)
  c1 = a^2 * u^2 + 1
  c2 = a^2 * u^2 + beta
  y[idx1] = sigma * beta/c1^nu/c2
  y[idx2] = sigma * 2 * beta/c1^nu/c2/gamma(nu) * (b/2 * h[idx2] * (c1/c2)^0.5)^nu * 
    besselK(b * h[idx2] * (c1/c2)^0.5, nu)
  return(y)
}

kernel = function(h, u, nu, par, sigma) {
  # another kernel.Only change the parameter of function
  Kernel(h, u, nu, par[1], par[2], par[3], sigma)
}


# empirical variogram generator
vv_p = function(data, location, start = c(1998, 1), end = c(2008, 5), keep = NULL) {
  location = SpatialPoints(location, proj4string = CRS("+proj=longlat +datum=WGS84"))
  time = seq(ISOdate(start[1], start[2], 28), ISOdate(end[1], end[2], 
                                                      28), by = "1 months")
  da = as.data.frame(as.vector(t(data)))
  colnames(da) = "inf"
  inflation = STFDF(location, time, da)
  summary(sort(spDists(inflation@sp)))
  plot(sort(spDists(inflation@sp)))
  vv = variogram(inf ~ 1, inflation, width = 400, cutoff = 2000, tlags = 0:30)
  vv[1, 2:3] = 0
  
  na_indi = which(is.na(vv), arr.ind = TRUE)
  
  
  for (i in 1:nrow(na_indi)) {
    
    vv[na_indi[i, 1], na_indi[i, 2]] = 0.5 * (vv[na_indi[i, 1] - 1, 
                                                 na_indi[i, 2]] + vv[na_indi[i, 1] + 1, na_indi[i, 2]])
    
  }
  vv
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
    gammaMod = sigma - sapply(1:nrow(data0), function(i) Kernel(data0$avgDist[i], 
                                                                data0$timelag[i], nu, a, b, beta, sigma))
    sum = sum(weightFun * (data0$gamma - gammaMod)^2)
    return(sum)
  }
  vv[1, 2:3] = 0
  k = optim(par = par0, fit.model, data0 = vv, nu = nu, sigma = sigma, 
            method = "L-BFGS-B", lower = c(1e-10, 1e-10, 0.1))
  k
}



# choice from multiple optimal points of WLS estimates
fit.count = function(vv, a, b, beta) {
  op = expand.grid(1:length(a), 1:length(b), 1:length(beta))
  fitv = function(vv, a, b, beta) {
    k = fit.vv(vv, a, b, beta, nu)
    return(c(k$par, k$value))
  }
  k = sapply(1:nrow(op), function(i) fitv(vv, a = a[op[i, 1]], b = b[op[i, 
                                                                        2]], beta = beta[op[i, 3]]))
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


# score calculation for JCEF estimation, suitable for setting in this
# paper but not used in this paper
bivariate.score.fast.wrap = function(dat, pair, a, b, beta, sigma, nu) {
  h = dd[pair[1], pair[2]]
  u = (pair[3] - pair[4]) * 30.42742
  mu1 = dat[pair[3], pair[1]]
  mu2 = dat[pair[4], pair[2]]
  corrcoef = Kernel(h, u, nu, a, b, beta, sigma)
  # to calculate the score of kernel
  dr1 = sigma * beta/((a^2 * u^2 + 1)^nu * (a^2 * u^2 + beta))
  dr2 = b * ((a^2 * u^2 + 1)/(a^2 * u^2 + beta))^0.5 * h
  d_dr2_beta = b * h * 0.5 * ((a^2 * u^2 + 1)/(a^2 * u^2 + beta))^(-0.5) * 
    (a^2 * u^2 + 1) * (-1)/(a^2 * u^2 + beta)^2
  d_dr2_a = b * h * 0.5 * ((a^2 * u^2 + 1)/(a^2 * u^2 + beta))^(-0.5) * 
    2 * a * u^2 * (beta - 1)/(a^2 * u^2 + beta)^2
  # the element of calculation
  d_beta0 = sigma/(a^2 * u^2 + 1)^nu * a^2 * u^2/(a^2 * u^2 + beta)^2
  d_a0 = dr1 * 2 * a * u^2 * (-nu/(a^2 * u^2 + 1) - 1/(a^2 * u^2 + beta))
  
  if (h == 0) {
    d_beta = d_beta0
    d_a = d_a0
    d_b = 0
  }
  if (h > 0) {
    d_beta = d_beta0/gamma(nu) * (dr2/2)^nu * besselK(dr2, nu) + dr1/gamma(nu) * 
      nu * (dr2/2)^(nu - 1) * d_dr2_beta/2 * besselK(dr2, nu) + dr1/gamma(nu) * 
      (dr2/2)^nu * (nu/dr2 * besselK(dr2, nu) - besselK(dr2, nu + 
                                                          1)) * d_dr2_beta
    d_a = d_a0/gamma(nu) * (dr2/2)^nu * besselK(dr2, nu) + dr1/gamma(nu) * 
      nu * (dr2/2)^(nu - 1) * d_dr2_a/2 * besselK(dr2, nu) + dr1/gamma(nu) * 
      (dr2/2)^nu * (nu/dr2 * besselK(dr2, nu) - besselK(dr2, nu + 
                                                          1)) * d_dr2_a
    d_b = dr1/gamma(nu) * nu * (dr2/2)^(nu - 1) * dr2/2/b * besselK(dr2, 
                                                                    nu) + dr1/gamma(nu) * (dr2/2)^nu * (nu/dr2 * besselK(dr2, nu) - 
                                                                                                          besselK(dr2, nu + 1)) * dr2/b
    d_beta = d_beta * 2
    d_a = d_a * 2
    d_b = d_b * 2
  }
  # calculate score in gaussian model
  rfun1 = 1/(1 - corrcoef^2)
  score_corrcoef = corrcoef * rfun1 - corrcoef * rfun1 * rfun1 * (mu1 * 
                                                                    mu1 + mu2 * mu2 - 2 * corrcoef * mu1 * mu2) + rfun1 * mu1 * mu2
  score = matrix(c(score_corrcoef * d_a, score_corrcoef * d_b, score_corrcoef * 
                     d_beta), nrow = 3, ncol = 1)
  return(score)
}


# sum of scores from relevant samples, this is a function from
# GeoCopula, but also can be used to WCL estimate. In WCL setting score
# means likelihood
bivariate.score.all = function(dat, pt, pf, pft, a, b, beta, sigma, nu) {
  score_spatial = sapply(1:nrow(pf), function(i) bivariate.score.fast.wrap(dat, 
                                                                           pf[i, ], a, b, beta, sigma, nu))
  score_temporal = sapply(1:nrow(pt), function(i) bivariate.score.fast.wrap(dat, 
                                                                            pt[i, ], a, b, beta, sigma, nu))
  score_s_and_p = sapply(1:nrow(pft), function(i) bivariate.score.fast.wrap(dat, 
                                                                            pft[i, ], a, b, beta, sigma, nu))
  meanscore_s = apply(score_spatial[-1, ], 1, mean)
  meanscore_t = apply(score_temporal[-2, ], 1, mean)
  meanscore_t_s = apply(score_s_and_p, 1, mean)
  return(list(meanscore_s = meanscore_s, meanscore_t = meanscore_t, meanscore_t_s = meanscore_t_s))
}


# for JCEF estimation, not used in this paper
Q.weight = function(x, nu, sigma, dat, pt, pf, pft, approach = "simple") {
  a = x[1]
  b = x[2]
  beta = x[3]
  scores = bivariate.score.all(dat, pt, pf, pft, a, b, beta, sigma, nu)
  allscore = c(scores$meanscore_s, scores$meanscore_t, scores$meanscore_t_s)
  allscore = as.matrix(allscore)
  Q.weight = allscore %*% t(allscore)
  return(Q.weight)
}
# weight function for JECF
weight_boot = function(par) {
  data = gen(par, 125, 8)
  par0 = c_p1(data)
  score = Q.weight(par0, nu, sigma, data, pt, pf, pft)
  as.vector(score)
}

Q.fun = function(x, nu, sigma, dat, pt, pf, pft, approach = "simple", weight) {
  a = x[1]
  b = x[2]
  beta = x[3]
  scores = bivariate.score.all(dat, pt, pf, pft, a, b, beta, sigma, nu)
  if (approach == "boot") {
    allscore = c(scores$meanscore_s, scores$meanscore_t, scores$meanscore_t_s)
    Ds = nrow(pf)
    Dt = nrow(pt)
    Dc = nrow(pft)
    allnum = c(Ds, Ds, Dt, Dt, Dc, Dc, Dc)
    N = diag(1/sqrt(allnum))
    Q.fun = allscore %*% N %*% solve(weight) %*% N %*% allscore
  } else {
    Q.fun = allscore %*% allscore
    return(Q.fun)
  }
}
# likelihood function for WCL

logh_1 = function(dat, pair, a, b, beta, sigma, nu) {
  h = dd[pair[1], pair[2]]
  u = (pair[3] - pair[4]) * 30.42742
  mu1 = dat[pair[3], pair[1]]
  mu2 = dat[pair[4], pair[2]]
  k = Kernel(h, u, nu, a, b, beta, sigma)
  
  value = -0.5 * log(1 - k * k) - 0.5 * (mu1 * mu1 - 2 * k * mu1 * mu2 + 
                                           mu2 * mu2)/(1 - k * k)
  value
}
logh_sum = function(x, nu, sigma, dat, pt, pf, pft) {
  a = x[1]
  b = x[2]
  beta = x[3]
  score_spatial = sapply(1:nrow(pf), function(i) logh_1(dat, pf[i, ], 
                                                        a, b, beta, sigma, nu))
  score_temporal = sapply(1:nrow(pt), function(i) logh_1(dat, pt[i, ], 
                                                         a, b, beta, sigma, nu))
  score_s_and_p = sapply(1:nrow(pft), function(i) logh_1(dat, pft[i, 
                                                                  ], a, b, beta, sigma, nu))
  l = sum(score_spatial) + sum(score_temporal) + sum(score_s_and_p)
  return(-l)
}

# three pair producer of tree kinds, time(pt), space(pf) and
# spatiotemporal(pft)
pft_f = function(data, tlag) {
  pf = which(dd <= cutoff & dd > 0, arr.in = TRUE)  #double counted
  pt = matrix(0, (nrow(data) - tlag) * tlag + tlag * (tlag - 1)/2, 2)
  m = 1
  for (i in 1:(nrow(data) - 1)) {
    if (i <= nrow(data) - tlag) {
      k = tlag
    } else {
      k = nrow(data) - i
    }
    pt[m:(m + k - 1), 1] = i
    pt[m:(m + k - 1), 2] = c((i + 1):(k + i))
    m = m + k
  }
  a = nrow(pf)
  b = nrow(pt)
  pft = expand.grid(1:a, 1:b)
  pft[, c(3:4)] = pf[pft[, 1], ]
  pft[, c(5:6)] = pt[pft[, 2], ]  #create D3 pair. 
  pft = pft[, -1]
  pft = pft[, -1]
  # pf
  t = nrow(data)
  for (i in c(1:t)) {
    pf0 = cbind(pf, matrix(i, nrow(pf), 2))
    if (i == 1) {
      pf1 = pf0
    } else {
      pf1 = rbind(pf1, pf0)
    }
  }
  # pt
  for (i in c(1:ncol(data))) {
    pf0 = cbind(matrix(i, nrow(pt), 2), pt)
    if (i == 1) {
      pf2 = pf0
    } else {
      pf2 = rbind(pf2, pf0)
    }
  }
  pf = pf1
  pt = pf2
  rm(pf1, pf2, pf0)
  pft = as.matrix(pft)
  k = list()
  k[[1]] = pf
  k[[2]] = pt
  k[[3]] = pft
  k
}

# optimum_function for WCL

optim_p = function(opti, dat0, pt, pf, pft, sigma = 1) {
  
  opt = function(par, dat0, nu, sigma, pt, pf, pft) {
    est = optim(par, logh_sum, dat = dat0, nu = nu, sigma = sigma, 
                 pt = pt, pf = pf, pft = pft, control = list(maxit = maxit, 
                                                             trace = trace, reltol = reltol))
    return(c(est$par, est$value))
  }
  
  k = sapply(1:nrow(opti), function(i) opt(par = c(opti[i, 1], opti[i, 2], opti[i, 3]), dat0, nu, sigma, pt, pf, pft))
  k = t(k)
  colnames(k) = c("a", "b", "beta", "value")
  k
}

# optimum_function for JECF
optim_s = function(opti, dat0, nu = 0.5, sigma = 1, weight) {
  
  opt = function(par, dat0, nu, sigma, weight) {
    est = optim(par, Q.fun, dat = dat0, nu = nu, sigma = sigma, pt = pt, 
                 pf = pf, pft = pft, approach = "boot", weight = weight, control = list(maxit = maxit, 
                                                                                        trace = trace, reltol = reltol))
    return(c(est$par, est$value))
  }
  k = sapply(1:nrow(opti), function(i) opt(par = c(opti[i, 1], opti[i, 2], opti[i, 3]), dat0, nu, sigma, weight))
  k = t(k)
  colnames(k) = c("a", "b", "beta", "value")
  k
}
# continuous esitimation for WCL
c_p = function(data, tp, start) {
  # continuous prediction
  d = data[1:tp, ]
  nc = ncol(d)
  nr = nrow(d)
  end = c(1:2)
  end[1] = start[1] + floor(tp/12)
  end[2] = start[2] - 1 + tp - floor(tp/12) * 12
  if (end[2] > 12) {
    end[1] = end[1] + 1
    end[2] = end[2] - 12
  }
  
  # empirical variogram
  vv = vv_p(d, location, start = start, end = end)
  vv[1, 2:3] = 0
  a = c(1e-10, 0.1)  # 
  b = c(1e-10, 1e-05, 0.001, 0.1)
  beta = c(1e-10, 0.001, 0.01, 0.1)
  # a=1e-10# b=1e-5# weile jiandan beta=1e-10#
  system.time(opti = fit.count(vv, a, b, beta))
  op = order(opti$value)
  opti_10 = opti[op, ][1:5, ]
  opti_10 = opti_10[4:6]
  PFT = pft_f(d, tlag)
  pf = PFT[[1]]
  pt = PFT[[2]]
  pft = PFT[[3]]
  reltol = 0.01  # gai 
  system.time(est_o = optim_p(as.matrix(opti_10[1, ]), dat0 = d, pt, 
                               pf, pft))
  k = which(est_o[, 4] == min(est_o[, 4]))
  par = as.matrix(est_o[k, ][1:3])
  x = c(1:4)
  x[1:3] = par
  x[4] = est_o[k, ][4]
  
  x
}
# continuous esitimation for WSL continuous prediction
c_w = function(data, tp, start) {
  d = data[1:tp, ]
  nc = ncol(d)
  nr = nrow(d)
  end = c(1:2)
  end = c(1:2)
  end[1] = start[1] + floor(tp/12)
  end[2] = start[2] - 1 + tp - floor(tp/12) * 12
  if (end[2] > 12) {
    end[1] = end[1] + 1
    end[2] = end[2] - 12
  }
  # empirical variogram
  vv = vv_p(d, location, start = start, end = end)
  vv[1, 2:3] = 0
  vv[is.na(vv)] = 0
  vv[vv == 0] = 0.1
  # optimum
  a = c(1e-10, 0.1)
  b = c(1e-10, 1e-05)
  beta = c(1e-10, 0.01)
  # a=1e-10# b=1e-5# weile jiandan beta=1e-10#
  system.time(opti = fit.count(vv, a, b, beta))
  op = order(opti$value)
  opti_1 = opti[op, ][1, ]
  opti_1 = opti_1[4:6]
  est_o = opti_1
  par = as.matrix(est_o[1:3])
  x = c(1:3)
  x[1:3] = par
  x
}

# estimation function for bootstrip in JECF continuous prediction
c_p1 = function(data) {
  tp = nrow(data)
  d = data[1:tp, ]
  nc = ncol(d)
  nr = nrow(d)
  end = c(1:2)
  if (tp%%12 == 0) {
    end[1] = 1998 + floor(tp/12) - 1
    end[2] = 12
  } else {
    end[1] = 1998 + floor(tp/12)
    end[2] = tp%%12
  }
  # empirical variogram
  vv = vv_p(d, location, end = end)
  vv[1, 2:3] = 0
  # optimum
  a = c(1e-10, 0.1)
  b = c(1e-10, 1e-05)
  beta = c(1e-10, 0.001)
  # a=1e-10# b=1e-5# weile jiandan beta=1e-10#
  system.time(opti = fit.count(vv, a, b, beta))
  op = order(opti$value)
  opti_10 = opti[op, ][1:5, ]
  opti_10 = opti_10[4:6]
  # reltol=0.01# gai
  system.time(est_o = optim_p(as.matrix(opti_10[1, ]), dat0 = d, pt, 
                               pf, pft))
  k = which(est_o[, 4] == min(est_o[, 4]))
  par = as.matrix(est_o[k, ][1:3])
  # prediction ??????????????????????
  # y=pre_single1(d,tlag=3,mis=country,par.fit=par)
  x = c(1:4)
  x[1:3] = par
  x[4] = est_o[k, ][4]
  # colnames(xx)=c('y','a','b','beta','value')
  x
}

# AR estimation
ar_p = function(mis, method = "in", tp, p = 3) {
  if (method == "in") {
    
    g4 = ts(d3[, mis], start = c(1998, 1), end = c(2008, 5), frequency = 12)
    fit4 = arima(g4, order = c(p, 0, 0), transform.pars = FALSE)
    d3_forcas = d3[, mis]
    y = matrix(0, length(tp), 2)
    y[, 2] = d3_forcas[tp]
    rownames(y) = rownames(d3[tp, ])
    for (i in tp) {
      y[i - tp[1] + 1, 1] = fit4$coef[p + 1]
      for (j in 1:p) {
        y[i - tp[1] + 1, 1] = y[i - tp[1] + 1, 1] + fit4$coef[j] * 
          (d3_forcas[i - j] - fit4$coef[p + 1])
      }
    }
    y = y[-1, ]
    y
  } else {
    ar_o = function(mis, tp) {
      end = c(1:2)
      if (tp%%12 == 0) {
        end[1] = 1998 + floor(tp/12) - 1
        end[2] = 12
      } else {
        end[1] = 1998 + floor(tp/12)
        end[2] = tp%%12
      }
      g4 = ts(d3[, mis][1:tp], start = c(1998, 1), end = end, frequency = 12)
      fit4 = arima(g4, order = c(p, 0, 0), transform.pars = FALSE)
      d3_forcas = d3[, mis][1:tp]
      y = fit4$coef[p + 1]
      for (j in 1:p) {
        y = y + fit4$coef[j] * (d3_forcas[tp - j + 1] - fit4$coef[p + 1])
      }
      y
    }
    y = matrix(0, length(tp), 2)
    y[, 2] = d3[, mis][tp]
    rownames(y) = rownames(d3[tp, ])
    y = y[-1, ]
    y[, 1] = sapply(1:(length(tp) - 1), function(i) ar_o(mis, tp[i]))
  }
  y
}

# best nu continuous prediction data = d4 tp = 62 startyear = start[1]
nu_p = function(data, tp, start, location) {
  d = data[1:tp, ]
  nc = ncol(d)
  nr = nrow(d)
  end = c(1:2)
  end[1] = start[1] + floor(tp/12)
  end[2] = start[2] - 1 + tp - floor(tp/12) * 12
  # empirical variogram
  vv = vv_p(d, location, start = start, end = end)
  a = c(1e-10)
  b = c(1e-10, 1e-05, 0.1)
  beta = c(1e-10, 0.1)
  # a=1e-10# b=1e-5# weile jiandan beta=1e-10#
  system.time(opti = fit.count(vv, a, b, beta))
  op = order(opti$value)
  opti_10 = opti[op, ][1:5, ]
  opti_10 = opti_10[4:6]
  PFT = pft_f(d, tlag)
  pf = PFT[[1]]
  pt = PFT[[2]]
  pft = PFT[[3]]
  reltol = 0.01  # gai 
  system.time(est_o = optim_p(as.matrix(opti_10[1, ]), dat0 = d, pt, 
                               pf, pft))
  k = which(est_o[, 4] == min(est_o[, 4]))
  par = as.matrix(est_o[k, ][1:3])
  x = c(1:5)
  x[1:3] = par
  x[4] = est_o[k, ][4]
  x[5] = nu
  x
}

# prediction normal prediction
pre_term1 = function(data, tlag, par) {
  nr = nrow(data)
  nc = ncol(data)
  par.fit = par
  tp = nr
  y = matrix(0, 1, nc)
  
  
  x = matrix(0, nc * tlag, 1)
  
  x = as.matrix(data[tp - 1, ])
  for (i in 2:tlag) {
    x = c(x, as.matrix(data[tp - i, ]))
    
  }
# this is for total prediction
cov_xy = matrix(0, nc, nc * tlag)
  for (i in 1:nc) {
    for (j in 1:tlag) {
      for (k in 1:nc) {
        cov_xy[i, (j - 1) * nc + k] = kernel(dd[i, k], j * 30.42742, 
                                             nu, par.fit, sigma)
      }
    }
  }
  
  cov_x = matrix(0, nc * tlag, nc * tlag)
  for (i in 1:nc) {
    for (j in 0:(tlag - 1)) {
      for (k in 1:nc) {
        for (l in 0:(tlag - 1)) {
          cov_x[j * nc + i, l * nc + k] = kernel(dd[i, k], (j - 
                                                              l) * 30.42742, nu, par.fit, sigma)
        }
      }
    }
  }
  # m_test=cov_xy%*%solve(cov_x)
  y = cov_xy %*% solve(cov_x) %*% x
  y
}
# missing value prediction
pre_single1 = function(data, tlag, mis, par.fit) {
  y = 0
  nr = nrow(data)
  tp1 = nr
  nc = ncol(data)
  index_c = mis
  keep = which(c(1:nc) != index_c)
   # constructing x
  x = matrix(0, nc * tlag + nc - 1, 1)
  
  x = as.matrix(data[tp1, ][keep])
  for (i in 1:tlag) {
    x = c(x, as.matrix(data[tp1 - i, ]))
    
  }
  cov_xy = matrix(0, 1, nc * tlag + nc - 1)
  
  for (i in 1:length(keep)) {
    cov_xy[1, i] = kernel(dd[keep[i], index_c], 0, nu, par.fit, sigma)
  }
  
  for (j in 1:tlag) {
    for (k in 1:nc) {
      cov_xy[1, j * nc + k - 1] = kernel(dd[index_c, k], j * 30.42742, 
                                         nu, par.fit, sigma)
    }
  }
  
  cov_x = matrix(0, nc * tlag + nc - 1, nc * tlag + nc - 1)
  for (i in 1:(nc - 1)) {
    for (j in 1:(nc - 1)) {
      cov_x[i, j] = kernel(dd[keep[i], keep[j]], 0, nu, par.fit, 
                           sigma)
    }
  }
  for (i in 1:(nc - 1)) {
    for (j in 1:(tlag)) {
      for (k in 1:nc) {
        cov_x[i, j * nc + k - 1] = kernel(dd[keep[i], k], j, nu, 
                                          par.fit, sigma)
      }
    }
  }
  for (i in 1:(nc - 1)) {
    for (j in 1:(tlag)) {
      for (k in 1:nc) {
        cov_x[j * nc + k - 1, i] = kernel(dd[keep[i], k], j, nu, 
                                          par.fit, sigma)
      }
    }
  }
  for (i in 1:nc) {
    for (j in 1:tlag) {
      for (k in 1:nc) {
        for (l in 1:tlag) {
          cov_x[j * nc + i - 1, l * nc + k - 1] = kernel(dd[i, k], (j - l) * 30.42742, nu, par.fit, sigma)
        }
      }
    }
  }
  y = cov_xy %*% solve(cov_x) %*% x
  y
}

# plot predicted series
plot_result = function(pre_data, mis = 4, m1 = "term", m2 = "n", m3 = "in", 
                       tp = c(55:125)) 
# ????if??????plot function????????????????
  
{
  if (m1 == "single") {
    yy = c(1:(length(tp) - 1))
    pa0 = as.matrix(pre_data[1, ][1:3])
    
    for (i in tp[2]:tp[length(tp)]) {
      k = i - tp[1] + 1
      if (m3 == "in") {
        pa = pa0
      } else {
        pa = as.matrix(pre_data[k - 1, ][1:3])
      }
      yy[k - 1] = pre_single1(d4[1:i, ], tlag = 3, mis = mis, par.fit = pa)
     }
    yy = as.matrix(yy)
    y2 = as.matrix(d4[tp[2]:tp[length(tp)], mis])
    yyy = cbind(yy, y2)
    colnames(yyy) = c("predicted", "real")
    y3 = yyy
    if (m2 == "m" | m2 == "r") {
      y3[, 1] = m_a_m(as.matrix(y3[, 1]), 1)
    }
    
  } else {
    yy_g = matrix(0, length(tp) - 1, 8)
    pa0 = as.matrix(pre_data[1, ][1:3])
    for (i in tp[2]:tp[length(tp)]) {
      k = i - tp[1] + 1
      if (m3 == "in") {
        pa = pa0
      } else {
        pa = as.matrix(pre_data[k - 1, ][1:3])
      }
      yy_g[k - 1, ] = pre_term1(d4[1:i, ], tlag = 3, par = pa)
    }
    yy = as.matrix(yy_g)
    y2 = as.matrix(d4[tp[2]:tp[length(tp)], mis])
    yyy = cbind(yy[, mis], y2)
    colnames(yyy) = c("predicted", "real")
    y3 = yyy
    if (m2 == "m" | m2 == "r") {
      y3[, 1] = m_a_m(as.matrix(y3[, 1]), 1)
    }
  }
  
  if (m2 == "m" | m2 == "n") {
    png(filename = paste(m3, m1, m2, mis, ".png", sep = "_"), width = 500, 
        height = 300)
    plot(y3[, 1], type = "l", ylim = c(-3, 3), xaxt = "n", xlab = "month", 
         ylab = "value", col = "blue")
    axis(1, at = seq(6, 71, 7), labels = rownames(y3)[seq(6, 71, 7)])  # Transformed Value in N(0,1) Inflation rate
    par(new = TRUE)
    plot(y3[, 2], type = "l", col = "red", xaxt = "n", axes = FALSE, 
         ylim = c(-3, 3), xlab = "", ylab = "")
    dev.off()
    re = y3
  } else {
    y4 = y3
    y4[, 1] = recov(y3[, 1], d3, mis)
    y4[, 2] = as.matrix(d3[tp[2]:tp[length(tp)], mis])
    png(filename = paste(m3, m1, m2, mis, ".png", sep = "_"), width = 500, 
        height = 300)
    plot(y4[, 1], type = "l", ylim = c(-1, 5), xaxt = "n", xlab = "month", 
         ylab = "inflation rate", col = "blue")
    axis(1, at = seq(6, 70, 7), labels = rownames(y3)[seq(6, 70, 7)])  # Transformed Value in N(0,1) Inflation rate
    par(new = TRUE)
    plot(y4[, 2], type = "l", col = "red", xaxt = "n", axes = FALSE, 
         ylim = c(-1, 5), xlab = "", ylab = "")
    dev.off()
    re = y4
  }
  re
}

# plot function for variogram
plot_v = function(data, vv, par) {
  k = list()
  k$par = par
  k$par[4] = 1  #sigma=1
  data0 = vv
  gammaMod = k$par[4] - sapply(1:nrow(data0), function(i) Kernel(data0$avgDist[i], 
                                                                 data0$timelag[i], nu, k$par[1], k$par[2], k$par[3], k$par[4]))
  vv1 = vv
  vv1$gamma = gammaMod
  plot(vv1, wireframe = T, xlab = list("distance (km)", rot = 30), ylab = list("time lag (days)", 
                                                                               rot = -35), scales = list(arrows = F, z = list(distance = 5)), 
       zlim = c(0, 1.2))
}

# moving average smoother
m_a_m = function(d, n) {
  # use 2n as detrend
  k = nrow(d)
  for (i in (n + 1):(k - n)) {
    sum = d[i - n, ]
    for (j in (i - n + 1):(i + n)) {
      sum = sum + d[j, ]
    }
    d[i, ] = 1/(2 * n + 1) * sum
  }
  return(d)
}
# generating function for bootstrip for JECF. Not used in this paper

gen = function(par.fit, tlag, nc) {
  
  cov_x = matrix(0, nc * tlag, nc * tlag)
  for (i in 1:nc) {
    for (j in 0:(tlag - 1)) {
      for (k in 1:nc) {
        for (l in 0:(tlag - 1)) {
          cov_x[j * nc + i, l * nc + k] = kernel(dd[i, k], (j - l) * 30.42742, nu, par.fit, sigma)
        }
      }
    }
  }
  
  simu_data = rmvnorm(1, sigma = cov_x)
  d6 = matrix(0, 125, 8)
  for (i in 1:125) {
    d6[126 - i, ] = simu_data[((i - 1) * 8 + 1):((i - 1) * 8 + 8)]
  }
  colnames(d6) = colnames(d4)
  d6
}

library(tidyverse)
library(furrr)
library(dtw)
library(DTWBI)

plan(multiprocess, workers = 10)


order = c(2,1,1)
ar = c(0.8, -0.5)
ma = (-0.2)
n = 100
noise_mean = 0
noise_sd = 1
alpha1 = 0
alpha2 = 0
stretch = 2

arima_ts = arima.sim(list(order = order, ar = ar, ma = ma), n = n)
noise = rnorm(n + 1, mean = noise_mean, sd = noise_sd)
noise[1] = 0

ts1 = alpha1 + arima_ts
ts2_short = alpha2 + arima_ts + noise
ts2 = approx(ts2_short, n = (n + 1) * stretch, method = "constant")$y
# ts2 = approx(ts1, n = (n + 1) * stretch, method = "constant")$y
# ts2 = alpha2 + approx(arima_ts + noise, n = (n + 1) * stretch)$y

dts1 = local.derivative.ddtw(ts1)
dts2 = local.derivative.ddtw(ts2)

align1 = dtw(ts1, ts2, keep = TRUE)
dtwPlotTwoWay(align1)
dtwPlotThreeWay(align1)

align2 = dtw(dts1, dts2, keep = TRUE)
dtwPlotTwoWay(align2)
dtwPlotThreeWay(align2)

align3 = dtw(ts1, ts2_short, keep = TRUE)
dtwPlotTwoWay(align3)
dtwPlotThreeWay(align3)

ts2_warpped = ts2[warp(align1, index.reference = TRUE)]

dist_short = mean((ts1 - ts2_short)^2)
dist = mean((ts1 - ts2_warpped)^2)






ts_simulation = function(order = c(1,1,0), ar = c(0.8), 
                         ma = NULL, n = 100, 
                         noise_mean = 0, noise_sd = 1,
                         stretch = 2, method = "constant",
                         random_state = NULL){
  if (!is.null(random_state)) {
    set.seed(random_state)
  }
  
  # ARIMA ts
  arima_ts = arima.sim(list(order = order, ar = ar, ma = ma), n = n)
  
  # white noise
  noise = rnorm(n + 1, mean = noise_mean, sd = noise_sd)
  noise[1] = 0
  
  # ts1
  ts1 = arima_ts
  
  # ts2
  ts2_short = arima_ts + noise
  ts2 = approx(ts2_short, n = (n + 1) * stretch, method = method)$y
  
  return(list(order = order, ar = ar,
              ma = ma, n = n, 
              noise_mean = noise_mean,
              noise_sd = noise_sd,
              stretch = stretch,
              method = method,
              random_state = random_state,
              ts1 = ts1,
              ts2_short = ts2_short,
              ts2 = ts2))
}

ts_sim = ts_simulation(order = c(2,1,1), ar = c(0.8, -0.5), ma = c(-0.2))


do_dtw = function(ts1, ts2, method = "dtw"){
  ts1_bk = ts1
  ts2_bk = ts2
  if (method == "ddtw") {
    ts1 = DTWBI::local.derivative.ddtw(ts1)
    ts2 = DTWBI::local.derivative.ddtw(ts2)
  }
  
  # dtw
  alignment = dtw::dtw(ts1, ts2, keep = TRUE)
  fig_TwoWay = dtw::dtwPlotTwoWay(alignment)
  fig_ThreeWay = dtw::dtwPlotThreeWay(alignment)
  
  # warp
  ts1_warpped = ts1[warp(alignment, index.reference = FALSE)]
  ts2_warpped = ts2[warp(alignment, index.reference = TRUE)]
  
  return(list(
    ts1 = ts1_bk,
    ts2 = ts2_bk,
    alignment = alignment,
    fig_TwoWay = fig_TwoWay,
    fig_ThreeWay = fig_ThreeWay,
    ts1_warpped = ts1_warpped,
    ts2_warpped = ts2_warpped
  ))
}

align = do_dtw(ts_sim$ts1, ts_sim$ts2)

comp_dist = function(ts1, ts2_short, ts2_warpped){
  dist_short = mean((ts1 - ts2_short)^2)
  dist = mean((ts1 - ts2_warpped)^2)
  prop = dist/dist_short
  
  return(list(dist_short = dist_short,
              dist = dist,
              prop = prop))
}

comp = comp_dist(align$ts1, ts_sim$ts2_short, align$ts2_warpped)


seed_list = runif(50, 0, 1000)


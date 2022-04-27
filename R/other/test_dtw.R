library(tidyverse)
library(furrr)
library(dtw)
library(DTWBI)

plan(multisession, workers = 6)


## Functions -----------------------------------------------------------------------------------------
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


do_dtw = function(ts1, ts2, method = "dtw"){
  ts1_bk = ts1
  ts2_bk = ts2
  if (method == "ddtw") {
    ts1 = DTWBI::local.derivative.ddtw(ts1)
    ts2 = DTWBI::local.derivative.ddtw(ts2)
    ts1[1] = ts1[2]
    ts1[length(ts1)] = ts1[length(ts1)-1]
    ts2[1] = ts2[2]
    ts2[length(ts2)] = ts2[length(ts2)-1]
  }
  
  # dtw
  alignment = dtw::dtw(ts1, ts2, keep = TRUE)
  fig_TwoWay = dtw::dtwPlotTwoWay(alignment)
  fig_ThreeWay = dtw::dtwPlotThreeWay(alignment)
  
  # warp
  ts1_warpped = suppressWarnings(ts1[warp(alignment, index.reference = FALSE)])
  ts2_warpped = suppressWarnings(ts2[warp(alignment, index.reference = TRUE)])
  
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


comp_dist = function(ts1, ts2_short, ts2_warpped){
  dist_short = mean((ts1 - ts2_short)^2)
  dist = mean((ts1 - ts2_warpped)^2)
  prop = dist/dist_short
  
  return(list(dist_short = dist_short,
              dist = dist,
              prop = prop))
}


do_test = function(order = c(2,1,1), ar = c(0.8, -0.5),
                   ma = c(-0.3), n = 100, n_seed = 50, 
                   noise_mean = 0,
                   noise_sd = 1, stretch = 2,
                   stretch_method = "constant",
                   dtw_method = "dtw"){
  # seeds
  seed_list = runif(n_seed, 0, 1000)
  
  # do test
  res = seed_list %>% 
    as.list %>% 
    future_map(
      ~{
        seed = .
        ts_sim = ts_simulation(order = order, ar = ar,
                               ma = ma, n = n,
                               noise_mean = noise_mean,
                               noise_sd = noise_sd,
                               stretch = stretch,
                               method = stretch_method,
                               random_state = seed)
        res_dtw = do_dtw(ts_sim$ts1, ts_sim$ts2, method = dtw_method)
        res_comp = comp_dist(res_dtw$ts1, ts_sim$ts2_short, res_dtw$ts2_warpped)
        res_comp$prop
      }
    ) %>% 
    do.call("rbind",.)
  
  return(res)
}

## test dtw ------------------------------------------------------------------------------------------
dtw_method = "dtw"
# noise mean
mean_list = c(0, 0.2, 0.5, 1, 2, 3, 4, 5)

for (noise_mean in mean_list) {
  res = do_test(noise_mean = noise_mean, dtw_method = dtw_method)
  res_mean = mean(res)
  print(paste0("noise_mean: ", noise_mean, ". ratio_mean: ", res_mean, "."))
}


# noise sd
sd_list = c(0.001, 0.2, 0.5, 1, 2, 3, 4, 5, 6)

for (noise_sd in sd_list) {
  res = do_test(noise_sd = noise_sd, dtw_method = dtw_method)
  res_mean = mean(res)
  print(paste0("noise_sd: ", noise_sd, ". ratio_mean: ", res_mean, "."))
}

# stretch method
stretch_method_list = c("constant", "linear")

for (stretch_method in stretch_method_list) {
  res = do_test(stretch_method = stretch_method, dtw_method = dtw_method)
  res_mean = mean(res)
  print(paste0("stretch_method: ", stretch_method, ". ratio_mean: ", res_mean, "."))
}

# stretch
stretch_list = c(1, 2, 3, 4, 5)

for (stretch in stretch_list) {
  res = do_test(stretch = stretch, dtw_method = dtw_method)
  res_mean = mean(res)
  print(paste0("stretch: ", stretch, ". ratio_mean: ", res_mean, "."))
}

# ar
ar_list = list(c(0.8, -0.5), c(0.7, -0.5), c(0.5, -0.5), c(0.3, -0.5))

for (ar in ar_list) {
  res = do_test(ar = ar, dtw_method = dtw_method)
  res_mean = mean(res)
  print(paste0("ar-1: ", ar[1], ". ratio_mean: ", res_mean, "."))
}

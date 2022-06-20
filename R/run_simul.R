## Setup -----------------------------------------------------------------------
library(checkpoint)
checkpoint("2022-04-01")

library(tidyverse)
library(furrr)
plan(multisession, workers = 7)
options(future.rng.onMisuse="ignore")
options(stringsAsFactors = FALSE)

source("./R/TwoStepDTW.R")
source("./R/synthetic_control.R")
source("./R/comp_methods.R")
set.seed(20220407)


## Functions -------------------------------------------------------------------
simulate_data = function(n = 1000,
                         burn_in = 100,
                         n_lag = 50,
                         beta1 = 0.9,
                         noise_mean = 0,
                         noise_sd = 0.1){
  # speed profiles
  phi1 = 1/(1 + exp(cumsum(rnorm(n + burn_in))))
  phi2 = 1/(1 + exp(cumsum(rnorm(n + burn_in))))
  phi3 = 1/(1 + exp(cumsum(rnorm(n + burn_in))))
  
  # initialize y and z
  y <- z <- yt <- zt <- u <- ut <- NULL
  # starting values of y and z
  ylag <- zlag <- ulag <- 1
  # Create a sequence of n values for y and z
  x = cumsum(rnorm(n = n + burn_in, mean = 0))
  
  # time series
  for (i in (1 + burn_in):(n + burn_in)) {
    yt = beta1*ylag + phi1[i]*x[i] +
      (1 - phi1[i])*x[i-n_lag] + 
      rnorm(n = 1, mean = noise_mean, sd = noise_sd)
    zt = beta1*zlag + phi3[i]*x[i] + 
      (1 - phi3[i])*x[i-n_lag] + 
      rnorm(n = 1, mean = noise_mean, sd = noise_sd)
    ut = beta1*ulag + phi2[i]*x[i] + 
      (1 - phi2[i])*x[i-n_lag] + 
      rnorm(n = 1, mean = noise_mean, sd = noise_sd)
    y <- c(y, yt)
    z <- c(z, zt)
    u <- c(u, ut)
    ylag = yt
    zlag = zt
    ulag = ut
  }
  
  # data frame
  data = data.frame(unit = c(rep('A', n), rep('B', n), rep('C',n)),
                    time = rep(1:n, 3),
                    value = c(u,y,z)) %>% 
    mutate(id = case_when(unit == "A" ~ 1,
                          unit == "B" ~ 2,
                          TRUE ~ 3),
           value_raw = value) %>% 
    select(c("id", "unit", "time", "value", "value_raw"))
  
  return(data)
}


run_simul = function(data, 
                     start_time = 1,
                     end_time = 200,
                     t_treat = 120,
                     width_range = (1:3)*2+3,
                     k_range = 4:6,
                     dtw1_range = 140:143,
                     n_mse = 20
                     ){
  # grid search
  grid_search = NULL
  i = 1
  for (width in width_range) {
    for (k in k_range) {
      for (dtw1_time in dtw1_range) {
        grid_search[[i]] = data.frame(width = width,
                            k = k,
                            dtw1_time = dtw1_time,
                            mse_original = NA_real_,
                            mse_new = NA_real_,
                            mse_ratio = NA_real_)
        i = i + 1
      }
    }
  }
  
  result = grid_search %>% 
    future_map(
      ~{
        search = .
        width = search$width
        k = search$k
        dtw1_time = search$dtw1_time
        
        data = preprocessing(data, filter_width = width)
        
        res = SimDesign::quiet(compare_methods(data = data,
                                               start_time = start_time,
                                               end_time = end_time,
                                               treat_time = t_treat,
                                               dtw1_time = dtw1_time,
                                               dependent = "A",
                                               dependent_id = 1,
                                               normalize_method = "t",
                                               k = k,
                                               synth_fun = "simulation",
                                               filter_width = width,
                                               plot_figures = FALSE,
                                               step.pattern = dtw::symmetricP2))
        
        synth_original = res$synth_origin$synthetic
        synth_new = res$synth_new$synthetic
        value_raw = res$synth_origin$value
        
        mse_original = mean((synth_original - value_raw)[t_treat:(t_treat + n_mse)]^2, rm.na = T)
        mse_new = mean((synth_new - value_raw)[t_treat:(t_treat + n_mse)]^2, rm.na = T)
        
        mse_ratio = mse_new/mse_original
        
        data.frame(width = width,
                   k = k,
                   dtw1_time = dtw1_time,
                   mse_original = mse_original,
                   mse_new = mse_new,
                   mse_ratio = mse_ratio)
      }
    )
  
  result = result %>%
    do.call("rbind", .)
  
  return(result)
}


## Data Simulation -------------------------------------------------------------
n_simulation = 1000
data_list = NULL

for (i in 1:n_simulation) {
  data_list[[i]] = simulate_data(n = 200,
                                 burn_in = 40,
                                 n_lag = 20,
                                 beta1 = 0.9,
                                 noise_mean = 0,
                                 noise_sd = 0.1)
}


## Run -------------------------------------------------------------------------
result = NULL

for (i in 1:length(data_list)) {
  cat(paste0("Simulation ", i, "..."))
  result[[i]] = run_simul(data_list[[i]],
                          start_time = 1,
                          end_time = 200,
                          t_treat = 120,
                          width_range = (1:10)*2+3,
                          k_range = 4:20,
                          dtw1_range = 120:140,
                          n_mse = 20)
  cat("Done.\n")
}


## Plot Result -----------------------------------------------------------------
# width = 21
# k = 7

width = 7
k = 5


data = preprocessing(data, filter_width = width)

res = compare_methods(data = data,
                      start_time = 1,
                      end_time = 200,
                      treat_time = 120,
                      dtw1_time = 140,
                      dependent = "A",
                      dependent_id = 1,
                      normalize_method = "t",
                      k = k,
                      synth_fun = "simulation",
                      filter_width = width,
                      plot_figures = TRUE,
                      step.pattern = dtw::symmetricP2)

df = rbind(res$df %>%
             filter(time <= 200) %>%
             select(c("unit", "time", "value_raw")) %>% 
             `colnames<-`(c("unit", "time", "value")),
           data.frame(unit = "w/o TSDTW",
                      time = 1:200,
                      value = res$synth_origin$synthetic[1:200]),
           data.frame(unit = "w/ TSDTW",
                      time = 1:200,
                      value = res$synth_new$synthetic[1:200])
           )

ggplot(df, aes(x = time, y = value, color = unit)) +
  geom_line(aes(linetype = unit)) + 
  geom_vline(xintercept = 120, linetype="dashed") +
  theme_bw() +
  scale_linetype_manual(values=c("solid", "dashed", "twodash", "solid", "solid")) +
  scale_color_manual(values=c("#4a4e4d", "grey70", "grey80", "#fe4a49","#2ab7ca"))


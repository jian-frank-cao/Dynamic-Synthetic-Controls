## Setup -----------------------------------------------------------------------
library(checkpoint)
checkpoint("2022-04-01")

library(tidyverse)
library(furrr)
plan(multisession, workers = 7)
options(future.rng.onMisuse="ignore")
options(stringsAsFactors = FALSE)
source("./R/misc.R")
source("./R/TFDTW.R")
source("./R/synth.R")
source("./R/compare.R")
source("./R/simulate.R")
set.seed(20220407)


## Data Simulation -------------------------------------------------------------
n_simulation = 1000
length = 1000
n = 5

# generate sobol sequence
sobol_seq = qrng::sobol(n_simulation*1, d = n - 1, randomize = "Owen",
                        seed = 20220401, skip = 100)
rnd_speeds = cbind(rep(0.5, n_simulation), sobol_seq*0.6 + 0.1)

# simulate
data_list = NULL
for (i in 1:n_simulation) {
  data_list[[i]] = SimData_ShapeSpeed(n = 5,
                                      length = 1000,
                                      rnd_speed = rnd_speeds[i,],
                                      n_SMA = 10,
                                      ar_x = 0.9,
                                      weight_speed = TRUE,
                                      speed_upper = 1,
                                      speed_lower = 0.5,
                                      t_treat = 800,
                                      shock = 50)
}

data_list[[21]] %>% ggplot(aes(x = time, y = value, color = unit)) +
  geom_line() +
  geom_vline(xintercept = 800, linetype="dashed")

saveRDS(data_list, "./data/simul_data_list_0811.Rds")


## Run -------------------------------------------------------------------------
data_list = readRDS("./data/simul_data_list_0811.Rds")
start_time = 1
end_time = 1000
treat_time = 800
dtw1_time = 900
n_mse = 100
k = 15
plot_figures = TRUE
step.pattern1 = dtw::symmetricP2
step.pattern2 = dtw::asymmetricP2
legend_position = c(0.3, 0.8)

result = data_list[1:100] %>% 
  future_map(
    ~{
      data = .
      res = SimDesign::quiet(compare_methods(data = data,
                            start_time = start_time,
                            end_time = end_time,
                            treat_time = treat_time,
                            dtw1_time = dtw1_time,
                            dependent = "A",
                            dependent_id = 1,
                            n_mse = n_mse,
                            k = k,
                            plot_figures = plot_figures,
                            step.pattern1 = step.pattern1,
                            step.pattern2 = step.pattern2,
                            predictors.origin = NULL,
                            special.predictors.origin = list(list("value_raw", 700:799, c("mean"))),
                            time.predictors.prior.origin = 1:799,
                            time.optimize.ssr.origin = 1:799,
                            predictors.new = NULL,
                            special.predictors.new = list(list("value_warped", 700:799, c("mean"))),
                            time.predictors.prior.new = 1:799,
                            time.optimize.ssr.new = 1:799,
                            legend_position = legend_position))
      
      synth_original = res$synth_origin$synthetic
      synth_new = res$synth_new$synthetic
      value_raw = res$synth_origin$value
      
      mse_original = mean((synth_original - value_raw)[1:(treat_time - 1)]^2, rm.na = T)
      mse_new = mean((synth_new - value_raw)[1:(treat_time - 1)]^2, rm.na = T)
      
      mse_ratio = mse_new/mse_original

      list(mse_original = mse_original,
           mse_new = mse_new,
           mse_ratio = mse_ratio,
           synth_original = synth_original,
           synth_new = synth_new,
           value_raw = value_raw)
    }
  )
  
  
  








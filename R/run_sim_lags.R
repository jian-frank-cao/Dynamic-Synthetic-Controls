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
length = 100
n = 5

# generate sobol sequence
sobol_seq = qrng::sobol(n_simulation*3, d = n - 1, randomize = "Owen",
                        seed = 20220401, skip = 100)
rnd_nCycles = cbind(rep(0.9, n_simulation),
                    sobol_seq[1:n_simulation,])
rnd_shift = cbind(rep(0.9, n_simulation),
                  sobol_seq[(n_simulation + 1):(2*n_simulation),])
rnd_lag = cbind(rep(0.9, n_simulation),
                sobol_seq[(2*n_simulation + 1):(3*n_simulation),])

# simulate
data_list = NULL
for (i in 1:n_simulation) {
  data_list[[i]] = SimData_Lags(n = n,
                                length = length,
                                rnd_nCycles = rnd_nCycles[i,],
                                rnd_shift = rnd_shift[i,],
                                rnd_lag = rnd_lag[i,],
                                nCycles_min = 6,
                                nCycles_max = 12,
                                noise_mean = 0,
                                noise_sd = 0.01,
                                n_lag_min = 5,
                                n_lag_max = 15,
                                extra_x = 20,
                                beta = 0.9,
                                ar_x = 0.9,
                                t_treat = 80,
                                shock = 50)
}

data_list[[169]] %>% ggplot(aes(x = time, y = value, color = unit)) +
  geom_line() +
  geom_vline(xintercept = 80, linetype="dashed")

saveRDS(data_list, "./data/simul_data_list_0906.Rds")


## Run -------------------------------------------------------------------------
data_list = readRDS("./data/simul_data_list_0906.Rds")
start_time = 1
end_time = 100
treat_time = 80
dtw1_time = 90
n_mse = 10
k = 15
dist_quantile = 0.95
n_IQR = 3
plot_figures = FALSE
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
                            n_IQR = n_IQR,
                            dist_quantile = dist_quantile,
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
  
saveRDS(result, "./data/res_sim_0831.Rds")
# result = readRDS("./data/res_sim_0819.Rds")

## Plot result -----------------------------------------------------------------
# placebo test figure
df = future_map2(
  result,
  as.list(1:length(result)),
  ~{
    item = .x
    id = .y
    gap_origin =  item[["value_raw"]] - item[["synth_original"]]
    gap_new = item[["value_raw"]] - item[["synth_new"]]
    data.frame(time = 1:length(gap_new),
               gap_origin = gap_origin,
               gap_new = gap_new,
               id = id)
  }
) %>% 
  do.call("rbind", .)


treatment = c(rep(0, 800),
              seq(0, 50, length.out = 100),
              rep(50, 100))

percent = df %>%
  group_by(time) %>% 
  summarise(ci_origin_upper = quantile(gap_origin, 0.975, na.rm = T),
            ci_origin_mean = mean(gap_origin, na.rm = T),
            ci_origin_lower = quantile(gap_origin, 0.025, na.rm = T),
            ci_new_upper = quantile(gap_new, 0.975, na.rm = T),
            ci_new_mean = mean(gap_new, na.rm = T),
            ci_new_lower = quantile(gap_new, 0.025, na.rm = T)) %>% 
  mutate(artifical_effect = treatment,
         id = 0)

fig = df %>% 
  ggplot(aes(x = time, group = id)) +
  geom_line(aes(y = gap_origin), col = "#4d648d", alpha=0.1) +
  geom_line(aes(y = gap_new), col = "#feb2a8", alpha=0.1) +
  geom_line(aes(x = time, y = ci_origin_upper), data = percent, col = "#2ab7ca", alpha=0.8) +
  geom_line(aes(x = time, y = ci_origin_lower), data = percent, col = "#2ab7ca", alpha=0.8) +
  geom_line(aes(x = time, y = ci_new_upper), data = percent, col = "#fe4a49", alpha=0.8) +
  geom_line(aes(x = time, y = ci_new_lower), data = percent, col = "#fe4a49", alpha=0.8) +
  geom_line(aes(x = time, y = ci_origin_mean), data = percent, col = "#2ab7ca", alpha=1) +
  geom_line(aes(x = time, y = ci_new_mean), data = percent, col = "#fe4a49", alpha=1) +
  geom_line(aes(x = time, y = artifical_effect), data = percent, col = "#008744", alpha=1) +
  geom_vline(xintercept = 800, linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  xlab("Time") +
  ylab("Synthetic Control - True Value") +
  theme_bw()

ggsave("./figures/placebo_sim_0813.pdf",
       fig, width = 6, height = 4,
       units = "in", limitsize = FALSE)



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
simulate_data_arima = function(n = 100,
                         burn_in = 20,
                         n_lag = 10,
                         ar_phi = 0.9,
                         ar_x = 0.9,
                         beta1 = 0.9,
                         noise_mean = 0,
                         noise_sd = 0.1){
  # speed profiles
  phi1 = arima.sim(list(order = c(1,0,0), ar = ar_phi), n = n + burn_in)
  phi2 = arima.sim(list(order = c(1,0,0), ar = ar_phi), n = n + burn_in)
  phi3 = arima.sim(list(order = c(1,0,0), ar = ar_phi), n = n + burn_in)
  
  phi1 = normalize(phi1, normalize_method = "minmax")
  phi2 = normalize(phi2, normalize_method = "minmax")
  phi3 = normalize(phi3, normalize_method = "minmax")
  
  # initialize y and z
  y <- z <- yt <- zt <- u <- ut <- NULL
  # starting values of y and z
  ylag <- zlag <- ulag <- 1
  # Create a sequence of n values for y and z
  x = arima.sim(list(order = c(1,0,0), ar = ar_x), n = n + burn_in)
  
  # n_lag
  if (length(n_lag) == 2) {
    n_lag = round(runif(n = 1, min = min(n_lag), max = max(n_lag)), 0)
  }
  
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


simulate_data_sin = function(n = 3,
                             length = 100,
                             rnd_nCycles = seq(0.1, 0.9, length.out = n),
                             rnd_shift = seq(0.9, 0.1, length.out = n),
                             rnd_trend = seq(0.2, 0.8, length.out = n),
                             nCycles_min = 6,
                             nCycles_max = 12,
                             trend_min = 0.01,
                             trend_max = 0.1,
                             length_divide = 50,
                             shock = 0.5){
  
  # prepare random numbers
  nCycles = rnd_nCycles * (nCycles_max - nCycles_min) + nCycles_min
  shifts = rnd_shift * 2 * pi
  trends = rnd_trend * (trend_max - trend_min) + trend_min
  
  # simulate
  data = NULL
  for (i in 1:n) {
    x = seq(0, nCycles[i] * pi, length.out = length)
    x = cumsum(sin(x + shifts[i])/2 + 1)/(length/length_divide)
    y = sin(x)
    if (i == 1) {
      trend = rep(trends[i], length) + c(rep(0, length*4/5),
                                   seq(0, shock, length.out = length/20),
                                   seq(shock, 0, length.out = length/20),
                                   rep(0, length/5-length/10))
    }else{
      trend = trends[i]
    }
    y = cumsum(y/2 + trend)
    data = rbind(data,
                 data.frame(id = i,
                            unit = LETTERS[i],
                            time = 1:length,
                            value = y,
                            value_raw = y))
  }
  return(data)
}


run_simul = function(data, 
                     start_time = 1,
                     end_time = 200,
                     t_treat = 120,
                     width_range = (1:3)*2+3,
                     k_range = 4:6,
                     dtw1_range = 140:145,
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
        
        
        mse = data.frame(width = width,
                   k = k,
                   dtw1_time = dtw1_time,
                   mse_original = mse_original,
                   mse_new = mse_new,
                   mse_ratio = mse_ratio)
        
        list(mse = mse,
             synth_original = synth_original,
             synth_new = synth_new,
             value_raw = value_raw)
      }
    )
  
  return(result)
}


## Data Simulation -------------------------------------------------------------
n_simulation = 1000
length = 100
n = 5
data_list = NULL

# generate sobol sequence
sobol_seq = qrng::sobol(n_simulation*3, d = n, randomize = "Owen",
                        seed = 20220401, skip = 100)
rnd_nCycles = sobol_seq[1:n_simulation,]
rnd_shift = sobol_seq[(n_simulation + 1):(2*n_simulation),]
rnd_trend = sobol_seq[(2*n_simulation + 1):(3*n_simulation),]


# for (i in 1:n_simulation) {
#   data_list[[i]] = simulate_data_arima(n = 200,
#                                  burn_in = 40,
#                                  n_lag = 20,
#                                  # n_lag = 2,
#                                  # n_lag = c(2,20),
#                                  beta1 = 0.9,
#                                  ar_phi = 0.9,
#                                  ar_x = 0.9,
#                                  noise_mean = 0,
#                                  noise_sd = 0.1)
# }

for (i in 1:n_simulation) {
  data_list[[i]] = simulate_data_sin(n = n,
                                     length = length,
                                     rnd_nCycles = rnd_nCycles[i,],
                                     rnd_shift = rnd_shift[i,],
                                     rnd_trend = rnd_trend[i,],
                                     nCycles_min = 6,
                                     nCycles_max = 12,
                                     trend_min = 0.01,
                                     trend_max = 0.1,
                                     length_divide = 50,
                                     shock = 0.5)
}


saveRDS(data_list, "./data/simul_data_list_0720.Rds")


## Run -------------------------------------------------------------------------
data_list = readRDS("./data/simul_data_list_0720.Rds")
result = NULL

for (i in 1:length(data_list)) {
  cat(paste0("Simulation ", i, "..."))
  result[[i]] = run_simul(data_list[[i]],
                          start_time = 1,
                          end_time = 200,
                          t_treat = 120,
                          width_range = (1:3)*2+3,
                          k_range = 4:6,
                          dtw1_range = 135:140,
                          n_mse = 20)
  cat("Done.\n")
}

saveRDS(result, "./data/res_simul_0711_v1.Rds")

min_ratio = result %>% 
  map(
    ~{
      res = lapply(., "[[", "mse") %>% do.call("rbind", .)
      min(res$mse_ratio)
    }
  ) %>% 
  do.call("c", .)

log_min_ratio = log(min_ratio)

t.test(log_min_ratio)

boxplot(log_min_ratio, outline = F, 
        # xlab = "Simulation Data",
        # ylab = latex2exp::TeX("$log(MSE_{w/ TFDTW}/MSE_{w/o TFDTW})$")
        ylab = "Log Ratio"
        )
abline(h = 0, lty = 5)
text(1,-0.1,"t test: P = 0")

# placebo test figure
df = future_map2(
  result,
  as.list(1:length(result)),
  ~{
    item = .x
    id = .y
    mse = lapply(., "[[", "mse") %>% do.call("rbind", .)
    n = which(mse$mse_ratio == min(mse$mse_ratio))[1]
    gap_origin = item[[n]][["synth_original"]] - item[[n]][["value_raw"]]
    gap_new = item[[n]][["synth_new"]] - item[[n]][["value_raw"]]
    data.frame(time = 1:200, 
               gap_origin = gap_origin,
               gap_new = gap_new,
               id = id)
  }
) %>% 
  do.call("rbind", .)

percent = df %>%
  group_by(time) %>% 
  summarise(ci_origin_upper = quantile(gap_origin, 0.975, na.rm = T),
            ci_origin_mean = mean(gap_origin, na.rm = T),
            ci_origin_lower = quantile(gap_origin, 0.025, na.rm = T),
            ci_new_upper = quantile(gap_new, 0.975, na.rm = T),
            ci_new_mean = mean(gap_new, na.rm = T),
            ci_new_lower = quantile(gap_new, 0.025, na.rm = T)) %>% 
  mutate(id = 0)



df %>% 
  # filter(unit %in% (mse %>% filter(mse1_pre < 2*10000) %>% .[["dependent"]])) %>% 
  ggplot(aes(x = time, group = id)) +
  geom_line(aes(y = gap_origin), col = "blue", alpha=0.1) +
  geom_line(aes(y = gap_new), col = "purple", alpha=0.1) +
  geom_line(aes(x = time, y = ci_origin_upper), data = percent, col = "black", alpha=1) +
  geom_line(aes(x = time, y = ci_origin_lower), data = percent, col = "black", alpha=1) +
  geom_line(aes(x = time, y = ci_new_upper), data = percent, col = "red", alpha=1) +
  geom_line(aes(x = time, y = ci_new_lower), data = percent, col = "red", alpha=1) +
  geom_vline(xintercept = 120, linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  # coord_cartesian(ylim=c(-32, 32)) +
  xlab("Time") +
  ylab("Synthetic Control - True Value") +
  theme_bw()


## Plot Result -----------------------------------------------------------------
# Plot example
# width = 21
# k = 7

data = data_list[[105]]
width = 9
k = 6
dtw1_time = 138
max_x = 150

data = preprocessing(data, filter_width = width)

res = compare_methods(data = data,
                      start_time = 1,
                      end_time = 200,
                      treat_time = 120,
                      dtw1_time = dtw1_time,
                      dependent = "A",
                      dependent_id = 1,
                      normalize_method = "t",
                      k = k,
                      synth_fun = "simulation",
                      filter_width = width,
                      plot_figures = FALSE,
                      step.pattern = dtw::symmetricP2)

df = rbind(res$df %>%
             filter(time <= max_x) %>%
             select(c("unit", "time", "value_raw")) %>% 
             `colnames<-`(c("unit", "time", "value")),
           data.frame(unit = "w/o TFDTW",
                      time = 1:max_x,
                      value = res$synth_origin$synthetic[1:max_x]),
           data.frame(unit = "w/ TFDTW",
                      time = 1:max_x,
                      value = res$synth_new$synthetic[1:max_x])
           )

fig = ggplot(df, aes(x = time, y = value, color = unit)) +
  geom_line(aes(linetype = unit)) + 
  geom_vline(xintercept = 120, linetype="dashed") +
  theme_bw() +
  scale_linetype_manual(values=c("solid", "dashed", "twodash", "solid", "solid")) +
  scale_color_manual(values=c("#4a4e4d", "grey70", "grey80", "#fe4a49","#2ab7ca"))


ggsave("./figures/simul_example.pdf",
       fig, width = 6, height = 4,
       units = "in", limitsize = FALSE)

# 

args = commandArgs(trailingOnly=TRUE)
i = as.integer(args[1])
job.start = Sys.time()

## Setup -----------------------------------------------------------------------
library(checkpoint)
checkpoint("2022-04-01")

library(tidyverse)
library(furrr)
plan(multisession, workers = 8)
options(future.rng.onMisuse="ignore")
options(stringsAsFactors = FALSE)
source("./R/misc.R")
source("./R/TFDTW.R")
source("./R/synth.R")
source("./R/implement.R")
source("./R/simulate.R")
source("./R/grid.search.R")
set.seed(20220407)


## Data Simulation -------------------------------------------------------------
# n_simulation = 1000
# length = 100
# n = 10
# rescale = 0.7
# 
# # generate sobol sequence
# sobol_seq = qrng::sobol(n_simulation*3, d = n - 1, randomize = "Owen",
#                         seed = 20220401, skip = 100)
# rnd_nCycles = cbind(rep(0.9, n_simulation),
#                     sobol_seq[1:n_simulation,]*rescale)
# rnd_shift = cbind(rep(0.9, n_simulation),
#                   sobol_seq[(n_simulation + 1):(2*n_simulation),]*rescale)
# rnd_lag = cbind(rep(0.9, n_simulation),
#                 sobol_seq[(2*n_simulation + 1):(3*n_simulation),]*rescale)
# 
# # simulate
# data_list = NULL
# for (i in 1:n_simulation) {
#   data_list[[i]] = SimData_Lags(n = n,
#                                 length = length,
#                                 rnd_nCycles = rnd_nCycles[i,],
#                                 rnd_shift = rnd_shift[i,],
#                                 rnd_lag = rnd_lag[i,],
#                                 nCycles_min = 5,
#                                 nCycles_max = 15,
#                                 noise_mean = 0,
#                                 noise_sd = 0.01,
#                                 n_lag_min = 5,
#                                 n_lag_max = 15,
#                                 extra_x = 20,
#                                 beta = 0.9,
#                                 ar_x = 0.9,
#                                 t_treat = 80,
#                                 shock = 50)
# }
# 
# data_list[[1]] %>% ggplot(aes(x = time, y = value, color = unit)) +
#   geom_line() +
#   geom_vline(xintercept = 80, linetype="dashed")
# 
# saveRDS(data_list, "./data/simul_data_list_0916.Rds")


## Run -------------------------------------------------------------------------
data_list = readRDS("./data/simul_data_list_0919.Rds")

# parameters
filter.width.range = (1:9)*2+3
k.range = 4:9
step.pattern.range = list(
  # symmetricP0 = dtw::symmetricP0, # too bumpy
  # symmetricP05 = dtw::symmetricP05,
  symmetricP1 = dtw::symmetricP1,
  symmetricP2 = dtw::symmetricP2,
  # asymmetricP0 = dtw::asymmetricP0, # too bumpy
  # asymmetricP05 = dtw::asymmetricP05,
  asymmetricP1 = dtw::asymmetricP1,
  asymmetricP2 = dtw::asymmetricP2,
  typeIc = dtw::typeIc,
  # typeIcs = dtw::typeIcs,
  # typeIIc = dtw::typeIIc,  # jumps
  # typeIIIc = dtw::typeIIIc, # jumps
  # typeIVc = dtw::typeIVc,  # jumps
  typeId = dtw::typeId,
  # typeIds = dtw::typeIds,
  # typeIId = dtw::typeIId, # jumps
  mori2006 = dtw::mori2006
)
grid.search.parallel = TRUE


args.TFDTW = list(buffer = 0, match.method = "fixed",
                  dist.quant = 0.95, 
                  window.type = "sakoechiba",
                  ## other
                  norm.method = "t",
                  step.pattern2 = dtw::asymmetricP2,
                  n.burn = 3, n.IQR = 3,
                  ma = 3, ma.na = "original",
                  default.margin = 3,
                  n.q = 1, n.r = 1)

args.synth = list(predictors = NULL,
                  special.predictors = 
                    expression(list(list(dep.var, 70:79, c("mean")),
                                    list(dep.var, 60:69, c("mean")),
                                    list(dep.var, 50:59, c("mean")))),
                  time.predictors.prior = 1:79,
                  time.optimize.ssr = 1:79)

args.TFDTW.synth = list(start.time = 1, end.time = 100, treat.time = 80,
                        args.TFDTW = args.TFDTW, args.synth = args.synth,
                        ## 2nd
                        n.mse = 10, 
                        ## other
                        plot.figures = FALSE,
                        plot.path = "./figures/",
                        legend.pos = c(0.3, 0.3))

args.TFDTW.synth.all.units = list(target = "A",
                                  # data = data, 
                                  args.TFDTW.synth = args.TFDTW.synth,
                                  ## 2nd
                                  all.units.parallel = FALSE)

cat("Simulation data set ", i, "...")
args.TFDTW.synth.all.units[["data"]] = data_list[[i]]
results = SimDesign::quiet(
  grid.search(filter.width.range = filter.width.range,
              k.range = k.range,
              step.pattern.range = step.pattern.range,
              args.TFDTW.synth.all.units = args.TFDTW.synth.all.units,
              grid.search.parallel = grid.search.parallel)
)
cat("Done.\n")

saveRDS(results, paste0("./data/res_sim_0921_", i, ".Rds"))
job.end = Sys.time()
print(job.end - job.start)

## Plot result -----------------------------------------------------------------
results = NULL
res.files = list.files("./data/res_sim")
for (res.file in res.files) {
  results = c(results, list(readRDS(paste0("./data/res_sim/", res.file))))
}

length = 100
shock = 50
treatment = c(rep(0, length*4/5),
              seq(0, shock, length.out = length/20),
              seq(shock, 0, length.out = length/20),
              rep(0, length/5-length/10))

causal_effect = NULL
causal_effect_lag = 0
beta = 0.9
treat_time = 80
n_mse = 10
for (j in 1:length) {
  temp = treatment[j] + beta*causal_effect_lag
  causal_effect <- c(causal_effect, temp)
  causal_effect_lag = temp
}

# placebo test figure
res = future_map2(
  results,
  as.list(1:length(results)),
  ~{
    item = .x
    id = .y
    neg.ratio = lapply(item, "[[", "neg.ratio") %>% do.call("c", .)
    max.neg.ratio = which(neg.ratio == max(neg.ratio, na.rm = T))[1]
    # min_ind = 82
    # gap_origin =  -item[[min_ind]][["gap_original"]]
    # gap_new = -item[[min_ind]][["gap_new"]]
    # data.frame(time = 1:length(gap_new),
    #            gap_origin = gap_origin,
    #            gap_new = gap_new,
    #            id = id)
    synth_original = item[[min_ind]][["res.synth.target.raw"]][["synthetic"]]
    synth_new = item[[min_ind]][["res.synth.target.TFDTW"]][["synthetic"]]
    value_raw = item[[min_ind]][["res.synth.target.raw"]][["value"]]
    
    gap_original = value_raw - synth_original
    gap_new = value_raw - synth_new

    diff_original = value_raw - synth_original - causal_effect
    diff_new = value_raw - synth_new - causal_effect

    mse_original = mean((diff_original)[treat_time:(treat_time + n_mse)]^2, rm.na = T)
    mse_new = mean((diff_new)[treat_time:(treat_time + n_mse)]^2, rm.na = T)
    mse = data.frame(mse_original = mse_original,
                     mse_new = mse_new,
                     log_ratio = log(mse_new/mse_original))
    list(df = data.frame(time = 1:length(value_raw),
                         id = id,
                         value_raw = value_raw,
                         synth_original = synth_original,
                         synth_new = synth_new,
                         gap_origin = gap_original,
                         gap_new = gap_new,
                         diff_original = diff_original,
                         diff_new = diff_new),
         mse = mse)
  }
)

mse = lapply(res, "[[", "mse") %>% do.call("rbind", .)
t.test(mse$log_ratio)

df = lapply(res, "[[", "df") %>% do.call("rbind", .)
  


percent = df %>%
  group_by(time) %>%
  summarise(ci_origin_upper = quantile(gap_origin, 0.975, na.rm = T),
            ci_origin_mean = mean(gap_origin, na.rm = T),
            ci_origin_lower = quantile(gap_origin, 0.025, na.rm = T),
            ci_new_upper = quantile(gap_new, 0.975, na.rm = T),
            ci_new_mean = mean(gap_new, na.rm = T),
            ci_new_lower = quantile(gap_new, 0.025, na.rm = T)) %>%
  mutate(artifical_effect = causal_effect,
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
  geom_vline(xintercept = 80, linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  xlab("Time") +
  ylab("Synthetic Control - True Value") +
  theme_bw()

ggsave("./figures/placebo_sim_0813.pdf",
       fig, width = 6, height = 4,
       units = "in", limitsize = FALSE)

# results = readRDS("./data/res_sim_0916.Rds")
# length = 100
# shock = 50
# treatment = c(rep(0, length*4/5),
#               seq(0, shock, length.out = length/20),
#               seq(shock, 0, length.out = length/20),
#               rep(0, length/5-length/10))
# 
# causal_effect = NULL
# causal_effect_lag = 0
# beta = 0.9
# treat_time = 80
# n_mse = 10
# for (j in 1:length) {
#   temp = treatment[j] + beta*causal_effect_lag
#   causal_effect <- c(causal_effect, temp)
#   causal_effect_lag = temp
# }
# 
# # placebo test figure
# df = future_map2(
#   result,
#   as.list(1:length(result)),
#   ~{
#     item = .x
#     id = .y
#     mse = lapply(item, "[[", "mse") %>% do.call("rbind", .)
#     min_ind = which(mse$mse_pre_new == min(mse$mse_pre_new, na.rm = T))[1]
#     # min_ind = 82
#     # gap_origin =  -item[[min_ind]][["gap_original"]]
#     # gap_new = -item[[min_ind]][["gap_new"]]
#     # data.frame(time = 1:length(gap_new),
#     #            gap_origin = gap_origin,
#     #            gap_new = gap_new,
#     #            id = id)
#     
#     temp = mse[min_ind,]
#     synth_original = item[[min_ind]][["synth_original"]]
#     synth_new = item[[min_ind]][["synth_new"]]
#     value_raw = item[[min_ind]][["value_raw"]]
#     
#     diff_original = value_raw - synth_original - causal_effect
#     diff_new = value_raw - synth_new - causal_effect
#     
#     mse_original = mean((diff_original)[treat_time:(treat_time + n_mse)]^2, rm.na = T)
#     mse_new = mean((diff_new)[treat_time:(treat_time + n_mse)]^2, rm.na = T)
#     temp$mse_original = mse_original
#     temp$mse_new = mse_new
#     temp
#   }
# ) %>% 
#   do.call("rbind", .)
# 
# df = df %>% mutate(log_ratio = log(mse_new/mse_original))
# t.test(df$log_ratio)
# 
# treatment = c(rep(0, 80),
#               seq(0, 50, length.out = 10),
#               rep(50, 10))
# 
# percent = df %>%
#   group_by(time) %>% 
#   summarise(ci_origin_upper = quantile(gap_origin, 0.975, na.rm = T),
#             ci_origin_mean = mean(gap_origin, na.rm = T),
#             ci_origin_lower = quantile(gap_origin, 0.025, na.rm = T),
#             ci_new_upper = quantile(gap_new, 0.975, na.rm = T),
#             ci_new_mean = mean(gap_new, na.rm = T),
#             ci_new_lower = quantile(gap_new, 0.025, na.rm = T)) %>% 
#   mutate(artifical_effect = treatment,
#          id = 0)
# 
# fig = df %>% 
#   ggplot(aes(x = time, group = id)) +
#   geom_line(aes(y = gap_origin), col = "#4d648d", alpha=0.1) +
#   geom_line(aes(y = gap_new), col = "#feb2a8", alpha=0.1) +
#   geom_line(aes(x = time, y = ci_origin_upper), data = percent, col = "#2ab7ca", alpha=0.8) +
#   geom_line(aes(x = time, y = ci_origin_lower), data = percent, col = "#2ab7ca", alpha=0.8) +
#   geom_line(aes(x = time, y = ci_new_upper), data = percent, col = "#fe4a49", alpha=0.8) +
#   geom_line(aes(x = time, y = ci_new_lower), data = percent, col = "#fe4a49", alpha=0.8) +
#   geom_line(aes(x = time, y = ci_origin_mean), data = percent, col = "#2ab7ca", alpha=1) +
#   geom_line(aes(x = time, y = ci_new_mean), data = percent, col = "#fe4a49", alpha=1) +
#   geom_line(aes(x = time, y = artifical_effect), data = percent, col = "#008744", alpha=1) +
#   geom_vline(xintercept = 80, linetype="dashed") +
#   geom_hline(yintercept = 0, linetype="dashed") +
#   xlab("Time") +
#   ylab("Synthetic Control - True Value") +
#   theme_bw()
# 
# ggsave("./figures/placebo_sim_0813.pdf",
#        fig, width = 6, height = 4,
#        units = "in", limitsize = FALSE)



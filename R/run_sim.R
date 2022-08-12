## Setup -----------------------------------------------------------------------
library(checkpoint)
checkpoint("2022-04-01")

library(tidyverse)
library(furrr)
plan(multisession, workers = 7)
options(future.rng.onMisuse="ignore")
options(stringsAsFactors = FALSE)

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
n_simulation = 1000
length = 1000
n = 5
result = NULL



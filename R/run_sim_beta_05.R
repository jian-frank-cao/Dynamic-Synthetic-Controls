## Setup -----------------------------------------------------------------------
library(checkpoint)
checkpoint("2022-04-01")

library(parallel)
n.cores = detectCores()
library(tidyverse)
library(furrr)
plan(multisession, workers = n.cores - 1)
options(future.rng.onMisuse="ignore")
options(stringsAsFactors = FALSE)
source("./R/utility/misc.R")
source("./R/utility/TFDTW.R")
source("./R/utility/synth.R")
source("./R/utility/implement.R")
source("./R/utility/simulate.R")
source("./R/utility/grid.search.R")
set.seed(20220407)


## Data Simulation -------------------------------------------------------------
n.simulation = 100
length = 100
n = 10
beta = 0.5
shock = 10

# simulate
data.list = NULL
for (i in 1:n.simulation) {
  data.list[[i]] = sim.data(n = n, length = length,
                            t.treat = 60, shock = shock, ar.x = 0.6,
                            n.SMA = 1, n.diff = 1,
                            speed.upper = 2,
                            speed.lower = 0.5,
                            reweight = TRUE,
                            rescale = TRUE,
                            rescale.multiplier = 20,
                            beta = beta)
}

saveRDS(data.list, paste0("./data/simul_data_beta_", beta, ".Rds"))


## Run -------------------------------------------------------------------------
data.list = readRDS(paste0("./data/simul_data_beta_", beta, ".Rds"))

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


args.TFDTW = list(buffer = 20, match.method = "open.end",
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
                    expression(list(list(dep.var, 50:59, c("mean")),
                                    list(dep.var, 40:49, c("mean")),
                                    list(dep.var, 30:39, c("mean")))),
                  time.predictors.prior = 1:59,
                  time.optimize.ssr = 1:59)

args.TFDTW.synth = list(start.time = 1, end.time = 100, treat.time = 60,
                        args.TFDTW = args.TFDTW, args.synth = args.synth,
                        ## 2nd
                        n.mse = 10, 
                        ## other
                        plot.figures = FALSE,
                        plot.path = "./figures/",
                        legend.pos = c(0.3, 0.7))

args.TFDTW.synth.all.units = list(target = "A",
                                  # data = data, 
                                  args.TFDTW.synth = args.TFDTW.synth,
                                  ## 2nd
                                  all.units.parallel = FALSE)

for (i in 1:length(data.list)) {
  cat("Simulation data set ", i, "...")
  args.TFDTW.synth.all.units[["data"]] = data.list[[i]]
  results = SimDesign::quiet(
    grid.search(filter.width.range = filter.width.range,
                k.range = k.range,
                step.pattern.range = step.pattern.range,
                args.TFDTW.synth.all.units = args.TFDTW.synth.all.units,
                grid.search.parallel = grid.search.parallel)
  )
  saveRDS(results, paste0("./data/res_sim/beta", beta,
                          "/res_sim_beta_", beta, "_", i, ".Rds"))
  cat("Done.\n")
}


## Test result -----------------------------------------------------------------
folder = paste0("./data/res_sim/beta", beta, "/")
file.list = as.list(list.files(folder))

length = 100
treat_time = 60
n_mse = 10
treatment = c(rep(0, treat_time),
              seq(0, shock, length.out = round(0.1*length)),
              rep(shock, round(0.9*length - treat_time)))

pre.start = 51
pre.end = 60
post.start = 61
post.end = 70

# de.mse
df.mse = future_map2(
  file.list,
  as.list(1:length(file.list)),
  ~{
    file.name = .x
    data.id = .y
    data.list = readRDS(paste0(folder, file.name))
    mse = future_map2(
      data.list,
      as.list(names(data.list)),
      ~{
        result = .x
        grid.id = .y
        value = result$res.synth.target.raw$value
        synth.raw = result$res.synth.target.raw$synthetic
        synth.TFDTW = result$res.synth.target.TFDTW$synthetic
        gap.raw = value - synth.raw - treatment
        gap.TFDTW = value - synth.TFDTW - treatment
        data.frame(grid.id = grid.id,
                   mse.preT.raw = mean(gap.raw[pre.start:pre.end]^2, na.rm = T),
                   mse.preT.TFDTW = mean(gap.TFDTW[pre.start:pre.end]^2, na.rm = T),
                   mse.postT.raw = mean(gap.raw[post.start:post.end]^2, na.rm = T),
                   mse.postT.TFDTW = mean(gap.TFDTW[post.start:post.end]^2, na.rm = T))
      }
    ) %>% do.call("rbind", .)
    mse$data.id = data.id
    
    mse %>% 
      top_n(-1, mse.preT.TFDTW) %>% 
      top_n(-1, grid.id)
  }
) %>% do.call("rbind", .)

saveRDS(df.mse, paste0("./data/df.mse_sim_beta_", beta, ".Rds"))

# t.test for log(MSEdsc/MSEsc)
df.mse = readRDS(paste0("./data/df.mse_sim_beta_", beta, ".Rds"))
df.mse = df.mse %>% 
  mutate(log.ratio = log(mse.postT.TFDTW/mse.postT.raw))
t.test(df.mse$log.ratio)

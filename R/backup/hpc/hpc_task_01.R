## Setup -----------------------------------------------------------------------
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
beta = 1
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


## Run -------------------------------------------------------------------------
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

args.TFDTW.synth.target.only = list(target = "A", id = 1,
                                    args.TFDTW.synth = args.TFDTW.synth)


for (i in 1:10) {
  job.start = Sys.time()
  cat("Simulation data set ", i, "...")
  set.seed(20220407)
  args.TFDTW.synth.target.only[["data"]] = data.list[[i]]
  results = SimDesign::quiet(
    grid.search.opt(filter.width.range = filter.width.range,
                    k.range = k.range,
                    step.pattern.range = step.pattern.range,
                    args.TFDTW.synth.target.only = args.TFDTW.synth.target.only,
                    grid.search.parallel = grid.search.parallel)
  )
  cat("Done.\n")
  
  saveRDS(results, paste0("./data/res_sim_rep/res_sim_", i, ".Rds"))
  job.end = Sys.time()
  print(job.end - job.start)
}


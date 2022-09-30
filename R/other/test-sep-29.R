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


data.list = readRDS("./data/simul_data_list_0926.Rds")
i = 200


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


args.TFDTW = list(buffer = 10, match.method = "open.end",
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
                        plot.figures = TRUE,
                        plot.path = "./figures/",
                        legend.pos = c(0.3, 0.3))

args.TFDTW.synth.all.units = list(target = "A",
                                  # data = data, 
                                  args.TFDTW.synth = args.TFDTW.synth,
                                  ## 2nd
                                  all.units.parallel = TRUE)

args.TFDTW.synth.all.units[["data"]] = data.list[[i]]

if (grid.search.parallel) {
  fun.map = furrr::future_map
}else{
  fun.map = purrr::map
}

# vanilla synthetic control
data = args.TFDTW.synth.all.units[["data"]]
units = data[c("id", "unit")] %>% distinct
units.list = units %>% split(., seq(nrow(units)))

args.synth = args.TFDTW.synth.all.units$args.TFDTW.synth$args.synth
args.synth[["df"]] = data
args.synth[["dep.var"]] = "value_raw"

res.synth.raw.list = units.list %>% 
  set_names(units$unit) %>% 
  fun.map(
    ~{
      item = .
      dependent.id = item$id
      args.synth[["dependent.id"]] = dependent.id
      res = do.call(do.synth, args.synth)
    }
  )

# grid search space
search.grid = expand.grid(filter.width.range, k.range,
                          names(step.pattern.range)) %>% 
  `colnames<-`(c("filter.width", "k", "step.pattern"))
search.grid.list = search.grid %>% split(., seq(nrow(search.grid)))

task = search.grid.list[[200]]

args.TFDTW.synth.all.units[["filter.width"]] = task$filter.width
args.TFDTW.synth.all.units$args.TFDTW.synth$args.TFDTW[["k"]] = task$k
args.TFDTW.synth.all.units$args.TFDTW.synth$args.TFDTW[["step.pattern1"]] =
  step.pattern.range[[task$step.pattern]]
args.TFDTW.synth.all.units[["res.synth.raw.list"]] = res.synth.raw.list
# res.TFDTW.synth.all.units = do.call(TFDTW.synth.all.units, args.TFDTW.synth.all.units)

target = args.TFDTW.synth.all.units$target
args.TFDTW.synth = args.TFDTW.synth.all.units$args.TFDTW.synth
all.units.parallel = args.TFDTW.synth.all.units$all.units.parallel
data = args.TFDTW.synth.all.units$data
filter.width = args.TFDTW.synth.all.units$filter.width
res.synth.raw.list = args.TFDTW.synth.all.units$res.synth.raw.list

if (!is.null(filter.width)) {
  data = preprocessing(data, filter.width)
}
args.TFDTW.synth[["data"]] = data
units = data[c("id", "unit")] %>% distinct
units.list = units %>% split(., seq(nrow(units)))

# run TFDTW.synth
if (all.units.parallel) {
  fun.map = furrr::future_map
}else{
  fun.map = purrr::map
}

item = units.list[[1]]

dependent = item$unit
dependent.id = item$id
args.TFDTW.synth[["dependent"]] = dependent
args.TFDTW.synth[["dependent.id"]] = dependent.id
args.TFDTW.synth[["res.synth.raw"]] = res.synth.raw.list[[dependent]]
# res.TFDTW.synth = do.call(TFDTW.synth, args.TFDTW.synth)

start.time = args.TFDTW.synth$start.time
end.time = args.TFDTW.synth$end.time
treat.time = args.TFDTW.synth$treat.time
args.TFDTW = args.TFDTW.synth$args.TFDTW
args.synth = args.TFDTW.synth$args.synth
n.mse = args.TFDTW.synth$n.mse
plot.figures = args.TFDTW.synth$plot.figures
args.TFDTW.synth$plot.path
args.TFDTW.synth$legend.pos
args.TFDTW.synth$data
args.TFDTW.synth$dependent
args.TFDTW.synth$dependent.id
args.TFDTW.synth$res.synth.raw


dtw::dtwPlotTwoWay(res.1stDTW$alignment, 
                   xts = data %>% filter(unit == "A") %>% .[["value_raw"]] %>% .[1:80],
                   yts = data %>% filter(unit == "E") %>% .[["value_raw"]] %>% .[1:90] + 50)
plot(ts(data %>% filter(unit == "A") %>% .[["value_raw"]] %>% .[1:80]))


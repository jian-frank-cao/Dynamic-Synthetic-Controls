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


## 1. Speed Difference ---------------------------------------------------------
# function - sim.data
sim.data = function(n = 3, nCycles = 6, length = 100, extra.x = round(0.2*length),
                    t.treat = 60, shock = 0, ar.x = 0.6, n.SMA = 1,
                    n.diff = 1, speed.upper = 2, speed.lower = 0.5,
                    reweight = TRUE, rescale = TRUE,
                    rescale.multiplier = 20, beta = 1){
  # common exogenous shocks
  x = cos(seq(0, nCycles * pi, length.out = length + extra.x + n.SMA + n.diff - 1))
  x = x/2 + 0.5
  
  xt = 0
  x.stack = NULL
  for (j in 1:length(x)) {
    xt = ar.x*xt + x[j]
    x.stack = c(x.stack, xt)
  }
  x = x.stack
  
  # smoothing
  x.SMA = ts(TTR::SMA(x, n = n.SMA)[-(1:(n.SMA - 1))])
  
  # difference
  x.diff = diff(x.SMA, difference = n.diff)
  pos.diff = x.diff > 0
  if (reweight) {
    pos.ratio = sum(pos.diff)/sum(!pos.diff)
  }
  
  # speeds
  log.speeds = seq(log(speed.lower), log(speed.upper), length.out = n)
  rnd.ind = 2
  log.speeds = c(log.speeds[rnd.ind], log.speeds[-rnd.ind])
  
  # simulate
  data = NULL
  for (i in 1:n) {
    # speed profile
    log.speed = log.speeds[i]
    if (reweight) {
      if (pos.ratio > 1) {
        pos.speed = exp(log.speed*(1/pos.ratio))
        neg.speed = exp(-log.speed)
      }else{
        pos.speed = exp(log.speed)
        neg.speed = exp(-log.speed*pos.ratio)
      }
    }else{
      pos.speed = exp(log.speed)
      neg.speed = exp(-log.speed)
    }
    
    phi.shape = rep(NA, length.out = length + extra.x)
    phi.shape[pos.diff] = pos.speed
    phi.shape[!pos.diff] = neg.speed
    
    log.phi.mean = mean(log(phi.shape), na.rm = T)
    log.phi.sd = sd(log(phi.shape), na.rm = T)
    
    phi.random = exp(rnorm(n = length + extra.x,
                           mean = log.phi.mean,
                           sd = log.phi.sd))
    
    # treatment
    if (i == 1) {
      treatment = c(rep(0, t.treat),
                    seq(0, shock, length.out = round(0.1*length)),
                    rep(shock, round(0.9*length - t.treat)))
    }else{
      treatment = 0
    }
    
    phi = beta*phi.shape + (1 - beta)*phi.random
    
    y = warpWITHweight(x[1:(length + extra.x)], phi)[1:length]
    
    if (rescale) {
      y = minmax.normalize(y, reference = y[1:t.treat])*rescale.multiplier
    }
    
    y = y + treatment
    y = y + rnorm(n = length(y), mean = 0, sd = 0.2)
    
    data = rbind(data,
                 data.frame(id = i,
                            unit = LETTERS[i],
                            time = 1:length,
                            value = y))
  }
  return(data)
}

# funcation TFDTW.synth.all.units
TFDTW.synth.all.units = function(data, target, 
                                 args.TFDTW.synth,
                                 filter.width = NULL,
                                 res.synth.raw.list = NULL,
                                 detailed.output = FALSE,
                                 all.units.parallel = FALSE){
  # prepare data
  if (!is.null(filter.width)) {
    data = preprocessing(data, filter.width)
  }
  args.TFDTW.synth[["data"]] = data
  units = data[c("id", "unit")] %>% distinct
  units.list = units %>% split(., seq(nrow(units))) %>% .[1]
  
  # run TFDTW.synth
  if (all.units.parallel) {
    fun.map = furrr::future_map
  }else{
    fun.map = purrr::map
  }
  results = units.list %>% 
    set_names(units$unit[1]) %>% 
    fun.map(
      ~{
        item = .
        dependent = item$unit
        dependent.id = item$id
        args.TFDTW.synth[["dependent"]] = dependent
        args.TFDTW.synth[["dependent.id"]] = dependent.id
        args.TFDTW.synth[["res.synth.raw"]] = res.synth.raw.list[[dependent]]
        do.call(TFDTW.synth, args.TFDTW.synth)
      }
    )
  
  # compute log ratio
  mse = lapply(results, '[[', "mse") %>% 
    do.call("rbind", .) %>%
    mutate(ratio = mse.postT.TFDTW/mse.postT.raw,
           log.ratio = log(ratio))
  
  
  # output
  res.synth.target.raw = results[[target]]$res.synth.raw
  res.synth.target.TFDTW = results[[target]]$res.synth.TFDTW
  if (!detailed.output) {
    args.TFDTW.synth = NULL
    results = NULL
  }
  
  return(list(target = target,
              filter.width = filter.width,
              args.TFDTW.synth = args.TFDTW.synth,
              results.TFDTW.synth = results,
              res.synth.target.raw = res.synth.target.raw,
              res.synth.target.TFDTW = res.synth.target.TFDTW,
              mse = mse))
}


# simulate data
data = sim.data(n = 3, nCycles = 8, length = 100,
                t.treat = 60, shock = 0, ar.x = 0.8,
                n.SMA = 1, n.diff = 1,
                speed.upper = 2,
                speed.lower = 0.5,
                reweight = TRUE,
                rescale = TRUE,
                rescale.multiplier = 20,
                beta = 1)


count = 1
for (item in c("B", "C")) {
  value = data %>% filter(unit == item) %>% .[["value"]]
  for (lag in 1:5) {
    value_lag = c(rep(value[1], lag),  value[1:(100-lag)])
    data = rbind(data,
                 data.frame(id = 3 + count,
                            unit = paste0(item, "_lag", lag),
                            time = 1:100,
                            value = minmax.normalize(value_lag,
                                                     reference = value_lag[1:60])*20))
    count = count + 1
  }
  for (power in 2:5) {
    value_power = value %>% .^power
    data = rbind(data,
                 data.frame(id = 3 + count,
                            unit = paste0(item, "_pow", power),
                            time = 1:100,
                            value = minmax.normalize(value_power,
                                                     reference = value_power[1:60])*20))
    count = count + 1
  }
}

data$time[1:90] = 11:100
data$time[91:100] = 1:10
data = data %>% filter(time %in% 16:100)
data$value_raw = data$value

data %>% ggplot(aes(x = time, y = value, color = unit)) + geom_line()
# data = data %>% 
#   filter(unit %in% c("A", "B", "C"))


## estimate synthetic control
# (with lags) width = 9, k = 7, step.pattern = asymmetricP1
# (without lags) width = 11, k = 5, step.pattern = asymmetricP1
# parameters
filter.width.range = 9
k.range = 7
step.pattern.range = list(
  # symmetricP0 = dtw::symmetricP0, # too bumpy
  # symmetricP05 = dtw::symmetricP05,
  # symmetricP1 = dtw::symmetricP1,
  # symmetricP2 = dtw::symmetricP2,
  # asymmetricP0 = dtw::asymmetricP0, # too bumpy
  # asymmetricP05 = dtw::asymmetricP05,
  asymmetricP1 = dtw::asymmetricP1#,
  # asymmetricP2 = dtw::asymmetricP2,
  # typeIc = dtw::typeIc,
  # typeIcs = dtw::typeIcs,
  # typeIIc = dtw::typeIIc,  # jumps
  # typeIIIc = dtw::typeIIIc, # jumps
  # typeIVc = dtw::typeIVc,  # jumps
  # typeId = dtw::typeId,
  # typeIds = dtw::typeIds,
  # typeIId = dtw::typeIId, # jumps
  # mori2006 = dtw::mori2006
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
                  time.predictors.prior = 16:59,
                  time.optimize.ssr = 16:59)

args.TFDTW.synth = list(start.time = 16, end.time = 100, treat.time = 60,
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

args.TFDTW.synth.all.units[["data"]] = data
results = SimDesign::quiet(
  grid.search(filter.width.range = filter.width.range,
              k.range = k.range,
              step.pattern.range = step.pattern.range,
              args.TFDTW.synth.all.units = args.TFDTW.synth.all.units,
              grid.search.parallel = grid.search.parallel)
)

mse = lapply(results, "[[", "mse") %>% 
  do.call("rbind", .)

ind.opt = which(mse$ratio == min(mse$ratio, na.rm = TRUE))

synthetic.original = results[[ind.opt]]$res.synth.target.raw$synthetic
synthetic.dsc = results[[ind.opt]]$res.synth.target.TFDTW$synthetic
n.na = sum(is.na(synthetic.dsc))


## plot figure
df = rbind(
  data.frame(id = 1, unit = "Unit T", time = 151:1000,
             value = approx(data[1:85,4], n = 850)$y),
  data.frame(id = 2, unit = "Unit C1", time = 151:1000,
             value = approx(data[86:170,4], n = 850)$y),
  data.frame(id = 3, unit = "Unit C2", time = 151:1000,
             value = approx(data[171:255,4], n = 850)$y),
  data.frame(id = 4, unit = "SC", time = 151:1000,
             value = approx(synthetic.original, n = 850)$y),
  data.frame(id = 5, unit = "DSC", time = 151:1000,
             value = c(approx(synthetic.dsc, n = (85-n.na)*10)$y, rep(NA, n.na*10)))
)

df$value = df$value + rnorm(4250, mean = 0, sd = 0.1)

fig = df %>% 
  ggplot(aes(x = time, y = value, color = unit, linetype = unit)) +
  geom_line(size = 0.7) + 
  scale_linetype_manual(name = NULL,
                        values = c("Unit T" = "solid", "Unit C1" = "dashed",
                                   "Unit C2" = "dotted", "SC" = "solid",
                                   "DSC" = "solid")) +
  scale_color_manual(name = NULL,
                     values = c("Unit T" = "#4a4e4d", "Unit C1" = "#aaaaaa",
                                "Unit C2" = "#aaaaaa", "SC" = "#2ab7ca",
                                "DSC" = "#fe4a49")) +
  geom_vline(xintercept = 600, linetype="dashed", col = "grey20") +
  annotate("text", x = 590, y = 18.5,
           label = "Treatment", col = "grey20",
           angle = 90) +
  xlim(350, 750) +
  xlab("Time") +
  ylab("Y") +
  theme_minimal() +
  theme(legend.position=c(0.23,0.2), 
        legend.box = "horizontal",
        legend.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())

ggsave("./figures/speed_problem.pdf",
       fig, width = 6, height = 4.5,
       units = "in", limitsize = FALSE)



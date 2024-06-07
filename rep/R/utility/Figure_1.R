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


## 1. Speed Difference ---------------------------------------------------------
sim.data = function(n = 3, nCycles = 6, length = 100, extra.x = round(0.2*length),
                    t.treat = 60, shock = 0, ar.x = 0.6, n.SMA = 1,
                    n.diff = 1, speed.upper = 2, speed.lower = 0.5,
                    reweight = TRUE, rescale = TRUE,
                    rescale.multiplier = 20, beta = 1){
  # common exogenous shocks
  x = cos(seq(0, nCycles * pi, length.out = length + extra.x + n.SMA + n.diff - 1))
  x = x/2 + 0.5
  
  xt = 0
  x2 = NULL
  for (j in 1:length(x)) {
    xt = 0.8*xt + x[j]
    x2 = c(x2, xt)
  }
  x = x2
  
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
                            value = y,
                            value_raw = y))
  }
  return(data)
}

data = sim.data(n = 3, nCycles = 8, length = 100,
                t.treat = 60, shock = 0, ar.x = 0.6,
                n.SMA = 1, n.diff = 1,
                speed.upper = 2,
                speed.lower = 0.5,
                reweight = TRUE,
                rescale = TRUE,
                rescale.multiplier = 20,
                beta = 1)

data %>% ggplot(aes(x = time, y = value, color = unit)) + geom_line()


# parameters
filter.width.range = 19
k.range = 7
step.pattern.range = list(
  # symmetricP0 = dtw::symmetricP0, # too bumpy
  # symmetricP05 = dtw::symmetricP05,
  # symmetricP1 = dtw::symmetricP1,
  # symmetricP2 = dtw::symmetricP2,
  # asymmetricP0 = dtw::asymmetricP0, # too bumpy
  # asymmetricP05 = dtw::asymmetricP05,
  # asymmetricP1 = dtw::asymmetricP1,
  # asymmetricP2 = dtw::asymmetricP2,
  typeIc = dtw::typeIc#,
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

args.TFDTW.synth.all.units[["data"]] = data
results = SimDesign::quiet(
  grid.search(filter.width.range = filter.width.range,
              k.range = k.range,
              step.pattern.range = step.pattern.range,
              args.TFDTW.synth.all.units = args.TFDTW.synth.all.units,
              grid.search.parallel = grid.search.parallel)
)


df = rbind(
  data.frame(id = 1, unit = "Unit T", time = 1:1000,
             value = approx(data[1:100,4], n = 1000)$y),
  data.frame(id = 2, unit = "Unit C1", time = 1:1000,
             value = approx(data[101:200,4], n = 1000)$y),
  data.frame(id = 3, unit = "Unit C2", time = 1:1000,
             value = approx(data[201:300,4], n = 1000)$y),
  data.frame(id = 4, unit = "SC", time = 1:1000,
             value = approx(results[[1]]$res.synth.target.raw$synthetic, n = 1000)$y),
  data.frame(id = 5, unit = "DSC", time = 1:1000,
             value = c(approx(results[[1]]$res.synth.target.TFDTW$synthetic, n = 770)$y, rep(NA, 230)))
)


df$unit = factor(df$unit, levels = c("Unit T", "Unit C1", "Unit C2",
                                     "SC","DSC"))
df$time[c(1:1000, 4001:5000)] = df$time[c(1:1000, 4001:5000)] + 100
df$value = df$value + rnorm(1000, mean = 0, sd = 0.1)

fig = df %>% 
  ggplot(aes(x = time, y = value, color = unit, linetype = unit)) +
  geom_line(size = 0.7) + 
  scale_linetype_manual(name = NULL,
                        labels = c("Unit T" = expression(paste("Unit ", bold(y)[1])),
                                   "Unit C1" = expression(paste("Unit ", bold(y)[2])),
                                   "Unit C2" = expression(paste("Unit ", bold(y)[3])),
                                   "SC" = "SC", "DSC" = "DSC"),
                        values = c("Unit T" = "solid", "Unit C1" = "dashed",
                                   "Unit C2" = "dotted", "SC" = "solid",
                                   "DSC" = "solid")) +
  scale_color_manual(name = NULL,
                     labels = c("Unit T" = expression(paste("Unit ", bold(y)[1])),
                                "Unit C1" = expression(paste("Unit ", bold(y)[2])),
                                "Unit C2" = expression(paste("Unit ", bold(y)[3])),
                                "SC" = "SC", "DSC" = "DSC"),
                     values = c("Unit T" = "#4a4e4d", "Unit C1" = "#aaaaaa",
                                "Unit C2" = "#aaaaaa", "SC" = "#3da4ab",
                                "DSC" = "#fe8a71")) +
  geom_vline(xintercept = 600, linetype="dashed", col = "grey20") +
  annotate("segment", x = 730, y = df$value[3730]+0.1,
           xend = 730, yend = df$value[630],
           arrow = arrow(ends = "both", length = unit(.2,"cm")),
           colour = "grey20", size = 0.8) +
  annotate("label", x = 680, y = 18, size = 3,
           label = "SC estimated\nTreatment Effect", col = "grey20") +
  annotate("text", x = 590, y = 18.5,
           label = "Treatment", col = "grey20",
           angle = 90) +
  xlim(350, 750) +
  xlab("Time") +
  ylab("Y") +
  theme_minimal() +
  theme(legend.position=c(0.23,0.2), 
        legend.box = "horizontal",
        legend.text.align = 0,
        legend.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())

ggsave("./figures/Figure_1.pdf",
       fig, width = 6, height = 4.5,
       units = "in", limitsize = FALSE)



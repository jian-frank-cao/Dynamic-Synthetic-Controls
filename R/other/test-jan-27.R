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

sim.data = function(n = 10, length = 100, extra.x = round(0.2*length),
                    t.treat = 60, shock = 10, nCycles = 2, n.SMA = 1,
                    n.diff = 1, speed.upper = 2, speed.lower = 0.5,
                    reweight = TRUE, rescale = TRUE, 
                    rescale.multiplier = 20, beta = 1){
  # common exogenous shocks
  x = sin(seq(0, ((extra.x+length)/length)*nCycles * pi, length.out = length + extra.x))
  
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
  rnd.ind = sample(c(1:round(0.3*n), round(0.7*n):n), size = 1)
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
    
    data = rbind(data,
                 data.frame(id = i,
                            unit = LETTERS[i],
                            time = 1:length,
                            value = y,
                            value_raw = y))
  }
  return(data)
}


n = 3
length = 20
extra.x = 0
t.treat = 12
shock = 0
nCycles = 2
n.SMA = 1
n.diff = 1
speed.upper = 3/2
speed.lower = 2/3
reweight = FALSE
rescale = FALSE
rescale.multiplier = 1
beta = 1

set.seed(1)
data = sim.data(n = n, length = length, extra.x = extra.x,
                t.treat = t.treat, shock = shock, nCycles = nCycles, n.SMA = n.SMA,
                n.diff = n.diff, speed.upper = speed.upper, speed.lower = speed.lower,
                reweight = reweight, rescale = rescale, 
                rescale.multiplier = rescale.multiplier, beta = beta)

data %>%
  ggplot(aes(x = time, y = value, color = unit)) +
  geom_line()

df = reshape2::dcast(data, time ~ unit)
df = df %>% 
  mutate(B.sq = B^2,
         C.sq = C^2)

model1 = lm(A ~ B + C, data = df %>% filter(time < 12))
lm.predict1 = predict(model1, df[,3:4])

model2 = lm(A ~ B + C + B.sq + C.sq, data = df %>% filter(time < 12))
lm.predict2 = predict(model2, df[,3:6])

plot(ts(df$A))
lines(1:20, df$B, col = "blue")
lines(1:20, df$C, col = "green")
# lines(1:20, lm.predict1, col = "red")
# lines(1:20, lm.predict2, col = "orange")
lines(1:20, results[[18]]$res.synth.target.raw$synthetic, col = "purple")
lines(1:20, results[[18]]$res.synth.target.TFDTW$synthetic, col = "pink")


# -------------------------------------
filter.width.range = (1:3)*2+3
k.range = 4:7
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


args.TFDTW = list(buffer = 5, match.method = "open.end",
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
                    expression(list(list(dep.var, 9:11, c("mean")),
                                    list(dep.var, 6:8, c("mean")),
                                    list(dep.var, 3:5, c("mean")))),
                  time.predictors.prior = 1:11,
                  time.optimize.ssr = 1:11)

args.TFDTW.synth = list(start.time = 1, end.time = 20, treat.time = 12,
                        args.TFDTW = args.TFDTW, args.synth = args.synth,
                        ## 2nd
                        n.mse = 5, 
                        ## other
                        plot.figures = FALSE,
                        plot.path = "./figures/",
                        legend.pos = c(0.3, 0.7))

args.TFDTW.synth.all.units = list(target = "A",
                                  data = data,
                                  args.TFDTW.synth = args.TFDTW.synth,
                                  ## 2nd
                                  all.units.parallel = FALSE)

filter.width.range = filter.width.range
k.range = k.range
step.pattern.range = step.pattern.range
args.TFDTW.synth.all.units = args.TFDTW.synth.all.units
grid.search.parallel = grid.search.parallel

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

results = NULL

for (i in 33:length(search.grid.list)) {
  print(i)
  task = search.grid.list[[i]]
  args.TFDTW.synth.all.units[["filter.width"]] = task$filter.width
  args.TFDTW.synth.all.units$args.TFDTW.synth$args.TFDTW[["k"]] = task$k
  args.TFDTW.synth.all.units$args.TFDTW.synth$args.TFDTW[["step.pattern1"]] =
    step.pattern.range[[task$step.pattern]]
  args.TFDTW.synth.all.units[["res.synth.raw.list"]] = res.synth.raw.list
  results[[i]] = do.call(TFDTW.synth.all.units, args.TFDTW.synth.all.units)
}

res = results %>% 
  future_map(
    ~{
      item = .
      value = item$res.synth.target.raw$value
      synth.sc = item$res.synth.target.raw$synthetic
      synth.dsc = item$res.synth.target.TFDTW$synthetic
      gap.sc = value - synth.sc
      gap.dsc = value - synth.dsc
      n.na = sum(is.na(gap.dsc))
      mse.dsc = mean((gap.dsc[12:(20-n.na)])^2)
      mse.sc = mean((gap.sc[12:(20-n.na)])^2)
      data.frame(n.na = n.na,
                 mse.sc = mse.sc,
                 mse.dsc = mse.dsc)
    }
  ) %>% do.call("rbind", .)


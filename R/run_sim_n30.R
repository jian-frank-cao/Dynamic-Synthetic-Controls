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


## Function --------------------------------------------------------------------
sim.data = function(n = 10, length = 100, extra.x = round(0.2*length),
                    t.treat = 60, shock = 10, arima.order = c(1,1,0),
                    ar.x = 0.6, ma.x = NULL, n.SMA = 1, n.diff = 1,
                    speed.upper = 2, speed.lower = 0.5,
                    treat.last = 0.1, reweight = TRUE, rescale = TRUE,
                    rescale.multiplier = 20, beta = 1){
  # common exogenous shocks
  x = arima.sim(list(order = arima.order, ar = ar.x, ma = ma.x),
                n = length + extra.x + n.SMA + n.diff - 2)
  
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
                    seq(0, shock, length.out = round(treat.last*length)),
                    rep(shock, round((1 - treat.last)*length - t.treat)))
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



## Data Simulation -------------------------------------------------------------
n.simulation = 150
length = 30
n = 10
beta = 1
shock = 10

# simulate
data.list = NULL
for (i in 1:n.simulation) {
  data.list[[i]] = sim.data(n = n, length = length,
                      t.treat = 20, shock = shock, ar.x = 0.6,
                      n.SMA = 1, n.diff = 1,
                      speed.upper = 3,
                      speed.lower = 1/3,
                      treat.last = 0.1,
                      reweight = TRUE,
                      rescale = TRUE,
                      rescale.multiplier = 10,
                      beta = beta)
  # data.list[[i]] = data.sim %>%
  #   filter(time %in% c((1:20)*5)) %>%
  #   group_by(unit) %>%
  #   mutate(time = 1:20) %>%
  #   ungroup %>% 
  #   data.frame
}

saveRDS(data.list, paste0("./data/simul_data_beta1_n30.Rds"))

# data.list[[5]] %>%
#   ggplot(aes(x = time, y = value, color = unit)) +
#   geom_line()

## Run -------------------------------------------------------------------------
data.list = readRDS("./data/simul_data_beta1_n30.Rds")

# parameters
filter.width.range = (1:3)*2+3
k.range = 4:6
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


args.TFDTW = list(buffer = 3, match.method = "open.end",  # check this
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
                    expression(list(list(dep.var, 16:19, c("mean")),
                                    list(dep.var, 12:15, c("mean")),
                                    list(dep.var, 8:11, c("mean")))),
                  time.predictors.prior = 1:19,
                  time.optimize.ssr = 1:19)

args.TFDTW.synth = list(start.time = 1, end.time = 30, treat.time = 20,
                        args.TFDTW = args.TFDTW, args.synth = args.synth,
                        ## 2nd
                        n.mse = 5, 
                        ## other
                        plot.figures = FALSE,
                        plot.path = "./figures/",
                        legend.pos = c(0.3, 0.7))

args.TFDTW.synth.all.units = list(target = "A",
                                  # data = data, 
                                  args.TFDTW.synth = args.TFDTW.synth,
                                  ## 2nd
                                  all.units.parallel = FALSE)

for (i in 87:length(data.list)) {
  cat("Simulation data set ", i, "...")
  args.TFDTW.synth.all.units[["data"]] = data.list[[i]]
  set.seed(20220407)
  results = SimDesign::quiet(
    grid.search(filter.width.range = filter.width.range,
                k.range = k.range,
                step.pattern.range = step.pattern.range,
                args.TFDTW.synth.all.units = args.TFDTW.synth.all.units,
                grid.search.parallel = grid.search.parallel)
  )
  saveRDS(results, paste0("./data/res_sim/n30/res_sim_n30_", i, ".Rds"))
  cat("Done.\n")
}


## Test result -----------------------------------------------------------------
folder = paste0("./data/res_sim/n30/")
file.list = as.list(list.files(folder))

length = 20
treat_time = 12
n_mse = 5
shock = 10
treatment = c(rep(0, treat_time),
              seq(0, shock, length.out = round(0.15*length)),
              rep(shock, round(0.85*length - treat_time)))

pre.start = 7
pre.end = 12
post.start = 13
post.end = 17

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

saveRDS(df.mse, paste0("./data/df.mse_sim_n20_2.Rds"))

# t.test for log(MSEdsc/MSEsc)
df.mse = readRDS(paste0("./data/df.mse_sim_n20_2.Rds"))
df.mse = df.mse %>%
  mutate(log.ratio = log(mse.postT.TFDTW/mse.postT.raw))
t.test(df.mse$log.ratio)
sum(df.mse$log.ratio<0)/nrow(df.mse)
wilcox.test(df.mse$log.ratio)

## Plot result -----------------------------------------------------------------
# df.gap
df.mse = readRDS(paste0("./data/df.mse_sim_n20_2.Rds"))
folder = paste0("./data/res_sim/n20_2/")
file.list = as.list(list.files(folder))
results = file.list %>%
  future_map(
    ~{
      file.name = .
      readRDS(paste0(folder, file.name))
    }
  )

df.gap = NULL
for (i in 1:nrow(df.mse)) {
  data.id = df.mse$data.id[i]
  grid.id = df.mse$grid.id[i]
  df.gap[[i]] = data.frame(
    time = 1:20,
    data.id = data.id,
    grid.id = grid.id,
    value = results[[data.id]][[grid.id]][["res.synth.target.raw"]][[1]],
    synth.sc = results[[data.id]][[grid.id]][["res.synth.target.raw"]][[3]],
    synth.dsc = results[[data.id]][[grid.id]][["res.synth.target.TFDTW"]][[3]]
  )
  print(i)
}

df.gap = df.gap %>%
  do.call("rbind", .) %>%
  mutate(
    gap.sc = value - synth.sc,
    gap.dsc = value - synth.dsc,
    group = paste0(data.id, "-", grid.id)
  )

saveRDS(df.gap, paste0("./data/df.gap_sim_n20_2.Rds"))

# plot
df.gap = readRDS(paste0("./data/df.gap_sim_n20_2.Rds"))

shock = 10
length = 20
treat_time = 12
n_mse = 5
treatment = c(rep(0, treat_time),
              seq(0, shock, length.out = round(0.15*length)),
              rep(shock, round(0.85*length - treat_time)))

df.quantile = df.gap %>%
  group_by(time) %>%
  summarise(mean.sc = mean(gap.sc, na.rm = T),
            mean.dsc = mean(gap.dsc, na.rm = T),
            quantile.sc.975 = quantile(gap.sc, 0.975, na.rm = T),
            quantile.sc.025 = quantile(gap.sc, 0.025, na.rm = T),
            quantile.dsc.975 = quantile(gap.dsc, 0.975, na.rm = T),
            quantile.dsc.025 = quantile(gap.dsc, 0.025, na.rm = T)) %>%
  mutate(group = "quantile",
         treatment = treatment)

# df.quantile[95:100, c("mean.dsc", "quantile.dsc.975", "quantile.dsc.025")] = NA


color.sc = "#2ab7ca"
color.dsc = "#fe4a49"
# color.sc = "grey70"
# color.dsc = "grey30"

colors = c("TE" = "grey20",
           "Avg. TE (SC)" = color.sc,
           "Avg. TE (DSC)" = color.dsc)

linetypes = c("TE" = "solid",
              "Avg. TE (SC)" = "dashed",
              "Avg. TE (DSC)" = "dashed")

fills = c("95% Quantile (SC)" = color.sc,
          "95% Quantile (DSC)" = color.dsc)

fig.big = df.gap %>%
  ggplot(aes(x = time, group = group)) +
  annotate("rect", xmin = 12, xmax = 17,
           ymin = -25, ymax = 35, alpha = .3) +
  geom_line(aes(y = gap.sc), col = color.sc, alpha=0.1) +
  geom_line(aes(y = gap.dsc), col = color.dsc, alpha=0.1) +
  geom_ribbon(aes(ymin = quantile.sc.025, ymax = quantile.sc.975,
                  fill = "95% Quantile (SC)"), data = df.quantile, alpha=0.6) +
  geom_ribbon(aes(ymin = quantile.dsc.025, ymax = quantile.dsc.975,
                  fill = "95% Quantile (DSC)"), data = df.quantile, alpha=0.6) +
  geom_line(aes(x = time, y = treatment, color = "TE", linetype = "TE"),
            data = df.quantile, size = 1) +
  geom_line(aes(x = time, y = mean.sc, color = "Avg. TE (SC)", linetype = "Avg. TE (SC)"),
            data = df.quantile, size = 0.6) +
  geom_line(aes(x = time, y = mean.dsc, color = "Avg. TE (DSC)", linetype = "Avg. TE (DSC)"),
            data = df.quantile, size = 0.6) +
  scale_color_manual(name = NULL, values = colors) +
  scale_fill_manual(name = NULL, values = fills) +
  scale_linetype_manual(name = NULL, values = linetypes) +
  geom_vline(xintercept = 12, linetype="dashed", col = "grey20") +
  geom_hline(yintercept = 0, linetype="dashed", col = "grey20") +
  annotate("text", x = 11, y = 25, label = "Treatment",
           col = "grey20", angle = 90) +
  coord_cartesian(ylim = c(-20, 30)) +
  xlab("Time") +
  ylab("Treatment Effect (TE)") +
  theme_bw() +
  theme(legend.position=c(0.3,0.15),
        legend.box = "horizontal",
        legend.background = element_rect(fill=NA))


# df.t.test = data.frame(Beta = c(0, 0.5, 1),
#                        t = c(-4.7862, -5.8137, -8.2442),
#                        P = c(0.0001, 0.0001, 0.0001))
#
# fig.small = df.t.test %>%
#   ggplot(aes(x = Beta, y = t)) +
#   annotate("rect", xmin = -0.2, xmax = 1.2,
#            ymin = -15, ymax = -2.63, alpha = .3) +
#   geom_line(size = 1) +
#   geom_point(size = 2, col = color.dsc) +
#   geom_hline(yintercept = 0, linetype="solid", col = "black") +
#   geom_hline(yintercept = -2.63, linetype="dashed", col = "grey20") +
#   annotate("text", x = 0.75, y = -4,
#            label = "P < 0.01", col = "grey20") +
#   scale_x_continuous(breaks = c(0, 0.5, 1),
#                      labels = c("0", "0.5", "1")) +
#   scale_y_continuous(breaks = c(-10, -5, 0)) +
#   coord_cartesian(ylim = c(-10, 0), xlim = c(0, 1)) +
#   xlab(expression(psi)) +
#   ylab("t") +
#   theme_bw() +
#   theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
#
#
# fig.sim = fig.big + annotation_custom(ggplotGrob(fig.small),
#                                       xmin = 5, xmax = 45,
#                                       ymin = 7, ymax = 30)

ggsave(paste0("./figures/sim_n20_2.pdf"),
       fig.big, width = 6, height = 4.5,
       units = "in", limitsize = FALSE)

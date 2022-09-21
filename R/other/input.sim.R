# data
data_list = readRDS("./data/simul_data_list_0916.Rds")
# data = data_list[[694]]

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


args.TFDTW = list(## 2nd
                  buffer = 0, match.method = "fixed",
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
                  special.predictors = expression(list(list(dep.var, 70:79, c("mean")),
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
                                  all.units.parallel = FALSE)

results = NULL
for (i in 1:5) {
  cat("Simulation data set ", i, "...")
  args.TFDTW.synth.all.units[["data"]] = data_list[[i]]
  results[[i]] = SimDesign::quiet(grid.search(filter.width.range = filter.width.range,
                                              k.range = k.range,
                                              step.pattern.range = step.pattern.range,
                                              args.TFDTW.synth.all.units = args.TFDTW.synth.all.units,
                                              grid.search.parallel = TRUE))
  cat("Done.\n")
}

saveRDS(results, "./data/res_sim_0916.Rds")

res = lapply(results[[1]], "[[", "neg.ratio") %>% do.call("c", .)

ind = as.numeric(which(res == max(res)))

p.values = lapply(results[[1]][ind], "[[", "p.value") %>% do.call("c", .)


# # # check res.synth.raw.list
# # res = res.synth.raw.list$E
# # plot(ts(res$value))
# # lines(res$average, col = "blue")
# # lines(res$synthetic, col = "red")
# 
# 
# # TFDTW.synth.all.units
# data = args.TFDTW.synth.all.units$data
# target = args.TFDTW.synth.all.units$target
# args.TFDTW.synth = args.TFDTW.synth.all.units$args.TFDTW.synth
# all.units.parallel = args.TFDTW.synth.all.units$all.units.parallel
# filter.width = args.TFDTW.synth.all.units$filter.width
# res.synth.raw.list = args.TFDTW.synth.all.units$res.synth.raw.list
# 
# 
# # TFDTW.synth
# start.time = args.TFDTW.synth$start.time
# end.time = args.TFDTW.synth$end.time
# treat.time = args.TFDTW.synth$treat.time
# args.TFDTW = args.TFDTW.synth$args.TFDTW
# args.synth = args.TFDTW.synth$args.synth
# n.mse = args.TFDTW.synth$n.mse
# plot.figures = args.TFDTW.synth$plot.figures
# plot.path = args.TFDTW.synth$plot.path
# legend.pos = args.TFDTW.synth$legend.pos
# data = args.TFDTW.synth$data
# dependent = args.TFDTW.synth$dependent
# dependent.id = args.TFDTW.synth$dependent.id
# res.synth.raw = args.TFDTW.synth$res.synth.raw
# 
# 
# # TFDTW
# buffer = args.TFDTW$buffer
# match.method = args.TFDTW$match.method
# dist.quant = args.TFDTW$dist.quant
# step.pattern2 = args.TFDTW$step.pattern2
# n.burn = args.TFDTW$n.burn
# n.IQR = args.TFDTW$n.IQR
# ma = args.TFDTW$ma
# ma.na = args.TFDTW$ma.na
# default.margin = args.TFDTW$default.margin
# n.q = args.TFDTW$n.q
# n.r = args.TFDTW$n.r
# k = args.TFDTW$k
# step.pattern1 = args.TFDTW$step.pattern1
# y = args.TFDTW$y
# t.treat = args.TFDTW$t.treat
# plot.figures = args.TFDTW$plot.figures
# x = args.TFDTW$x
# 
# 
# # synth
# predictors = args.synth$predictors
# special.predictors = args.synth$special.predictors
# time.predictors.prior = args.synth$time.predictors.prior
# time.optimize.ssr = args.synth$time.optimize.ssr
# df = args.synth$df
# dependent.id = args.synth$dependent.id
# dep.var = args.synth$dep.var

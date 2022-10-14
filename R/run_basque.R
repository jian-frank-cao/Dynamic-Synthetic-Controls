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
source("./R/grid.search.R")
set.seed(20220407)


## Basque Terrorism Data -------------------------------------------------------
data(basque, package = "Synth")
data = basque
colnames(data)[1:4] = c("id", "unit", "time", "value")
data = data %>% mutate(invest_ratio = invest/value,
                       value_raw = value)

# rescale
df.rescale = data %>% 
  filter(time <= 1990) %>% 
  group_by(unit) %>% 
  summarise(value.min = min(value),
            value.max = max(value),
            multiplier = 5/(value.max - value.min)) %>% 
  ungroup()

data = left_join(data, df.rescale, by = "unit")
data = data %>% mutate(value.bak = value_raw,
                       value_raw = (value_raw - value.min)*multiplier,
                       value = value_raw)


## Grid Search Basque ----------------------------------------------------------
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
                  window.type = "none",
                  ## other
                  norm.method = "t",
                  step.pattern2 = dtw::asymmetricP2,
                  n.burn = 3, n.IQR = 3,
                  ma = 3, ma.na = "original",
                  default.margin = 3,
                  n.q = 1, n.r = 1)

args.synth = list(predictors = NULL,
                  special.predictors = 
                    expression(list(
                      list(dep.var, 1960:1969, c("mean")),
                      list("invest_ratio", 1964:1969, c("mean")),
                      list("popdens", 1969, c("mean")),
                      list("sec.agriculture", 1961:1969, c("mean")),
                      list("sec.energy", 1961:1969, c("mean")),
                      list("sec.industry", 1961:1969, c("mean")),
                      list("sec.construction", 1961:1969, c("mean")),
                      list("sec.services.venta", 1961:1969, c("mean")),
                      list("sec.services.nonventa", 1961:1969, c("mean")),
                      list("school.illit", 1964:1969, c("mean")),
                      list("school.prim", 1964:1969, c("mean")),
                      list("school.med", 1964:1969, c("mean")),
                      list("school.high", 1964:1969, c("mean")),
                      list("school.post.high", 1964:1969, c("mean"))
                    )),
                  time.predictors.prior = 1955:1969,
                  time.optimize.ssr = 1955:1969)

args.TFDTW.synth = list(start.time = 1955, end.time = 1990, treat.time = 1970,
                        args.TFDTW = args.TFDTW, args.synth = args.synth,
                        ## 2nd
                        n.mse = 10, 
                        ## other
                        plot.figures = FALSE,
                        plot.path = "./figures/",
                        legend.pos = c(0.3, 0.3))

args.TFDTW.synth.all.units = list(target = "Basque Country (Pais Vasco)",
                                  # data = data, 
                                  args.TFDTW.synth = args.TFDTW.synth,
                                  ## 2nd
                                  all.units.parallel = FALSE)

cat("Start...")
args.TFDTW.synth.all.units[["data"]] = data
results = SimDesign::quiet(
  grid.search(filter.width.range = filter.width.range,
              k.range = k.range,
              step.pattern.range = step.pattern.range,
              args.TFDTW.synth.all.units = args.TFDTW.synth.all.units,
              grid.search.parallel = grid.search.parallel)
)
cat("Done.\n")

saveRDS(results, paste0("./data/res_basque_1013", i, ".Rds"))
job.end = Sys.time()
print(job.end - job.start)


## Optimal Run Tobacco ---------------------------------------------------------
# # prepare data
# start_time = 1955
# end_time = 1990
# treat_time = 1970
# dtw1_time = 1970
# filter_width = 5
# k = 5
# type_dtw1 = "fixed"
# n_mse = 10
# n_IQR = 3
# dist_quantile = 0.95
# plot_figures = TRUE
# step.pattern1 = dtw::mori2006
# step.pattern2 = dtw::asymmetricP2
# predictors.origin = NULL
# special.predictors.origin = list(
#   list("value_raw", 1960:1969, c("mean")),
#   list("invest_ratio", 1964:1969, c("mean")),
#   list("popdens", 1969, c("mean")),
#   list("sec.agriculture", 1961:1969, c("mean")),
#   list("sec.energy", 1961:1969, c("mean")),
#   list("sec.industry", 1961:1969, c("mean")),
#   list("sec.construction", 1961:1969, c("mean")),
#   list("sec.services.venta", 1961:1969, c("mean")),
#   list("sec.services.nonventa", 1961:1969, c("mean")),
#   list("school.illit", 1964:1969, c("mean")),
#   list("school.prim", 1964:1969, c("mean")),
#   list("school.med", 1964:1969, c("mean")),
#   list("school.high", 1964:1969, c("mean")),
#   list("school.post.high", 1964:1969, c("mean"))
# )
# time.predictors.prior.origin = 1955:1969
# time.optimize.ssr.origin = 1955:1969
# predictors.new = NULL
# special.predictors.new = list(
#   list("value_warped", 1960:1969, c("mean")),
#   list("invest_ratio", 1964:1969, c("mean")),
#   list("popdens", 1969, c("mean")),
#   list("sec.agriculture", 1961:1969, c("mean")),
#   list("sec.energy", 1961:1969, c("mean")),
#   list("sec.industry", 1961:1969, c("mean")),
#   list("sec.construction", 1961:1969, c("mean")),
#   list("sec.services.venta", 1961:1969, c("mean")),
#   list("sec.services.nonventa", 1961:1969, c("mean")),
#   list("school.illit", 1964:1969, c("mean")),
#   list("school.prim", 1964:1969, c("mean")),
#   list("school.med", 1964:1969, c("mean")),
#   list("school.high", 1964:1969, c("mean")),
#   list("school.post.high", 1964:1969, c("mean"))
# )
# time.predictors.prior.new = 1955:1969
# time.optimize.ssr.new = 1955:1969
# legend_position = c(0.3, 0.3)
# ... = NULL
# normalize_method = "t"
# ma = 3
# ma_na = "original"
# n_q = 1
# n_r = 1
# 
# 
# data = preprocessing(data, filter_width)
# units = data[c("id", "unit")] %>% distinct
# 
# # run
# result = as.list(1:nrow(units)) %>% 
#   future_map(
#     ~{
#       i = .
#       dependent = units$unit[i]
#       dependent_id = units$id[i]
#       # print(paste0(dependent, ":", i, "-", k, " start..."))
#       res = SimDesign::quiet(compare_methods(data = data,
#                                              start_time = start_time,
#                                              end_time = end_time,
#                                              treat_time = treat_time,
#                                              dtw1_time = dtw1_time,
#                                              dependent = dependent,
#                                              dependent_id = dependent_id,
#                                              n_mse = n_mse,
#                                              k = k,
#                                              type_dtw1 = type_dtw1,
#                                              n_IQR = n_IQR,
#                                              dist_quantile = dist_quantile,
#                                              plot_figures = plot_figures,
#                                              step.pattern1 = step.pattern1,
#                                              step.pattern2 = step.pattern2,
#                                              predictors.origin = predictors.origin,
#                                              special.predictors.origin = special.predictors.origin,
#                                              time.predictors.prior.origin = time.predictors.prior.origin,
#                                              time.optimize.ssr.origin = time.optimize.ssr.origin,
#                                              predictors.new = predictors.new,
#                                              special.predictors.new = special.predictors.new,
#                                              time.predictors.prior.new = time.predictors.prior.new,
#                                              time.optimize.ssr.new = time.optimize.ssr.new,
#                                              legend_position = legend_position))
#       # print(paste0(dependent, ":", i, "-", k, " start...Done."))
#       res$mse = res$mse %>% mutate(dependent = dependent, k = k)
#       res
#     }
#   )
# 
# saveRDS(result, "./data/grid_search_v6/result_basque_70_fixed.Rds")
# 
# 
# # mse
# mse = result %>% 
#   lapply(., "[[", "mse") %>% 
#   do.call("rbind", .) %>% 
#   mutate(ratio = mse2_post/mse1_post,
#          log_ratio = log(ratio))
# mse = mse %>% filter(dependent != "Basque Country (Pais Vasco)")
# length(which(mse$log_ratio < 0))/nrow(mse)
# boxplot(mse$log_ratio, outline = FALSE)
# abline(h = 0, lty = 5)
# 
# t.test(mse$log_ratio)
# 
# saveRDS(mse, "./data/grid_search_v6/mse_basque_70.Rds")
# 
# 
# # plot time series figure
# df = rbind(data.frame(unit = "Basque Country",
#                       time = 1955:1990,
#                       value = result[[4]]$synth_origin$value),
#            data.frame(unit = "Synthetic Control w/o TFDTW",
#                       time = 1955:1990,
#                       value = result[[4]]$synth_origin$synthetic),
#            data.frame(unit = "Synthetic Control w/ TFDTW",
#                       time = 1955:1990,
#                       value = result[[4]]$synth_new$synthetic))
# 
# fig = ggplot(df, aes(x = time, y = value, color = unit)) +
#   geom_line() + 
#   geom_vline(xintercept = 1970, linetype="dashed") +
#   theme_bw() +
#   scale_color_manual(values=c("grey40", "#fe4a49","#2ab7ca")) +
#   xlab("Time") +
#   ylab("Per Capita GDP (1986 USD Thousand)") + 
#   theme(legend.title=element_blank(),
#         legend.position = c(0.8, 0.2))
# 
# saveRDS(fig, "./data/grid_search_v6/fig_comp_method_basque_70.Rds")
# 
# 
# # plot placebo figure
# df = result %>% 
#   map(
#     ~{
#       data.frame(unit = .[["mse"]][["dependent"]],
#                  time = 1955:1990,
#                  value = .[["synth_origin"]][["value"]],
#                  synth_origin = .[["synth_origin"]][["synthetic"]],
#                  synth_new = .[["synth_new"]][["synthetic"]])
#     }
#   ) %>% 
#   do.call("rbind", .) %>% 
#   mutate(
#     color = case_when(unit == "Basque Country (Pais Vasco)" ~ "black",
#                       TRUE ~ "grey 70"),
#     gap_origin = value - synth_origin,
#     gap_new = value - synth_new
#   )
# 
# df %>% 
#   filter(unit %in% (mse %>% filter(mse1_pre < 2*0.001) %>% .[["dependent"]])) %>% 
#   ggplot(aes(x = time, group = unit)) +
#   geom_line(aes(y = gap_origin), col = "#adcbe3") +
#   geom_line(aes(y = gap_new), col = "#fec8c1") +
#   geom_line(aes(y = gap_origin), data = df %>% filter(unit == "Basque Country (Pais Vasco)"), col = "#2ab7ca", size = 1) +
#   geom_line(aes(y = gap_new), data = df %>% filter(unit == "Basque Country (Pais Vasco)"), col = "#fe4a49", size = 1) +
#   geom_vline(xintercept = 1970, linetype="dashed") +
#   geom_hline(yintercept = 0, linetype="dashed") +
#   coord_cartesian(ylim=c(-1.2, 1.2)) +
#   xlab("year") +
#   ylab("Per Capita GDP") +
#   theme_bw()









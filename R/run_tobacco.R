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


## California Tobacco Data -----------------------------------------------------
load("./data/smoking.rda")
prop99 = read.csv("./data/prop99.csv")

exclude_states = c("Massachusetts", "Arizona", "Oregon", "Florida",
                   "Alaska", "Hawaii", "Maryland", "Michigan",
                   "New Jersey", "New York",
                   "Washington", "District of Columbia")
include_states = sort(setdiff(unique(prop99$LocationDesc),
                              exclude_states))
states = data.frame(id = 1:length(include_states),
                    unit = include_states)
smoking = smoking %>% mutate_all(as.numeric)
colnames(smoking)[1:3] = c("id", "time", "value")
smoking = right_join(states, smoking, by = "id")
smoking = smoking %>%
  mutate(value_raw = value,
         age15to24 = age15to24*100)

data = smoking

# rescale
df.rescale = data %>% 
  filter(time <= 1989) %>% 
  group_by(unit) %>% 
  summarise(value.min = min(value),
            value.max = max(value),
            multiplier = 45/(value.max - value.min)) %>% 
  ungroup()

data = left_join(data, df.rescale, by = "unit")
data = data %>% mutate(value.bak = value_raw,
                       value_raw = (value_raw - value.min)*multiplier,
                       value = value_raw)


## Grid Search Tobacco ---------------------------------------------------------
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
                      list(dep.var, 1988, c("mean")),
                      list(dep.var, 1980, c("mean")),
                      list(dep.var, 1975, c("mean")),
                      list("beer", 1984:1988, c("mean")),
                      list("lnincome", 1980:1988, c("mean")),
                      list("age15to24", 1980:1988, c("mean")),
                      list("retprice", 1980:1988, c("mean"))
                    )),
                  time.predictors.prior = 1970:1988,
                  time.optimize.ssr = 1970:1988)

args.TFDTW.synth = list(start.time = 1970, end.time = 2000, treat.time = 1989,
                        args.TFDTW = args.TFDTW, args.synth = args.synth,
                        ## 2nd
                        n.mse = 10, 
                        ## other
                        plot.figures = FALSE,
                        plot.path = "./figures/",
                        legend.pos = c(0.3, 0.3))

args.TFDTW.synth.all.units = list(target = "California",
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

saveRDS(results, "./data/res_tobacco_1013.Rds")
job.end = Sys.time()
print(job.end - job.start)

## Result ----------------------------------------------------------------------
results = readRDS("./data/res_tobacco_1013.Rds")
neg.ratio = lapply(results, "[[", "neg.ratio") %>% 
  do.call("c", .)
p.value = lapply(results, "[[", "p.value") %>% 
  do.call("c", .)
max.neg = which(neg.ratio == max(neg.ratio))
min.p = which(p.value == min(p.value))

neg.ratio[max.neg]
p.value[max.neg]

neg.ratio[min.p]
p.value[min.p]

## Optimal Run Tobacco ---------------------------------------------------------
# parameters
filter.width.range = 13
k.range = 8
step.pattern.range = list(
  # symmetricP0 = dtw::symmetricP0, # too bumpy
  # symmetricP05 = dtw::symmetricP05,
  # symmetricP1 = dtw::symmetricP1#,
  # symmetricP2 = dtw::symmetricP2,
  # asymmetricP0 = dtw::asymmetricP0, # too bumpy
  # asymmetricP05 = dtw::asymmetricP05,
  # asymmetricP1 = dtw::asymmetricP1,
  # asymmetricP2 = dtw::asymmetricP2,
  # typeIc = dtw::typeIc,
  # typeIcs = dtw::typeIcs,
  # typeIIc = dtw::typeIIc,  # jumps
  # typeIIIc = dtw::typeIIIc, # jumps
  # typeIVc = dtw::typeIVc,  # jumps
  typeId = dtw::typeId#,
  # typeIds = dtw::typeIds,
  # typeIId = dtw::typeIId, # jumps
  # mori2006 = dtw::mori2006
)
grid.search.parallel = FALSE


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
                      list(dep.var, 1988, c("mean")),
                      list(dep.var, 1980, c("mean")),
                      list(dep.var, 1975, c("mean")),
                      list("beer", 1984:1988, c("mean")),
                      list("lnincome", 1980:1988, c("mean")),
                      list("age15to24", 1980:1988, c("mean")),
                      list("retprice", 1980:1988, c("mean"))
                    )),
                  time.predictors.prior = 1970:1988,
                  time.optimize.ssr = 1970:1988)

args.TFDTW.synth = list(start.time = 1970, end.time = 2000, treat.time = 1989,
                        args.TFDTW = args.TFDTW, args.synth = args.synth,
                        ## 2nd
                        n.mse = 10, 
                        ## other
                        plot.figures = TRUE,
                        plot.path = "./figures/",
                        legend.pos = c(0.3, 0.3))

args.TFDTW.synth.all.units = list(target = "California",
                                  # data = data, 
                                  args.TFDTW.synth = args.TFDTW.synth,
                                  detailed.output = TRUE,
                                  ## 2nd
                                  all.units.parallel = TRUE)

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

saveRDS(results, "./data/res_tobacco_optimal_1018.Rds")

# plot placebo figure
df = results[[1]][["results.TFDTW.synth"]] %>%
  map(
    ~{
      item = .
      data.frame(unit = item[["dependent"]],
                 time = 1970:2000,
                 value = item$res.synth.raw$value,
                 synth_origin = item$res.synth.raw$synthetic,
                 synth_new = item$res.synth.TFDTW$synthetic)
    }
  ) %>%
  do.call("rbind", .) %>%
  mutate(
    color = case_when(unit == "California" ~ "black",
                      TRUE ~ "grey 70"),
    gap_origin = value - synth_origin,
    gap_new = value - synth_new
  )

df = left_join(df, df.rescale, by = "unit")
df = df %>% 
  mutate(gap_origin = gap_origin/multiplier,
         gap_new = gap_new/multiplier)

mse = results[["1"]][["mse"]]

fig = df %>%
  filter(unit %in% (mse %>% filter(mse.preT.raw < 5*9 & unit != "California") %>% .[["unit"]])) %>%
  ggplot(aes(x = time, group = unit)) +
  geom_line(aes(y = gap_origin), col = "#adcbe3") +
  geom_line(aes(y = gap_new), col = "#fec8c1") +
  geom_line(aes(y = gap_origin), data = df %>% filter(unit == "California"), col = "#2ab7ca", size = 1) +
  geom_line(aes(y = gap_new), data = df %>% filter(unit == "California"), col = "#fe4a49", size = 1) +
  geom_vline(xintercept = 1990, linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  # coord_cartesian(ylim=c(-32, 32)) +
  xlab("Year") +
  ylab("True Data - Synthetic Control") +
  ggtitle(paste0("Percent=", round(results[[1]]$neg.ratio, 4), ", P.value=", round(results[[1]]$p.value, 4))) +
  theme_bw()

ggsave("./figures/placebo_tobacco_1018.pdf",
       fig, width = 6, height = 4,
       units = "in", limitsize = FALSE)




# # prepare data
# start_time = 1970
# end_time = 2000
# treat_time = 1989
# dtw1_time = 1989
# filter_width = 13
# k = 5
# dtw1_method = "fixed"
# n_mse = 10
# n_burn = 3
# n_IQR = 3
# dist_quantile = 0.95
# plot_figures = TRUE
# step.pattern1 = dtw::mori2006
# step.pattern2 = dtw::asymmetricP2
# predictors.origin = NULL
# special.predictors.origin = list(
#   list("value_raw", 1988, c("mean")),
#   list("value_raw", 1980, c("mean")),
#   list("value_raw", 1975, c("mean")),
#   list("beer", 1984:1988, c("mean")),
#   list("lnincome", 1980:1988, c("mean")),
#   list("age15to24", 1980:1988, c("mean")),
#   list("retprice", 1980:1988, c("mean"))
# )
# time.predictors.prior.origin = 1970:1988
# time.optimize.ssr.origin = 1970:1988
# predictors.new = NULL
# special.predictors.new = list(
#   list("value_warped", 1988, c("mean")),
#   list("value_warped", 1980, c("mean")),
#   list("value_warped", 1975, c("mean")),
#   list("beer", 1984:1988, c("mean")),
#   list("lnincome", 1980:1988, c("mean")),
#   list("age15to24", 1980:1988, c("mean")),
#   list("retprice", 1980:1988, c("mean"))
# )
# time.predictors.prior.new = 1970:1988
# time.optimize.ssr.new = 1970:1988
# legend_position = c(0.3, 0.3)
# ... = NULL
# normalize_method = "t"
# ma = 3
# ma_na = "original"
# n_q = 1
# n_r = 1
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
#                                              n_burn = n_burn,
#                                              k = k,
#                                              dtw1_method = dtw1_method,
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
# saveRDS(result, "./data/grid_search_v6/result_tobacco_89_fixed.Rds")
# 
# 
# # mse
# mse = result %>% 
#   lapply(., "[[", "mse") %>% 
#   do.call("rbind", .) %>% 
#   mutate(ratio = mse2_post/mse1_post,
#          log_ratio = log(ratio))
# mse = mse %>% filter(dependent != "California")
# # mse = mse %>% filter(mse1_pre < 5*3)
# length(which(mse$log_ratio < 0))/nrow(mse)
# boxplot(mse$log_ratio, outline = FALSE)
# abline(h = 0, lty = 5)
# 
# t.test(mse$log_ratio)
# 
# saveRDS(mse, "./data/grid_search_v6/mse_tobacco_89.Rds")
# 
# 
# # plot time series figure
# df = rbind(data.frame(unit = "California",
#                       time = 1970:2000,
#                       value = result[[3]]$synth_origin$value),
#            data.frame(unit = "Synthetic Control w/o TFDTW",
#                       time = 1970:2000,
#                       value = result[[3]]$synth_origin$synthetic),
#            data.frame(unit = "Synthetic Control w/ TFDTW",
#                       time = 1970:2000,
#                       value = result[[3]]$synth_new$synthetic))
# 
# fig = ggplot(df, aes(x = time, y = value, color = unit)) +
#   geom_line() + 
#   geom_vline(xintercept = 1988, linetype="dashed") +
#   theme_bw() +
#   scale_color_manual(values=c("grey40", "#fe4a49","#2ab7ca")) +
#   xlab("Time") +
#   ylab("Cigarette Sale Per Capita (in Packs)") + 
#   theme(legend.title=element_blank(),
#         legend.position = c(0.2, 0.2))
# 
# saveRDS(fig, "./data/grid_search_v6/fig_comp_method_tobacco_89.Rds")
# 
# # plot placebo figure
# df = result %>% 
#   map(
#     ~{
#       data.frame(unit = .[["mse"]][["dependent"]],
#                  time = 1970:2000,
#                  value = .[["synth_origin"]][["value"]],
#                  synth_origin = .[["synth_origin"]][["synthetic"]],
#                  synth_new = .[["synth_new"]][["synthetic"]])
#     }
#   ) %>% 
#   do.call("rbind", .) %>% 
#   mutate(
#     color = case_when(unit == "California" ~ "black",
#                       TRUE ~ "grey 70"),
#     gap_origin = value - synth_origin,
#     gap_new = value - synth_new
#   )
# 
# df %>% 
#   filter(unit %in% (mse %>% filter(mse1_pre < 2*3) %>% .[["dependent"]])) %>%
#   ggplot(aes(x = time, group = unit)) +
#   geom_line(aes(y = gap_origin), col = "#adcbe3") +
#   geom_line(aes(y = gap_new), col = "#fec8c1") +
#   geom_line(aes(y = gap_origin), data = df %>% filter(unit == "California"), col = "#2ab7ca", size = 1) +
#   geom_line(aes(y = gap_new), data = df %>% filter(unit == "California"), col = "#fe4a49", size = 1) +
#   geom_vline(xintercept = 1988, linetype="dashed") +
#   geom_hline(yintercept = 0, linetype="dashed") +
#   coord_cartesian(ylim=c(-32, 32)) +
#   xlab("year") +
#   ylab("gap in per-capita cigarette sales (in packs)") +
#   theme_bw()









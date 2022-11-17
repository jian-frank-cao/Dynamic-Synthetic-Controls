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
            value.max = max(value)) %>% 
  ungroup()

mean.diff = mean(df.rescale$value.max - df.rescale$value.min)

df.rescale = df.rescale %>% 
  mutate(
    multiplier = mean.diff/(value.max - value.min)
  )

data = left_join(data, df.rescale, by = "unit")
data = data %>% 
  mutate(
    value.bak = value_raw,
    value_raw = (value_raw - value.min)*multiplier,
    value = value_raw
  )


# ## Grid Search Tobacco ---------------------------------------------------------
# # parameters
# filter.width.range = (1:9)*2+3
# k.range = 4:9
# step.pattern.range = list(
#   # symmetricP0 = dtw::symmetricP0, # too bumpy
#   # symmetricP05 = dtw::symmetricP05,
#   symmetricP1 = dtw::symmetricP1,
#   symmetricP2 = dtw::symmetricP2,
#   # asymmetricP0 = dtw::asymmetricP0, # too bumpy
#   # asymmetricP05 = dtw::asymmetricP05,
#   asymmetricP1 = dtw::asymmetricP1,
#   asymmetricP2 = dtw::asymmetricP2,
#   typeIc = dtw::typeIc,
#   # typeIcs = dtw::typeIcs,
#   # typeIIc = dtw::typeIIc,  # jumps
#   # typeIIIc = dtw::typeIIIc, # jumps
#   # typeIVc = dtw::typeIVc,  # jumps
#   typeId = dtw::typeId,
#   # typeIds = dtw::typeIds,
#   # typeIId = dtw::typeIId, # jumps
#   mori2006 = dtw::mori2006
# )
# grid.search.parallel = TRUE
# 
# 
# args.TFDTW = list(buffer = 0, match.method = "fixed",
#                   dist.quant = 0.95, 
#                   window.type = "none",
#                   ## other
#                   norm.method = "t",
#                   step.pattern2 = dtw::asymmetricP2,
#                   n.burn = 3, n.IQR = 3,
#                   ma = 3, ma.na = "original",
#                   default.margin = 3,
#                   n.q = 1, n.r = 1)
# 
# args.synth = list(predictors = NULL,
#                   special.predictors = 
#                     expression(list(
#                       list(dep.var, 1988, c("mean")),
#                       list(dep.var, 1980, c("mean")),
#                       list(dep.var, 1975, c("mean")),
#                       list("beer", 1984:1988, c("mean")),
#                       list("lnincome", 1980:1988, c("mean")),
#                       list("age15to24", 1980:1988, c("mean")),
#                       list("retprice", 1980:1988, c("mean"))
#                     )),
#                   time.predictors.prior = 1970:1988,
#                   time.optimize.ssr = 1970:1988)
# 
# args.TFDTW.synth = list(start.time = 1970, end.time = 2000, treat.time = 1989,
#                         args.TFDTW = args.TFDTW, args.synth = args.synth,
#                         ## 2nd
#                         n.mse = 10, 
#                         ## other
#                         plot.figures = FALSE,
#                         plot.path = "./figures/",
#                         legend.pos = c(0.3, 0.3))
# 
# args.TFDTW.synth.all.units = list(target = "California",
#                                   # data = data, 
#                                   args.TFDTW.synth = args.TFDTW.synth,
#                                   ## 2nd
#                                   all.units.parallel = FALSE)
# 
# args.TFDTW.synth.all.units[["data"]] = data
# results = SimDesign::quiet(
#   grid.search(filter.width.range = filter.width.range,
#               k.range = k.range,
#               step.pattern.range = step.pattern.range,
#               args.TFDTW.synth.all.units = args.TFDTW.synth.all.units,
#               grid.search.parallel = grid.search.parallel)
# )
# 
# saveRDS(results, "./data/res_tobacco_1019.Rds")
# job.end = Sys.time()
# print(job.end - job.start)
# 
# ## Result ----------------------------------------------------------------------
# results = readRDS("./data/res_tobacco_1019.Rds")
# neg.ratio = lapply(results, "[[", "neg.ratio") %>% 
#   do.call("c", .)
# p.value = lapply(results, "[[", "p.value") %>% 
#   do.call("c", .)
# max.neg = which(neg.ratio == max(neg.ratio))
# min.p = which(p.value == min(p.value))
# 
# neg.ratio[max.neg]
# p.value[max.neg]
# 
# neg.ratio[min.p]
# p.value[min.p]
# 
# ## Optimal Run Tobacco ---------------------------------------------------------
# # parameters
# filter.width.range = 15
# k.range = 7
# step.pattern.range = list(
#   # symmetricP0 = dtw::symmetricP0, # too bumpy
#   # symmetricP05 = dtw::symmetricP05,
#   # symmetricP1 = dtw::symmetricP1#,
#   # symmetricP2 = dtw::symmetricP2,
#   # asymmetricP0 = dtw::asymmetricP0, # too bumpy
#   # asymmetricP05 = dtw::asymmetricP05,
#   # asymmetricP1 = dtw::asymmetricP1,
#   # asymmetricP2 = dtw::asymmetricP2,
#   typeIc = dtw::typeIc#,
#   # typeIcs = dtw::typeIcs,
#   # typeIIc = dtw::typeIIc,  # jumps
#   # typeIIIc = dtw::typeIIIc, # jumps
#   # typeIVc = dtw::typeIVc,  # jumps
#   # typeId = dtw::typeId,
#   # typeIds = dtw::typeIds,
#   # typeIId = dtw::typeIId, # jumps
#   # mori2006 = dtw::mori2006
# )
# grid.search.parallel = FALSE
# 
# 
# args.TFDTW = list(buffer = 0, match.method = "fixed",
#                   dist.quant = 0.95,
#                   window.type = "none",
#                   ## other
#                   norm.method = "t",
#                   step.pattern2 = dtw::asymmetricP2,
#                   n.burn = 3, n.IQR = 3,
#                   ma = 3, ma.na = "original",
#                   default.margin = 3,
#                   n.q = 1, n.r = 1)
# 
# args.synth = list(predictors = NULL,
#                   special.predictors =
#                     expression(list(
#                       list(dep.var, 1988, c("mean")),
#                       list(dep.var, 1980, c("mean")),
#                       list(dep.var, 1975, c("mean")),
#                       list("beer", 1984:1988, c("mean")),
#                       list("lnincome", 1980:1988, c("mean")),
#                       list("age15to24", 1980:1988, c("mean")),
#                       list("retprice", 1980:1988, c("mean"))
#                     )),
#                   time.predictors.prior = 1970:1988,
#                   time.optimize.ssr = 1970:1988)
# 
# args.TFDTW.synth = list(start.time = 1970, end.time = 2000, treat.time = 1989,
#                         args.TFDTW = args.TFDTW, args.synth = args.synth,
#                         ## 2nd
#                         n.mse = 10,
#                         ## other
#                         plot.figures = TRUE,
#                         plot.path = "./figures/",
#                         legend.pos = c(0.3, 0.3))
# 
# args.TFDTW.synth.all.units = list(target = "California",
#                                   # data = data,
#                                   args.TFDTW.synth = args.TFDTW.synth,
#                                   detailed.output = TRUE,
#                                   ## 2nd
#                                   all.units.parallel = TRUE)
# 
# args.TFDTW.synth.all.units[["data"]] = data
# results = SimDesign::quiet(
#   grid.search(filter.width.range = filter.width.range,
#               k.range = k.range,
#               step.pattern.range = step.pattern.range,
#               args.TFDTW.synth.all.units = args.TFDTW.synth.all.units,
#               grid.search.parallel = grid.search.parallel)
# )

# saveRDS(results, "./data/res_tobacco_optimal_1020.Rds")
results = readRDS("./data/res_tobacco_optimal_1020.Rds")

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

df.interval = df %>% 
  filter(unit != "California") %>% 
  group_by(time) %>% 
  summarise(mean_origin = mean(gap_origin, na.rm = T),
            sd_origin = sd(gap_origin, na.rm = T),
            mean_new = mean(gap_new, na.rm = T),
            sd_new = sd(gap_new, na.rm = T),
            # ci_origin_upper = mean_origin + qt(0.975, df = 39-1)*sd_origin,
            # ci_origin_lower = mean_origin - qt(0.975, df = 39-1)*sd_origin,
            # ci_new_upper = mean_new + qt(0.975, df = 39-1)*sd_new,
            # ci_new_lower = mean_new - qt(0.975, df = 39-1)*sd_new,
            ci_origin_upper = quantile(gap_origin, 0.975, na.rm = T),
            ci_origin_lower = quantile(gap_origin, 0.025, na.rm = T),
            ci_new_upper = quantile(gap_new, 0.975, na.rm = T),
            ci_new_lower = quantile(gap_new, 0.025, na.rm = T)) %>% 
  mutate(unit = "area")

df.interval[28:31, 8:9] = NA

color_original = "#2ab7ca"
color_new = "#fe4a49"
# color_original = "grey70"
# color_new = "grey30"

colors = c("Gap (Original)" = color_original,
           "Gap (TFDTW)" = color_new)

fills = c("95% CI (Original)" = color_original,
          "95% CI (TFDTW)" = color_new)

fig_tobacco = df %>%
  filter(unit %in% (mse %>% filter(unit != "California") %>% .[["unit"]])) %>%
  ggplot(aes(x = time, group = unit)) +
  annotate("rect", xmin = 1989, xmax = 1999, ymin = -60, ymax = 70, 
           alpha = .3) +
  geom_line(aes(y = gap_origin), col = color_original, alpha = 0.4) +
  geom_line(aes(y = gap_new), col = color_new, alpha = 0.4) +
  geom_ribbon(aes(ymin = ci_origin_lower, ymax = ci_origin_upper, fill="95% CI (Original)"),
              data = df.interval, alpha=0.5) +
  geom_ribbon(aes(ymin = ci_new_lower, ymax = ci_new_upper, fill="95% CI (TFDTW)"),
              data = df.interval, alpha=0.5) +
  geom_line(aes(y = gap_origin, color = "Gap (Original)"), data = df %>% filter(unit == "California"), size = 1) +
  geom_line(aes(y = gap_new, color = "Gap (TFDTW)"), data = df %>% filter(unit == "California"), size = 1) +
  scale_color_manual(name = NULL, values = colors) +
  scale_fill_manual(name = NULL, values = fills) +
  geom_vline(xintercept = 1989, linetype="dashed", col = "grey20") +
  geom_hline(yintercept = 0, linetype="dashed", col = "grey20") +
  annotate("text", x = 1994, y = 33,
           label = "P = 0.023", col = "grey20") +
  annotate("text", x = 1988, y = 25, angle = 90,
           label = "Treatment", col = "grey20") +
  coord_cartesian(xlim=c(1969, 2009), ylim=c(-40,40)) +
  xlab("Year") +
  ylab("y - Synthetic Control") +
  theme_bw() + 
  theme(legend.position = "none", 
        legend.box = "horizontal",
        legend.background = element_rect(fill=NA))

# ggsave("./figures/placebo_tobacco_1103.pdf",
#        fig_tobacco, width = 8, height = 6,
#        units = "in", limitsize = FALSE)





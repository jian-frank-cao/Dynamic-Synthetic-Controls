args = commandArgs(trailingOnly=TRUE)
index = as.integer(args[1])
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
  filter(time <= 1970) %>%
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

# data list
data.list = list(list(target = "Basque Country (Pais Vasco)",
                 data = data))
ids = data$id %>%
  unique
for (i in ids) {
  data.temp = data %>% filter(!(id %in% c(i, 17)))
  data.list = c(data.list, 
                list(list(target = data.temp$unit[1],
                     data = data.temp)))
}


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

args.TFDTW.synth = list(start.time = 1955, end.time = 1997, treat.time = 1970,
                        args.TFDTW = args.TFDTW, args.synth = args.synth,
                        ## 2nd
                        n.mse = 10,
                        ## other
                        plot.figures = FALSE,
                        plot.path = "./figures/",
                        legend.pos = c(0.3, 0.3))

args.TFDTW.synth.all.units = list(target = data.list[[index]]$target,
                                  # data = data,
                                  args.TFDTW.synth = args.TFDTW.synth,
                                  ## 2nd
                                  detailed.output = TRUE,
                                  all.units.parallel = FALSE)

args.TFDTW.synth.all.units[["data"]] = data.list[[index]]$data
results = SimDesign::quiet(
  grid.search(filter.width.range = filter.width.range,
              k.range = k.range,
              step.pattern.range = step.pattern.range,
              args.TFDTW.synth.all.units = args.TFDTW.synth.all.units,
              grid.search.parallel = grid.search.parallel)
)

saveRDS(results, paste0("./data/res_basque_1204_", index, ".Rds"))
job.end = Sys.time()
print(job.end - job.start)

# ## Result ----------------------------------------------------------------------
# results = readRDS("./data/res_basque_1019.Rds")
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
# 
# ## Optimal Run -----------------------------------------------------------------
# # parameters
# filter.width.range = 21
# k.range = 5
# step.pattern.range = list(
#   # symmetricP0 = dtw::symmetricP0, # too bumpy
#   # symmetricP05 = dtw::symmetricP05,
#   # symmetricP1 = dtw::symmetricP1,
#   symmetricP2 = dtw::symmetricP2#,
#   # asymmetricP0 = dtw::asymmetricP0, # too bumpy
#   # asymmetricP05 = dtw::asymmetricP05,
#   # asymmetricP1 = dtw::asymmetricP1,
#   # asymmetricP2 = dtw::asymmetricP2,
#   # typeIc = dtw::typeIc,
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
#                       list(dep.var, 1960:1969, c("mean")),
#                       list("invest_ratio", 1964:1969, c("mean")),
#                       list("popdens", 1969, c("mean")),
#                       list("sec.agriculture", 1961:1969, c("mean")),
#                       list("sec.energy", 1961:1969, c("mean")),
#                       list("sec.industry", 1961:1969, c("mean")),
#                       list("sec.construction", 1961:1969, c("mean")),
#                       list("sec.services.venta", 1961:1969, c("mean")),
#                       list("sec.services.nonventa", 1961:1969, c("mean")),
#                       list("school.illit", 1964:1969, c("mean")),
#                       list("school.prim", 1964:1969, c("mean")),
#                       list("school.med", 1964:1969, c("mean")),
#                       list("school.high", 1964:1969, c("mean")),
#                       list("school.post.high", 1964:1969, c("mean"))
#                     )),
#                   time.predictors.prior = 1955:1969,
#                   time.optimize.ssr = 1955:1969)
# 
# args.TFDTW.synth = list(start.time = 1955, end.time = 1997, treat.time = 1970,
#                         args.TFDTW = args.TFDTW, args.synth = args.synth,
#                         ## 2nd
#                         n.mse = 10,
#                         ## other
#                         plot.figures = TRUE,
#                         plot.path = "./figures/",
#                         legend.pos = c(0.3, 0.3))
# 
# args.TFDTW.synth.all.units = list(target = "Basque Country (Pais Vasco)",
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
# 
# # saveRDS(results, "./data/res_basque_optimal_1020.Rds")
# results = readRDS("./data/res_basque_optimal_1020.Rds")
# 
# # plot placebo figure
# df = results[[1]][["results.TFDTW.synth"]] %>%
#   map(
#     ~{
#       item = .
#       data.frame(unit = item[["dependent"]],
#                  time = 1955:1997,
#                  value = item$res.synth.raw$value,
#                  synth_origin = item$res.synth.raw$synthetic,
#                  synth_new = item$res.synth.TFDTW$synthetic)
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
# df = left_join(df, df.rescale, by = "unit")
# df = df %>% 
#   mutate(gap_origin = gap_origin/multiplier,
#          gap_new = gap_new/multiplier)
# 
# mse = results[["1"]][["mse"]]
# 
# df.interval = df %>% 
#   filter(unit != "Basque Country (Pais Vasco)") %>% 
#   group_by(time) %>% 
#   summarise(mean_origin = mean(gap_origin, na.rm = T),
#             sd_origin = sd(gap_origin, na.rm = T),
#             mean_new = mean(gap_new, na.rm = T),
#             sd_new = sd(gap_new, na.rm = T),
#             # ci_origin_upper = mean_origin + qt(0.975, df = 18-1)*sd_origin,
#             # ci_origin_lower = mean_origin - qt(0.975, df = 18-1)*sd_origin,
#             # ci_new_upper = mean_new + qt(0.975, df = 18-1)*sd_new,
#             # ci_new_lower = mean_new - qt(0.975, df = 18-1)*sd_new,
#             ci_origin_upper = quantile(gap_origin, 0.975, na.rm = T),
#             ci_origin_lower = quantile(gap_origin, 0.025, na.rm = T),
#             ci_new_upper = quantile(gap_new, 0.975, na.rm = T),
#             ci_new_lower = quantile(gap_new, 0.025, na.rm = T)) %>% 
#   mutate(unit = "area")
# 
# color_original = "#2ab7ca"
# color_new = "#fe4a49"
# # color_original = "grey70"
# # color_new = "grey30"
# 
# colors = c("Gap (Original)" = color_original,
#            "Gap (TFDTW)" = color_new)
# 
# fills = c("95% CI (Original)" = color_original,
#           "95% CI (TFDTW)" = color_new)
# 
# fig_basque = df %>%
#   filter(unit %in% (mse %>% filter(unit != "Basque Country (Pais Vasco)") %>% .[["unit"]])) %>%
#   ggplot(aes(x = time, group = unit)) +
#   annotate("rect", xmin = 1970, xmax = 1980, ymin = -3, ymax = 3, 
#            alpha = .3) +
#   geom_line(aes(y = gap_origin), col = color_original, alpha = 0.4) +
#   geom_line(aes(y = gap_new), col = color_new, alpha = 0.4) +
#   geom_ribbon(aes(ymin = ci_origin_lower, ymax = ci_origin_upper, fill="95% CI (Original)"),
#               data = df.interval, alpha=0.5) +
#   geom_ribbon(aes(ymin = ci_new_lower, ymax = ci_new_upper, fill="95% CI (TFDTW)"),
#               data = df.interval, alpha=0.5) +
#   geom_line(aes(y = gap_origin, color = "Gap (Original)"), data = df %>% filter(unit == "Basque Country (Pais Vasco)"), size = 1) +
#   geom_line(aes(y = gap_new, color = "Gap (TFDTW)"), data = df %>% filter(unit == "Basque Country (Pais Vasco)"), size = 1) +
#   scale_color_manual(name = NULL, values = colors) +
#   scale_fill_manual(name = NULL, values = fills) +
#   geom_vline(xintercept = 1970, linetype="dashed", col = "grey20") +
#   geom_hline(yintercept = 0, linetype="dashed", col = "grey20") +
#   annotate("text", x = 1975, y = 0.85,
#            label = "P = 0.026", col = "grey20") +
#   annotate("text", x = 1969, y = 0.6, angle = 90,
#            label = "Treatment", col = "grey20") +
#   coord_cartesian(xlim=c(1950, 1990), ylim = c(-1, 1)) +
#   xlab("Year") +
#   ylab("y - Synthetic Control") +
#   theme_bw() + 
#   theme(legend.position = "none", 
#         legend.box = "horizontal",
#         legend.background = element_rect(fill=NA))
# 
# # ggsave("./figures/placebo_basque_1103.pdf",
# #        fig_basque, width = 8, height = 6,
# #        units = "in", limitsize = FALSE)
# 
# 
## Placebo v2 ------------------------------------------------------------------
results = NULL
folder = "./data/placebo/basque/"
res.files = list.files(folder)
for (res.file in res.files) {
  results = c(results, list(readRDS(paste0(folder, res.file))))
}

df = future_map2(
  results,
  as.list(1:length(results)),
  ~{
    item = .x
    index = .y
    mse = future_map2(
      item,
      names(item),
      ~{
        item = .x
        id = .y
        item$mse %>% mutate(id = id)
      }
    ) %>%
      do.call("rbind", .)

    units = unique(mse$unit)

    df.gap.list = NULL
    for (i in 1:length(units)) {
      target = units[[i]]
      scores = mse %>%
        filter(unit != target) %>%
        group_by(id) %>%
        summarise(percent = mean(log.ratio < 0),
                  p.value = t.test(log.ratio)$p.value)
      max.percent = which(scores$percent == max(scores$percent))
      min.p = which(scores$p.value[max.percent] == min(scores$p.value[max.percent])[1])[1]
      opt.ind = as.numeric(scores$id[max.percent[min.p]])
      df.gap.list[[i]] = data.frame(
        unit = paste0("d", index, "-", target),
        time = 1955:1997,
        value = item[[opt.ind]]$results.TFDTW.synth[[target]]$res.synth.raw$value,
        gap_original = item[[opt.ind]]$results.TFDTW.synth[[target]]$gap.raw,
        gap_new = item[[opt.ind]]$results.TFDTW.synth[[target]]$gap.TFDTW
      )
      # df.gap.list[[i]] = item[[opt.ind]]$mse %>% filter(unit == target) #%>%
        #mutate(unit = paste0("d", index, "-", target))
    }
    df.gap.list %>% do.call("rbind", .)
  }
) %>%
  do.call("rbind", .)

# ICC::ICCest(unit, log.ratio, data = df, CI.type = "S")

t.interval = 1971:1980
df = df %>% filter(time %in% t.interval)
n.t = length(t.interval)
n.datasets = nrow(df)/length(t.interval)


# log(diff/diff) t test ==========================================
df = df %>% 
  mutate(log.ratio = log(abs(gap_new)/abs(gap_original)),
         id = factor(str_split(unit, "-", simplify = TRUE)[,1]),
         target = factor(str_split(unit, "-", simplify = TRUE)[,2]))

res.aov = aov(log.ratio ~ id*target, df)
summary.aov = summary(res.aov)

BMS = summary.aov[[1]]$`Mean Sq`[1]
JMS = summary.aov[[1]]$`Mean Sq`[2]
IMS = summary.aov[[1]]$`Mean Sq`[3]
EMS = summary.aov[[1]]$`Mean Sq`[4]
k = 17
n = 18
l = 10
target = (BMS-IMS)/(l*k) + (JMS-IMS)/(l*n) + (IMS-EMS)/(l)
res.icc = target/(target + EMS)
res.vif = 1 + (l - 1)*res.icc
DF = nrow(df)/res.vif

t.value = t.test(df$log.ratio)$statistic
p.value = pt(t.value, df = DF, lower.tail = TRUE)*2
# ===============================================================


df_original = reshape2::dcast(df[c("unit", "time", "gap_original")],
                      time ~ unit, value.var = "gap_original")
value.icc.sc = irr::icc(
  df_original[,-1], model = "twoway",
  type = "agreement", unit = "single"
)$value
vif.sc = (nrow(df_original) - 1)*value.icc.sc + 1
DF.sc = (dim(df_original)[1]*dim(df_original)[2])/vif.sc

df_new = reshape2::dcast(df[c("unit", "time", "gap_new")],
                              time ~ unit, value.var = "gap_new")
value.icc.dsc = irr::icc(
  df_new[,-1], model = "twoway",
  type = "agreement", unit = "single"
)$value
vif.dsc = (nrow(df_new) - 1)*value.icc.dsc + 1
DF.dsc = (dim(df_new)[1]*dim(df_new)[2])/vif.dsc



var.sc = df %>% group_by(unit) %>%
  summarise(variance = var(gap_original)*(n.t - 1)) %>%
  ungroup %>%
  .[["variance"]] %>%
  sum(., na.rm = T)/(n.datasets*(n.t - 1))

var.dsc = df %>% group_by(unit) %>%
  summarise(variance = var(gap_new)*(n.t - 1)) %>%
  ungroup %>%
  .[["variance"]] %>%
  sum(., na.rm = T)/(n.datasets*(n.t - 1))

f.value = var.dsc/var.sc
f.value = round(f.value, 4)
p.value = pf(f.value, DF.dsc,
             DF.sc, lower.tail = TRUE)
p.value = round(p.value, 4)

# results = readRDS("./data/res_basque_1019.Rds")
# 
# mse = future_map2(
#   results,
#   names(results),
#   ~{
#     item = .x
#     id = .y
#     item$mse %>% mutate(id = id)
#   }
# ) %>%
#   do.call("rbind", .) %>%
#   filter(unit != "Basque Country (Pais Vasco)")
# 
# units = unique(mse$unit)
# opt.grid = data.frame(unit = units, id = 0)
# for (i in 1:nrow(opt.grid)) {
#   target = opt.grid$unit[i]
#   scores = mse %>%
#     filter(unit != target) %>%
#     group_by(id) %>%
#     summarise(percent = mean(log.ratio < 0),
#               p.value = t.test(log.ratio)$p.value)
#   max.percent = which(scores$percent == max(scores$percent))
#   min.p = which(scores$p.value[max.percent] == min(scores$p.value[max.percent])[1])[1]
#   opt.grid[i,2] = as.numeric(scores$id[max.percent[min.p]])
# }
# 
# tasks = unique(opt.grid$id)
# 
# 
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
# search.grid = expand.grid(filter.width.range, k.range,
#                           names(step.pattern.range)) %>%
#   `colnames<-`(c("filter.width", "k", "step.pattern"))
# tasks = cbind(data.frame(id = tasks),
#               search.grid[tasks,])
# grid.search.parallel = FALSE
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
#                       list(dep.var, 1960:1969, c("mean")),
#                       list("invest_ratio", 1964:1969, c("mean")),
#                       list("popdens", 1969, c("mean")),
#                       list("sec.agriculture", 1961:1969, c("mean")),
#                       list("sec.energy", 1961:1969, c("mean")),
#                       list("sec.industry", 1961:1969, c("mean")),
#                       list("sec.construction", 1961:1969, c("mean")),
#                       list("sec.services.venta", 1961:1969, c("mean")),
#                       list("sec.services.nonventa", 1961:1969, c("mean")),
#                       list("school.illit", 1964:1969, c("mean")),
#                       list("school.prim", 1964:1969, c("mean")),
#                       list("school.med", 1964:1969, c("mean")),
#                       list("school.high", 1964:1969, c("mean")),
#                       list("school.post.high", 1964:1969, c("mean"))
#                     )),
#                   time.predictors.prior = 1955:1969,
#                   time.optimize.ssr = 1955:1969)
# 
# args.TFDTW.synth = list(start.time = 1955, end.time = 1997, treat.time = 1970,
#                         args.TFDTW = args.TFDTW, args.synth = args.synth,
#                         ## 2nd
#                         n.mse = 10,
#                         ## other
#                         plot.figures = TRUE,
#                         plot.path = "./figures/",
#                         legend.pos = c(0.3, 0.3))
# 
# args.TFDTW.synth.all.units = list(target = "Basque Country (Pais Vasco)",
#                                   # data = data,
#                                   args.TFDTW.synth = args.TFDTW.synth,
#                                   detailed.output = TRUE,
#                                   ## 2nd
#                                   all.units.parallel = TRUE)
# 
# 
# results = tasks %>%
#   split(., seq(nrow(tasks))) %>%
#   set_names(tasks$id) %>%
#   map(
#     ~{
#       search = .
# 
#       args.TFDTW.synth.all.units[["data"]] = data
#       results = SimDesign::quiet(
#         grid.search(filter.width.range = search$filter.width,
#                     k.range = search$k,
#                     step.pattern.range = step.pattern.range[search$step.pattern],
#                     args.TFDTW.synth.all.units = args.TFDTW.synth.all.units,
#                     grid.search.parallel = grid.search.parallel)
#       )
#       results
#     }
#   )
# 
# saveRDS(results, "./data/placebo_v2_basque.Rds")
# results = readRDS("./data/placebo_v2_basque.Rds")
# 
# # plot figures
# df = future_map2(
#   results[as.character(opt.grid$id)],
#   units %>% as.list,
#   ~{
#     item = .x[[1]]
#     unit = .y
#     value = item$results.TFDTW.synth[[unit]]$res.synth.raw$value
#     gap_original = item$results.TFDTW.synth[[unit]]$gap.raw
#     gap_new = item$results.TFDTW.synth[[unit]]$gap.TFDTW
#     data.frame(unit = unit,
#                time = 1955:1997,
#                value = value,
#                gap_original = gap_original,
#                gap_new = gap_new)
#   }
# ) %>%
#   do.call("rbind", .)
# 
# t.interval = 1971:1980
# df = df %>% filter(time %in% t.interval)
# n.t = length(t.interval)
# n.datasets = nrow(df2)/length(t.interval)
# 
# 
# df_original = reshape2::dcast(df[c("unit", "time", "gap_original")],
#                       time ~ unit, value.var = "gap_original")
# value.icc.sc = irr::icc(
#   df_original[,-1], model = "twoway",
#   type = "agreement", unit = "single"
# )$value
# vif.sc = (nrow(df_original) - 1)*value.icc.sc + 1
# DF.sc = (dim(df_original)[1]*dim(df_original)[2])/vif.sc
# 
# df_new = reshape2::dcast(df[c("unit", "time", "gap_new")],
#                               time ~ unit, value.var = "gap_new")
# value.icc.dsc = irr::icc(
#   df_new[,-1], model = "twoway",
#   type = "agreement", unit = "single"
# )$value
# vif.dsc = (nrow(df_new) - 1)*value.icc.dsc + 1
# DF.dsc = (dim(df_new)[1]*dim(df_new)[2])/vif.dsc
# 
# 
# var.original = df2 %>% group_by(unit) %>%
#   summarise(variance = var(gap_original)*(n.t - 1)) %>%
#   ungroup %>%
#   .[["variance"]] %>%
#   sum(., na.rm = T)/(n.datasets*(n.t - 1))
# 
# var.new = df2 %>% group_by(unit) %>%
#   summarise(variance = var(gap_new)*(n.t - 1)) %>%
#   ungroup %>%
#   .[["variance"]] %>%
#   sum(., na.rm = T)/(n.datasets*(n.t - 1))
# 
# f.value = var.new/var.original
# f.value = round(f.value, 4)
# p.value = pf(f.value, DF.dsc,
#              DF.sc, lower.tail = TRUE)
# p.value = round(p.value, 4)




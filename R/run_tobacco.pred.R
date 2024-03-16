# args = commandArgs(trailingOnly=TRUE)
# task_id = as.integer(args[1])
# job.start = Sys.time()

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
source("./R/utility/grid.search.R")
set.seed(20220407)


## Function --------------------------------------------------------------------
impute_values <- function(column) {
  # Use last observation carried forward (LOCF) imputation
  column <- zoo::na.locf(column, na.rm = FALSE)
  
  # Use next observation carried backward (NOCB) for leading NAs
  if (is.na(column[1])) {
    column <- zoo::na.locf(column, fromLast = TRUE)
  }
  return(column)
}

TFDTW.synth = function(data, start.time, end.time, treat.time,
                       dependent, dependent.id,
                       args.TFDTW, args.synth,
                       pred.vars = NULL,
                       res.synth.raw = NULL,
                       n.mse = 10, plot.figures = FALSE,
                       plot.path = "./figures/",
                       legend.pos = c(0.3, 0.3)){
  # prepare data for TFDTW
  t.treat = (treat.time - start.time) + 1
  
  y.raw = data %>% 
    filter(unit == dependent & time <= end.time) %>%
    .[["value_raw"]]
  
  y.processed = data %>% 
    filter(unit == dependent & time <= end.time) %>%
    .[["value"]]
  
  x.list = data %>% 
    filter(unit != dependent & time <= end.time) %>% 
    # select(c("unit", "time", "value", "value_raw")) %>%
    group_by(unit) %>% 
    group_split(.keep = TRUE)
  
  args.TFDTW[["y"]] = y.processed
  args.TFDTW[["t.treat"]] = t.treat
  args.TFDTW[["plot.figures"]] = plot.figures
  
  # TFDTW
  results = NULL
  for (item in x.list) {
    unit = as.character(item$unit[1])
    x.raw = item$value_raw
    args.TFDTW[["x"]] = item$value
    
    res = do.call(TFDTW, args.TFDTW)
    
    x.warped = c(
      warpWITHweight(x.raw[1:res$cutoff], res$weight.a)[1:t.treat],
      warpWITHweight(x.raw[res$cutoff:length(x.raw)], res$avg.weight)[-1]
    )
    
    item$value_warped = x.warped[1:length(y.raw)]
    
    if (!is.null(pred.vars)) {
      for (pred.var in pred.vars) {
        pred.target = item[[pred.var]]
        pred.target.warped = c(
          warpWITHweight(pred.target[1:res$cutoff], res$weight.a)[1:t.treat],
          warpWITHweight(pred.target[res$cutoff:length(pred.target)],
                         res$avg.weight)[-1]
        )
        item[[pred.var]] = pred.target.warped[1:length(y.raw)]
      }
    }
    
    df.synth = item %>% 
      mutate(time = 1:length(y.raw) + start.time - 1,
             unit = unit)
    
    fig.warp = NULL
    if (plot.figures) {
      fig.warp = plot.warped(unit = unit, dependent = dependent,
                             t.treat = t.treat, y.raw = y.raw,
                             x.raw = x.raw, x.warped = x.warped,
                             legend.pos = legend.pos)
    }
    
    results[[unit]] = list(unit = unit, 
                           x.warped = x.warped,
                           df.synth = df.synth,
                           fig.warp = fig.warp)
  }
  
  # plot warped data
  if (plot.figures) {
    file.name = paste0(plot.path, "warped-", dependent, ".pdf")
    stack.warped(fig.list = lapply(results,"[[","fig.warp"),
                 file.name = file.name, ncol = 3)
  }
  
  # prepare data for Synth
  df.synth = lapply(results, "[[", "df.synth") %>% 
    do.call("rbind", .) %>% 
    `row.names<-`(NULL)
  
  df.dependent = data %>% 
    filter(unit == dependent & time <= end.time) %>% 
    mutate(time = 1:length(y.raw) + start.time - 1,
           unit = dependent,
           value_warped = y.raw)
  df.synth = rbind(df.synth, df.dependent)
  df.synth = data.frame(df.synth)
  
  args.synth[["df"]] = df.synth
  args.synth[["dependent.id"]] = dependent.id
  
  # Synth
  if (is.null(res.synth.raw)) {
    args.synth[["dep.var"]] = "value_raw"
    res.synth.raw = do.call(do.synth, args.synth)
  }
  
  args.synth[["dep.var"]] = "value_warped"
  res.synth.TFDTW = do.call(do.synth, args.synth)
  
  # plot synthetic control
  if (plot.figures) {
    file.name = paste0(plot.path, "synth_", dependent, ".pdf")
    plot.synth(res.synth.raw = res.synth.raw,
               res.synth.TFDTW = res.synth.TFDTW,
               dependent = dependent, start.time = start.time,
               end.time = end.time, treat.time = treat.time,
               file.name = file.name, legend.pos = legend.pos)
  }
  
  # output
  gap.raw = res.synth.raw$value - res.synth.raw$synthetic
  gap.TFDTW = res.synth.TFDTW$value - res.synth.TFDTW$synthetic
  
  mse = data.frame(
    unit = dependent,
    mse.preT.raw = mean(gap.raw[1:(t.treat-1)]^2, na.rm = T),
    mse.preT.TFDTW = mean(gap.TFDTW[1:(t.treat-1)]^2, na.rm = T),
    mse.postT.raw = mean(gap.raw[(t.treat+1):(t.treat+n.mse)]^2, na.rm = T),
    mse.postT.TFDTW = mean(gap.TFDTW[(t.treat+1):(t.treat+n.mse)]^2, na.rm = T)
  )
  
  return(list(dependent = dependent, dependent.id = dependent.id,
              # args.TFDTW = args.TFDTW, args.synth = args.synth,
              # results.TFDTW = results, df.synth = df.synth,
              res.synth.raw = res.synth.raw,
              res.synth.TFDTW = res.synth.TFDTW,
              gap.raw = gap.raw, gap.TFDTW = gap.TFDTW,
              mse = mse))
}


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
pred.vars = c("lnincome", "beer", "age15to24", "retprice")
data = data %>% 
  group_by(unit) %>% 
  mutate(across(pred.vars, ~impute_values(.))) %>% 
  ungroup() %>% 
  data.frame

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

# data list
data.list = list(list(target = "California",
                      data = data))
ids = data$id %>%
  unique
for (i in ids) {
  data.temp = data %>% filter(!(id %in% c(i, 3)))
  data.list = c(data.list, 
                list(list(target = data.temp$unit[1],
                          data = data.temp)))
}

select2 = combn(setdiff(ids, 3), 2, simplify = TRUE)[,1:100]

for (i in 1:ncol(select2)) {
  data.temp = data %>% filter(!(id %in% c(select2[, i], 3)))
  data.list = c(data.list, 
                list(list(target = data.temp$unit[1],
                          data = data.temp)))
}


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
                        n.mse = 10, pred.vars = pred.vars,
                        ## other
                        plot.figures = FALSE,
                        plot.path = "./figures/",
                        legend.pos = c(0.3, 0.3))

# args.TFDTW.synth.all.units = list(target = data.list[[index]]$target,
#                                   # data = data,
#                                   args.TFDTW.synth = args.TFDTW.synth,
#                                   ## 2nd
#                                   detailed.output = TRUE,
#                                   all.units.parallel = FALSE)
# 
# for (index in 1:length(data.list)) {
#   args.TFDTW.synth.all.units[["data"]] = data.list[[index]]$data
#   results = SimDesign::quiet(
#     grid.search(filter.width.range = filter.width.range,
#                 k.range = k.range,
#                 step.pattern.range = step.pattern.range,
#                 args.TFDTW.synth.all.units = args.TFDTW.synth.all.units,
#                 grid.search.parallel = grid.search.parallel)
#   )
#   
#   saveRDS(results, paste0("./data/placebo/tobacco/res_tobacco_", index, ".Rds"))
#   print(index)
# }

original_vector = setdiff(1:100, c((1:9)*10-9, (1:9)*10-8))

# Create an index to group every four elements
group_index <- ceiling(seq_along(original_vector) / 4)

# Split the original vector into a list of vectors, each containing 4 elements
task_list <- split(original_vector, group_index)


for (index in task_list[[task_id]]) {
  args.TFDTW.synth.all.units = list(target = data.list[[index]]$target,
                                    # data = data,
                                    args.TFDTW.synth = args.TFDTW.synth,
                                    ## 2nd
                                    detailed.output = TRUE,
                                    all.units.parallel = FALSE)
  args.TFDTW.synth.all.units[["data"]] = data.list[[index]]$data
  cat(paste0("Tobacco data set ", index, "..."))
  set.seed(20220407)
  results = SimDesign::quiet(
    grid.search(filter.width.range = filter.width.range,
                k.range = k.range,
                step.pattern.range = step.pattern.range,
                args.TFDTW.synth.all.units = args.TFDTW.synth.all.units,
                grid.search.parallel = grid.search.parallel)
  )
  
  saveRDS(results, paste0("./data/pred/tobacco/res_tobacco_", index, ".Rds"))
  cat("Done.\n")
}

job.end = Sys.time()
print(job.end - job.start)


## Placebo ---------------------------------------------------------------------
folder = "./data/pred/tobacco/"
file.list = as.list(list.files(folder))
get_data_id = function(file_name){
  return(strsplit(strsplit(file_name, "_")[[1]][3], "\\.")[[1]][1])
}
all.ids = sapply(list.files(folder), get_data_id)

pre.start = 11
pre.end = 20
post.start = 21
post.end = 28

# df.mse
df.mse = future_map2(
  file.list,
  as.list(all.ids),
  ~{
    file.name = .x
    data.id = .y
    data.list = readRDS(paste0(folder, file.name))
    mse = future_map2(
      data.list,
      as.list(names(data.list)),
      ~{
        result.synth = .x[["results.TFDTW.synth"]]
        grid.id = .y
        mse = result.synth %>% 
          map(
            ~{
              task = .
              unit = task$dependent
              ###
              scales = df.rescale %>% filter(unit == task$dependent)
              # value.min = scales$value.min
              multiplier = scales$multiplier
              gap.raw = task$gap.raw/multiplier
              gap.TFDTW = task$gap.TFDTW/multiplier
              ###
              data.frame(unit = unit,
                         mse.preT.raw = mean(gap.raw[pre.start:pre.end]^2, na.rm = T),
                         mse.preT.TFDTW = mean(gap.TFDTW[pre.start:pre.end]^2, na.rm = T),
                         mse.postT.raw = mean(gap.raw[post.start:post.end]^2, na.rm = T),
                         mse.postT.TFDTW = mean(gap.TFDTW[post.start:post.end]^2, na.rm = T))
            }
          ) %>% do.call("rbind", .)
        mse %>% mutate(grid.id = grid.id)
      }
    ) %>% do.call("rbind", .)
    mse$data.id = data.id
    
    mse %>% 
      group_by(unit) %>% 
      top_n(-1, mse.preT.TFDTW) %>% 
      top_n(-1, grid.id)
  }
) %>% do.call("rbind", .)

saveRDS(df.mse, "./data/df.mse_tobacco_pred2.Rds")

# t.test for log(MSEdsc/MSEsc)
df.mse = readRDS("./data/df.mse_tobacco_pred2.Rds")
df.mse = df.mse %>% 
  mutate(log.ratio = log(mse.postT.TFDTW/mse.postT.raw),
         data.id = as.character(data.id))

res.aov = aov(log.ratio ~ data.id*unit, df.mse)
summary.aov = summary(res.aov)

BMS = summary.aov[[1]]$`Mean Sq`[1]
JMS = summary.aov[[1]]$`Mean Sq`[2]
EMS = summary.aov[[1]]$`Mean Sq`[3]
n = summary.aov[[1]]$`Df`[1] + 1
k = summary.aov[[1]]$`Df`[2] + 1
res.icc = (BMS - EMS)/(BMS + (k - 1)*EMS + k*(JMS - EMS)/n)
res.vif = 1 + (k - 1)*res.icc
DF = nrow(df.mse)/res.vif

t.value = t.test(df.mse$log.ratio)$statistic
p.value = pt(t.value, df = DF, lower.tail = TRUE)*2


## Plot results ----------------------------------------------------------------
# df.target
results.target = readRDS("./data/pred/tobacco/res_tobacco_1.Rds")
target = "California"
scales = df.rescale %>% filter(unit == target)
multiplier = scales$multiplier

pre.start = 11
pre.end = 20
post.start = 21
post.end = 28

mse = future_map2(
  results.target,
  as.list(names(results.target)),
  ~{
    result.synth = .x[["results.TFDTW.synth"]][[target]]
    grid.id = .y
    gap.raw = result.synth$gap.raw/multiplier
    gap.TFDTW = result.synth$gap.TFDTW/multiplier
    data.frame(grid.id = grid.id,
               mse.preT.raw = mean(gap.raw[pre.start:pre.end]^2, na.rm = T),
               mse.preT.TFDTW = mean(gap.TFDTW[pre.start:pre.end]^2, na.rm = T),
               mse.postT.raw = mean(gap.raw[post.start:post.end]^2, na.rm = T),
               mse.postT.TFDTW = mean(gap.TFDTW[post.start:post.end]^2, na.rm = T))
  }
) %>% do.call("rbind", .)

opt.grid.id = mse %>% 
  top_n(-1, mse.preT.TFDTW) %>% 
  top_n(-1, grid.id) %>% 
  .[["grid.id"]]

df.target = data.frame(
  time = 1970:2000,
  unit = target,
  data.id = 0,
  grid.id = opt.grid.id,
  value = results.target[[opt.grid.id]][[4]][[target]][[3]][["value"]],
  synth.sc = results.target[[opt.grid.id]][[4]][[target]][[3]][["synthetic"]],
  synth.dsc = results.target[[opt.grid.id]][[4]][[target]][[4]][["synthetic"]]
)

df.target = df.target %>% 
  mutate(
    gap.sc = (value - synth.sc)/multiplier,
    gap.dsc = (value - synth.dsc)/multiplier,
    group = paste0(data.id, "-", grid.id, "-", unit)
  )

saveRDS(df.target, "./data/df.target_tobacco_pred2.Rds")


# df.gap
df.mse = readRDS("./data/df.mse_tobacco_pred2.Rds")
folder = "./data/pred/tobacco/"
file.list = as.list(list.files(folder))
all.ids = sapply(list.files(folder), get_data_id)
results = file.list %>%
  set_names(all.ids) %>% 
  future_map(
    ~{
      file.name = .
      readRDS(paste0(folder, file.name))
    }
  )


df.gap = NULL
for (i in 1:nrow(df.mse)) {
  unit = df.mse$unit[i]
  data.id = df.mse$data.id[i]
  grid.id = df.mse$grid.id[i]
  scales = df.rescale %>% filter(unit == df.mse$unit[i])
  value.min = scales$value.min
  multiplier = scales$multiplier
  value = results[[data.id]][[grid.id]][[4]][[unit]][[3]][["value"]]
  synth.sc = results[[data.id]][[grid.id]][[4]][[unit]][[3]][["synthetic"]]
  synth.dsc = results[[data.id]][[grid.id]][[4]][[unit]][[4]][["synthetic"]]
  
  df.gap[[i]] = data.frame(
    time = 1970:2000,
    unit = unit,
    data.id = data.id,
    grid.id = grid.id,
    value = value/multiplier + value.min,
    synth.sc = synth.sc/multiplier + value.min,
    synth.dsc = synth.dsc/multiplier + value.min
  )
  print(i)
}

df.gap = df.gap %>%
  do.call("rbind", .) %>%
  mutate(
    gap.sc = value - synth.sc,
    gap.dsc = value - synth.dsc,
    group = paste0(data.id, "-", grid.id, "-", unit)
  )

saveRDS(df.gap, "./data/df.gap_tobacco_pred2.Rds")

avg.log.mse.sc = mean(log(df.mse$mse.postT.raw), na.rm = TRUE)
avg.log.mse.dsc = mean(log(df.mse$mse.postT.TFDTW), na.rm = TRUE)

# plot
df.target = readRDS("./data/df.target_tobacco_pred2.Rds")
df.gap = readRDS("./data/df.gap_tobacco_pred2.Rds")

df.quantile = df.gap %>%
  group_by(time) %>%
  summarise(quantile.sc.975 = quantile(gap.sc, 0.975, na.rm = T),
            quantile.sc.025 = quantile(gap.sc, 0.025, na.rm = T),
            quantile.dsc.975 = quantile(gap.dsc, 0.975, na.rm = T),
            quantile.dsc.025 = quantile(gap.dsc, 0.025, na.rm = T)) %>% 
  mutate(group = "quantile")

df.quantile[29:31,4:5] = NA

color.sc = "#2ab7ca"
color.dsc = "#fe4a49"
# color.sc = "grey70"
# color.dsc = "grey30"

colors = c("TE (SC)" = color.sc,
           "TE (DSC)" = color.dsc)

fills = c("95% Quantile (SC)" = color.sc,
          "95% Quantile (DSC)" = color.dsc)

set.seed(20230812)
group.sample = sample(unique(df.gap$group), 100)

fig_tobacco = df.gap %>%
  filter(group %in% group.sample) %>% 
  ggplot(aes(x = time, group = group)) +
  annotate("rect", xmin = 1989, xmax = 1999, 
           ymin = -60, ymax = 70, alpha = .3) +
  geom_line(aes(y = gap.sc), col = color.sc, alpha = 0.4) +
  geom_line(aes(y = gap.dsc), col = color.dsc, alpha = 0.4) +
  geom_ribbon(aes(ymin = quantile.sc.025,
                  ymax = quantile.sc.975,
                  fill="95% Quantile (SC)"),
              data = df.quantile, alpha=0.5) +
  geom_ribbon(aes(ymin = quantile.dsc.025,
                  ymax = quantile.dsc.975,
                  fill="95% Quantile (DSC)"),
              data = df.quantile, alpha=0.5) +
  geom_line(aes(y = gap.sc, color = "TE (SC)"),
            data = df.target, size = 1) +
  geom_line(aes(y = gap.dsc, color = "TE (DSC)"), 
            data = df.target, size = 1) +
  scale_color_manual(name = NULL, values = colors) +
  scale_fill_manual(name = NULL, values = fills) +
  geom_vline(xintercept = 1989, linetype="dashed", col = "grey20") +
  geom_hline(yintercept = 0, linetype="dashed", col = "grey20") +
  annotate("text", x = 1994, y = 33,
           label = "t = -9.093\nP < 0.0001", col = "grey20") +
  annotate("text", x = 1988, y = 25, angle = 90,
           label = "Treatment", col = "grey20") +
  annotate("text", x = 2005, y = -20, label = "bar(log(MSE))[SC]==4.08", parse = TRUE,
           col = "grey20", size = 4, fontface = "bold") +
  annotate("text", x = 2005, y = -30, label = "bar(log(MSE))[DSC]==3.92", parse = TRUE,
           col = "grey20", size = 4, fontface = "bold") +
  coord_cartesian(xlim=c(1969, 2009), ylim=c(-40,40)) +
  xlab("Year") +
  ylab("TE(y - Synthetic Control)") +
  theme_bw() +
  theme(legend.position = "none",
        legend.box = "horizontal",
        legend.background = element_rect(fill=NA))


saveRDS(fig_tobacco, "./data/placebo_tobacco_pred.Rds")

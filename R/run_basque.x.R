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
        item[[pred.var]] = pred.target[1:length(y.raw)]
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



## Basque Terrorism Data -------------------------------------------------------
data(basque, package = "Synth")
data = basque
colnames(data)[1:4] = c("id", "unit", "time", "value")
pred.vars = c("sec.agriculture", "sec.energy", "sec.industry",
              "sec.construction", "sec.services.venta",
              "sec.services.nonventa", "school.illit",
              "school.prim", "school.med", "school.high",
              "school.post.high", "popdens", "invest", "invest_ratio")
data = data %>% mutate(invest_ratio = invest/value,
                       value_raw = value) %>% 
  group_by(unit) %>% 
  mutate(across(pred.vars, ~impute_values(.))) %>% 
  ungroup()

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

select2 = combn(setdiff(ids, 17), 2, simplify = TRUE)[,1:100]

for (i in 1:ncol(select2)) {
  data.temp = data %>% filter(!(id %in% c(select2[, i], 17)))
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
                        n.mse = 10, pred.vars = pred.vars,
                        ## other
                        plot.figures = FALSE,
                        plot.path = "./figures/",
                        legend.pos = c(0.3, 0.3))

for (index in 1:length(data.list)) {
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
  
  saveRDS(results, paste0("./data/placebo/basque/res_basque_", index, ".Rds"))
  print(index)
}


# ## Placebo ---------------------------------------------------------------------
# folder = "./data/placebo/basque/"
# file.list = as.list(list.files(folder))
# 
# pre.start = 7
# pre.end = 16
# post.start = 17
# post.end = 26
# 
# # df.mse
# df.mse = file.list %>% 
#   future_map(
#   ~{
#     file.name = .
#     data.id = strsplit(file.name, "_")[[1]][4]
#     data.id = as.numeric(strsplit(data.id, ".R")[[1]][1])
#     data.list = readRDS(paste0(folder, file.name))
#     mse = future_map2(
#       data.list,
#       as.list(names(data.list)),
#       ~{
#         result.synth = .x[["results.TFDTW.synth"]]
#         grid.id = .y
#         mse = result.synth %>% 
#           map(
#             ~{
#               task = .
#               unit = task$dependent
#               gap.raw = task$gap.raw
#               gap.TFDTW = task$gap.TFDTW
#               data.frame(unit = unit,
#                          mse.preT.raw = mean(gap.raw[pre.start:pre.end]^2, na.rm = T),
#                          mse.preT.TFDTW = mean(gap.TFDTW[pre.start:pre.end]^2, na.rm = T),
#                          mse.postT.raw = mean(gap.raw[post.start:post.end]^2, na.rm = T),
#                          mse.postT.TFDTW = mean(gap.TFDTW[post.start:post.end]^2, na.rm = T))
#             }
#           ) %>% do.call("rbind", .)
#         mse %>% mutate(grid.id = grid.id)
#       }
#     ) %>% do.call("rbind", .)
#     mse$data.id = data.id
#     
#     mse %>% 
#       group_by(unit) %>% 
#       top_n(-1, mse.preT.TFDTW) %>% 
#       top_n(-1, grid.id)
#   }
# ) %>% do.call("rbind", .)
# 
# saveRDS(df.mse, "./data/df.mse_basque.x.Rds")
# 
# 
# ## Preparation -----------------------------------------------------------------
# df.mse = readRDS("./data/df.mse_basque.x.Rds")
# 
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
# 
# # grid search space
# search.grid = expand.grid(filter.width.range, k.range,
#                           names(step.pattern.range)) %>% 
#   `colnames<-`(c("filter.width", "k", "step.pattern"))
# search.grid.list = search.grid %>% split(., seq(nrow(search.grid)))


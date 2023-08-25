# args = commandArgs(trailingOnly=TRUE)
# treat_time  = as.integer(args[1])
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


## Data ------------------------------------------------------------------------
# symbols = c("TIP", "INXG.L", "XRB.TO", "IBCI.DE") #,"GTIP" ,"ILB.XA","XSTH.TO"
# 
# data_list = NULL
# 
# # Loop over the symbols and download the data
# for(symbol in symbols) {
#   response = quantmod::getSymbols(symbol, auto.assign = FALSE,
#                              from = "2010-02-01", to = "2022-12-31")
#   
#   data_list[[symbol]] = response %>% 
#     data.frame(.) %>% 
#     mutate(date = as.Date(rownames(.))) %>% 
#     `rownames<-`(NULL) %>% 
#     .[,c(7:6)] %>% 
#     `colnames<-`(c("date", symbol))
# }
# 
# data = Reduce(function(...) merge(..., by = "date", all = TRUE), data_list)
# cols_to_fill <- names(data)[-1]
# data = data %>% fill(all_of(cols_to_fill))
# 
# saveRDS(data, "./data/bond.Rds")

data = readRDS("./data/bond.Rds")
data = reshape2::melt(data, id.vars = 1)
data_daily = data %>% 
  mutate(date = as.Date(date),
         month = format(date, "%Y-%m"),
         unit = as.character(variable),
         id = as.numeric(variable)) %>% 
  group_by(id, unit, date) %>% 
  # summarise(value = mean(value, na.rm = T)) %>% 
  # summarise(value = first(value)) %>% 
  summarise(value = last(value)) %>%
  mutate(time = 1:3339,
         value_raw = value) %>% 
  ungroup() %>% 
  select(c("id", "unit", "time", "value", "value_raw", "date")) %>% 
  data.frame(.) %>% 
  filter(date > as.Date("2018-01-01")) %>% 
  mutate(time = rep(1:1292,4))

## Functions -------------------------------------------------------------------
# grid.search = function(filter.width.range, k.range, step.pattern.range,
#                        args.TFDTW.synth.all.units,
#                        grid.search.parallel = TRUE){
#   if (grid.search.parallel) {
#     fun.map = furrr::future_map
#   }else{
#     fun.map = purrr::map
#   }
#   
#   # vanilla synthetic control
#   target = args.TFDTW.synth.all.units[["target"]]
#   data = args.TFDTW.synth.all.units[["data"]]
#   units = data[c("id", "unit")] %>% distinct %>% filter(unit == target)
#   units.list = units %>% split(., seq(nrow(units)))
#   
#   args.synth = args.TFDTW.synth.all.units$args.TFDTW.synth$args.synth
#   args.synth[["df"]] = data
#   args.synth[["dep.var"]] = "value_raw"
#   
#   res.synth.raw.list = units.list %>% 
#     set_names(units$unit) %>% 
#     fun.map(
#       ~{
#         item = .
#         dependent.id = item$id
#         args.synth[["dependent.id"]] = dependent.id
#         res = do.call(do.synth, args.synth)
#       }
#     )
#   
#   # grid search space
#   search.grid = expand.grid(filter.width.range, k.range,
#                             names(step.pattern.range)) %>% 
#     `colnames<-`(c("filter.width", "k", "step.pattern"))
#   search.grid.list = search.grid %>% split(., seq(nrow(search.grid)))
#   
#   # search start
#   results = search.grid.list %>% 
#     fun.map(
#       ~{
#         task = .
#         args.TFDTW.synth.all.units[["filter.width"]] = task$filter.width
#         args.TFDTW.synth.all.units$args.TFDTW.synth$args.TFDTW[["k"]] = task$k
#         args.TFDTW.synth.all.units$args.TFDTW.synth$args.TFDTW[["step.pattern1"]] =
#           step.pattern.range[[task$step.pattern]]
#         args.TFDTW.synth.all.units[["res.synth.raw.list"]] = res.synth.raw.list
#         do.call(TFDTW.synth.target.only, args.TFDTW.synth.all.units)
#       }
#     )
#   
#   return(results)
# }






## Grid Search -----------------------------------------------------------------
target = "TIP"
treat_time = 516

# parameters
filter.width.range = 3*30+1
k.range = 3*30
step.pattern.range = list(
  # symmetricP0 = dtw::symmetricP0, # too bumpy
  # symmetricP05 = dtw::symmetricP05,
  # symmetricP1 = dtw::symmetricP1,
  symmetricP2 = dtw::symmetricP2,
  # asymmetricP0 = dtw::asymmetricP0, # too bumpy
  # asymmetricP05 = dtw::asymmetricP05,
  # asymmetricP1 = dtw::asymmetricP1,
  # asymmetricP2 = dtw::asymmetricP2,
  # typeIc = dtw::typeIc,
  # typeIcs = dtw::typeIcs,
  # typeIIc = dtw::typeIIc,  # jumps
  # typeIIIc = dtw::typeIIIc, # jumps
  # typeIVc = dtw::typeIVc,  # jumps
  # typeId = dtw::typeId,
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
                    list(
                      list("value_raw",(treat_time - 12*30):(treat_time - 1),c("mean"))
                    ),
                  time.predictors.prior = 1:(treat_time - 1),
                  time.optimize.ssr = 1:(treat_time - 1))

args.TFDTW.synth = list(start.time = 1, end.time = 1292, treat.time = treat_time,
                        args.TFDTW = args.TFDTW, args.synth = args.synth,
                        ## 2nd
                        n.mse = 24*20,
                        ## other
                        plot.figures = FALSE,
                        plot.path = "./figures/",
                        legend.pos = c(0.3, 0.3))

args.TFDTW.synth.all.units = list(target = target,
                                  # data = data,
                                  args.TFDTW.synth = args.TFDTW.synth,
                                  ## 2nd
                                  detailed.output = TRUE,
                                  all.units.parallel = FALSE)

args.TFDTW.synth.all.units[["data"]] = data_daily
results = SimDesign::quiet(
  grid.search(filter.width.range = filter.width.range,
              k.range = k.range,
              step.pattern.range = step.pattern.range,
              args.TFDTW.synth.all.units = args.TFDTW.synth.all.units,
              grid.search.parallel = grid.search.parallel)
)

saveRDS(results, paste0("./data/bond_daily_", treat_time, ".Rds"))
job.end = Sys.time()
print(job.end - job.start)

# plot
df = results[[2]]
target = "TIP"
treat_time = 516
rbind(data.frame(time = 1:1292,
                 value = df$results.TFDTW.synth[[target]]$res.synth.raw$value - 
                   df$results.TFDTW.synth[[target]]$res.synth.raw$synthetic,
                 unit = "sc"),
      data.frame(time = 1:1292,
                 value = df$results.TFDTW.synth[[target]]$res.synth.raw$value - 
                   df$results.TFDTW.synth[[target]]$res.synth.TFDTW$synthetic,
                 unit = "dsc")) %>% 
  ggplot(aes(x = time, y = value, color = unit)) +
  geom_line()  +
  ylab("Observed - Synthetic Control") +
  labs(color = "") +
  geom_vline(xintercept = treat_t)




## Results ---------------------------------------------------------------------
# folder = "./data/bond/"
# file.list = list.files(folder)
# 
# results = NULL
# for (file in file.list) {
#   print(file)
#   df = readRDS(paste0(folder, file))
#   treat_time = as.numeric(strsplit(strsplit(file, "[.]")[[1]][1], "_")[[1]][2])
#   mse = future_map2(
#     as.list(1:length(df)),
#     df,
#     ~{
#       id = .x
#       item = .y
#       item$mse$id = id
#       item$mse$diff = 0
#       item$mse$error20 = 0
#       for (i in 1:4) {
#         unit = item$mse$unit[i]
#         error = (item$results.TFDTW.synth[[unit]]$res.synth.raw$value -
#                    item$results.TFDTW.synth[[unit]]$res.synth.TFDTW$synthetic)
#         difference = (item$results.TFDTW.synth[[unit]]$res.synth.raw$synthetic -
#           item$results.TFDTW.synth[[unit]]$res.synth.TFDTW$synthetic)
#         msd = mean((difference[treat_time:(treat_time + 23)])^2, na.rm = TRUE)
#         item$mse$diff[i] = msd
#         error20 = mean((error[(treat_time-20):(treat_time-1)])^2, na.rm = TRUE)
#         item$mse$error20[i] = error20
#       }
#       item$mse
#     }
#   ) %>%
#     do.call("rbind", .)
#   mse$file = file
#   results[[file]] = mse
# }
# 
# results = results %>%
#   do.call("rbind", .) %>%
#   `rownames<-`(NULL)
# 
# saveRDS(results, "./data/bond_results.Rds")
results = readRDS("./data/bond_results.Rds")

target = "TIP"
mse = results %>% 
  filter(unit == target) %>% 
  group_by(file) %>% 
  filter(error20 == min(error20)) %>% 
  filter(id == min(id))

data_monthly %>% 
  ggplot(aes(x = time, y = value, color = unit)) +
  geom_line() +
  geom_vline(xintercept = 123) +
  scale_x_continuous(breaks = c(0, 48, 96, 144),
                     labels = c("2010", "2014", "2018", "2022"))

treat_t = 119
id = 145
#TIP 110 75, 119 4, 111 5, 
#TIP last 119 145
df = readRDS(paste0("./data/bond3/bond_last_", treat_t, ".Rds"))[[id]] 

rbind(data.frame(time = 1:155,
                 value = df$results.TFDTW.synth[[target]]$res.synth.raw$value - 
                   df$results.TFDTW.synth[[target]]$res.synth.raw$synthetic,
                 unit = "sc"),
      data.frame(time = 1:155,
                 value = df$results.TFDTW.synth[[target]]$res.synth.raw$value - 
                   df$results.TFDTW.synth[[target]]$res.synth.TFDTW$synthetic,
                 unit = "dsc")) %>% 
  ggplot(aes(x = time, y = value, color = unit)) +
  geom_line()  +
  ylab("Observed - Synthetic Control") +
  labs(color = "") +
  geom_vline(xintercept = treat_t) +
  scale_x_continuous(breaks = c(0, 48, 96, 144),
                     labels = c("2010", "2014", "2018", "2022"))

# rbind(data.frame(time = 1:155,
#                  value = df$results.TFDTW.synth$TIP$res.synth.raw$value,
#                  unit = "value"),
#       data.frame(time = 1:155,
#                  value = df$results.TFDTW.synth$TIP$res.synth.raw$synthetic,
#                  unit = "sc"),
#       data.frame(time = 1:155,
#                  value = df$results.TFDTW.synth$TIP$res.synth.TFDTW$synthetic,
#                  unit = "dsc")) %>% 
#   ggplot(aes(x = time, y = value, color = unit)) +
#   geom_line()  +
#   geom_vline(xintercept = 111)
# 
# rbind(data.frame(time = 1:155,
#                  value = df$res.synth.target.raw$value,
#                  unit = "value"),
#       data.frame(time = 1:155,
#                  value = df$res.synth.target.raw$synthetic,
#                  unit = "sc"),
#       data.frame(time = 1:155,
#                  value = df$res.synth.target.TFDTW$synthetic,
#                  unit = "dsc")) %>% 
#   ggplot(aes(x = time, y = value, color = unit)) +
#   geom_line()  +
#   geom_vline(xintercept = 29)

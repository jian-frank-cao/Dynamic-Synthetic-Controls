## Setup -----------------------------------------------------------------------
library(checkpoint)
checkpoint("2022-04-01")

library(tidyverse)
library(furrr)
plan(multisession, workers = 7)
options(future.rng.onMisuse="ignore")
options(stringsAsFactors = FALSE)

source("./R/misc.R")
source("./R/TFDTW.R")
source("./R/synth.R")
source("./R/compare.R")
set.seed(20220407)


## Basque Terrorism Data --------------------------------------------------
data(basque, package = "Synth")
data = basque
colnames(data)[1:4] = c("id", "unit", "time", "value")
data = data %>% mutate(invest_ratio = invest/value,
                       value_raw = value)
# data = data %>% filter(unit != "Basque Country (Pais Vasco)")


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
         age15to24 = age15to24*100) #%>% 
#filter(unit != "Rhode Island")

data = smoking
# data = data %>% filter(unit != "California")



## Germany Reunification Data --------------------------------------------------
data = foreign::read.dta("./data/repgermany.dta")
colnames(data)[1:4] = c("id", "unit", "time", "value")
data = data %>% mutate(value_raw = value)
# data = data %>% filter(unit != "West Germany")


## Grid Search Tobacco ---------------------------------------------------------
# search space
width_range = (1:9)*2+3
k_range = 4:9
step_pattern_range = list(
  # symmetricP0 = dtw::symmetricP0, # too bumpy
  symmetricP05 = dtw::symmetricP05,
  symmetricP1 = dtw::symmetricP1,
  symmetricP2 = dtw::symmetricP2,
  # asymmetricP0 = dtw::asymmetricP0, # too bumpy
  asymmetricP05 = dtw::asymmetricP05,
  asymmetricP1 = dtw::asymmetricP1,
  asymmetricP2 = dtw::asymmetricP2,
  typeIc = dtw::typeIc,
  typeIcs = dtw::typeIcs,
  # typeIIc = dtw::typeIIc,  # jumps
  # typeIIIc = dtw::typeIIIc, # jumps
  # typeIVc = dtw::typeIVc,  # jumps
  typeId = dtw::typeId,
  typeIds = dtw::typeIds,
  # typeIId = dtw::typeIId, # jumps
  mori2006 = dtw::mori2006
)

res_grid = expand.grid(width_range, k_range,
                       names(step_pattern_range)) %>% 
  `colnames<-`(c("width", "k", "step_pattern")) %>% 
  mutate(mse_pre_original = NA_real_,
         mse_pre_new = NA_real_,
         pos_ratio = NA_real_,
         t_test = NA_real_)

res_grid_filename = "./data/grid_search_v6/res_grid_tobacco_89.Rds"
# saveRDS(res_grid, res_grid_filename)
res_grid = readRDS(res_grid_filename)

# search
start_time = 1970
treat_time = 2000
end_time = 1989
dtw1_time = 1994
n_mse = 10
n_IQR = 3
dist_quantile = 0.95
plot_figures = FALSE
step.pattern2 = dtw::asymmetricP2
predictors.origin = NULL
special.predictors.origin = list(
  list("value_raw", 1988, c("mean")),
  list("value_raw", 1980, c("mean")),
  list("value_raw", 1975, c("mean")),
  list("beer", 1984:1988, c("mean")),
  list("lnincome", 1980:1988, c("mean")),
  list("age15to24", 1980:1988, c("mean")),
  list("retprice", 1980:1988, c("mean"))
)
time.predictors.prior.origin = 1970:1988
time.optimize.ssr.origin = 1970:1988
predictors.new = NULL
special.predictors.new = list(
  list("value_warped", 1988, c("mean")),
  list("value_warped", 1980, c("mean")),
  list("value_warped", 1975, c("mean")),
  list("beer", 1984:1988, c("mean")),
  list("lnincome", 1980:1988, c("mean")),
  list("age15to24", 1980:1988, c("mean")),
  list("retprice", 1980:1988, c("mean"))
)
time.predictors.prior.new = 1970:1988
time.optimize.ssr.new = 1970:1988
legend_position = c(0.8, 0.8)

for (i in which(is.na(res_grid$pos_ratio))) {
  width = res_grid$width[i]
  k = res_grid$k[i]
  pattern_name = res_grid$step_pattern[i]
  step.pattern1 = step_pattern_range[[pattern_name]]
  
  cat(paste0("[Task-", i, "]: width-", width, ", k-", k, 
             ", ", pattern_name, "......"))
  
  res = SimDesign::quiet(run_all_units(data = data,
                                       start_time = start_time,
                                       end_time = end_time,
                                       treat_time = treat_time,
                                       dtw1_time = dtw1_time,
                                       n_mse = n_mse,
                                       k = k,
                                       n_IQR = n_IQR,
                                       dist_quantile = dist_quantile,
                                       plot_figures = plot_figures,
                                       step.pattern1 = step.pattern1,
                                       step.pattern2 = step.pattern2,
                                       predictors.origin = predictors.origin,
                                       special.predictors.origin = special.predictors.origin,
                                       time.predictors.prior.origin = time.predictors.prior.origin,
                                       time.optimize.ssr.origin = time.optimize.ssr.origin,
                                       predictors.new = predictors.new,
                                       special.predictors.new = special.predictors.new,
                                       time.predictors.prior.new = time.predictors.prior.new,
                                       time.optimize.ssr.new = time.optimize.ssr.new,
                                       legend_position = legend_position))
  
  res_grid$mse_pre_original[i] = res$mse_pre_original
  res_grid$mse_pre_new[i] = res$mse_pre_new
  res_grid$pos_ratio[i] = res$pos_ratio
  res_grid$t_test[i] = res$t_test
  
  cat("Done.\n")
  gc()
}

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

res_grid = expand.grid(width_range, k_range,
                       names(step_pattern_range)) %>% 
  `colnames<-`(c("width", "k", "step_pattern")) %>% 
  mutate(mse_pre_original = NA_real_,
         mse_pre_new = NA_real_,
         mse_post_original = NA_real_,
         mse_post_new = NA_real_,
         pos_ratio = NA_real_,
         t_test = NA_real_)

res_grid_filename = "./data/grid_search_v6/res_grid_tobacco_89_fixed.Rds"
# saveRDS(res_grid, res_grid_filename)
res_grid = readRDS(res_grid_filename)

# search
start_time = 1970
end_time = 2000
treat_time = 1989
dtw1_time = 1989
TSDTW_type = "fixed"
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
legend_position = c(0.3, 0.3)

mse_list = NULL
for (i in which(is.na(res_grid$pos_ratio))) {  # which(is.na(res_grid$pos_ratio))
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
                                       TSDTW_type = TSDTW_type,
                                       filter_width = width,
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
  
  res_grid$mse_pre_original[i] = median(res$mse$mse1_pre)
  res_grid$mse_pre_new[i] = median(res$mse$mse2_pre)
  res_grid$mse_post_original[i] = median(res$mse$mse1_post)
  res_grid$mse_post_new[i] = median(res$mse$mse2_post)
  res_grid$pos_ratio[i] = res$pos_ratio
  res_grid$t_test[i] = res$t_test
  mse_list[[i]] = res$mse
  
  cat("Done.\n")
  gc()
}

## Optimal Run Tobacco ---------------------------------------------------------
# prepare data
start_time = 1970
end_time = 2000
treat_time = 1989
dtw1_time = 1989
filter_width = 9
k = 6
TSDTW_type = "fixed"
n_mse = 10
n_IQR = 3
dist_quantile = 0.95
plot_figures = TRUE
step.pattern1 = dtw::mori2006
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
legend_position = c(0.3, 0.3)
... = NULL
normalize_method = "t"
ma = 3
ma_na = "original"
n_q = 1
n_r = 1


data = preprocessing(data, filter_width)
units = data[c("id", "unit")] %>% distinct

# run
result = as.list(1:nrow(units)) %>% 
  future_map(
    ~{
      i = .
      dependent = units$unit[i]
      dependent_id = units$id[i]
      # print(paste0(dependent, ":", i, "-", k, " start..."))
      res = SimDesign::quiet(compare_methods(data = data,
                                             start_time = start_time,
                                             end_time = end_time,
                                             treat_time = treat_time,
                                             dtw1_time = dtw1_time,
                                             dependent = dependent,
                                             dependent_id = dependent_id,
                                             n_mse = n_mse,
                                             k = k,
                                             TSDTW_type = TSDTW_type,
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
      # print(paste0(dependent, ":", i, "-", k, " start...Done."))
      res$mse = res$mse %>% mutate(dependent = dependent, k = k)
      res
    }
  )

saveRDS(result, "./data/grid_search_v6/result_tobacco_89.Rds")


# mse
mse = result %>% 
  lapply(., "[[", "mse") %>% 
  do.call("rbind", .) %>% 
  mutate(ratio = mse2_post/mse1_post,
         log_ratio = log(ratio))
mse = mse %>% filter(dependent != "California")
# mse = mse %>% filter(mse1_pre < 5*3)
length(which(mse$log_ratio < 0))/nrow(mse)
boxplot(mse$log_ratio, outline = FALSE)
abline(h = 0, lty = 5)

t.test(mse$log_ratio)

saveRDS(mse, "./data/grid_search_v6/mse_tobacco_89.Rds")


# plot time series figure
df = rbind(data.frame(unit = "California",
                      time = 1970:2000,
                      value = result[[3]]$synth_origin$value),
           data.frame(unit = "Synthetic Control w/o TFDTW",
                      time = 1970:2000,
                      value = result[[3]]$synth_origin$synthetic),
           data.frame(unit = "Synthetic Control w/ TFDTW",
                      time = 1970:2000,
                      value = result[[3]]$synth_new$synthetic))

fig = ggplot(df, aes(x = time, y = value, color = unit)) +
  geom_line() + 
  geom_vline(xintercept = 1988, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values=c("grey40", "#fe4a49","#2ab7ca")) +
  xlab("Time") +
  ylab("Cigarette Sale Per Capita (in Packs)") + 
  theme(legend.title=element_blank(),
        legend.position = c(0.2, 0.2))

saveRDS(fig, "./data/grid_search_v2/fig_comp_method_tobacco_89.Rds")

# plot placebo figure
df = result %>% 
  map(
    ~{
      data.frame(unit = .[["mse"]][["dependent"]],
                 time = 1970:2000,
                 value = .[["synth_origin"]][["value"]],
                 synth_origin = .[["synth_origin"]][["synthetic"]],
                 synth_new = .[["synth_new"]][["synthetic"]])
    }
  ) %>% 
  do.call("rbind", .) %>% 
  mutate(
    color = case_when(unit == "California" ~ "black",
                      TRUE ~ "grey 70"),
    gap_origin = value - synth_origin,
    gap_new = value - synth_new
  )

df %>% 
  # filter(unit %in% (mse %>% filter(mse1_pre < 2*3) %>% .[["dependent"]])) %>% 
  ggplot(aes(x = time, group = unit)) +
  geom_line(aes(y = gap_origin), col = "#adcbe3") +
  geom_line(aes(y = gap_new), col = "#fec8c1") +
  geom_line(aes(y = gap_origin), data = df %>% filter(unit == "California"), col = "#2ab7ca", size = 1) +
  geom_line(aes(y = gap_new), data = df %>% filter(unit == "California"), col = "#fe4a49", size = 1) +
  geom_vline(xintercept = 1988, linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  coord_cartesian(ylim=c(-32, 32)) +
  xlab("year") +
  ylab("gap in per-capita cigarette sales (in packs)") +
  theme_bw()









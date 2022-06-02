## Setup -----------------------------------------------------------------------
library(checkpoint)
checkpoint("2022-04-01")

library(tidyverse)
library(furrr)
plan(multisession, workers = 11)
options(future.rng.onMisuse="ignore")
options(stringsAsFactors = FALSE)

source("./R/TwoStepDTW.R")
source("./R/synthetic_control.R")
source("./R/comp_methods.R")
set.seed(20220407)


## Basque Terrorism Data --------------------------------------------------
data(basque, package = "Synth")
data = basque
colnames(data)[1:4] = c("id", "unit", "time", "value")
data = data %>% mutate(invest_ratio = invest/value,
                       value_raw = value)
data = data %>% filter(unit != "Basque Country (Pais Vasco)")


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
         age15to24 = age15to24*100) %>% 
  filter(unit != "Rhode Island")

data = smoking
data = data %>% filter(unit != "California")



## Germany Reunification Data --------------------------------------------------
data = foreign::read.dta("./data/repgermany.dta")
colnames(data)[1:4] = c("id", "unit", "time", "value")
data = data %>% mutate(value_raw = value)


## Grid Search -----------------------------------------------------------------
# search space
width_range = (1:7)*2+3
k_range = 4:9
dtw1_time_range = 1989:1996
start_time = 1970
end_time = 1999
treat_time = 1989

res_grid = NULL      
for (width in width_range) {
  for (k in k_range) {
    for (dtw1_time in dtw1_time_range) {
      temp = data.frame(width = width,
                        k = k,
                        start_time = start_time,
                        end_time = end_time,
                        treat_time = treat_time,
                        dtw1_time = dtw1_time,
                        pos_ratio = NA_real_,
                        t_test = NA_real_)
      res_grid = rbind(res_grid, temp)
    }
  }
}

res_grid_filename = "./data/res_grid_tobacco_89.Rds"
# saveRDS(res_grid, res_grid_filename)
res_grid = readRDS(res_grid_filename)

# search
synth_fun = "tobacco-89"

for (i in which(is.na(res_grid$pos_ratio))) {
  width = res_grid$width[i]
  k = res_grid$k[i]
  dtw1_time = res_grid$dtw1_time[i]
  start_time = res_grid$start_time[i]
  end_time = res_grid$end_time[i]
  treat_time = res_grid$treat_time[i]
  
  cat(paste0("[Task-", i, "]: width-", width, ", k-", k,
               ", dtw1_time-", dtw1_time, "......"))
  
  data = preprocessing(data, filter_width = width)
  
  res = SimDesign::quiet(run_all_units(data = data,
                                       start_time = start_time,
                                       end_time = end_time,
                                       treat_time = treat_time,
                                       dtw1_time = dtw1_time,
                                       # plot_figures = FALSE, 
                                       # normalize_method = "t",
                                       # step.pattern = dtw::symmetricP2,
                                       # legend_position = c(0.3, 0.3),
                                       filter_width = width,
                                       k = k,
                                       synth_fun = synth_fun))
  
  res_grid$pos_ratio[i] = res$pos_ratio
  res_grid$t_test[i] = res$t_test
  
  cat("Done.\n")
  gc()
}



## Optimal Run -----------------------------------------------------------------
# prepare data
start_time = 1955
end_time = 1980
treat_time = 1970
dtw1_time = 1971
plot_figures = FALSE
normalize_method = "t"
dtw_method = "dtw"
step.pattern = dtw::symmetricP2
legend_position = c(0.3, 0.3)
filter_width = 5
k = 6
synth_fun = "basque-70"

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
      res = compare_methods(data = data,
                            start_time = start_time,
                            end_time = end_time,
                            treat_time = treat_time,
                            dtw1_time = dtw1_time,
                            dependent = dependent,
                            dependent_id = dependent_id,
                            normalize_method = "t",
                            filter_width = width,
                            k = k,
                            plot_figures = F,
                            synth_fun = synth_fun,
                            step.pattern = dtw::symmetricP2)
      # print(paste0(dependent, ":", i, "-", k, " start...Done."))
      res$mse %>% mutate(dependent = dependent, k = k)
    }
  )

result = result %>% 
  do.call("rbind", .) %>% 
  mutate(ratio = mse2_post/mse1_post,
         log_ratio = log(ratio))
length(which(result$log_ratio < 0))/nrow(result)
boxplot(result$log_ratio, outline = FALSE)
abline(h = 0, lty = 5)

t.test(result$log_ratio)

saveRDS(result, "./data/result_tobacco_1986.Rds")

## Results ---------------------------------------------------------------------
result_1985 = readRDS("./data/result_tobacco_1985.Rds")
result_1986 = readRDS("./data/result_tobacco_1986.Rds")
result_1987 = readRDS("./data/result_tobacco_1987.Rds")
result_1989 = readRDS("./data/result_tobacco_1989.Rds")
result_1990 = readRDS("./data/result_gdp_1990.Rds")


result = rbind(result_1985 %>% mutate(treatment = "Tobacco_1985"),
               result_1986 %>% mutate(treatment = "Tobacco_1986"),
               result_1987 %>% mutate(treatment = "Tobacco_1987"),
               result_1989 %>% mutate(treatment = "Tobacco_1989"),
               result_1990 %>% mutate(treatment = "GDP_1990"))

result = result %>% filter(dependent != "Rhode Island")
result = result %>% mutate(ratio = mse2_post/mse1_post,
                           log_ratio = log(ratio))

ggplot(result, aes(x=treatment, y=log_ratio)) + 
  geom_boxplot() +
  theme_bw() +
  # coord_cartesian(ylim = c(-40, 70)) +
  geom_hline(yintercept=0, linetype="dashed")

t_test = result %>% group_by(treatment) %>% 
  summarise(t_test = t.test(log_ratio)$p.value)





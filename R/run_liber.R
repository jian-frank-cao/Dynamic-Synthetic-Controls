## Setup -----------------------------------------------------------------------
library(checkpoint)
checkpoint("2022-04-01")

library(tidyverse)
library(furrr)
plan(multisession, workers = 7)
options(future.rng.onMisuse="ignore")
options(stringsAsFactors = FALSE)

source("./R/TwoStepDTW.R")
source("./R/synthetic_control.R")
source("./R/comp_methods.R")
set.seed(20220407)


## Liberalization Data ---------------------------------------------------------
data = foreign::read.dta("./data/PTpanelRESTAT.dta") %>%
  filter(ctycode %in% c(104,22,121,162,6,68,69,114,124,174)) %>% 
  mutate(countryname = case_when(countryname == "" ~ NA_character_,
                                 TRUE ~ countryname),
         inv_ratio = inv_ratio*100,
         unit = zoo::na.locf(countryname),
         id = ctycode,
         time = year,
         value = rgdppp,
         value_raw = value) %>% 
  select(c("id","unit","time","value","school2","pop_growth",
           "inflation","democracy","inv_ratio","value_raw")) %>% 
  filter(time >= 1964)

# mexico = data %>% filter(unit == "Mexico" & time %in% c(1960:1985))
# colMeans(mexico[c("value","school2","pop_growth",
#                   "inflation","democracy","inv_ratio")],na.rm = T)

## Synth Function --------------------------------------------------------------
do_synth_mexico_86 = function(df, dep_var, dependent_id,
                          start_time = 1964, end_time = 2005,
                          t_treat = 1986){
  # find v
  dataprep.out <-
    Synth::dataprep(
      foo = df,
      predictors    = dep_var,
      dependent     = dep_var,
      unit.variable = 1,
      time.variable = 3,
      special.predictors = list(
        list("school2", start_time:(t_treat - 1), c("mean")),
        list("pop_growth", start_time:(t_treat - 1), c("mean")),
        list("inflation", start_time:(t_treat - 1), c("mean")),
        list("democracy", start_time:(t_treat - 1), c("mean")),
        list("inv_ratio", start_time:(t_treat - 1), c("mean"))
      ),
      treatment.identifier = dependent_id,
      controls.identifier = setdiff(unique(df$id), dependent_id),
      time.predictors.prior = start_time:(t_treat - 1),
      time.optimize.ssr = start_time:(t_treat - 1), 
      unit.names.variable = 2,
      time.plot = start_time:end_time
    )
  
  # fit training model
  synth.out <- 
    Synth::synth(
      data.prep.obj=dataprep.out,
      Margin.ipop=.005,Sigf.ipop=7,Bound.ipop=6
    )
  
  value = df %>% filter(id == dependent_id) %>% `$`(value_raw)
  average = df %>% filter(id != dependent_id) %>% group_by(time) %>% 
    summarise(average = mean(value_raw, na.rm = TRUE)) %>% `$`(average)
  synthetic = dataprep.out$Y0%*%synth.out$solution.w %>% as.numeric
  
  return(list(value = value,
              average = average,
              synthetic = synthetic))
}

## Grid Search -----------------------------------------------------------------
# search space
width_range = (1:7)*2+3
k_range = 4:9
start_time = 1964
treat_time = 1986
end_time = 2005
dtw1_time_range = treat_time:(treat_time + 7)



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

res_grid_filename = "./data/res_grid_mexico_86_P2.Rds"
# saveRDS(res_grid, res_grid_filename)
res_grid = readRDS(res_grid_filename)

# search
synth_fun = "mexico-86"

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
                                       plot_figures = F,
                                       # normalize_method = "t",
                                       step.pattern = dtw::symmetric2,
                                       # legend_position = c(0.3, 0.3),
                                       filter_width = width,
                                       n_mse = 10,
                                       k = k,
                                       synth_fun = synth_fun))
  
  res_grid$pos_ratio[i] = res$pos_ratio
  res_grid$t_test[i] = res$t_test
  
  cat("Done.\n")
  gc()
}



## Optimal Run Basque ----------------------------------------------------------
# prepare data
start_time = 1964
end_time = 2005
treat_time = 1986
dtw1_time = 1990
plot_figures = FALSE
normalize_method = "t"
dtw_method = "dtw"
step.pattern = dtw::symmetricP2
legend_position = c(0.3, 0.3)
filter_width = 13
k = 9
synth_fun = "mexico-86"

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
                            n_mse = 10,
                            k = k,
                            plot_figures = F,
                            synth_fun = synth_fun,
                            step.pattern = step.pattern)
      # print(paste0(dependent, ":", i, "-", k, " start...Done."))
      res$mse = res$mse %>% mutate(dependent = dependent, k = k)
      res
    }
  )

saveRDS(result, "./data/grid_search_v2/result_basque_70.Rds")


# mse
mse = result %>% 
  lapply(., "[[", "mse") %>% 
  do.call("rbind", .) %>% 
  mutate(ratio = mse2_post/mse1_post,
         log_ratio = log(ratio))
mse = mse %>% filter(!(dependent %in% c("Panama", "Mexico")))
length(which(mse$log_ratio < 0))/nrow(mse)
boxplot(mse$log_ratio, outline = FALSE)
abline(h = 0, lty = 5)

t.test(mse$log_ratio)

saveRDS(mse, "./data/grid_search_v2/mse_basque_70.Rds")


# plot figure
df = rbind(data.frame(unit = "Mexico",
                      time = 1964:2005,
                      value = result[[7]]$synth_origin$value),
           data.frame(unit = "Synthetic Control w/o TFDTW",
                      time = 1964:2005,
                      value = result[[7]]$synth_origin$synthetic),
           data.frame(unit = "Synthetic Control w/ TFDTW",
                      time = 1964:2005,
                      value = result[[7]]$synth_new$synthetic))

fig = ggplot(df, aes(x = time, y = value, color = unit)) +
  geom_line() + 
  geom_vline(xintercept = 1986, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values=c("grey40", "#fe4a49","#2ab7ca")) +
  xlab("Time") +
  ylab("Per Capita GDP (1986 USD Thousand)") + 
  theme(legend.title=element_blank(),
        legend.position = c(0.8, 0.2))

saveRDS(fig, "./data/grid_search_v2/fig_comp_method_basque_70.Rds")



## Optimal Run Tobacco ---------------------------------------------------------
# prepare data
start_time = 1970
end_time = 1999
treat_time = 1989
dtw1_time = 1995
plot_figures = FALSE
normalize_method = "t"
dtw_method = "dtw"
step.pattern = dtw::symmetricP2
legend_position = c(0.3, 0.3)
filter_width = 7
k = 8
synth_fun = "tobacco-89"

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
      res$mse = res$mse %>% mutate(dependent = dependent, k = k)
      res
    }
  )

saveRDS(result, "./data/grid_search_v2/result_tobacco_89.Rds")


# mse
mse = result %>% 
  lapply(., "[[", "mse") %>% 
  do.call("rbind", .) %>% 
  mutate(ratio = mse2_post/mse1_post,
         log_ratio = log(ratio))
mse = mse %>% filter(dependent != "California")
length(which(mse$log_ratio < 0))/nrow(mse)
boxplot(mse$log_ratio, outline = FALSE)
abline(h = 0, lty = 5)

t.test(mse$log_ratio)

saveRDS(mse, "./data/grid_search_v2/mse_tobacco_89.Rds")


# plot figure
df = rbind(data.frame(unit = "California",
                      time = 1970:1999,
                      value = result[[3]]$synth_origin$value),
           data.frame(unit = "Synthetic Control w/o TFDTW",
                      time = 1970:1999,
                      value = result[[3]]$synth_origin$synthetic),
           data.frame(unit = "Synthetic Control w/ TFDTW",
                      time = 1970:1999,
                      value = result[[3]]$synth_new$synthetic))

fig = ggplot(df, aes(x = time, y = value, color = unit)) +
  geom_line() + 
  geom_vline(xintercept = 1989, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values=c("grey40", "#fe4a49","#2ab7ca")) +
  xlab("Time") +
  ylab("Cigarette Sale Per Capita (in Packs)") + 
  theme(legend.title=element_blank(),
        legend.position = c(0.2, 0.2))

saveRDS(fig, "./data/grid_search_v2/fig_comp_method_tobacco_89.Rds")


## Optimal Run Germany ---------------------------------------------------------
# prepare data
start_time = 1960
end_time = 2003
treat_time = 1990
dtw1_time = 1997
plot_figures = FALSE
normalize_method = "t"
dtw_method = "dtw"
step.pattern = dtw::symmetricP2
legend_position = c(0.3, 0.3)
filter_width = 5
k = 6
synth_fun = "germany-90"

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
      res$mse = res$mse %>% mutate(dependent = dependent, k = k)
      res
    }
  )

saveRDS(result, "./data/grid_search_v2/result_germany_90.Rds")


# mse
mse = result %>% 
  lapply(., "[[", "mse") %>% 
  do.call("rbind", .) %>% 
  mutate(ratio = mse2_post/mse1_post,
         log_ratio = log(ratio))
mse = mse %>% filter(dependent != "West Germany")
length(which(mse$log_ratio < 0))/nrow(mse)
boxplot(mse$log_ratio, outline = FALSE)
abline(h = 0, lty = 5)

t.test(mse$log_ratio)

saveRDS(mse, "./data/grid_search_v2/mse_germany_90.Rds")


# plot figure
df = rbind(data.frame(unit = "West Germany",
                      time = 1960:2003,
                      value = result[[17]]$synth_origin$value),
           data.frame(unit = "Synthetic Control w/o TFDTW",
                      time = 1960:2003,
                      value = result[[17]]$synth_origin$synthetic),
           data.frame(unit = "Synthetic Control w/ TFDTW",
                      time = 1960:2003,
                      value = result[[17]]$synth_new$synthetic))

df$unit = factor(df$unit, levels = c("West Germany",
                                     "Synthetic Control w/o TFDTW",
                                     "Synthetic Control w/ TFDTW"))

fig = ggplot(df, aes(x = time, y = value, color = unit)) +
  geom_line() + 
  geom_vline(xintercept = 1990, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values=c("grey40", "#fe4a49","#2ab7ca")) +
  xlab("Time") +
  ylab("Per Capita GDP (PPP, 2002 USD)") + 
  theme(legend.title=element_blank(),
        legend.position = c(0.3, 0.8))

saveRDS(fig, "./data/grid_search_v2/fig_comp_method_germany_90.Rds")


## Results ---------------------------------------------------------------------
# boxplot
mse_basque = readRDS("./data/grid_search_v2/mse_basque_70.Rds")
mse_tobacco = readRDS("./data/grid_search_v2/mse_tobacco_89.Rds")
mse_germany = readRDS("./data/grid_search_v2/mse_germany_90.Rds")

mse = rbind(mse_basque %>% mutate(data = "Basque"),
            mse_tobacco %>% mutate(data = "Tobacco"),
            mse_germany %>% mutate(data = "Germany"))

t_test = mse %>% group_by(data) %>% 
  summarise(t_test = t.test(log_ratio)$p.value)

fig = ggplot(mse, aes(x=data, y=log_ratio)) + 
  geom_boxplot() +
  theme_bw() +
  # coord_cartesian(ylim = c(-40, 70)) +
  geom_hline(yintercept=0, linetype="dashed") +
  xlab("") +
  ylab("Log Ratio") +
  annotate("text", x = 1:3, y = -0.08, label = c("P=0.05", "P=0.05", "P=0.01"))

ggsave("./figures/abadie_boxplot.pdf",
       fig, width = 5, height = 4,
       units = "in", limitsize = FALSE)

# plot time series
fig_basque = readRDS("./data/grid_search_v2/fig_comp_method_basque_70.Rds")
fig_tobacco = readRDS("./data/grid_search_v2/fig_comp_method_tobacco_89.Rds")
fig_germany = readRDS("./data/grid_search_v2/fig_comp_method_germany_90.Rds")

ggsave("./figures/comp_method_basque_70.pdf",
       fig_basque, width = 6, height = 5,
       units = "in", limitsize = FALSE)

ggsave("./figures/comp_method_tobacco_89.pdf",
       fig_tobacco, width = 6, height = 5,
       units = "in", limitsize = FALSE)

ggsave("./figures/comp_method_germany_90.pdf",
       fig_germany, width = 6, height = 5,
       units = "in", limitsize = FALSE)


fig_comp_method = gridExtra::marrangeGrob(list(basque = fig_basque,
                                   tobacco = fig_tobacco,
                                   germany = fig_germany),
                              ncol = 3,
                              nrow = 1)
ggsave("./figures/comp_method_20220627.pdf",
       fig_comp_method, width = 3*4, height = 1*4,
       units = "in", limitsize = FALSE)


# result_1985 = readRDS("./data/result_tobacco_1985.Rds")
# result_1986 = readRDS("./data/result_tobacco_1986.Rds")
# result_1987 = readRDS("./data/result_tobacco_1987.Rds")
# result_1989 = readRDS("./data/result_tobacco_1989.Rds")
# result_1990 = readRDS("./data/result_gdp_1990.Rds")
# 
# 
# result = rbind(result_1985 %>% mutate(treatment = "Tobacco_1985"),
#                result_1986 %>% mutate(treatment = "Tobacco_1986"),
#                result_1987 %>% mutate(treatment = "Tobacco_1987"),
#                result_1989 %>% mutate(treatment = "Tobacco_1989"),
#                result_1990 %>% mutate(treatment = "GDP_1990"))
# 
# result = result %>% filter(dependent != "Rhode Island")
# result = result %>% mutate(ratio = mse2_post/mse1_post,
#                            log_ratio = log(ratio))
# 
# ggplot(result, aes(x=treatment, y=log_ratio)) + 
#   geom_boxplot() +
#   theme_bw() +
#   # coord_cartesian(ylim = c(-40, 70)) +
#   geom_hline(yintercept=0, linetype="dashed")
# 
# t_test = result %>% group_by(treatment) %>% 
#   summarise(t_test = t.test(log_ratio)$p.value)
# 




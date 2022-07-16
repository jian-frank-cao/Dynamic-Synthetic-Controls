## Setup -----------------------------------------------------------------------
library(checkpoint)
checkpoint("2022-04-01")

library(tidyverse)
library(furrr)
plan(multisession, workers = 7)
options(future.rng.onMisuse="ignore")
options(stringsAsFactors = FALSE)

source("./R/TwoStepDTW_OpenEnd.R")
source("./R/synthetic_control.R")
source("./R/comp_methods.R")
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


## Grid Search -----------------------------------------------------------------
# search space
width_range = (1:9)*2+3
k_range = 4:9
start_time = 1970
treat_time = 1989
end_time = 2000
dtw1_time_range = 1994
step_pattern_range = list(
  # symmetricP0 = dtw::symmetricP0, # too bumpy
  symmetricP05 = dtw::symmetricP05,
  symmetricP1 = dtw::symmetricP1,
  symmetricP2 = dtw::symmetricP2,
  asymmetricP0 = dtw::asymmetricP0,
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


res_grid = NULL      
for (width in width_range) {
  for (k in k_range) {
    for (dtw1_time in dtw1_time_range) {
      for (pattern_name in names(step_pattern_range)) {
        temp = data.frame(width = width,
                          k = k,
                          step_pattern = pattern_name,
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
}

res_grid_filename = "./data/grid_search_v4/res_grid_tobacco_89.Rds"
# saveRDS(res_grid, res_grid_filename)
res_grid = readRDS(res_grid_filename)

# search
synth_fun = "tobacco-89"

for (i in which(is.na(res_grid$pos_ratio))) {
  width = res_grid$width[i]
  k = res_grid$k[i]
  pattern_name = res_grid$step_pattern[i]
  step.pattern = step_pattern_range[[pattern_name]]
  dtw1_time = res_grid$dtw1_time[i]
  start_time = res_grid$start_time[i]
  end_time = res_grid$end_time[i]
  treat_time = res_grid$treat_time[i]
  
  cat(paste0("[Task-", i, "]: width-", width, ", k-", k,
               ", dtw1-", dtw1_time, ", ", pattern_name, "......"))
  
  data = preprocessing(data, filter_width = width)
  
  res = SimDesign::quiet(run_all_units(data = data,
                                       start_time = start_time,
                                       end_time = end_time,
                                       treat_time = treat_time,
                                       dtw1_time = dtw1_time,
                                       # plot_figures = FALSE, 
                                       # normalize_method = "t",
                                       step.pattern = step.pattern,
                                       # legend_position = c(0.3, 0.3),
                                       filter_width = width,
                                       k = k,
                                       synth_fun = synth_fun))
  
  res_grid$pos_ratio[i] = res$pos_ratio
  res_grid$t_test[i] = res$t_test
  
  cat("Done.\n")
  gc()
}



## Optimal Run Basque ----------------------------------------------------------
# prepare data
start_time = 1955
end_time = 1980
treat_time = 1970
dtw1_time = 1970
plot_figures = FALSE
normalize_method = "t"
dtw_method = "dtw"
step.pattern = dtw::symmetricP2
legend_position = c(0.3, 0.3)
filter_width = 5
k = 5
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
                            filter_width = filter_width,
                            k = k,
                            plot_figures = F,
                            synth_fun = synth_fun,
                            step.pattern = dtw::symmetricP2)
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
mse = mse %>% filter(dependent != "Basque Country (Pais Vasco)")
length(which(mse$log_ratio < 0))/nrow(mse)
boxplot(mse$log_ratio, outline = FALSE)
abline(h = 0, lty = 5)

t.test(mse$log_ratio)

saveRDS(mse, "./data/grid_search_v2/mse_basque_70.Rds")


# plot time series figure
df = rbind(data.frame(unit = "Basque Country",
                      time = 1955:1980,
                      value = result[[4]]$synth_origin$value),
           data.frame(unit = "Synthetic Control w/o TFDTW",
                      time = 1955:1980,
                      value = result[[4]]$synth_origin$synthetic),
           data.frame(unit = "Synthetic Control w/ TFDTW",
                      time = 1955:1980,
                      value = result[[4]]$synth_new$synthetic))

fig = ggplot(df, aes(x = time, y = value, color = unit)) +
  geom_line() + 
  geom_vline(xintercept = 1970, linetype="dashed") +
  theme_bw() +
  scale_color_manual(values=c("grey40", "#fe4a49","#2ab7ca")) +
  xlab("Time") +
  ylab("Per Capita GDP (1986 USD Thousand)") + 
  theme(legend.title=element_blank(),
        legend.position = c(0.8, 0.2))

saveRDS(fig, "./data/grid_search_v2/fig_comp_method_basque_70.Rds")


# plot placebo figure
df = result %>% 
  map(
    ~{
      data.frame(unit = .[["mse"]][["dependent"]],
                 time = 1955:1980,
                 value = .[["synth_origin"]][["value"]],
                 synth_origin = .[["synth_origin"]][["synthetic"]],
                 synth_new = .[["synth_new"]][["synthetic"]])
    }
  ) %>% 
  do.call("rbind", .) %>% 
  mutate(
    color = case_when(unit == "Basque Country (Pais Vasco)" ~ "black",
                      TRUE ~ "grey 70"),
    gap_origin = value - synth_origin,
    gap_new = value - synth_new
  )

df %>% 
  filter(unit %in% (mse %>% filter(mse1_pre < 2*0.001) %>% .[["dependent"]])) %>% 
  ggplot(aes(x = time, group = unit)) +
  geom_line(aes(y = gap_origin), col = "#adcbe3") +
  geom_line(aes(y = gap_new), col = "#fec8c1") +
  geom_line(aes(y = gap_origin), data = df %>% filter(unit == "Basque Country (Pais Vasco)"), col = "#2ab7ca", size = 1) +
  geom_line(aes(y = gap_new), data = df %>% filter(unit == "Basque Country (Pais Vasco)"), col = "#fe4a49", size = 1) +
  geom_vline(xintercept = 1970, linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  coord_cartesian(ylim=c(-0.75, 0.75)) +
  xlab("year") +
  ylab("gap in per-capita cigarette sales (in packs)") +
  theme_bw()

## Optimal Run Tobacco ---------------------------------------------------------
# prepare data
start_time = 1970
end_time = 2000
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
                            filter_width = filter_width,
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
# mse = mse %>% filter(mse1_pre < 5*3)
length(which(mse$log_ratio < 0))/nrow(mse)
boxplot(mse$log_ratio, outline = FALSE)
abline(h = 0, lty = 5)

t.test(mse$log_ratio)

saveRDS(mse, "./data/grid_search_v2/mse_tobacco_89.Rds")


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
  filter(unit %in% (mse %>% filter(mse1_pre < 2*3) %>% .[["dependent"]])) %>% 
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


## Optimal Run Germany ---------------------------------------------------------
# prepare data
start_time = 1960
end_time = 2003
treat_time = 1990
dtw1_time = 1995
plot_figures = FALSE
normalize_method = "t"
dtw_method = "dtw"
step.pattern = dtw::symmetric2
legend_position = c(0.3, 0.3)
filter_width = 5
k = 5
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
                            filter_width = filter_width,
                            k = k,
                            plot_figures = F,
                            synth_fun = synth_fun,
                            step.pattern = step.pattern)
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


# plot time series figure
df = rbind(data.frame(unit = "West Germany",
                      time = 1960:2000,
                      value = result[[17]]$synth_origin$value[-c(42:44)]),
           data.frame(unit = "Synthetic Control w/o TFDTW",
                      time = 1960:2000,
                      value = result[[17]]$synth_origin$synthetic[-c(42:44)]),
           data.frame(unit = "Synthetic Control w/ TFDTW",
                      time = 1960:2000,
                      value = result[[17]]$synth_new$synthetic[-c(42:44)]))

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


# plot placebo figure
df = result %>% 
  map(
    ~{
      data.frame(unit = .[["mse"]][["dependent"]],
                 time = 1960:2000,
                 value = .[["synth_origin"]][["value"]][-c(42:44)],
                 synth_origin = .[["synth_origin"]][["synthetic"]][-c(42:44)],
                 synth_new = .[["synth_new"]][["synthetic"]][-c(42:44)])
    }
  ) %>% 
  do.call("rbind", .) %>% 
  mutate(
    color = case_when(unit == "West Germany" ~ "black",
                      TRUE ~ "grey 70"),
    gap_origin = value - synth_origin,
    gap_new = value - synth_new
  )

df %>% 
  filter(unit %in% (mse %>% filter(mse1_pre < 2*10000) %>% .[["dependent"]])) %>% 
  ggplot(aes(x = time, group = unit)) +
  geom_line(aes(y = gap_origin), col = "#adcbe3") +
  geom_line(aes(y = gap_new), col = "#fec8c1") +
  geom_line(aes(y = gap_origin), data = df %>% filter(unit == "West Germany"), col = "#2ab7ca", size = 1) +
  geom_line(aes(y = gap_new), data = df %>% filter(unit == "West Germany"), col = "#fe4a49", size = 1) +
  geom_vline(xintercept = 1990, linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  # coord_cartesian(ylim=c(-32, 32)) +
  xlab("year") +
  ylab("gap in per-capita cigarette sales (in packs)") +
  theme_bw()


## Results ---------------------------------------------------------------------
# Abadie boxplot
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


# Abadie + simulation boxplot
mse_basque = readRDS("./data/grid_search_v2/mse_basque_70.Rds")
mse_tobacco = readRDS("./data/grid_search_v2/mse_tobacco_89.Rds")
mse_germany = readRDS("./data/grid_search_v2/mse_germany_90.Rds")

mse_basque = mse_basque %>% filter(dependent != "Basque Country (Pais Vasco)")
mse_tobacco = mse_basque %>% filter(dependent != "California")
mse_germany = mse_basque %>% filter(dependent != "West Germany")

result_simul = readRDS("./data/res_simul_0626_v1.Rds")

mse_simul = result %>% 
  map(
    ~{
      res = .
      min_ratio = min(res$mse_ratio)
      res %>% filter(mse_ratio == min_ratio) %>% .[1,]
    }
  ) %>% 
  do.call("rbind", .)

df = rbind(mse_basque[c("mse1_post", "mse2_post")] %>% mutate(group = "Basque"),
           mse_tobacco[c("mse1_post", "mse2_post")] %>% mutate(group = "Tobacco"),
           mse_germany[c("mse1_post", "mse2_post")] %>% mutate(group = "Germany"),
           mse_simul %>% mutate(mse1_post = mse_original,
                                mse2_post = mse_new,
                                group = "Simulation") %>% select(c("mse1_post", "mse2_post", "group")))

df = df %>% mutate(mse1 = log(mse1_post),
                   mse2 = log(mse2_post),
                   log_ratio = log(mse2_post/mse1_post),
                   group = factor(group, levels = c("Basque", "Tobacco", "Germany", "Simulation")))

fig = ggplot(df, aes(x=group, y=log_ratio)) + 
  geom_boxplot() +
  theme_bw() +
  # coord_cartesian(ylim = c(-40, 70)) +
  geom_hline(yintercept=0, linetype="dashed") +
  xlab("") +
  ylab("Log Ratio") +
  annotate("text", x = 1:4 + 0.2, y = 0.1, label = c("P=0.05", "P=0.05", "P=0.01", "P=0"))

ggsave("./figures/comp_method_boxplot.pdf",
       fig, width = 6, height = 5,
       units = "in", limitsize = FALSE)

# cloud plot
fig = df %>% 
  ggplot(aes(x = mse1, y = mse2, color = group)) +
  geom_mark_hull(aes(fill = group, label = group), concavity = 3) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  xlab("Log MSE (Abadie)") +
  ylab("Log MSE (TFDTW)") +
  theme_bw()

ggsave("./figures/comp_method_cloud.pdf",
       fig, width = 6, height = 5,
       units = "in", limitsize = FALSE)

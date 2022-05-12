# setwd("F:/Trinity/Dynamic-Synthetic-Control")

library(checkpoint)
checkpoint("2022-04-01")

library(tidyverse)
library(furrr)
plan(multisession, workers = 6)
options(future.rng.onMisuse="ignore")


source("./R/TwoStepDTW.R")
source("./R/synthetic_control.R")
set.seed(20220407)

## Data ------------------------------------------------------------------------
load("./data/smoking.rda")
prop99 = read.csv("./data/prop99.csv")

exclude_states = c("Massachusetts", "Arizona", "Oregon", "Florida",
                   "Alaska", "Hawaii", "Maryland", "Michigan",
                   "New Jersey", "New York",
                   "Washington", "District of Columbia")
states = data.frame(id = 1:39,
                    state = sort(setdiff(unique(prop99$LocationDesc),
                                   exclude_states)))
smoking = smoking %>% mutate_all(as.numeric)
colnames(smoking)[1] = "id"
smoking = right_join(states, smoking, by = "id")
colnames(smoking)[2:4] = c("unit", "time", "value")
smoking = smoking %>% mutate(value_raw = value)

## Pre-processing --------------------------------------------------------------
values = reshape2::dcast(smoking %>% select(c("unit", "time", "value_raw")),
                         time ~ unit, value.var = "value_raw")

# transform
values = values %>% mutate_at(setdiff(colnames(values), "time"),
                              ~normalize(., "t"))

# adding buffer
add_buffer = function(TS, n){
  model_right = forecast::auto.arima(TS)
  right <- as.numeric(forecast::forecast(model_right, h = n)$mean)
  model_left = forecast::auto.arima(rev(TS))
  left <- rev(as.numeric(forecast::forecast(model_left, h = n)$mean))
  
  return(c(left, TS, right))
}

width = 11
values2 = sapply(values %>% select(-time),
                 add_buffer, n = (width - 1)/2) %>% 
  data.frame(.)

# derivative
values2 = values2 %>%
  mutate_all(~signal::sgolayfilt(., 3, width, 2)) %>%
  .[((width - 1)/2 + 1):((width - 1)/2 + nrow(values)),]
values[-1] = values2

# join data 
df <- reshape2::melt(values ,  id.vars = 'time',
                     variable.name = 'unit')

smoking = right_join(df, smoking %>% select(-value), by = c("time", "unit"))
smoking$age15to24 = smoking$age15to24*100
smoking = smoking[c("id", "unit", "time", "value", 
                    "lnincome", "beer", "age15to24",
                    "retprice", "value_raw")]

## Function --------------------------------------------------------------------
compare_methods = function(data,
                           start_time,
                           end_time,
                           treat_time,
                           dtw1_time,
                           dependent,
                           dependent_id,
                           t_treat = (treat_time - start_time) + 1,
                           n = (end_time - start_time) + 1,
                           n_dtw1 = (dtw1_time - start_time) + 1,
                           k = 6,
                           n_q = 1,
                           n_r = 1,
                           plot_figures = FALSE,
                           normalize_method = "t",
                           dtw_method = "dtw",
                           margin = 10,
                           step.pattern = dtw::symmetric2,
                           legend_position = c(0.3, 0.3), ...){
  
  
  # prepare data
  y_raw = data %>% 
    filter(unit == dependent &
             time <= end_time) %>%
    .[["value_raw"]]
  
  y_processed = data %>% 
    filter(unit == dependent &
             time <= end_time) %>%
    .[["value"]]
  
  x_list = data %>% 
    filter(unit != dependent &
             time <= end_time) %>% 
    select(c("unit", "time", "value", "value_raw")) %>% 
    group_by(unit) %>% 
    group_split(.keep = TRUE)
  
  # TSDTW
  results = NULL
  for (z in 1:length(x_list)) {
    item = x_list[[z]]
    unit = item$unit[1]
    x_processed = item$value
    x_raw = item$value_raw
    
    x = x_processed
    y = y_processed
    
    res = TwoStepDTW(x, y, t_treat, k, n_dtw1, dtw_method = dtw_method,
                     normalize_method = normalize_method,
                     step.pattern = step.pattern, ...)
    
    x_warped = c(warp_ts(res$W_a, x_raw[1:res$cutoff]),
                 warp_using_weight(x_raw[-(1:(res$cutoff - 1))],
                                   res$avg_weight)[-1])
    
    
    df = data.frame(time = 1:length(x_raw),
                    y = y_raw[1:length(x_raw)],
                    x = x_raw[1:length(x_raw)],
                    warped = x_warped[1:length(x_raw)]) %>%
      `colnames<-`(c("time", dependent, unit, paste0(unit, "-Warped"))) %>%
      reshape2::melt(., id.vars = "time") %>%
      `colnames<-`(c("time", "unit", "value"))
    
    fig = df %>%
      ggplot(aes(x = time, y = value, color = unit)) +
      geom_line() +
      scale_color_manual(values = c("#2a4d69", "#ee4035", "#7bc043")) +
      geom_vline(xintercept = t_treat, linetype="dashed",
                 color = "grey30", size = 0.3) +
      theme_bw() +
      ggtitle(unit) +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = legend_position)

    results[[unit]] = list(unit = unit,
                           x = x_raw,
                           res = res,
                           x_warped = x_warped,
                           df = df,
                           fig = fig)
  }
  
  # plot warped data
  if (plot_figures) {
    plot_warped(lapply(results,"[[","fig"), dependent, k, ncol = 3)
  }
  
  # prepare data for synthetic control
  df = results %>% 
    future_map(
      ~{
        item = .
        unit = item$unit
        x_warped = item$x_warped
        
        res = data.frame(
          time = 1:length(y_raw) + start_time - 1,
          unit = unit,
          value_warped = NA
        )
        res$value_warped = x_warped[1:((length(y_raw)-1) + 1)]  
        
        res
      }
    ) %>% 
    do.call("rbind", .) %>% 
    `row.names<-`(NULL)
  
  df = rbind(df,
             data.frame(time = 1:length(y_raw) + start_time - 1,
                        unit = dependent,
                        value_warped = y_raw))
  
  df = right_join(data, df, by = c("unit", "time"))
  df = data.frame(df)
  
  # w/o TSDTW
  synth_origin = do_synth_tobacco_87(df, "value_raw", dependent_id, start_time, n)
  if (plot_figures) {
    plot_synth_tobacco(synth_origin, "without_TSDTW", dependent, treat_time, k,
                       start_time, end_time)
  }
  
  # w/ TSDTW
  synth_new = do_synth_tobacco_87(df, "value_warped", dependent_id, start_time, n)
  if (plot_figures) {
    plot_synth_tobacco(synth_new, "TSDTW", dependent, treat_time, k,
                       start_time, end_time)
  }

  # mse
  diff1 = (synth_origin$synthetic - synth_origin$value)
  diff2 = (synth_new$synthetic - synth_new$value)
  
  mse1_1 = mean(diff1[1:(t_treat - 1)]^2, na.rm = T)
  mse2_1 = mean(diff2[1:(t_treat - 1)]^2, na.rm = T)
  
  mse1_2 = mean(diff1[t_treat:((end_time - start_time - 3))]^2, na.rm = T)
  mse2_2 = mean(diff2[t_treat:((end_time - start_time - 3))]^2, na.rm = T)
  
  
  return(list(
    dtw_results = results,
    df = df,
    synth_origin = synth_origin,
    synth_new = synth_new,
    mse = data.frame(mse1_pre = mse1_1,
                     mse2_pre = mse2_1,
                     mse1_post = mse1_2,
                     mse2_post = mse2_2)
  ))
}


## Run -------------------------------------------------------------------------
# data = smoking
# start_time = 1970
# end_time = 2000
# treat_time = 1985
# dtw1_time = 1990
# dependent = "California"
# dependent_id = 3
# t_treat = (treat_time - start_time) + 1
# n = (end_time - start_time) + 1
# n_dtw1 = (dtw1_time - start_time) + 1
# k = 6
# n_q = 1
# n_r = 1
# normalize_method = "t"
# dtw_method = "dtw"
# margin = 10
# step.pattern = dtw::symmetricP2
# ... = NULL
# legend_position = c(0.3, 0.3)



units = smoking[c("id", "unit")] %>% distinct
k = 5
result = as.list(1:nrow(units)) %>% 
  future_map(
    ~{
      i = .
      dependent = units$unit[i]
      dependent_id = units$id[i]
      print(paste0(dependent, ":", i, "-", k, " start..."))
      res = compare_methods(data = smoking,
                            start_time = 1970,
                            end_time = 1997,
                            treat_time = 1987,
                            dtw1_time = 1993,
                            dependent = dependent,
                            dependent_id = dependent_id,
                            normalize_method = "t",
                            k = k,
                            plot_figures = TRUE,
                            step.pattern = dtw::symmetricP2)
      print(paste0(dependent, ":", i, "-", k, " start...Done."))
      res$mse %>% mutate(dependent = dependent, k = k)
    }
  )

result = result %>% 
  do.call("rbind", .) %>% 
  mutate(improve = mse1_post - mse2_post)
length(which(result$improve>0))/39
boxplot(result$improve, outline = FALSE)
abline(h = 0, lty = 5)

t.test(result$improve)



result$ratio = result$improve/result$mse1_post
boxplot(result$ratio, outline = FALSE)
abline(h = 0, lty = 5)

t.test(result$ratio)



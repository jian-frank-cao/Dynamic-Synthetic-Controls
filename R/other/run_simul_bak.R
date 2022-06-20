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


## Simulation Data -------------------------------------------------------------
data = read.csv("./data/simulData_v3.csv", stringsAsFactors = FALSE)
data = data %>% 
  mutate(id = case_when(unit == "A" ~ 1,
                        unit == "B" ~ 2,
                        TRUE ~ 3),
         value = Y,
         value_raw = value) %>% 
  select(c("id", "unit", "time", "value", "value_raw"))



## Without TSDTW ---------------------------------------------------------------
# res_synth_old = do_synth_simul(data, "value_raw", 1)
# 
# plot_synth(res_synth_old,
#            "A", 50, 1, 80,
#            "./figures/synth_simul_without_TSDTW.pdf")


## With TSDTW ------------------------------------------------------------------
# grid search
for (width in (1:9)*2+3) {
  for (k in 4:9) {
    for (dtw1_time in 70:70) {
      data = preprocessing(data, filter_width = width)
      
      res = SimDesign::quiet(compare_methods(data = data,
                            start_time = 1,
                            end_time = 80,
                            treat_time = 50,
                            dtw1_time = dtw1_time,
                            dependent = "A",
                            dependent_id = 1,
                            normalize_method = "t",
                            k = k,
                            synth_fun = "simulation",
                            filter_width = width,
                            plot_figures = FALSE,
                            step.pattern = dtw::symmetricP2))
      mse = mean((res$synth_new$synthetic - res$synth_new$value)[50:60]^2, rm.na = T)
      print(paste0("width: ", width, ", k: ", k, ", dtw1: ", dtw1_time, ", MSE: ", mse))
    }
  }
}






# width = 21
# k = 7

width = 19
k = 5


data = preprocessing(data, filter_width = width)

res = compare_methods(data = data,
                      start_time = 1,
                      end_time = 80,
                      treat_time = 50,
                      dtw1_time = 70,
                      dependent = "A",
                      dependent_id = 1,
                      normalize_method = "t",
                      k = k,
                      synth_fun = "simulation",
                      filter_width = width,
                      plot_figures = TRUE,
                      step.pattern = dtw::symmetricP2)

df = rbind(res$df %>%
             filter(time <= 80) %>%
             select(c("unit", "time", "value_raw")) %>% 
             `colnames<-`(c("unit", "time", "value")),
           data.frame(unit = "w/o TSDTW",
                      time = 1:80,
                      value = res$synth_origin$synthetic[1:80]),
           data.frame(unit = "w/ TSDTW",
                      time = 1:80,
                      value = res$synth_new$synthetic[1:80])
           )

ggplot(df, aes(x = time, y = value, color = unit)) +
  geom_line(aes(linetype = unit)) + 
  geom_vline(xintercept = 50, linetype="dashed") +
  theme_bw() +
  scale_linetype_manual(values=c("solid", "dashed", "twodash", "solid", "solid")) +
  scale_color_manual(values=c("#4a4e4d", "grey70", "grey80", "#fe4a49","#2ab7ca"))


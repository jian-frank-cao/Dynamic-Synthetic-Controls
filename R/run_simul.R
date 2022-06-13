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
data = read.csv("./data/simulData.csv", stringsAsFactors = FALSE)
data = data %>% 
  mutate(id = case_when(unit == "A" ~ 1,
                        unit == "B" ~ 2,
                        TRUE ~ 3),
         value = Y,
         value_raw = value) %>% 
  select(c("id", "unit", "time", "value", "value_raw"))



## Without TSDTW ---------------------------------------------------------------
res_synth_old = do_synth_simul(data, "value_raw", 1)

plot_synth(res_synth_old,
           "A", 60, 1, 150,
           "./figures/synth_simul_without_TSDTW.pdf")


## With TSDTW ------------------------------------------------------------------
width = 31
k = 6

data = preprocessing(data, filter_width = width)

res = compare_methods(data = data,
                      start_time = 1,
                      end_time = 150,
                      treat_time = 60,
                      dtw1_time = 80,
                      dependent = "A",
                      dependent_id = 1,
                      normalize_method = "t",
                      k = k,
                      synth_fun = "simulation",
                      filter_width = width,
                      plot_figures = TRUE,
                      step.pattern = dtw::symmetricP2)

df = rbind(res$df %>%
             filter(time <= 100) %>%
             select(c("unit", "time", "value_raw")) %>% 
             `colnames<-`(c("unit", "time", "value")),
           data.frame(unit = "w/o TSDTW",
                      time = 1:100,
                      value = res$synth_origin$synthetic[1:100]),
           data.frame(unit = "w/ TSDTW",
                      time = 1:100,
                      value = res$synth_new$synthetic[1:100])
           )

ggplot(df, aes(x = time, y = value, color = unit)) +
  geom_line(aes(linetype = unit)) + 
  geom_vline(xintercept = 60, linetype="dashed") +
  theme_bw() +
  scale_linetype_manual(values=c("solid", "dashed", "twodash", "solid", "solid")) +
  scale_color_manual(values=c("#4a4e4d", "grey70", "grey80", "#fe4a49","#2ab7ca"))

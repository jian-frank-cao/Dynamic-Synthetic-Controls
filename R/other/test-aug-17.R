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
source("./R/simulate.R")
set.seed(20220407)


data_list = readRDS("./data/simul_data_list_0818.Rds")
start_time = 1
end_time = 1000
treat_time = 800
dtw1_time = 900
n_mse = 100
k = 15
plot_figures = TRUE
step.pattern1 = dtw::symmetricP2
step.pattern2 = dtw::asymmetricP2
legend_position = c(0.3, 0.8)
... = NULL
normalize_method = "t"
n_q = 1
n_r = 1
default_margin = 3
ma = 3
ma_na = "original"

data = data_list[[26]]
dependent = "A"
dependent_id = 1
predictors.origin = NULL
special.predictors.origin = list(list("value_raw", 700:799, c("mean")))
time.predictors.prior.origin = 1:799
time.optimize.ssr.origin = 1:799
predictors.new = NULL
special.predictors.new = list(list("value_warped", 700:799, c("mean")))
time.predictors.prior.new = 1:799
time.optimize.ssr.new = 1:799

t_treat = (treat_time - start_time) + 1
n = (end_time - start_time) + 1
n_dtw1 = (dtw1_time - start_time) + 1



x_list = data %>% 
  filter(unit != dependent &
           time <= end_time) %>% 
  select(c("unit", "time", "value", "value_raw")) %>% 
  group_by(unit) %>% 
  group_split(.keep = TRUE)

y_raw = data %>% 
  filter(unit == dependent &
           time <= end_time) %>%
  .[["value_raw"]]
y_processed = data %>% 
  filter(unit == dependent &
           time <= end_time) %>%
  .[["value"]]
# y_raw = x_list[[4]]$value_raw
# y_processed = x_list[[4]]$value


item = x_list[[1]]
unit = as.character(item$unit[1])
x = item$value
x_raw = item$value_raw
y = y_processed

## Run -------------------------------------------------------------------------
res_1stDTW = first_dtw(x, y, k, n_dtw1, t_treat,
                       normalize_method,
                       step.pattern1, plot_figures, ...)
x_pre = res_1stDTW$x_pre
x_post = res_1stDTW$x_post
W_a = res_1stDTW$W_a
cutoff = res_1stDTW$cutoff

# compute weight a
weight_a_o = warping2weight(W_a)
weight_a = as.numeric(stats::filter(weight_a_o, rep(1/ma, ma)))
weight_a = zoo::na.locf(weight_a, na.rm = FALSE)
if (ma_na == "one") {
  weight_a[is.na(weight_a)] = 1
}else if(ma_na == "first-available") {
  weight_a[is.na(weight_a)] = weight_a[!is.na(weight_a)][1]
}else if (ma_na == "original") {
  weight_a[is.na(weight_a)] = weight_a_o[is.na(weight_a)]
}

# 2nd dtw
res_2ndDTW = second_dtw(x_post, x_pre, 
                        weight_a, k, normalize_method,
                        dist_quantile = dist_quantile,
                        n_q, n_r, step.pattern = step.pattern2, ...)
# avg_weight = res_2ndDTW$avg_weight[-(1:(k - 3))]
avg_weight = res_2ndDTW$avg_weight

w_2fdtw = cumsum(avg_weight)
w_dtw = dtw::warp(dtw::dtw(x_post, y[t_treat:1000], step.pattern = step.pattern2), index.reference = FALSE)

plot(ts(w_dtw))
lines(w_2fdtw, 1:length(w_2fdtw))

x_warp1 = c(warp_ts(W_a, x_raw[1:cutoff]), x_post[w_dtw][-1])
x_warp2 = c(warp_using_weight(x_pre, weight_a)[1:800], warp_using_weight(x_raw[-(1:(cutoff - 1))],
                                                             avg_weight)[-1])

plot(ts(y_raw))
lines(x_raw, col = "blue")
lines(x_warp1, col = "green")
lines(x_warp2, col = "red")

## Test ------------------------------------------------------------------------
res_1stDTW = first_dtw(x, y, k, n_dtw1, t_treat,
                       normalize_method, 
                       step.pattern1, plot_figures, ...)
x_pre = res_1stDTW$x_pre
x_post = res_1stDTW$x_post
W_a = res_1stDTW$W_a
cutoff = res_1stDTW$cutoff

# compute weight a
weight_a_o = warping2weight(W_a)
weight_a = as.numeric(stats::filter(weight_a_o, rep(1/ma, ma)))
weight_a = zoo::na.locf(weight_a, na.rm = FALSE)
if (ma_na == "one") {
  weight_a[is.na(weight_a)] = 1
}else if(ma_na == "first-available") {
  weight_a[is.na(weight_a)] = weight_a[!is.na(weight_a)][1]
}else if (ma_na == "original") {
  weight_a[is.na(weight_a)] = weight_a_o[is.na(weight_a)]
}

n_pre = length(x_pre)
n_post = length(x_post)
# slide target window
i = 1
weight = NULL


Q = x_post[i:(i + k - 1)]
Q = normalize(Q, normalize_method)
costs_qr = NULL

# slide reference window
continue = TRUE
margin = default_margin
j = 1
while (continue) {  # j <= n_pre - k + 1
  # check if the search is finished
  if (j > n_pre - k + 3) {
    continue = FALSE
    next
  }
  # define R
  R = x_pre[j:min(j + k + margin - 1, n_pre)]
  R = normalize(R, normalize_method)
  # check if R is too short
  R_too_short = ref_too_short(Q, R, step.pattern = step.pattern2)
  if (R_too_short) {
    # check if the R can be extended
    if (j < n_pre - k - margin + 1) {
      margin = margin + 1
      next
    }else{
      margin = default_margin
      j = j + n_r
      next
    }
  }
  # match Q and R
  alignment_qr = dtw::dtw(Q, R, open.end = TRUE,
                          step.pattern = step.pattern2,
                          distance.only = TRUE)
  costs_qr = rbind(costs_qr,
                   data.frame(cost = alignment_qr$distance,
                              j = j,
                              margin = margin))
  j = j + n_r
  margin = default_margin
}
# find the minimum cost
min_cost = which(costs_qr$cost == min(costs_qr$cost))[1]
j_opt = costs_qr$j[min_cost]
margin_opt = costs_qr$margin[min_cost]

# obtain warping path W_pp_i: x_post -> n_pre
Rs = x_pre[j_opt:min(j_opt + k + margin_opt - 1, n_pre)]
Rs = normalize(Rs, normalize_method)
alignment_qrs = dtw::dtw(Q, Rs, open.end = TRUE,
                         step.pattern = step.pattern2, ...)
W_pp_i = Matrix::sparseMatrix(alignment_qrs$index1,
                              alignment_qrs$index2)

# obtain weight_b
weight_a_Rs = weight_a[j_opt:(j_opt + ncol(W_pp_i) - 1)]
weight_b = as.numeric((W_pp_i %*% weight_a_Rs)/rowSums(as.matrix(W_pp_i)))

# convert warping path to weight
weight_i = matrix(rep(NaN, n_post), nrow = 1)
weight_i[1, i:(i + k - 1)] = weight_b

# stack weight
weight = rbind(weight, weight_i)

# next
i = i + n_q

plot(ts(Q))
lines(Rs, col = "red")



## -------------------------
i = 1
treatment = c(rep(0, 800),
              seq(0, 50, length.out = 100),
              rep(50, 100))
item = result[[i]]
gap_original = item$value_raw - item$synth_original
gap_new = item$value_raw - item$synth_new
plot(ts(gap_original), col = "blue")
lines(treatment, col = "black")
lines(gap_new, col = "red")
i = i+1

data_list[[26]] %>% 
  ggplot(aes(x = time, y = value, color = unit)) +
  geom_line()


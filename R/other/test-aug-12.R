data = data_list[[211]]
start_time = 1
end_time = 1000
treat_time = 800
dtw1_time = 900
dependent = "A"
dependent_id = 1
normalize_method = "t"
k = k
synth_fun = "simulation"
filter_width = width
plot_figures = T
step.pattern = dtw::typeIIc
... = NULL
dtw_method = "dtw"


plot(ts(y_raw))
lines(x_raw, col = "blue")
lines(x_warped, col = "red")

alignment = dtw::dtw(y, x, keep = TRUE,
                     step.pattern = step.pattern,
                     open.end = TRUE, ...)
dtw::dtwPlotTwoWay(res_1stDTW$alignment, xts = y_raw[1:800], yts = x_raw[1:900] + 200)
dtw::dtwPlotTwoWay(alignment, xts = y, yts = x + 10)

i = 1
Q = x_post[i:(i + k - 1)]
Q = normalize(Q, normalize_method)

# partial match Q -> x_pre
x_pre_normalized = normalize(x_pre, normalize_method = normalize_method)
alignment_qx = dtw::dtw(Q, x_pre_normalized,
                        open.begin = TRUE,
                        open.end = TRUE,
                        keep = TRUE,
                        step.pattern = step.pattern, ...)

# obtain warping path W_pp_i: x_post -> x_pre
begin_x_pre = min(alignment_qx$index2)
end_x_pre = max(alignment_qx$index2)
Q = x_post[i:(i + k - 1)]
plot(ts(x_pre))
lines(begin_x_pre:(begin_x_pre + k - 1), Q, col = "red")
i = i + 1


i = 1
n_r = 1
n_q = 1
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
  R_too_short = ref_too_short(Q, R, step.pattern = step.pattern)
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
                          step.pattern = step.pattern,
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
Q = x_post[i:(i + k - 1)]
plot(ts(x_pre))
lines(j_opt:(j_opt + k - 1), Q, col = "red")
i = i + 1

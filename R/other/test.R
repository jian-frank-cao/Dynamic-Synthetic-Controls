last_year = 2003
dependent = "West Germany"
dependent_id = 7
stretch = 1
t_treat = (21 - 1) * stretch + 1
n = (44 - 1) * stretch + 1
n_1st = (30 - 1) * stretch + 1
k = 6 * stretch
n_q = 1
n_r = 1
normalize_method = "minmax"
dtw_method = "dtw"
margin = 10
... = NULL
step.pattern = dtw::symmetricP1

# last_year = 2003
# dependent = "Germany"
# dependent_id = 7
# stretch = 1
# t_treat = (81 - 1) * stretch + 1
# n = (176 - 1) * stretch + 1
# n_1st = (120 - 1) * stretch + 1
# k = 24 * stretch
# n_q = 1
# n_r = 1
# normalize_method = "minmax"
# dtw_method = "dtw"
# margin = 10
# ... = NULL
# step.pattern = dtw::symmetric2

y_raw = data %>% 
  filter(country == dependent & year <= last_year) %>%
  .[["gdp"]]

x_list = data %>% 
  filter(country != dependent & year <= last_year) %>% 
  select(c("country", "gdp")) %>% 
  group_by(country) %>% 
  group_split(.keep = TRUE)

item = x_list[[8]]
country = item[["country"]][1]
x = item[["gdp"]]

y = approx(y_raw, n = (length(x) - 1) * stretch + 1, method = "linear")$y
x = approx(x, n = (length(x) - 1) * stretch + 1, method = "linear")$y

############################################# 1st
# normalize
y_bak = y
x_bak = x
y = normalize(y[1:n_1st], normalize_method)
x = normalize(x[1:n_1st], normalize_method)

# ddtw
if (dtw_method == "ddtw") {
  y = derivatives(y)
  x = derivatives(x)
}

# plot(ts(y+pattern))
# lines(x+pattern, col = "red")
# 
# idx<-seq(0,6.28,len=25);
# pattern<-sin(idx)+runif(25)/10;


idx2 = seq(0, 12.56*2, len = 25)
pattern2 = cos(idx2) + runif(25)/10

pattern3 = pattern + 0.5*pattern2
plot(ts(pattern3))
lines(pattern, col = "red")
lines(pattern2, col = "blue")

A = sin(idx2)
B = cos(idx2)
alpha = 0.5
beta = 0.7
pattern = alpha*A + beta*B
plot(ts(pattern))
lines(alpha*A, col = "red")
lines(beta*B, col = "blue")

# dtw
alignment = dtw::dtw(x, y, keep = TRUE, step.pattern = step.pattern, ...)
fig_ThreeWay = dtw::dtwPlotThreeWay(alignment)
wq = suppressWarnings(dtw::warp(alignment, index.reference = FALSE))
wt = suppressWarnings(dtw::warp(alignment, index.reference = TRUE))
W = Matrix::sparseMatrix(alignment$index1, alignment$index2)
cutoff = round(wq[t_treat])

# partition warping path W
W_a = W[1:cutoff, 1:t_treat]

# cut x
x_pre = x_bak[1:cutoff]
x_post = x_bak[-(1:(cutoff - k + 2))]

################################################## 2nd

n_pre = length(x_pre)
n_post = length(x_post)

# slide target window
i = 1
weight = NULL
while (i <= n_post - k + 1) {
  Q = x_post[i:(i + k - 1)]
  Q = normalize(Q, normalize_method)
  j = 1
  costs_qr = c()
  
  # slide reference window
  while (j <= n_pre - k + 1) {
    R = x_pre[j:(j + k - 1)]
    R = normalize(R, normalize_method)
    alignment_qr = dtw::dtw(Q, R, open.end = FALSE,
                            step.pattern = step.pattern,
                            distance.only = TRUE)
    costs_qr = c(costs_qr, alignment_qr$distance)
    j = j + n_r
  }
  j_opt = which(costs_qr == min(costs_qr))
  
  # adjust for optimal R
  p = 0
  kp = k + p
  # p = p_min
  # costs_qrp = c()
  # while (p <= p_max) {
  #   Rp = x_pre[j_opt:(j_opt + k + p)]
  #   alignment_qrp = dtw::dtw(Q, Rp, distance.only = TRUE)
  #   costs_qrp = c(costs_qrp, alignment_qrp$distance)
  #   p = p + 1
  # }
  # p_opt = which(costs_qrp == min(costs_qrp))
  
  # obtain warping path W_pp_i: x_post -> n_pre
  Rs = x_pre[j_opt:(j_opt + kp - 1)]
  Rs = normalize(Rs, normalize_method)
  alignment_qrs = dtw::dtw(Q, Rs, step.pattern = step.pattern, ...)
  W_pp_i = Matrix::sparseMatrix(alignment_qrs$index1,
                                alignment_qrs$index2)
  
  # obtain warping path W_b_i: x_post -> y
  W_a_Rs = W_a[j_opt:(j_opt + kp - 1),]
  col_sums = colSums(as.matrix(W_a_Rs))
  ind_nonzero = which(col_sums > 0)
  n_ind = length(ind_nonzero)
  min_ind = min(ind_nonzero)
  max_ind = max(ind_nonzero)
  ind_left = min_ind - j_opt
  ind_right = max_ind - (j_opt + kp - 1)
  W_a_Rs = W_a_Rs[, min_ind:max_ind]
  W_b_i = W_pp_i %*% W_a_Rs
  
  # convert warping path to weight
  weight_i = matrix(rep(NaN, n_post), nrow = 1)
  weight_i[1, i:(i + k - 1)] = warping2weight(W_b_i)
  
  # stack weight
  weight = rbind(weight, weight_i)
  
  # next
  i = i + n_q
}

# average weight
avg_weight = colMeans(weight, na.rm = TRUE)
######################################################## 
x_w = c(warp_ts(W_a, x_bak[1:cutoff]),
        warp_using_weight(x_bak[-(1:(cutoff - 1))],
                          avg_weight[-(1:(k - 3))])[-1])

plot(ts(x_bak))
lines(x_w, col = "red")

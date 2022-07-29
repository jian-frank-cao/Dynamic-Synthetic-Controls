

Q = x_post[i:(i + k - 1)]
Q = normalize(Q, normalize_method)

# partial match Q -> x_pre
x_pre_normalized = normalize(x_pre, normalize_method = normalize_method)
alignment_qx = dtw::dtw(Q, x_pre_normalized,
                        open.begin = TRUE,
                        open.end = TRUE,
                        keep = TRUE,
                        step.pattern = dtw::rigid, ...)

# obtain warping path W_pp_i: x_post -> x_pre
begin_x_pre = min(alignment_qx$index2)
end_x_pre = max(alignment_qx$index2)

plot(ts(x_pre_normalized))

lines(begin_x_pre:(begin_x_pre + length(Q) - 1), Q, col = "red")

i = i + 1





Q = x_post[i:(i + k - 1)]
Q = normalize(Q, normalize_method)
j = 1
costs_qr = c()

# slide reference window
while (j <= n_pre - k + 1) {
  R = x_pre[j:(j + k - 1)]
  R = normalize(R, normalize_method)
  alignment_qr = dtw::dtw(Q, R, open.end = TRUE,
                          step.pattern = dtw::rigid,
                          distance.only = TRUE)
  costs_qr = c(costs_qr, alignment_qr$distance)
  j = j + n_r
}
j_opt = which(costs_qr == min(costs_qr))
plot(ts(x_pre_normalized))
lines(j_opt:(j_opt + length(Q) - 1), Q, col = "red")
i = i + 1
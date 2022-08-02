i = 1

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


Rs = x_pre[j_opt:(j_opt + kp - 1)]
Rs = normalize(Rs, normalize_method)
alignment_qrs = dtw::dtw(Q, Rs, step.pattern = dtw::rigid, ...)
W_pp_i = Matrix::sparseMatrix(alignment_qrs$index1,
                              alignment_qrs$index2)


plot(ts(x_pre))
lines(j_opt:(j_opt + k - 1), Rs, col = "red")

i = i + 1


plot(ts(x_pre))
lines(16:55, x_post, col = "red")

## functions -------------------------------------------------------------------
# normalization
t_normalize = function(data){
  return((data - mean(data))/sd(data))
}

minmax_normalize = function(data){
  return((data - min(data))/(max(data) - min(data)))
}

normalize = function(data, normalize_method){
  if (normalize_method == "minmax") {
    data = minmax_normalize(data)
  }else if (normalize_method == "t") {
    data = t_normalize(data)
  }
  return(data)
}


# compute derivatives
derivatives = function(data){
  data = DTWBI::local.derivative.ddtw(data)
  data[1] = data[2]
  data[length(data)] = data[length(data) - 1]
  return(data)
}


# 1st dtw
first_dtw = function(x, y, k, n_dtw1, t_treat,
                     normalize_method = "t",
                     dtw_method = "dtw", 
                     step.pattern = dtw::symmetricP2,
                     plot_figures = FALSE, ...){
  # normalize
  y_bak = y
  x_bak = x
  y = normalize(y[1:n_dtw1], normalize_method)
  x = normalize(x[1:n_dtw1], normalize_method)
  
  # ddtw
  if (dtw_method == "ddtw") {
    y = derivatives(y)
    x = derivatives(x)
  }
  
  # dtw
  alignment = dtw::dtw(x, y, keep = TRUE, step.pattern = step.pattern, ...)
  if (plot_figures) {
    fig_ThreeWay = dtw::dtwPlotThreeWay(alignment)
  }
  wq = suppressWarnings(dtw::warp(alignment, index.reference = FALSE))
  W = Matrix::sparseMatrix(alignment$index1, alignment$index2)
  cutoff = round(wq[t_treat])
  
  # partition warping path W
  W_a = W[1:cutoff, 1:t_treat]
  
  # cut x
  x_pre = x_bak[1:cutoff]
  x_post = x_bak[-(1:(cutoff - k + 2))]
  
  return(list(x = x_bak,
              y = y_bak,
              k = k,
              n = n,
              t_treat = t_treat,
              alignment = alignment,
              wq = wq,
              W_a = W_a,
              cutoff = cutoff,
              x_pre = x_pre,
              x_post = x_post))
}


# find maximum path
maximum_path_naive = function(W_b_bar){
  index1 = 1
  index2 = NULL
  finished = FALSE
  indices1 = c()
  indices2 = c()
  while (!finished) {
    if (is.null(index2)) {
      index2 = which(W_b_bar[1,] == max(W_b_bar[1,]))
    }else if (index1 == n_post) {
      now = W_b_bar[index1, index2]
      right = W_b_bar[index1, index2 + 1]
      if (right >= now) {
        index2 = index2 + 1
      }else{
        finished = TRUE
      }
    }else{
      right = W_b_bar[index1, index2 + 1]
      down = W_b_bar[index1 + 1, index2]
      right_down = W_b_bar[index1 + 1, index2 + 1]
      if (right == max(c(right, down, right_down))) {
        index2 = index2 + 1
      }else if (down == max(c(right, down, right_down))) {
        index1 = index1 + 1
      }else{
        index1 = index1 + 1
        index2 = index2 + 1
      }
    }
    indices1 = c(indices1, index1)
    indices2 = c(indices2, index2)
  }
  indices2 = indices2 - min(indices2) + 1
  
  return(list(indices1 = indices1,
              indices2 = indices2))
}

maximum_path = function(W, margin){
  n_row = nrow(W)
  W = W[, (margin + 1):(n_row + margin)]
  # col_sums = colSums(as.matrix(W))
  # ind_nonzero = which(col_sums > 0)
  # W = W[, min(ind_nonzero):max(ind_nonzero)]
  n_col = ncol(W)
  W[1, 1] = 1
  W[n_row, n_col] = 1
  
  # compute accumulative weight
  weight = as.matrix(W) * 0
  for (i in 1:n_row) {
    for (j in 1:n_col) {
      if (i == 1 & j == 1) {
        weight[i, j] = W[i, j]
        next
      }
      if (i == 1) {
        weight[i, j] = weight[i, j - 1] + W[i, j]
        next
      }
      if (j == 1) {
        weight[i, j] = weight[i - 1, j] + W[i, j]
        next
      }
      weight[i, j] = max(weight[i - 1, j - 1],
                         weight[i - 1, j],
                         weight[i, j - 1]) + W[i, j]
    }
  }
  
  # find the maximum path
  index1 = n_row
  index2 = n_col
  finished = FALSE
  indices1 = c(index1)
  indices2 = c(index2)
  while (!finished) {
    if (index1 == 1 & index2 == 1) {
      finished = TRUE
    }else if (index1 == 1) {
      index2 = index2 - 1
    }else if (index2 == 1) {
      index1 = index1 - 1
    }else{
      left = weight[index1, index2 - 1]
      up = weight[index1 - 1, index2]
      up_left = weight[index1 -1, index2 - 1]
      maxx = max(c(left, up, up_left))
      if (up_left == maxx) {
        index1 = index1 - 1
        index2 = index2 - 1
      }else if (left > up) {
        index2 = index2 - 1
      }else if (up > left) {
        index1 = index1 - 1
      }else if (up == left) {
        if (index1 > 2 & index2 > 2) {
          up_up = weight[index1 - 2, index2]
          left_left = weight[index1, index2 - 2]
          if (up_up >= left_left) {
            index1 = index1 - 1
          }else{
            index2 = index2 - 1
          }
        }else{
          index1 = index1 - 1
        }
      }
    }

    indices1 = c(index1, indices1)
    indices2 = c(index2, indices2)
  }
  
  return(list(indices1 = indices1,
              indices2 = indices2))
}


# warping path W to weight
warping2weight = function(W){
  w = as.matrix(W)
  count = rep(1/colSums(w), nrow(w)) %>% 
    matrix(.,
           nrow = ncol(w),
           ncol = nrow(w)) %>% 
    t(.)
  weight = rowSums(w * count)
  
  return(weight)
}

# use weight to warp target time series
warp_using_weight = function(ts, weight){
  n_ts = length(ts)
  n_weight = length(weight)
  if (n_ts != n_weight) {
    print("Lengths of ts and weight are different.")
    return()
  }
  
  # form tranformation matrix
  count = 1
  residual = 0
  w = matrix(0, nrow = n_weight, ncol = n_weight + 100)
  for (i in 1:n_weight) {
    z = weight[i]
    while (z > 0) {
      if (z + residual >= 1) {
        now = 1 - residual
        residual = 0
        z = z - now
        w[i, count] = now
        count = count + 1
      }else{
        residual = residual + z
        now = z
        z = 0
        w[i, count] = now
      }
    }
  }
  
  col_sum = colSums(w)
  ind_max = which(col_sum > 0) %>% 
    max
  w = w[, 1:ind_max]
  w[n_weight, ind_max] = 1 - sum(w[-n_weight, ind_max])
  
  # warp time series
  ts_warped = t(w) %*% matrix(ts, ncol = 1)
  
  return(ts_warped[,1])
}


# 2nd dtw
second_dtw = function(x_post, x_pre,
                      W_a, k, normalize_method = "t",
                      n_q = 1, n_r = 1,
                      margin = 10,
                      p_min = -5, p_max = 5,
                      step.pattern = dtw::symmetricP2, ...){
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
      alignment_qr = dtw::dtw(Q, R, open.end = TRUE,
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
    # n_ind = length(ind_nonzero)
    min_ind = min(ind_nonzero)
    max_ind = max(ind_nonzero)
    # ind_left = min_ind - j_opt
    # ind_right = max_ind - (j_opt + kp - 1)
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
  
  return(list(weight = weight,
              avg_weight = avg_weight))
}


# use W to get warp_ind
warp_ind = function(W){
  warp_ind = (W * matrix(rep(1:nrow(W), ncol(W)),
                         nrow = nrow(W))) %>% 
    as.matrix
  warp_ind[warp_ind == 0] = NaN
  warp_ind = colMeans(warp_ind, na.rm = TRUE)
  
  return(warp_ind)
}


# use W to warp ts
warp_ts = function(w, ts){
  w = as.matrix(w)
  ts_warped = (ts %*% w)/colSums(w)
  
  return(ts_warped)
}


# Two Step DTW
TwoStepDTW = function(x, y, t_treat, k, n_dtw1,
                      normalize_method = "t",
                      dtw_method = "dtw",
                      n_q = 1, n_r = 1, margin = 10,
                      p_min = -5, p_max = 5, 
                      step.pattern = dtw::symmetricP2, 
                      plot_figures = FALSE, ...){
  # 1st dtw
  res_1stDTW = first_dtw(x, y, k, n_dtw1, t_treat,
                         normalize_method, dtw_method,
                         step.pattern, plot_figures, ...)
  x_pre = res_1stDTW$x_pre
  x_post = res_1stDTW$x_post
  W_a = res_1stDTW$W_a
  cutoff = res_1stDTW$cutoff
  
  # 2nd dtw
  res_2ndDTW = second_dtw(x_post, x_pre, 
                          W_a, k, normalize_method,
                          n_q, n_r, margin,
                          step.pattern = step.pattern, ...)
  avg_weight = res_2ndDTW$avg_weight[-(1:(k - 3))]
  
  # warp x
  # warp_ind_a = warp_ind(W_a)
  # warp_ind_bs = warp_ind(W_bs)
  # x_w = c(x_original[1:cutoff][warp_ind_a],
  #         x_original[-(1:(cutoff - 1))][warp_ind_bs[-1]])
  # x_w = c(warp_ts(W_a, x_original[1:cutoff]),
  #         warp_using_weight(x_original[-(1:(cutoff - 1))],
  #                           avg_weight[-(1:(k - 3))])[-1])
  
  return(list(y = y,
              x = x,
              W_a = W_a,
              weight = res_2ndDTW$weight,
              avg_weight = avg_weight,
              t_treat = t_treat,
              cutoff = cutoff))
}

plot_warped = function(fig_list, ncol, file_name){
  nrow = ceiling(length(fig_list)/ncol)
  fig = gridExtra::marrangeGrob(fig_list, ncol = ncol,
                                nrow = nrow)
  ggsave(file_name,
         fig, width = ncol*4, height = nrow*4,
         units = "in", limitsize = FALSE)
}

# paste0("./figures/warped_",
#        paste0(c(dependent, k), collapse = "_"),
#        ".pdf")




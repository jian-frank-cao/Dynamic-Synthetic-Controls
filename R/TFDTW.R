## Functions -------------------------------------------------------------------
# 1st dtw
first_dtw = function(x, y, k, n_dtw1, t_treat,
                     normalize_method = "t",
                     step.pattern = dtw::symmetricP2,
                     plot_figures = FALSE, ...){
  # backup
  y_bak = y
  x_bak = x
  
  # normalize
  y = normalize(y_bak[1:t_treat], normalize_method)
  x = normalize(x_bak[1:n_dtw1], normalize_method, x_bak[1:t_treat])
  
  # check if x is too short
  x_too_short = ref_too_short(y, x, step.pattern = step.pattern)
  while(x_too_short & n_dtw1 < length(x_bak)){
    n_dtw1 = n_dtw1 + 1
    x = normalize(x_bak[1:n_dtw1], normalize_method, x_bak[1:t_treat])
    x_too_short = ref_too_short(y, x, step.pattern = step.pattern)
  }
  
  # dtw
  alignment = dtw::dtw(y, x, keep = TRUE,
                       step.pattern = step.pattern,
                       open.end = TRUE, ...)
  if (plot_figures) {
    fig_ThreeWay = dtw::dtwPlotThreeWay(alignment)
  }
  wr = suppressWarnings(dtw::warp(alignment, index.reference = TRUE))
  W = Matrix::sparseMatrix(alignment$index2, alignment$index1)
  cutoff = round(wr[t_treat])
  
  # partition warping path W
  W_a = W[1:cutoff, 1:t_treat]
  
  # cut x
  x_pre = x_bak[1:cutoff]
  # x_post = x_bak[-(1:(cutoff - k + 2))]
  x_post = x_bak[-(1:(cutoff - 1))]
  
  return(list(x = x_bak,
              y = y_bak,
              k = k,
              n = n,
              t_treat = t_treat,
              alignment = alignment,
              wr = wr,
              W_a = W_a,
              cutoff = cutoff,
              x_pre = x_pre,
              x_post = x_post))
}


# 2nd dtw
second_dtw = function(x_post, x_pre,
                      weight_a, k, normalize_method = "t",
                      dist_quantile = 1,
                      n_q = 1, n_r = 1,
                      default_margin = 3,
                      step.pattern = dtw::asymmetricP2, ...){
  n_pre = length(x_pre)
  n_post = length(x_post)
 
  # slide target window
  i = 1
  weight = NULL
  distance = NULL
  while (i <= n_post - k + 1) {
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
    margin_opt = costs_qr$margin[min_cost]
    
    # obtain warping path W_pp_i: x_post -> n_pre
    Rs = x_pre[j_opt:min(j_opt + k + margin_opt - 1, n_pre)]
    Rs = normalize(Rs, normalize_method)
    alignment_qrs = dtw::dtw(Q, Rs, open.end = TRUE,
                             step.pattern = step.pattern, ...)
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
    
    # distance
    distance = c(distance, alignment_qrs$distance)
    
    # next
    i = i + n_q
  }
  
  # handle misfits
  misfits = which(distance > quantile(distance, dist_quantile))
  weight = weight[-misfits,]
  
  # handle outliers
  
  
  # average weight
  avg_weight = colMeans(weight, na.rm = TRUE)
  avg_weight[is.na(avg_weight)] = 1
  
  return(list(weight = weight,
              avg_weight = avg_weight))
}


# Two Step DTW
TwoStepDTW = function(x, y, t_treat, k, n_dtw1,
                      normalize_method = "t",
                      ma = 3, ma_na = "original",
                      dist_quantile = 1,
                      n_q = 1, n_r = 1, 
                      step.pattern1 = dtw::symmetricP2,
                      step.pattern2 = dtw::asymmetricP2,
                      plot_figures = FALSE, ...){
  # 1st dtw
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
  
  return(list(y = y,
              x = x,
              W_a = W_a,
              weight_a = weight_a,
              weight = res_2ndDTW$weight,
              avg_weight = avg_weight,
              t_treat = t_treat,
              cutoff = cutoff))
}

## Functions -------------------------------------------------------------------
# normalization
t_normalize = function(data, reference = NULL){
  if (is.null(reference)) {
    mu = mean(data)
    sigma = sd(data)
  }else{
    mu = mean(reference)
    sigma = sd(reference)
  }
  res = (data - mu)/sigma
  return(res)
}

minmax_normalize = function(data, reference = NULL){
  if (is.null(reference)) {
    minimum = min(data)
    maximum = max(data)
  }else{
    minimum = min(reference)
    maximum = max(reference)
  }
  res = (data - minimum)/(maximum - minimum)
  return(res)
}

normalize = function(data, normalize_method, reference = NULL){
  if (normalize_method == "minmax") {
    data = minmax_normalize(data, reference)
  }else if (normalize_method == "t") {
    data = t_normalize(data, reference)
  }
  return(data)
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
  w = matrix(0, nrow = n_weight, ncol = 2*n_weight)
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


# use W to warp ts
warp_ts = function(w, ts){
  w = as.matrix(w)
  ts_warped = (ts %*% w)/colSums(w)
  
  return(ts_warped)
}


# check if reference is too short in dtw
ref_too_short = function(query, reference, 
                         step.pattern = dtw::symmetricP2){
  alignment = tryCatch(dtw::dtw(reference, query,
                                step.pattern = step.pattern,
                                open.end = TRUE),
                       error = function(e) return(NULL))
  if (is.null(alignment)) {
    return(FALSE)
  }
  if (length(unique(alignment$index2)) == 1) {
    return(TRUE)
  }
  wq = suppressWarnings(dtw::warp(alignment, index.reference = FALSE))
  return(round(max(wq)) < length(query))
}


# plot warped
plot_warped = function(fig_list, ncol, file_name){
  nrow = ceiling(length(fig_list)/ncol)
  fig = gridExtra::marrangeGrob(fig_list, ncol = ncol,
                                nrow = nrow)
  ggsave(file_name,
         fig, width = ncol*4, height = nrow*4,
         units = "in", limitsize = FALSE)
}


# add buffer
add_buffer = function(TS, n){
  model_right = forecast::auto.arima(TS)
  right <- as.numeric(forecast::forecast(model_right, h = n)$mean)
  model_left = forecast::auto.arima(rev(TS))
  left <- rev(as.numeric(forecast::forecast(model_left, h = n)$mean))
  
  return(c(left, TS, right))
}


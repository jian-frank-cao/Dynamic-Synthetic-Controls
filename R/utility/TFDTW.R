## Functions -------------------------------------------------------------------
# 1st dtw
first.dtw = function(x, y, t.treat, buffer = 10, 
                     norm.method = "t", match.method = "fixed",
                     step.pattern = dtw::symmetricP2,
                     window.type = "none",
                     window.size = NULL,
                     plot.figures = FALSE, ...){
  # backup
  y.bak = y
  x.bak = x
  
  # normalize
  y = normalize(y.bak[1:t.treat], norm.method)
  x = normalize(x.bak[1:(t.treat + buffer)],
                norm.method, x.bak[1:t.treat])
  
  if (match.method == "fixed") {
    alignment = dtw::dtw(y, x, keep = TRUE,
                         step.pattern = step.pattern,
                         window.type = window.type,
                         window.size = window.size,
                         open.end = FALSE, ...)
  }else if(match.method == "open.end"){
    # check if x is too short
    x.too.short = RefTooShort(y, x, step.pattern = step.pattern)
    while(x.too.short & (t.treat + buffer) < length(x.bak)){
      buffer = buffer + 1
      x = normalize(x.bak[1:(t.treat + buffer)],
                    norm.method, x.bak[1:t.treat])
      x.too.short = RefTooShort(y, x, step.pattern = step.pattern)
    }
    
    # dtw
    alignment = dtw::dtw(y, x, keep = TRUE,
                         step.pattern = step.pattern,
                         window.type = window.type,
                         window.size = window.size,
                         open.end = TRUE, ...)
  }
  
  if (plot.figures) {
    dtw::dtwPlotThreeWay(alignment)
  }
  wr = suppressWarnings(dtw::warp(alignment, index.reference = TRUE))
  W = Matrix::sparseMatrix(alignment$index2, alignment$index1)
  cutoff = round(wr[t.treat])
  
  # partition warping path W
  W.a = W[1:cutoff, 1:t.treat]
  
  return(list(x = x.bak, y = y.bak,
              t.treat = t.treat,
              buffer = buffer,
              step.pattern = step.pattern,
              alignment = alignment,
              wr = wr, W.a = W.a,
              cutoff = cutoff))
}


# 2nd dtw
second.dtw = function(x.post, x.pre, k, weight.a, 
                      norm.method = "t",
                      default.margin = 3,
                      n.q = 1, n.r = 1,
                      step.pattern = dtw::asymmetricP2,
                      window.type = "none",
                      window.size = NULL,
                      dist.quant = 1, n.IQR = 3, ...){
  n.pre = length(x.pre)
  n.post = length(x.post)
 
  # slide target window
  i = 1
  weight.stacked = NULL
  distance = NULL
  while (i <= n.post - k + 1) {
    Q = x.post[i:(i + k - 1)]
    Q = normalize(Q, norm.method)
    costs.qr = NULL
    
    # slide reference window
    continue = TRUE
    margin = default.margin
    j = 1
    while (continue) {  # j <= n.pre - k + 1
      # check if the search is finished
      if (j > n.pre - k + 3) {
        continue = FALSE
        next
      }
      # define R
      R = x.pre[j:min(j + k + margin - 1, n.pre)]
      R = normalize(R, norm.method)
      # check if R is too short
      R.too.short = RefTooShort(Q, R, step.pattern = step.pattern)
      if (R.too.short) {
        # check if the R can be extended
        if (j < n.pre - k - margin + 1) {
          margin = margin + 1
          next
        }else{
          margin = default.margin
          j = j + n.r
          next
        }
      }
      # match Q and R
      alignment.qr = dtw::dtw(Q, R, open.end = TRUE,
                              step.pattern = step.pattern,
                              window.type = window.type,
                              window.size = window.size,
                              distance.only = TRUE, ...)
      costs.qr = rbind(costs.qr,
                       data.frame(cost = alignment.qr$distance,
                                  j = j,
                                  margin = margin))
      j = j + n.r
      margin = default.margin
    }
    # find the minimum cost
    min.cost = which(costs.qr$cost == min(costs.qr$cost))[1]
    j.opt = costs.qr$j[min.cost]
    margin.opt = costs.qr$margin[min.cost]
    
    # obtain warping path W.pp.i: x.post -> n.pre
    Rs = x.pre[j.opt:min(j.opt + k + margin.opt - 1, n.pre)]
    Rs = normalize(Rs, norm.method)
    alignment.qrs = dtw::dtw(Q, Rs, open.end = TRUE,
                             step.pattern = step.pattern,
                             window.type = window.type,
                             window.size = window.size, ...)
    W.pp.i = Matrix::sparseMatrix(alignment.qrs$index1,
                                  alignment.qrs$index2)
    
    # obtain weight.b
    weight.a.Rs = weight.a[j.opt:(j.opt + ncol(W.pp.i) - 1)]
    weight.b = (W.pp.i %*% weight.a.Rs)/rowSums(as.matrix(W.pp.i))
    weight.b = as.numeric(weight.b)
    
    # convert warping path to weight
    weight.i = matrix(rep(NaN, n.post), nrow = 1)
    weight.i[1, i:(i + k - 1)] = weight.b
    
    # stack weight
    weight.stacked = rbind(weight.stacked, weight.i)
    
    # distance
    distance = c(distance, alignment.qrs$distance)
    
    # next
    i = i + n.q
  }
  
  # handle misfits
  misfits = which(distance > quantile(distance, dist.quant))
  if (length(misfits) > 0) {
    weight.stacked = weight.stacked[-misfits,]
  }

  # handle outliers
  weight.stacked = data.frame(weight.stacked) %>% 
    mutate_all(RemoveOutliers, n.IQR = n.IQR)
  
  # average weight
  avg.weight = colMeans(weight.stacked, na.rm = TRUE)
  avg.weight[is.na(avg.weight)] = 1
  
  return(list(x.post = x.post, x.pre =x.pre,
              k = k, weight.a = weight.a, 
              step.pattern = step.pattern,
              weight.stacked = weight.stacked,
              avg.weight = avg.weight,
              dist.quant = dist.quant, n.IQR = n.IQR))
}


# Two Step DTW
TFDTW = function(x, y, k, t.treat, buffer, 
                 norm.method = "t",
                 match.method = "fixed",
                 step.pattern1 = dtw::symmetricP2,
                 step.pattern2 = dtw::asymmetricP2,
                 plot.figures = FALSE, n.burn = 3,
                 ma = 3, ma.na = "original",
                 dist.quant = 1, n.IQR = 3,
                 window.type = "none",
                 default.margin = 3,
                 n.q = 1, n.r = 1, ...){
  # window size
  if (window.type == "sakoechiba") {
    window.size1 = as.integer(t.treat/2)
    window.size2 = as.integer(k/2)
  }else{
    window.size1 = NULL
    window.size2 = NULL
  }
  
  # 1st dtw
  res.1stDTW = first.dtw(x = x, y = y,
                         t.treat = t.treat, 
                         buffer = buffer, 
                         norm.method = norm.method, 
                         match.method = match.method,
                         step.pattern = step.pattern1,
                         window.type = window.type,
                         window.size = window.size1,
                         plot.figures = plot.figures, ...)
  cutoff = res.1stDTW$cutoff
  x.pre = x[1:cutoff]
  x.post = x[(cutoff - n.burn):length(x)]
  W.a = res.1stDTW$W.a
  
  # compute weight a
  weight.a.o = warp2weight(W.a)
  weight.a = as.numeric(stats::filter(weight.a.o, rep(1/ma, ma)))
  weight.a = zoo::na.locf(weight.a, na.rm = FALSE)
  if (ma.na == "one") {
    weight.a[is.na(weight.a)] = 1
  }else if(ma.na == "first.available") {
    weight.a[is.na(weight.a)] = weight.a[!is.na(weight.a)][1]
  }else if (ma.na == "original") {
    weight.a[is.na(weight.a)] = weight.a.o[is.na(weight.a)]
  }
  
  # 2nd dtw
  res.2ndDTW = second.dtw(x.post = x.post, x.pre = x.pre,
                          k = k, weight.a = weight.a, 
                          norm.method = norm.method,
                          step.pattern = step.pattern2,
                          dist.quant = dist.quant,
                          n.IQR = n.IQR, 
                          window.type = window.type,
                          window.size = window.size2,
                          default.margin = default.margin,
                          n.q = n.q, n.r = n.r, ...)
  avg.weight = res.2ndDTW$avg.weight
  avg.weight = avg.weight[(n.burn + 1):length(avg.weight)]
  
  return(list(y = y, x = x, k = k,
              t.treat = t.treat,
              cutoff = cutoff,
              n.burn = n.burn,
              W.a = W.a, weight.a = weight.a,
              weight.stacked = res.2ndDTW$weight.stacked,
              avg.weight = avg.weight))
}

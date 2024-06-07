## Functions -------------------------------------------------------------------
# normalization
t.normalize = function(data, reference = NULL){
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

minmax.normalize = function(data, reference = NULL){
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

normalize = function(data, norm.method, reference = NULL){
  if (norm.method == "minmax") {
    data = minmax.normalize(data, reference)
  }else if (norm.method == "t") {
    data = t.normalize(data, reference)
  }
  return(data)
}


# add buffer
add.buffer = function(TS, n){
  model.right = forecast::auto.arima(TS)
  right <- as.numeric(forecast::forecast(model.right, h = n)$mean)
  model.left = forecast::auto.arima(rev(TS))
  left <- rev(as.numeric(forecast::forecast(model.left, h = n)$mean))
  
  return(c(left, TS, right))
}


# pre-processing
preprocessing = function(data, filter.width = 5, norm.method = "t",
                         n.poly = 3, n.deri = 2, plot.data = FALSE){
  # transform
  values = reshape2::dcast(data %>% select(c("unit", "time", "value_raw")),
                           time ~ unit, value.var = "value_raw")
  
  # normalize
  values = values %>% mutate_at(setdiff(colnames(values), "time"),
                                ~normalize(., norm.method))
  
  # add buffer
  n.buffer = (filter.width - 1)/2
  values.w.buffer = sapply(values %>% select(-time),
                           add.buffer, n = n.buffer) %>% 
    data.frame(.)
  
  # derivative
  values.w.buffer = values.w.buffer %>%
    mutate_all(~signal::sgolayfilt(., n.poly, filter.width, n.deri)) 
  values[-1] = values.w.buffer[(n.buffer + 1):(n.buffer + nrow(values)),]
  
  # join
  df <- reshape2::melt(values, id.vars = 'time', variable.name = 'unit')
  data = right_join(df, data %>% select(-value), by = c("time", "unit"))
  data = data[c("id", "unit", "time", "value", colnames(data)[-(1:4)])]
  
  # plot
  if (plot.data) {
    ggplot(data, aes(x = time, y = value, color = unit)) +
      geom_line() +
      theme_bw()
  }
  
  return(data)
}


# transform warping path W to weight
warp2weight = function(W){
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
warpWITHweight = function(ts, weight){
  n.ts = length(ts)
  n.weight = length(weight)
  if (n.ts != n.weight) {
    print("Lengths of ts and weight are different.")
    return()
  }
  
  # form tranformation matrix
  count = 1
  residual = 0
  w = matrix(0, nrow = n.weight, ncol = 2*n.weight)
  for (i in 1:n.weight) {
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
  
  col.sum = colSums(w)
  ind.max = which(col.sum > 0) %>% 
    max
  w = w[, 1:ind.max]
  w[n.weight, ind.max] = 1 - sum(w[-n.weight, ind.max])
  
  # warp time series
  ts.warped = t(w) %*% matrix(ts, ncol = 1)
  
  return(ts.warped[,1])
}


# use W to warp ts
warpWITHpath = function(w, ts){
  w = as.matrix(w)
  ts.warped = (ts %*% w)/colSums(w)
  
  return(ts.warped)
}


# check if reference is too short in dtw
RefTooShort = function(query, reference, 
                       step.pattern = dtw::symmetricP2,
                       window.type = "none",
                       window.size = NULL){
  alignment = tryCatch(dtw::dtw(reference, query,
                                step.pattern = step.pattern,
                                window.type = window.type,
                                window.size = window.size,
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


# remove outliers in weight matrix
RemoveOutliers = function(data, n.IQR = 3){
  Q1 = quantile(data, 0.25, na.rm = TRUE)
  Q3 = quantile(data, 0.75, na.rm = TRUE)
  IQR = Q3 - Q1
  upper = Q3 + n.IQR*IQR
  lower = Q1 - n.IQR*IQR
  data[data > upper] = NaN
  data[data < lower] = NaN
  return(data)
}


# plot warped data
plot.warped = function(unit, dependent, t.treat,
                       y.raw, x.raw, x.warped,
                       legend.pos = c(0.3, 0.3)){
  unit.warped = paste0(unit, "-Warped")
  df.warp = data.frame(time = 1:length(x.raw),
                       y = y.raw[1:length(x.raw)],
                       x = x.raw[1:length(x.raw)],
                       warped = x.warped[1:length(x.raw)]) %>%
    `colnames<-`(c("time", dependent, unit, unit.warped)) %>%
    reshape2::melt(., id.vars = "time") %>%
    `colnames<-`(c("time", "unit", "value"))
  df.warp$unit = factor(df.warp$unit, 
                        levels = c(dependent, unit, unit.warped))
  
  fig.warp = df.warp %>%
    ggplot(aes(x = time, y = value, color = unit)) +
    geom_line() +
    scale_color_manual(values = c("#2a4d69", "#ee4035", "#7bc043")) +
    geom_vline(xintercept = t.treat, linetype="dashed",
               color = "grey30", size = 0.3) +
    theme_bw() +
    ggtitle(unit) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.pos = legend.pos)
  
  return(fig.warp)
}


# stack warped data figures
stack.warped = function(fig.list, file.name, ncol = 3){
  nrow = ceiling(length(fig.list)/ncol)
  fig = gridExtra::marrangeGrob(fig.list, ncol = ncol,
                                nrow = nrow)
  ggsave(file.name, fig, width = ncol*4, height = nrow*4,
         units = "in", limitsize = FALSE)
}


# plot synthetic control
plot.synth = function(res.synth.raw, res.synth.TFDTW,
                      dependent, start.time, end.time,
                      treat.time, file.name,
                      legend.pos = c(0.3, 0.3)){
  value = res.synth.raw$value
  avg.raw = res.synth.raw$average
  synth.raw = res.synth.raw$synthetic
  avg.TFDTW = res.synth.TFDTW$average
  synth.TFDTW = res.synth.TFDTW$synthetic
  
  df = data.frame(time = rep(start.time:end.time, 5),
                  unit = rep(c(dependent, "avg.raw", "avg.TFDTW",
                           "synth.raw", "synth.TFDTW"),
                           each = length(value)),
                  value = c(value, avg.raw, avg.TFDTW,
                            synth.raw, synth.TFDTW))
  df$unit = factor(df$unit, levels = c(dependent, "avg.raw", "avg.TFDTW",
                                       "synth.raw", "synth.TFDTW"))
  
  fig = ggplot(df, aes(x = time, y = value, color = unit)) +
    geom_line() +
    geom_vline(xintercept = treat.time, linetype="dashed") +
    scale_color_manual(values=c("#4a4e4d", "#0e9aa7", "#4b86b4",
                                "#f6cd61", "#fe8a71")) +
    theme_bw() +
    ggtitle(dependent) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.pos = legend.pos)
  
  ggsave(file.name, fig, width = 8, height = 6,
         units = "in", limitsize = FALSE)
}



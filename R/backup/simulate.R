## Functions -------------------------------------------------------------------
Diff2PhiUnif = function(xDiff,
                        speed_upper = 1,
                        speed_lower = 0.5,
                        weight_speed = TRUE,
                        rnd = 0.3
){
  pos_deriv = xDiff >= 0
  pos_ratio = sum(pos_deriv)/length(pos_deriv)
  speed = (rnd - 0.5)*2*(speed_upper - speed_lower) +
    sign(rnd - 0.5)*speed_lower
  if (weight_speed) {
    pos_speed = 1 + speed*(1 - pos_ratio)
    neg_speed = 1 - speed*(pos_ratio)
  }else{
    pos_speed = 1 + speed
    neg_speed = 1 - speed
    
  }
  phi = rep(0, length(xDiff))
  phi[pos_deriv] = pos_speed
  phi[!pos_deriv] = neg_speed
  phi = cumsum(phi)
  return(phi)
}


ApplyPhi = function(x, phi){
  output = NULL
  for (i in 1:length(phi)) {
    ind = floor(phi[i])
    weight = phi[i] - ind
    value = weight*(x[ind + 2] - x[ind + 1]) + x[ind + 1]
    output = c(output, value)
  }
  return(output)
}


SimData_ShapeSpeed = function(n = 5,
                              length = 1000,
                              rnd_speed = seq(0.1, 0.9, length.out = n),
                              n_SMA = 10,
                              ar_x = 0.9,
                              weight_speed = TRUE,
                              speed_upper = 1,
                              speed_lower = 0.5,
                              extra_x = round(0.2*length),
                              t_treat = 800,
                              shock = 50){
  # common exogenous shocks
  x = arima.sim(list(order = c(1,1,0), ar = ar_x), n = length + extra_x + n_SMA)
  xSMA10 = ts(TTR::SMA(x, n = n_SMA)[-(1:(n_SMA - 1))])
  xDiff = diff(xSMA10, difference = 1)
  
  # simulate
  data = NULL
  for (i in 1:n) {
    # speed profile
    phi = Diff2PhiUnif(xDiff[1:(length + extra_x)], 
                       weight_speed = weight_speed,
                       speed_upper = speed_upper,
                       speed_lower = speed_lower,
                       rnd = rnd_speed[i])
    
    # treatment
    if (i == 1) {
      treatment = c(rep(0, t_treat),
                    seq(0, shock, length.out = round(0.1*length)),
                    rep(shock, round(0.9*length + extra_x - t_treat)))
    }else{
      treatment = 0
    }
    
    y = ApplyPhi(x[-(1:(n_SMA - 1))], phi)
    y = y + treatment
    
    data = rbind(data,
                 data.frame(id = i,
                            unit = LETTERS[i],
                            time = 1:length,
                            value = y[1:length],
                            value_raw = y[1:length]))
  }
  return(data)
}


SimData_Lags = function(n = 3,
                        length = 100,
                        rnd_nCycles = seq(0.1, 0.9, length.out = n),
                        rnd_shift = seq(0.9, 0.1, length.out = n),
                        rnd_lag = seq(0.1, 0.9, length.out = n),
                        nCycles_min = 6,
                        nCycles_max = 12,
                        noise_mean = 0,
                        noise_sd = 0.01,
                        n_lag_min = 5,
                        n_lag_max = 15,
                        extra_x = 20,
                        beta = 0.9,
                        ar_x = 0.9,
                        t_treat = 80,
                        shock = 5){
  
  # prepare random numbers
  nCycles = rnd_nCycles * (nCycles_max - nCycles_min) + nCycles_min
  shifts = rnd_shift * 2 * pi
  n_lags = round(rnd_lag * (n_lag_max - n_lag_min) + n_lag_min, 0)
  
  # common exogenous shocks
  x = arima.sim(list(order = c(1,1,0), ar = ar_x), n = length + extra_x)
  # x = cumsum(sin(seq(0, 5*pi, length.out = length + n_lag))/2+0.5)
  
  # simulate
  data = NULL
  
  for (i in 1:n) {
    # speed profile
    phi = sin(seq(0, nCycles[i] * pi, length.out = length) + shifts[i])
    # trend
    if (i == 1) {
      trend = rep(0, length) + c(rep(0, length*4/5),
                                 seq(0, shock, length.out = length/20),
                                 seq(shock, 0, length.out = length/20),
                                 rep(0, length/5-length/10))
    }else{
      trend = rep(0, length)
    }
    y = NULL
    ylag = 1
    for (j in 1:length) {
      yt = trend[j] + beta*ylag + phi[j]*x[j + n_lags[i]] +
        (1 - phi[j])*x[j] + 
        rnorm(n = 1, mean = noise_mean, sd = noise_sd)
      y <- c(y, yt)
      ylag = yt
    }
    
    data = rbind(data,
                 data.frame(id = i,
                            unit = LETTERS[i],
                            time = 1:length,
                            value = y,
                            value_raw = y))
  }
  return(data)
}






run_simul = function(data, 
                     start_time = 1,
                     end_time = 100,
                     treat_time = 80,
                     width_range = (1:9)*2+3,
                     k_range = 4:12,
                     dtw1_time = 90,
                     dependent = "A",
                     dependent_id = 1,
                     dtw1_method = "open-end",
                     n_mse = 10,
                     n_IQR = 3,
                     n_burn = 3,
                     normalize_method = "t",
                     step.pattern2 = dtw::asymmetricP2,
                     dist_quantile = 1,
                     plot_figures = FALSE,
                     legend_position = c(0.3, 0.3),
                     ma = 3,
                     ma_na = "original",
                     step_pattern_range = list(
                       # symmetricP0 = dtw::symmetricP0, # too bumpy
                       # symmetricP05 = dtw::symmetricP05,
                       symmetricP1 = dtw::symmetricP1,
                       symmetricP2 = dtw::symmetricP2,
                       # asymmetricP0 = dtw::asymmetricP0, # too bumpy
                       # asymmetricP05 = dtw::asymmetricP05,
                       asymmetricP1 = dtw::asymmetricP1,
                       asymmetricP2 = dtw::asymmetricP2,
                       # typeIc = dtw::typeIc,
                       typeIcs = dtw::typeIcs,
                       # typeIIc = dtw::typeIIc,  # jumps
                       # typeIIIc = dtw::typeIIIc, # jumps
                       # typeIVc = dtw::typeIVc,  # jumps
                       # typeId = dtw::typeId,
                       typeIds = dtw::typeIds,
                       # typeIId = dtw::typeIId, # jumps
                       mori2006 = dtw::mori2006
                     ),
                     predictors.origin = NULL,
                     special.predictors.origin = list(list("value_raw", 70:79, c("mean")),
                                                      list("value_raw", 60:69, c("mean")),
                                                      list("value_raw", 50:59, c("mean"))),
                     time.predictors.prior.origin = 1:79,
                     time.optimize.ssr.origin = 1:79,
                     predictors.new = NULL,
                     special.predictors.new = list(list("value_warped", 70:79, c("mean")),
                                                   list("value_warped", 60:69, c("mean")),
                                                   list("value_warped", 50:59, c("mean"))),
                     time.predictors.prior.new = 1:79,
                     time.optimize.ssr.new = 1:79
){
  # grid
  grid_search = expand.grid(width_range, k_range,
                         names(step_pattern_range)) %>% 
    `colnames<-`(c("width", "k", "step_pattern")) %>% 
    mutate(mse_pre_original = NA_real_,
           mse_pre_new = NA_real_,
           mse_post_original = NA_real_,
           mse_post_new = NA_real_,
           pos_ratio = NA_real_,
           t_test = NA_real_)
  grid_search = grid_search %>% 
    split(., seq(nrow(grid_search)))
  
  # search
  result = grid_search %>% 
    future_map(
      ~{
        search = .
        width = search$width
        k = search$k
        step.pattern1 = step_pattern_range[[search$step_pattern]]
        
        data = preprocessing(data, filter_width = width)

        res = SimDesign::quiet(compare_methods(data = data,
                                               start_time = start_time,
                                               end_time = end_time,
                                               treat_time = treat_time,
                                               dtw1_time = dtw1_time,
                                               dependent = dependent,
                                               dependent_id = dependent_id,
                                               dtw1_method = dtw1_method,
                                               n_mse = n_mse,
                                               k = k,
                                               n_IQR = n_IQR,
                                               n_burn = n_burn,
                                               ma = ma,
                                               ma_na = ma_na,
                                               dist_quantile = dist_quantile,
                                               plot_figures = plot_figures,
                                               step.pattern1 = step.pattern1,
                                               step.pattern2 = step.pattern2,
                                               predictors.origin = predictors.origin,
                                               special.predictors.origin = special.predictors.origin,
                                               time.predictors.prior.origin = time.predictors.prior.origin,
                                               time.optimize.ssr.origin = time.optimize.ssr.origin,
                                               predictors.new = predictors.new,
                                               special.predictors.new = special.predictors.new,
                                               time.predictors.prior.new = time.predictors.prior.new,
                                               time.optimize.ssr.new = time.optimize.ssr.new,
                                               legend_position = legend_position))
        
        synth_original = res$synth_origin$synthetic
        synth_new = res$synth_new$synthetic
        value_raw = res$synth_origin$value
        
        gap_original = synth_original - value_raw
        gap_new = synth_new - value_raw
        
        mse_original = mean((gap_original)[treat_time:(treat_time + n_mse)]^2, rm.na = T)
        mse_new = mean((gap_new)[treat_time:(treat_time + n_mse)]^2, rm.na = T)
        
        mse_ratio = mse_new/mse_original
        
        
        mse = data.frame(width = width,
                         k = k,
                         step_pattern = search$step_pattern,
                         dtw1_time = dtw1_time,
                         mse_pre_original = res$mse$mse1_pre,
                         mse_pre_new = res$mse$mse2_pre,
                         mse_post_original = mse_original,
                         mse_post_new = mse_new,
                         mse_ratio = mse_ratio)
        
        list(mse = mse,
             synth_original = synth_original,
             synth_new = synth_new,
             value_raw = value_raw,
             gap_original = gap_original,
             gap_new = gap_new)
      }
    )
  
  return(result)
}

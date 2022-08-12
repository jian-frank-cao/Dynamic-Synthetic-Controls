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


run_simul = function(data, 
                     start_time = 1,
                     end_time = 100,
                     t_treat = 80,
                     width_range = (1:8)*2+3,
                     k_range = 4:12,
                     dtw1_range = 90,
                     step_pattern_range = list(
                       # symmetricP0 = dtw::symmetricP0, # too bumpy
                       symmetricP05 = dtw::symmetricP05,
                       symmetricP1 = dtw::symmetricP1,
                       symmetricP2 = dtw::symmetricP2,
                       # asymmetricP0 = dtw::asymmetricP0, # too bumpy
                       asymmetricP05 = dtw::asymmetricP05,
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
                     n_mse = 20
){
  # grid search
  grid_search = expand.grid(width_range, k_range,
                            names(step_pattern_range)) %>% 
    `colnames<-`(c("width", "k", "step_pattern")) %>% 
    mutate(mse_original = NA_real_,
           mse_new = NA_real_,
           mse_ratio = NA_real_)
  grid_search = grid_search %>% 
    split(., seq(nrow(grid_search)))
  
  result = grid_search %>% 
    future_map(
      ~{
        search = .
        width = search$width
        k = search$k
        step.pattern = step_pattern_range[[search$step_pattern]]
        
        data = preprocessing(data, filter_width = width)
        
        res = SimDesign::quiet(compare_methods(data = data,
                                               start_time = start_time,
                                               end_time = end_time,
                                               treat_time = t_treat,
                                               dtw1_time = dtw1_time,
                                               dependent = "A",
                                               dependent_id = 1,
                                               normalize_method = "t",
                                               k = k,
                                               synth_fun = "simulation",
                                               filter_width = width,
                                               plot_figures = FALSE,
                                               step.pattern = step.pattern))
        
        synth_original = res$synth_origin$synthetic
        synth_new = res$synth_new$synthetic
        value_raw = res$synth_origin$value
        
        mse_original = mean((synth_original - value_raw)[t_treat:(t_treat + n_mse)]^2, rm.na = T)
        mse_new = mean((synth_new - value_raw)[t_treat:(t_treat + n_mse)]^2, rm.na = T)
        
        mse_ratio = mse_new/mse_original
        
        
        mse = data.frame(width = width,
                         k = k,
                         step_pattern = search$step_pattern,
                         dtw1_time = dtw1_time,
                         mse_original = mse_original,
                         mse_new = mse_new,
                         mse_ratio = mse_ratio)
        
        list(mse = mse,
             synth_original = synth_original,
             synth_new = synth_new,
             value_raw = value_raw)
      }
    )
  
  return(result)
}

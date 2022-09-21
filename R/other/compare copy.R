## Functions -------------------------------------------------------------------
compare.methods = function(data,
                           start.time,
                           end.time,
                           treat.time,
                           # dtw1_time,
                           dependent,
                           dependent.id,
                           # dtw1_method = "fixed",
                           # filter_width = 5,
                           # k = 5,
                           args.TFDTW,
                           args.synth,
                           n.mse = 10,
                           # n_IQR = 3,
                           # n_burn = 3,
                           # dist_quantile = 1,
                           # ma = 3,
                           # ma_na = "original",
                           # n_q = 1,
                           # n_r = 1,
                           plot.figures = FALSE,
                           # normalize_method = "t",
                           # margin = 10,
                           # step.pattern1 = dtw::symmetricP2,
                           # step.pattern2 = dtw::asymmetricP2,
                           # predictors.origin,
                           # special.predictors.origin,
                           # time.predictors.prior.origin,
                           # time.optimize.ssr.origin,
                           # predictors.new,
                           # special.predictors.new,
                           # time.predictors.prior.new,
                           # time.optimize.ssr.new,
                           legend.position = c(0.3, 0.3), ...){
  
  
  # prepare data
  t.treat = (treat.time - start.time) + 1
  # n = (end.time - start.time) + 1
  # n_dtw1 = (dtw1_time - start.time) + 1
  
  y.raw = data %>% 
    filter(unit == dependent & time <= end.time) %>%
    .[["value_raw"]]
  
  y.processed = data %>% 
    filter(unit == dependent & time <= end.time) %>%
    .[["value"]]
  
  x.list = data %>% 
    filter(unit != dependent & time <= end.time) %>% 
    select(c("unit", "time", "value", "value_raw")) %>% 
    group_by(unit) %>% 
    group_split(.keep = TRUE)
  
  # TSDTW
  results = NULL
  for (item in x.list) {
    unit = as.character(item$unit[1])
    # x = item$value
    # x.raw = item$value_raw
    # y = y.processed
    args.TFDTW[["x"]] = item$value
    args.TFDTW[["y"]] = y.processed
    args.TFDTW[["t.treat"]] = t.treat
    args.TFDTW[["plot.figures"]] = plot.figures
    
    # res = TwoStepDTW(x = x, y = y, k = k,
    #                  n_dtw1 = n_dtw1, t.treat = t.treat, 
    #                  normalize_method = normalize_method,
    #                  dtw1_method = dtw1_method,
    #                  step.pattern1 = step.pattern1, 
    #                  plot.figures = plot.figures,
    #                  n_burn = n_burn,
    #                  ma = ma, ma_na = ma_na,
    #                  step.pattern2 = step.pattern2, 
    #                  dist_quantile = dist_quantile,
    #                  n_IQR = n_IQR, ...)
    res = do.call(TFDTW, args.TFDTW)
    
    x_warped = c(
      warpWITHweight(x.raw[1:res$cutoff], res$weight_a)[1:t.treat],
      warpWITHweight(x.raw[res$cutoff:length(x.raw)], res$avg_weight)[-1]
    )
    
    if (plot.figures) {
      df.warp = data.frame(time = 1:length(x.raw),
                      y = y.raw[1:length(x.raw)],
                      x = x.raw[1:length(x.raw)],
                      warped = x_warped[1:length(x.raw)]) %>%
        `colnames<-`(c("time", dependent, unit, paste0(unit, "-Warped"))) %>%
        reshape2::melt(., id.vars = "time") %>%
        `colnames<-`(c("time", "unit", "value"))
      
      fig.warp = df %>%
        ggplot(aes(x = time, y = value, color = unit)) +
        geom_line() +
        scale_color_manual(values = c("#2a4d69", "#ee4035", "#7bc043")) +
        geom_vline(xintercept = t.treat, linetype="dashed",
                   color = "grey30", size = 0.3) +
        theme_bw() +
        ggtitle(unit) +
        theme(plot.title = element_text(hjust = 0.5),
              legend.position = legend.position)
    }else{
      df.warp = NULL
      fig.warp = NULL
    }
    
    results[[unit]] = list(unit = unit,
                           x = x.raw,
                           res = res,
                           x_warped = x_warped,
                           df.warp = df.warp,
                           fig.warp = fig.warp)
  }
  
  # plot warped data
  if (plot.figures) {
    file_name = paste0("./figures/warped-", 
                       paste0(c(dependent, filter_width,
                                k, end.time, treat.time,
                                dtw1_time), collapse = "-"),
                       ".pdf")
    plot.warped(lapply(results,"[[","fig.warp"),
                ncol = 3,
                file_name)
  }
  
  # prepare data for synthetic control
  df.synth = results %>% 
    future_map(
      ~{
        item = .
        unit = item$unit
        x_warped = item$x_warped
        
        res = data.frame(
          time = 1:length(y.raw) + start.time - 1,
          unit = unit,
          value_warped = NA,
          stringsAsFactors = FALSE
        )
        res$value_warped = x_warped[1:((length(y.raw)-1) + 1)]  
        
        res
      }
    ) %>% 
    do.call("rbind", .) %>% 
    `row.names<-`(NULL)
  
  df.synth = rbind(df.synth,
             data.frame(time = 1:length(y.raw) + start.time - 1,
                        unit = dependent,
                        value_warped = y.raw))
  
  df.synth = right_join(data, df.synth, by = c("unit", "time"))
  df.synth = data.frame(df.synth)
  
  # synthetic control
  synth_origin = do_synth(df.synth, "value_raw", dependent.id,
                          predictors = predictors.origin,
                          special.predictors = special.predictors.origin,
                          time.predictors.prior = time.predictors.prior.origin,
                          time.optimize.ssr = time.optimize.ssr.origin)
  synth_new = do_synth(df.synth, "value_warped", dependent.id,
                       predictors = predictors.new,
                       special.predictors = special.predictors.new,
                       time.predictors.prior = time.predictors.prior.new,
                       time.optimize.ssr = time.optimize.ssr.new)
  
  # plot synthetic control
  if (plot.figures) {
    plot_synth(synth_origin, 
               dependent, treat.time,
               start.time, end.time,
               paste0("./figures/synth_",
                      paste0(c(dependent, "without_TSDTW", 
                               filter_width,k, end.time, treat.time,
                               dtw1_time), collapse = "_"),
                      ".pdf"))
    plot_synth(synth_new, 
               dependent, treat.time,
               start.time, end.time,
               paste0("./figures/synth_",
                      paste0(c(dependent, "TSDTW", 
                               filter_width,k, end.time, treat.time,
                               dtw1_time), collapse = "_"),
                      ".pdf"))
  }
  
  
  # MSE
  diff1 = (synth_origin$synthetic - synth_origin$value)
  diff2 = (synth_new$synthetic - synth_new$value)
  
  mse1_1 = mean(diff1[1:(t.treat - 1)]^2, na.rm = T)
  mse2_1 = mean(diff2[1:(t.treat - 1)]^2, na.rm = T)
  
  mse1_2 = mean(diff1[t.treat:(t.treat + n.mse - 1)]^2, na.rm = T)
  mse2_2 = mean(diff2[t.treat:(t.treat + n.mse - 1)]^2, na.rm = T)
  
  
  return(list(
    dtw_results = results,
    df.synth = df.synth,
    synth_origin = synth_origin,
    synth_new = synth_new,
    mse = data.frame(mse1_pre = mse1_1,
                     mse2_pre = mse2_1,
                     mse1_post = mse1_2,
                     mse2_post = mse2_2)
  ))
}


# pre-processing
preprocessing = function(data,
                         filter_width = 5,
                         normalize_method = "t",
                         n_poly = 3,
                         n_derivative = 2,
                         plot_data = FALSE
){
  
  # transform
  values = reshape2::dcast(data %>% select(c("unit", "time", "value_raw")),
                           time ~ unit, value.var = "value_raw")
  
  # normalize
  values = values %>% mutate_at(setdiff(colnames(values), "time"),
                                ~normalize(., normalize_method))
  
  # add buffer
  n_buffer = (filter_width - 1)/2
  values_w_buffer = sapply(values %>% select(-time),
                           add.buffer, n = n_buffer) %>% 
    data.frame(.)
  
  # derivative
  values_w_buffer = values_w_buffer %>%
    mutate_all(~signal::sgolayfilt(., n_poly, filter_width, n_derivative)) 
  values[-1] = values_w_buffer[(n_buffer + 1):(n_buffer + nrow(values)),]
  
  # join
  df <- reshape2::melt(values ,  id.vars = 'time',
                       variable.name = 'unit')
  
  data = right_join(df, data %>% select(-value), by = c("time", "unit"))
  data = data[c("id", "unit", "time", "value", colnames(data)[-(1:4)])]
  
  # plot
  if (plot_data) {
    ggplot(data, aes(x = time, y = value, color = unit)) +
      geom_line() +
      theme_bw()
  }
  
  return(data)
}




# compare methods for all units
run_all_units = function(data,
                         start.time,
                         end.time,
                         treat.time,
                         dtw1_time,
                         predictors.origin,
                         special.predictors.origin,
                         time.predictors.prior.origin,
                         time.optimize.ssr.origin,
                         predictors.new,
                         special.predictors.new,
                         time.predictors.prior.new,
                         time.optimize.ssr.new,
                         filter_width = NULL,
                         k = 6,
                         dtw1_method = "fixed",
                         plot.figures = FALSE, 
                         normalize_method = "t",
                         step.pattern1 = dtw::symmetricP2,
                         step.pattern2 = dtw::asymmetricP2,
                         legend.position = c(0.3, 0.3),
                         n.mse = 10,
                         n_IQR = 3,
                         n_burn = 3,
                         dist_quantile = 1,
                         ma = 3,
                         ma_na = "original",
                         # margin = 10,
                         detail = FALSE
){
  # prepare data
  if (!is.null(filter_width)) {
    data = preprocessing(data, filter_width)
  }
  units = data[c("id", "unit")] %>% distinct
  
  # run
  result = as.list(1:nrow(units)) %>% 
    future_map(
      ~{
        i = .
        dependent = units$unit[i]
        dependent.id = units$id[i]
        res = SimDesign::quiet(compare_methods(data = data,
                              start.time = start.time,
                              end.time = end.time,
                              treat.time = treat.time,
                              dtw1_time = dtw1_time,
                              dependent = dependent,
                              dependent.id = dependent.id,
                              dtw1_method = dtw1_method,
                              filter_width = filter_width,
                              k = k,
                              n.mse = n.mse,
                              n_IQR = n_IQR,
                              n_burn = n_burn,
                              dist_quantile = dist_quantile,
                              ma = ma,
                              ma_na = ma_na,
                              normalize_method = normalize_method,
                              plot.figures = plot.figures,
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
                              legend.position = legend.position))
        res$mse = res$mse %>% mutate(unit = dependent)
        if (detail == FALSE) {
          list(mse = res$mse)
        }else{
          res
        }
      }
    )
  
  # compute log ratio
  mse = lapply(result, '[[', "mse") %>% 
    do.call("rbind", .) %>% 
    mutate(ratio = mse2_post/mse1_post,
           log_ratio = log(ratio))
  
  pos_ratio = length(which(mse$log_ratio < 0))/nrow(mse)
  t_test = t.test(mse$log_ratio)
  
  return(list(result = result,
              width = filter_width,
              k = k,
              start.time = start.time,
              end.time = end.time,
              treat.time = treat.time,
              dtw1_time = dtw1_time,
              mse = mse,
              pos_ratio = pos_ratio,
              t_test = t_test$p.value))
}


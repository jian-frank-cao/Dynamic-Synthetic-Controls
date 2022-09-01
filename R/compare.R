## Functions -------------------------------------------------------------------
compare_methods = function(data,
                           start_time,
                           end_time,
                           treat_time,
                           dtw1_time,
                           dependent,
                           dependent_id,
                           filter_width = 5,
                           n_mse = 100,
                           k = 15,
                           dist_quantile = 1,
                           ma = 3,
                           ma_na = "original",
                           n_q = 1,
                           n_r = 1,
                           plot_figures = FALSE,
                           normalize_method = "t",
                           margin = 10,
                           step.pattern1 = dtw::symmetricP2,
                           step.pattern2 = dtw::asymmetricP2,
                           predictors.origin,
                           special.predictors.origin,
                           time.predictors.prior.origin,
                           time.optimize.ssr.origin,
                           predictors.new,
                           special.predictors.new,
                           time.predictors.prior.new,
                           time.optimize.ssr.new,
                           legend_position = c(0.3, 0.3), ...){
  
  
  # prepare data
  t_treat = (treat_time - start_time) + 1
  n = (end_time - start_time) + 1
  n_dtw1 = (dtw1_time - start_time) + 1
  
  y_raw = data %>% 
    filter(unit == dependent &
             time <= end_time) %>%
    .[["value_raw"]]
  
  y_processed = data %>% 
    filter(unit == dependent &
             time <= end_time) %>%
    .[["value"]]
  
  x_list = data %>% 
    filter(unit != dependent &
             time <= end_time) %>% 
    select(c("unit", "time", "value", "value_raw")) %>% 
    group_by(unit) %>% 
    group_split(.keep = TRUE)
  
  # TSDTW
  results = NULL
  for (item in x_list) {
    unit = as.character(item$unit[1])
    x = item$value
    x_raw = item$value_raw
    y = y_processed
    
    res = TwoStepDTW(x, y, t_treat, k, n_dtw1,
                     normalize_method = normalize_method,
                     ma = ma, ma_na = ma_na,
                     step.pattern1 = step.pattern1, 
                     step.pattern2 = step.pattern2, 
                     plot_figures = plot_figures, ...)
    
    # x_warped = c(warp_ts(res$W_a, x_raw[1:res$cutoff]),
    #              warp_using_weight(x_raw[-(1:(res$cutoff - 1))],
    #                                res$avg_weight)[-1])
    x_warped = c(warp_using_weight(x_raw[1:res$cutoff], res$weight_a)[1:t_treat],
                warp_using_weight(x_raw[-(1:(res$cutoff - 1))], res$avg_weight)[-1])
    
    
    if (plot_figures) {
      df = data.frame(time = 1:length(x_raw),
                      y = y_raw[1:length(x_raw)],
                      x = x_raw[1:length(x_raw)],
                      warped = x_warped[1:length(x_raw)]) %>%
        `colnames<-`(c("time", dependent, unit, paste0(unit, "-Warped"))) %>%
        reshape2::melt(., id.vars = "time") %>%
        `colnames<-`(c("time", "unit", "value"))
      
      fig = df %>%
        ggplot(aes(x = time, y = value, color = unit)) +
        geom_line() +
        scale_color_manual(values = c("#2a4d69", "#ee4035", "#7bc043")) +
        geom_vline(xintercept = t_treat, linetype="dashed",
                   color = "grey30", size = 0.3) +
        theme_bw() +
        ggtitle(unit) +
        theme(plot.title = element_text(hjust = 0.5),
              legend.position = legend_position)
    }else{
      df = NULL
      fig = NULL
    }
    
    results[[unit]] = list(unit = unit,
                           x = x_raw,
                           res = res,
                           x_warped = x_warped,
                           df = df,
                           fig = fig)
  }
  
  # plot warped data
  if (plot_figures) {
    file_name = paste0("./figures/warped-", 
                       paste0(c(dependent, filter_width,
                                k, end_time, treat_time,
                                dtw1_time), collapse = "-"),
                       ".pdf")
    plot_warped(lapply(results,"[[","fig"),
                ncol = 3,
                file_name)
  }
  
  # prepare data for synthetic control
  df = results %>% 
    future_map(
      ~{
        item = .
        unit = item$unit
        x_warped = item$x_warped
        
        res = data.frame(
          time = 1:length(y_raw) + start_time - 1,
          unit = unit,
          value_warped = NA,
          stringsAsFactors = FALSE
        )
        res$value_warped = x_warped[1:((length(y_raw)-1) + 1)]  
        
        res
      }
    ) %>% 
    do.call("rbind", .) %>% 
    `row.names<-`(NULL)
  
  df = rbind(df,
             data.frame(time = 1:length(y_raw) + start_time - 1,
                        unit = dependent,
                        value_warped = y_raw))
  
  df = right_join(data, df, by = c("unit", "time"))
  df = data.frame(df)
  
  # synthetic control
  synth_origin = do_synth(df, "value_raw", dependent_id,
                          predictors = predictors.origin,
                          special.predictors = special.predictors.origin,
                          time.predictors.prior = time.predictors.prior.origin,
                          time.optimize.ssr = time.optimize.ssr.origin)
  synth_new = do_synth(df, "value_warped", dependent_id,
                       predictors = predictors.new,
                       special.predictors = special.predictors.new,
                       time.predictors.prior = time.predictors.prior.new,
                       time.optimize.ssr = time.optimize.ssr.new)
  
  # plot synthetic control
  if (plot_figures) {
    plot_synth(synth_origin, 
               dependent, treat_time,
               start_time, end_time,
               paste0("./figures/synth_",
                      paste0(c(dependent, "without_TSDTW", 
                               filter_width,k, end_time, treat_time,
                               dtw1_time), collapse = "_"),
                      ".pdf"))
    plot_synth(synth_new, 
               dependent, treat_time,
               start_time, end_time,
               paste0("./figures/synth_",
                      paste0(c(dependent, "TSDTW", 
                               filter_width,k, end_time, treat_time,
                               dtw1_time), collapse = "_"),
                      ".pdf"))
  }
  
  
  # MSE
  diff1 = (synth_origin$synthetic - synth_origin$value)
  diff2 = (synth_new$synthetic - synth_new$value)
  
  mse1_1 = mean(diff1[1:(t_treat - 1)]^2, na.rm = T)
  mse2_1 = mean(diff2[1:(t_treat - 1)]^2, na.rm = T)
  
  mse1_2 = mean(diff1[t_treat:(t_treat + n_mse - 1)]^2, na.rm = T)
  mse2_2 = mean(diff2[t_treat:(t_treat + n_mse - 1)]^2, na.rm = T)
  
  
  return(list(
    dtw_results = results,
    df = df,
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
                           add_buffer, n = n_buffer) %>% 
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
                         start_time,
                         end_time,
                         treat_time,
                         dtw1_time,
                         predictors.origin,
                         special.predictors.origin,
                         time.predictors.prior.origin,
                         time.optimize.ssr.origin,
                         predictors.new,
                         special.predictors.new,
                         time.predictors.prior.new,
                         time.optimize.ssr.new,
                         plot_figures = FALSE, 
                         normalize_method = "t",
                         step.pattern1 = dtw::symmetricP2,
                         step.pattern2 = dtw::asymmetricP2,
                         legend_position = c(0.3, 0.3),
                         filter_width = NULL,
                         n_mse = 10,
                         k = 6,
                         n_IQR = 3,
                         dist_quantile = 0.95,
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
        dependent_id = units$id[i]
        res = SimDesign::quiet(compare_methods(data = data,
                              start_time = start_time,
                              end_time = end_time,
                              treat_time = treat_time,
                              dtw1_time = dtw1_time,
                              dependent = dependent,
                              dependent_id = dependent_id,
                              n_mse = n_mse,
                              k = k,
                              n_IQR = n_IQR,
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
              start_time = start_time,
              end_time = end_time,
              treat_time = treat_time,
              dtw1_time = dtw1_time,
              mse = mse,
              pos_ratio = pos_ratio,
              t_test = t_test$p.value))
}


## Functions -------------------------------------------------------------------
# do synth
do_synth = function(df, dep_var, dependent_id, predictors,
                    special.predictors, time.predictors.prior,
                    time.optimize.ssr){
  # find v
  dataprep.out <-
    Synth::dataprep(
      foo = df,
      predictors    = predictors,
      dependent     = dep_var,
      unit.variable = 1,
      time.variable = 3,
      special.predictors = special.predictors,
      treatment.identifier = dependent_id,
      controls.identifier = setdiff(unique(df$id), dependent_id),
      time.predictors.prior = time.predictors.prior,
      time.optimize.ssr = time.optimize.ssr, 
      unit.names.variable = 2,
      time.plot = sort(unique(df$time))
    )
  
  # fit training model
  synth.out <- 
    Synth::synth(
      data.prep.obj=dataprep.out,
      Margin.ipop=.005,Sigf.ipop=7,Bound.ipop=6
    )
  
  value = df %>% filter(id == dependent_id) %>% `$`(value_raw)
  average = df %>% filter(id != dependent_id) %>% group_by(time) %>% 
    summarise(average = mean(value_raw, na.rm = TRUE)) %>% `$`(average)
  synthetic = dataprep.out$Y0%*%synth.out$solution.w %>% as.numeric
  
  return(list(value = value,
              average = average,
              synthetic = synthetic))
}


# plot synth
plot_synth = function(res_synth, dependent, treat_time,
                      start_time, end_time, file_name){
  value = res_synth$value
  average = res_synth$average
  synthetic = res_synth$synthetic
  
  df = rbind(data.frame(time = start_time:end_time,
                        unit = dependent,
                        value = value),
             data.frame(time = start_time:end_time,
                        unit = "Average",
                        value = average),
             data.frame(time = start_time:end_time,
                        unit = "Synthetic",
                        value = synthetic))
  
  fig = ggplot(df, aes(time,value)) +
    geom_line(aes(colour = unit)) +
    geom_vline(xintercept = treat_time, linetype="dashed") +
    scale_color_manual(values=c("#ee4035", "#2a4d69", "#7bc043")) +
    theme_bw()
  
  ggsave(file_name,
         fig, width = 8, height = 6,
         units = "in", limitsize = FALSE)
}

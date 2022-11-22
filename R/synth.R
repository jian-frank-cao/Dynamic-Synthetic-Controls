## Functions -------------------------------------------------------------------
# do synth
do.synth = function(df, dep.var, dependent.id, predictors,
                    special.predictors, time.predictors.prior,
                    time.optimize.ssr){
  # find v
  dataprep.out <-
    Synth::dataprep(
      foo = df,
      predictors    = eval(predictors),
      dependent     = dep.var,
      unit.variable = 1,
      time.variable = 3,
      special.predictors = eval(special.predictors),
      treatment.identifier = dependent.id,
      controls.identifier = setdiff(unique(df$id), dependent.id),
      time.predictors.prior = time.predictors.prior,
      time.optimize.ssr = time.optimize.ssr, 
      unit.names.variable = 2,
      time.plot = df.synth %>%
        filter(id == dependent.id & !is.na(value_warped)) %>% 
        .[["time"]] %>% unique %>% sort
    )
  
  # fit training model
  synth.out <- 
    Synth::synth(
      data.prep.obj=dataprep.out,
      Margin.ipop=.005,Sigf.ipop=7,Bound.ipop=6
    )
  
  value = df %>% filter(id == dependent.id) %>% `$`(value_warped)
  average = df %>% filter(id != dependent.id) %>% group_by(time) %>% 
    summarise(average = mean(!!sym(dep.var), na.rm = TRUE)) %>% `$`(average)
  synthetic = dataprep.out$Y0%*%synth.out$solution.w %>% as.numeric
  
  return(list(value = value,
              average = average,
              synthetic = synthetic))
}


# plot synth
# plot.synth = function(res.synth, dependent, treat.time,
#                       start.time, end.time, file.name){
#   value = res.synth$value
#   average = res.synth$average
#   synthetic = res.synth$synthetic
#   
#   df = rbind(data.frame(time = start.time:end.time,
#                         unit = dependent,
#                         value = value),
#              data.frame(time = start.time:end.time,
#                         unit = "Average",
#                         value = average),
#              data.frame(time = start.time:end.time,
#                         unit = "Synthetic",
#                         value = synthetic))
#   
#   fig = ggplot(df, aes(time,value)) +
#     geom_line(aes(colour = unit)) +
#     geom_vline(xintercept = treat.time, linetype="dashed") +
#     scale_color_manual(values=c("#ee4035", "#2a4d69", "#7bc043")) +
#     theme_bw()
#   
#   ggsave(file.name,
#          fig, width = 8, height = 6,
#          units = "in", limitsize = FALSE)
# }



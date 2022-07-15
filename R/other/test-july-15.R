result = NULL

for (i in 2:length(step_pattern_range)) {
  step.pattern = step_pattern_range[[i]]
  result[[i]] = as.list(1:nrow(units)) %>% 
    future_map(
      ~{
        i = .
        dependent = units$unit[i]
        dependent_id = units$id[i]
        # print(paste0(dependent, ":", i, "-", k, " start..."))
        res = compare_methods(data = data,
                              start_time = start_time,
                              end_time = end_time,
                              treat_time = treat_time,
                              dtw1_time = dtw1_time,
                              dependent = dependent,
                              dependent_id = dependent_id,
                              normalize_method = "t",
                              filter_width = filter_width,
                              k = k,
                              plot_figures = F,
                              synth_fun = synth_fun,
                              step.pattern = step.pattern)
        # print(paste0(dependent, ":", i, "-", k, " start...Done."))
        res$mse = res$mse %>% mutate(dependent = dependent, k = k)
        res
      }
    )
  print(i)
}





for (i in 1:17) {
  dependent = units$unit[i]
  dependent_id = units$id[i]
  # print(paste0(dependent, ":", i, "-", k, " start..."))
  res = compare_methods(data = data,
                        start_time = start_time,
                        end_time = end_time,
                        treat_time = treat_time,
                        dtw1_time = dtw1_time,
                        dependent = dependent,
                        dependent_id = dependent_id,
                        normalize_method = "t",
                        filter_width = filter_width,
                        k = k,
                        plot_figures = F,
                        synth_fun = synth_fun,
                        step.pattern = step.pattern)
}


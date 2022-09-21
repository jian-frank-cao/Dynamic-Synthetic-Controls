## Functions -------------------------------------------------------------------
TFDTW.synth = function(data, start.time, end.time, treat.time,
                       dependent, dependent.id,
                       args.TFDTW, args.synth,
                       res.synth.raw = NULL,
                       n.mse = 10, plot.figures = FALSE,
                       plot.path = "./figures/",
                       legend.pos = c(0.3, 0.3)){
  # prepare data for TFDTW
  t.treat = (treat.time - start.time) + 1

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
  
  args.TFDTW[["y"]] = y.processed
  args.TFDTW[["t.treat"]] = t.treat
  args.TFDTW[["plot.figures"]] = plot.figures
  
  # TFDTW
  results = NULL
  for (item in x.list) {
    unit = as.character(item$unit[1])
    x.raw = item$value_raw
    args.TFDTW[["x"]] = item$value
    
    res = do.call(TFDTW, args.TFDTW)
    
    x.warped = c(
      warpWITHweight(x.raw[1:res$cutoff], res$weight.a)[1:t.treat],
      warpWITHweight(x.raw[res$cutoff:length(x.raw)], res$avg.weight)[-1]
    )
    
    df.synth = data.frame(
      time = 1:length(y.raw) + start.time - 1,
      unit = unit,
      value_warped = NA,
      stringsAsFactors = FALSE
    )
    df.synth$value_warped = x.warped[1:length(y.raw)]
    
    fig.warp = NULL
    if (plot.figures) {
      fig.warp = plot.warped(unit = unit, dependent = dependent,
                             t.treat = t.treat, y.raw = y.raw,
                             x.raw = x.raw, x.warped = x.warped,
                             legend.pos = legend.pos)
    }
    
    results[[unit]] = list(unit = unit, 
                           x.warped = x.warped,
                           df.synth = df.synth,
                           fig.warp = fig.warp)
  }
  
  # plot warped data
  if (plot.figures) {
    file.name = paste0(plot.path, "warped-", dependent, ".pdf")
    stack.warped(fig.list = lapply(results,"[[","fig.warp"),
                 file.name = file.name, ncol = 3)
  }
  
  # prepare data for Synth
  df.synth = lapply(results, "[[", "df.synth") %>% 
    do.call("rbind", .) %>% 
    `row.names<-`(NULL)
  
  df.synth = rbind(df.synth,
                   data.frame(time = 1:length(y.raw) + start.time - 1,
                              unit = dependent,
                              value_warped = y.raw))
  
  df.synth = right_join(data, df.synth, by = c("unit", "time"))
  df.synth = data.frame(df.synth)
  
  args.synth[["df"]] = df.synth
  args.synth[["dependent.id"]] = dependent.id
  
  # Synth
  if (is.null(res.synth.raw)) {
    args.synth[["dep.var"]] = "value_raw"
    res.synth.raw = do.call(do.synth, args.synth)
  }
  
  args.synth[["dep.var"]] = "value_warped"
  res.synth.TFDTW = do.call(do.synth, args.synth)
  
  # plot synthetic control
  if (plot.figures) {
    file.name = paste0(plot.path, "synth_", dependent, ".pdf")
    plot.synth(res.synth.raw = res.synth.raw,
               res.synth.TFDTW = res.synth.TFDTW,
               dependent = dependent, start.time = start.time,
               end.time = end.time, treat.time = treat.time,
               file.name = file.name, legend.pos = legend.pos)
  }
  
  # output
  gap.raw = res.synth.raw$value - res.synth.raw$synthetic
  gap.TFDTW = res.synth.TFDTW$value - res.synth.TFDTW$synthetic
  
  mse = data.frame(
    unit = dependent,
    mse.preT.raw = mean(gap.raw[1:(t.treat-1)]^2, na.rm = T),
    mse.preT.TFDTW = mean(gap.TFDTW[1:(t.treat-1)]^2, na.rm = T),
    mse.postT.raw = mean(gap.raw[(t.treat+1):(t.treat+n.mse)]^2, na.rm = T),
    mse.postT.TFDTW = mean(gap.TFDTW[(t.treat+1):(t.treat+n.mse)]^2, na.rm = T)
  )
  
  return(list(dependent = dependent, dependent.id = dependent.id,
              # args.TFDTW = args.TFDTW, args.synth = args.synth,
              # results.TFDTW = results, df.synth = df.synth,
              res.synth.raw = res.synth.raw,
              res.synth.TFDTW = res.synth.TFDTW,
              gap.raw = gap.raw, gap.TFDTW = gap.TFDTW,
              mse = mse))
}


TFDTW.synth.all.units = function(data, target, 
                                 args.TFDTW.synth,
                                 filter.width = NULL,
                                 res.synth.raw.list = NULL,
                                 all.units.parallel = FALSE){
  # prepare data
  if (!is.null(filter.width)) {
    data = preprocessing(data, filter.width)
  }
  args.TFDTW.synth[["data"]] = data
  units = data[c("id", "unit")] %>% distinct
  units.list = units %>% split(., seq(nrow(units)))
  
  # run TFDTW.synth
  if (all.units.parallel) {
    fun.map = furrr::future_map
  }else{
    fun.map = purrr::map
  }
  results = units.list %>% 
    set_names(units$unit) %>% 
    fun.map(
      ~{
        item = .
        dependent = item$unit
        dependent.id = item$id
        args.TFDTW.synth[["dependent"]] = dependent
        args.TFDTW.synth[["dependent.id"]] = dependent.id
        args.TFDTW.synth[["res.synth.raw"]] = res.synth.raw.list[[dependent]]
        do.call(TFDTW.synth, args.TFDTW.synth)
      }
    )
  
  # compute log ratio
  mse = lapply(results, '[[', "mse") %>% 
    do.call("rbind", .) %>%
    mutate(ratio = mse.postT.TFDTW/mse.postT.raw,
           log.ratio = log(ratio)) %>% 
    filter(unit != target)
  
  # tests
  neg.ratio = sum(mse$log.ratio < 0)/nrow(mse)
  p.value = t.test(mse$log.ratio)$p.value
  
  return(list(target = target,
              filter.width = filter.width,
              # args.TFDTW.synth = args.TFDTW.synth,
              # results.TFDTW.synth = results,
              res.synth.target.raw = results[[target]]$res.synth.raw,  #
              res.synth.target.TFDTW = results[[target]]$res.synth.TFDTW,  #
              mse = mse,
              neg.ratio = neg.ratio,
              p.value = p.value))
}



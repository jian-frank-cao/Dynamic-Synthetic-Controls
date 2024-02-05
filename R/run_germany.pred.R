args = commandArgs(trailingOnly=TRUE)
task_id = as.integer(args[1])

if (task_id >= 12 & task_id <= 21) {
  task_id = task_id - 12
  cat(paste0("basque-", task_id))
  source("./R/run_basque.pred.R")
}else if (task_id >= 22) {
  task_id = task_id - 22
  cat(paste0("tobacco-", task_id))
  source("./R/run_tobacco.pred.R")
}else{
  job.start = Sys.time()
  
  ## Setup -----------------------------------------------------------------------
  library(checkpoint)
  checkpoint("2022-04-01")
  
  library(parallel)
  n.cores = detectCores()
  library(tidyverse)
  library(furrr)
  plan(multisession, workers = n.cores)
  options(future.rng.onMisuse="ignore")
  options(stringsAsFactors = FALSE)
  
  source("./R/utility/misc.R")
  source("./R/utility/TFDTW.R")
  source("./R/utility/synth.R")
  source("./R/utility/implement.R")
  source("./R/utility/grid.search.R")
  set.seed(20220407)
  
  
  ## Function --------------------------------------------------------------------
  impute_values <- function(column) {
    # Use last observation carried forward (LOCF) imputation
    column <- zoo::na.locf(column, na.rm = FALSE)
    
    # Use next observation carried backward (NOCB) for leading NAs
    if (is.na(column[1])) {
      column <- zoo::na.locf(column, fromLast = TRUE)
    }
    return(column)
  }
  
  TFDTW.synth = function(data, start.time, end.time, treat.time,
                         dependent, dependent.id,
                         args.TFDTW, args.synth,
                         pred.vars = NULL,
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
      # select(c("unit", "time", "value", "value_raw")) %>%
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
      
      item$value_warped = x.warped[1:length(y.raw)]
      
      if (!is.null(pred.vars)) {
        for (pred.var in pred.vars) {
          pred.target = item[[pred.var]]
          pred.target.warped = c(
            warpWITHweight(pred.target[1:res$cutoff], res$weight.a)[1:t.treat],
            warpWITHweight(pred.target[res$cutoff:length(pred.target)],
                           res$avg.weight)[-1]
          )
          item[[pred.var]] = pred.target.warped[1:length(y.raw)]
        }
      }
      
      df.synth = item %>% 
        mutate(time = 1:length(y.raw) + start.time - 1,
               unit = unit)
      
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
    
    df.dependent = data %>% 
      filter(unit == dependent & time <= end.time) %>% 
      mutate(time = 1:length(y.raw) + start.time - 1,
             unit = dependent,
             value_warped = y.raw)
    df.synth = rbind(df.synth, df.dependent)
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
  
  
  ## Germany Reunification Data --------------------------------------------------
  data = readRDS("./data/repgermany.Rds")
  colnames(data)[1:4] = c("id", "unit", "time", "value")
  pred.vars = c("infrate", "trade", "schooling", "invest60", "invest70",
                "invest80", "industry")
  data = data %>% mutate(value_raw = value) %>% 
    group_by(unit) %>% 
    mutate(across(pred.vars, ~impute_values(.))) %>% 
    ungroup() %>% 
    data.frame
  
  # rescale
  df.rescale = data %>% 
    filter(time <= 1990) %>% 
    group_by(unit) %>% 
    summarise(value.min = min(value),
              value.max = max(value)) %>% 
    ungroup()
  
  mean.diff = mean(df.rescale$value.max - df.rescale$value.min)
  
  df.rescale = df.rescale %>% 
    mutate(
      multiplier = mean.diff/(value.max - value.min)
    )
  
  data = left_join(data, df.rescale, by = "unit")
  data = data %>% 
    mutate(
      value.bak = value_raw,
      value_raw = (value_raw - value.min)*multiplier,
      value = value_raw
    )
  
  # data list
  data.list = list(list(target = "West Germany",
                        data = data))
  ids = data$id %>%
    unique
  for (i in ids) {
    data.temp = data %>% filter(!(id %in% c(i, 7)))
    data.list = c(data.list, 
                  list(list(target = data.temp$unit[1],
                            data = data.temp)))
  }
  
  select2 = combn(setdiff(ids, 7), 2, simplify = TRUE)[,1:100]
  
  for (i in 1:ncol(select2)) {
    data.temp = data %>% filter(!(id %in% c(select2[, i], 7)))
    data.list = c(data.list, 
                  list(list(target = data.temp$unit[1],
                            data = data.temp)))
  }
  
  
  ## Grid Search Germany ---------------------------------------------------------
  # parameters
  filter.width.range = (1:9)*2+3
  k.range = 4:9
  step.pattern.range = list(
    # symmetricP0 = dtw::symmetricP0, # too bumpy
    # symmetricP05 = dtw::symmetricP05,
    symmetricP1 = dtw::symmetricP1,
    symmetricP2 = dtw::symmetricP2,
    # asymmetricP0 = dtw::asymmetricP0, # too bumpy
    # asymmetricP05 = dtw::asymmetricP05,
    asymmetricP1 = dtw::asymmetricP1,
    asymmetricP2 = dtw::asymmetricP2,
    typeIc = dtw::typeIc,
    # typeIcs = dtw::typeIcs,
    # typeIIc = dtw::typeIIc,  # jumps
    # typeIIIc = dtw::typeIIIc, # jumps
    # typeIVc = dtw::typeIVc,  # jumps
    typeId = dtw::typeId,
    # typeIds = dtw::typeIds,
    # typeIId = dtw::typeIId, # jumps
    mori2006 = dtw::mori2006
  )
  grid.search.parallel = TRUE
  
  
  args.TFDTW = list(buffer = 0, match.method = "fixed",
                    dist.quant = 0.95,
                    window.type = "none",
                    ## other
                    norm.method = "t",
                    step.pattern2 = dtw::asymmetricP2,
                    n.burn = 3, n.IQR = 3,
                    ma = 3, ma.na = "original",
                    default.margin = 3,
                    n.q = 1, n.r = 1)
  
  args.synth = list(predictors = expression(c(dep.var,"trade","infrate")),
                    special.predictors =
                      expression(list(list("industry", 1981:1989, c("mean")),
                                      list("schooling",c(1980,1985), c("mean")),
                                      list("invest80" ,1980, c("mean")))),
                    time.predictors.prior = 1981:1989,
                    time.optimize.ssr = 1960:1989)
  
  args.TFDTW.synth = list(start.time = 1960, end.time = 2003, treat.time = 1990,
                          args.TFDTW = args.TFDTW, args.synth = args.synth,
                          ## 2nd
                          n.mse = 10, pred.vars = pred.vars,
                          ## other
                          plot.figures = FALSE,
                          plot.path = "./figures/",
                          legend.pos = c(0.3, 0.8))
  
  # args.TFDTW.synth.all.units = list(target = data.list[[index]]$target,
  #                                   # data = data,
  #                                   args.TFDTW.synth = args.TFDTW.synth,
  #                                   ## 2nd
  #                                   detailed.output = TRUE,
  #                                   all.units.parallel = FALSE)
  # args.TFDTW.synth.all.units[["data"]] = data.list[[index]]$data
  # cat(paste0("Germany data set ", index, "..."))
  # results = SimDesign::quiet(
  #   grid.search(filter.width.range = filter.width.range,
  #               k.range = k.range,
  #               step.pattern.range = step.pattern.range,
  #               args.TFDTW.synth.all.units = args.TFDTW.synth.all.units,
  #               grid.search.parallel = grid.search.parallel)
  # )
  # saveRDS(results, paste0("./data/pred/germany/res_germany_", index, ".Rds"))
  # cat("Done.\n")
  # 
  # job.end = Sys.time()
  # print(job.end - job.start)
  
  for (index in (task_id*10+1):((task_id+1)*10)) {
    args.TFDTW.synth.all.units = list(target = data.list[[index]]$target,
                                      # data = data,
                                      args.TFDTW.synth = args.TFDTW.synth,
                                      ## 2nd
                                      detailed.output = TRUE,
                                      all.units.parallel = FALSE)
    args.TFDTW.synth.all.units[["data"]] = data.list[[index]]$data
    cat(paste0("Germany data set ", index, "..."))
    set.seed(20220407)
    results = SimDesign::quiet(
      grid.search(filter.width.range = filter.width.range,
                  k.range = k.range,
                  step.pattern.range = step.pattern.range,
                  args.TFDTW.synth.all.units = args.TFDTW.synth.all.units,
                  grid.search.parallel = grid.search.parallel)
    )
    
    saveRDS(results, paste0("./data/pred/germany/res_germany_", index, ".Rds"))
    cat("Done.\n")
  }
}









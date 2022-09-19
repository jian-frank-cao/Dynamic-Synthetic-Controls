## Functions -------------------------------------------------------------------
grid.search = function(filter.width.range, k.range, step.pattern.range,
                       args.TFDTW.synth.all.units,
                       grid.search.parallel = TRUE){
  # vanilla synthetic control
  data = args.TFDTW.synth.all.units[["data"]]
  units = data[c("id", "unit")] %>% distinct
  units.list = units %>% split(., seq(nrow(units)))
  
  if (grid.search.parallel) {
    fun.map = furrr::future_map
  }else{
    fun.map = purrr::map
  }
  
  args.synth = args.TFDTW.synth.all.units$args.TFDTW.synth$args.synth
  args.synth[["df"]] = data
  args.synth[["dep.var"]] = "value_raw"
  
  res.synth.raw.list = units.list %>% 
    set_names(units$unit) %>% 
    fun.map(
      ~{
        item = .
        dependent = item$unit
        dependent.id = item$id
        args.synth[["dependent.id"]] = dependent.id
        res = do.call(do.synth, args.synth)
      }
    )
  
  # grid search space
  search.grid = expand.grid(filter.width.range, k.range,
                             names(step.pattern.range)) %>% 
    `colnames<-`(c("filter.width", "k", "step.pattern"))
  search.grid.list = search.grid %>% split(., seq(nrow(search.grid)))
  
  results = search.grid.list %>% 
    fun.map(
      ~{
        task = .
        args.TFDTW.synth.all.units[["filter.width"]] = task$filter.width
        args.TFDTW.synth.all.units$args.TFDTW.synth$args.TFDTW[["k"]] = task$k
        args.TFDTW.synth.all.units$args.TFDTW.synth$args.TFDTW[["step.pattern1"]] =
          step.pattern.range[[task$step.pattern]]
        args.TFDTW.synth.all.units[["res.synth.raw.list"]] = res.synth.raw.list
        do.call(TFDTW.synth.all.units, args.TFDTW.synth.all.units)
      }
    )
  
  return(results)
}



results = NULL
folder = "./data/placebo/basque/"
res.files = list.files(folder)
for (res.file in res.files) {
  results = c(results, list(readRDS(paste0(folder, res.file))))
}

df = future_map2(
  results,
  as.list(1:length(results)),
  ~{
    item = .x
    index = .y
    
    mse = future_map2(
      item,
      as.list(names(item)),
      ~{
        data = .x[["results.TFDTW.synth"]]
        id = .y
        mse = data %>% 
          map(
            ~{
              task = .
              data.frame(unit = task$dependent,
                         gap.raw = task$gap.raw,
                         gap.TFDTW = task$gap.TFDTW)
            }
          ) %>% do.call("rbind", .)
        mse %>% mutate(id = id)
      }
    ) %>% do.call("rbind", .)
    
    mse = mse %>% 
      group_by(unit)

  }
) %>% do.call("rbind", .)
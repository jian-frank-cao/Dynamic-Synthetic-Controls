## Setup -----------------------------------------------------------------------
library(checkpoint)
checkpoint("2022-04-01")

library(parallel)
n.cores = detectCores()
library(tidyverse)
library(furrr)
plan(multisession, workers = n.cores - 1)
options(future.rng.onMisuse="ignore")
options(stringsAsFactors = FALSE)

source("./R/utility/misc.R")
source("./R/utility/TFDTW.R")
source("./R/utility/synth.R")
source("./R/utility/implement.R")
source("./R/utility/grid.search.R")
set.seed(20220407)


## Basque Terrorism Data -------------------------------------------------------
data(basque, package = "Synth")
data = basque
colnames(data)[1:4] = c("id", "unit", "time", "value")
data = data %>% mutate(invest_ratio = invest/value,
                       value_raw = value)

# rescale
df.rescale = data %>%
  filter(time <= 1970) %>%
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
data.list = list(list(target = "Basque Country (Pais Vasco)",
                      data = data))
ids = data$id %>%
  unique
for (i in ids) {
  data.temp = data %>% filter(!(id %in% c(i, 17)))
  data.list = c(data.list, 
                list(list(target = data.temp$unit[1],
                          data = data.temp)))
}

select2 = combn(setdiff(ids, 17), 2, simplify = TRUE)[,1:100]

for (i in 1:ncol(select2)) {
  data.temp = data %>% filter(!(id %in% c(select2[, i], 17)))
  data.list = c(data.list, 
                list(list(target = data.temp$unit[1],
                          data = data.temp)))
}


## Placebo ---------------------------------------------------------------------
folder = "./data/placebo/basque/"
file.list = as.list(list.files(folder))

pre.start = 7
pre.end = 16
post.start = 17
post.end = 26

# df.mse
df.mse = file.list %>% 
  future_map(
  ~{
    file.name = .
    data.id = strsplit(file.name, "_")[[1]][4]
    data.id = as.numeric(strsplit(data.id, ".R")[[1]][1])
    data.list = readRDS(paste0(folder, file.name))
    mse = future_map2(
      data.list,
      as.list(names(data.list)),
      ~{
        result.synth = .x[["results.TFDTW.synth"]]
        grid.id = .y
        mse = result.synth %>% 
          map(
            ~{
              task = .
              unit = task$dependent
              gap.raw = task$gap.raw
              gap.TFDTW = task$gap.TFDTW
              data.frame(unit = unit,
                         mse.preT.raw = mean(gap.raw[pre.start:pre.end]^2, na.rm = T),
                         mse.preT.TFDTW = mean(gap.TFDTW[pre.start:pre.end]^2, na.rm = T),
                         mse.postT.raw = mean(gap.raw[post.start:post.end]^2, na.rm = T),
                         mse.postT.TFDTW = mean(gap.TFDTW[post.start:post.end]^2, na.rm = T))
            }
          ) %>% do.call("rbind", .)
        mse %>% mutate(grid.id = grid.id)
      }
    ) %>% do.call("rbind", .)
    mse$data.id = data.id
    
    mse %>% 
      group_by(unit) %>% 
      top_n(-1, mse.preT.TFDTW) %>% 
      top_n(-1, grid.id)
  }
) %>% do.call("rbind", .)

saveRDS(df.mse, "./data/df.mse_basque.x.Rds")


## Placebo ---------------------------------------------------------------------
df.mse = readRDS("./data/df.mse_basque.x.Rds")
df.grid = 


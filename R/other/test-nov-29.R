## Setup -----------------------------------------------------------------------
library(checkpoint)
checkpoint("2022-04-01")

library(tidyverse)
library(furrr)
plan(multisession, workers = 8)
options(future.rng.onMisuse="ignore")
options(stringsAsFactors = FALSE)
source("./R/misc.R")
source("./R/TFDTW.R")
source("./R/synth.R")
source("./R/implement.R")
source("./R/simulate.R")
source("./R/grid.search.R")
set.seed(20220407)


# Data Simulation -------------------------------------------------------------
n.simulation = 100
length = 100
n = 10

# generate sobol sequence
# sobol.seq = qrng::sobol(n.simulation*1, d = n - 1, randomize = "Owen",
#                         seed = 20220401, skip = 100)
# rnd.speeds = cbind(rep(0.5, n.simulation), sobol.seq*0.5 + 0.1)


# simulate
data.list = NULL
for (i in 1:n.simulation) {
  data.list[[i]] = sim.data3(n = n, length = length,
                            t.treat = 60, shock = 0, ar.x = 0.6,
                            n.SMA = 1, n.diff = 1,
                            speed.upper = 2,
                            speed.lower = 0.5,
                            reweight = TRUE,
                            rescale = TRUE,
                            rescale.multiplier = 20,
                            beta = 0.5)
}

data.list[[12]] %>% ggplot(aes(x = time, y = value, color = unit)) +
  geom_line() +
  geom_vline(xintercept = 60, linetype="dashed")
saveRDS(data.list, "./data/simul_data_beta_103.Rds")

## Run -------------------------------------------------------------------------
data.list = readRDS(paste0("./data/simul_data_beta_103.Rds"))

args.synth = list(predictors = NULL,
                  special.predictors = 
                    expression(list(list(dep.var, 50:59, c("mean")),
                                    list(dep.var, 40:49, c("mean")),
                                    list(dep.var, 30:39, c("mean")))),
                  time.predictors.prior = 1:59,
                  time.optimize.ssr = 1:59)

results = data.list %>% 
  future_map(
    ~{
      data = .
      units = data[c("id", "unit")] %>% distinct
      units.list = units %>% split(., seq(nrow(units)))
      args.synth[["df"]] = data
      args.synth[["dep.var"]] = "value_raw"
      res.synth.raw.list = units.list %>% 
        set_names(units$unit) %>% 
        map(
          ~{
            item = .
            dependent.id = item$id
            args.synth[["dependent.id"]] = dependent.id
            res = do.call(do.synth, args.synth)
          }
        )
      res.synth.raw.list
    }
  )


## Results ---------------------------------------------------------------------
units = LETTERS[1:10]
for (UNIT in units) {
  res = results %>% 
    future_map(
      ~{
        item = .
        res.synth = map2(
          item,
          names(item),
          ~{
            result = .x
            unit = .y
            data.frame(unit = unit,
                       time = 1:100,
                       value = result$value,
                       original = result$synthetic)
          }
        ) %>% 
          do.call("rbind", .)
        res.synth = res.synth %>% 
          mutate(gap.original = value - original)
        gap.original = res.synth %>% filter(unit == UNIT) %>% .[["gap.original"]]
        res.synth = res.synth %>% 
          filter(unit != UNIT)
        df.quant = res.synth %>% 
          group_by(time) %>% 
          summarise(
            ci_origin_upper = quantile(gap.original, 0.975, na.rm = T),
            ci_origin_lower = quantile(gap.original, 0.025, na.rm = T),
            ci_new_upper = quantile(gap.TFDTW, 0.975, na.rm = T),
            ci_new_lower = quantile(gap.TFDTW, 0.025, na.rm = T)
          )
        df.quant$gap.original = gap.original
        df.quant = df.quant %>% 
          mutate(sig.original = ifelse(gap.original > ci_origin_upper | gap.original < ci_origin_lower, 1, 0),
          ) %>% 
          data.frame(.) %>% 
          `rownames<-`(1:100)
        df.quant
      }
    ) %>% 
    do.call("rbind", .)
  
  res = res %>% filter(time > 60 & time < 71)
  res = res[complete.cases(res),7]
  
  res.boot = NULL
  for (i in 1:1000) {
    ind = sample(1:length(res), 1000, replace = TRUE)
    boot = res[ind]
    res.boot = c(res.boot, mean(boot, na.rm = TRUE))
  }
  
  print(UNIT)
  print(mean(res.boot))
}


item = results[[14]]

res.synth = map2(
  item,
  names(item),
  ~{
    result = .x
    unit = .y
    data.frame(unit = unit,
               time = 1:100,
               value = result$value,
               original = result$synthetic)
  }
) %>% 
  do.call("rbind", .)
res.synth = res.synth %>% 
  mutate(gap.original = value - original)
gap.original = res.synth %>% filter(unit == UNIT) %>% .[["gap.original"]]
res.synth2 = res.synth
res.synth = res.synth %>% 
  filter(unit != UNIT)
df.quant = res.synth %>% 
  group_by(time) %>% 
  summarise(
    ci_origin_upper = quantile(gap.original, 0.975, na.rm = T),
    ci_origin_lower = quantile(gap.original, 0.025, na.rm = T),
    ci_new_upper = quantile(gap.TFDTW, 0.975, na.rm = T),
    ci_new_lower = quantile(gap.TFDTW, 0.025, na.rm = T)
  )
df.quant$gap.original = gap.original
df.quant = df.quant %>% 
  mutate(sig.original = ifelse(gap.original > ci_origin_upper | gap.original < ci_origin_lower, 1, 0),
  ) %>% 
  data.frame(.) %>% 
  `rownames<-`(1:100)

UNIT = "C"
colors = c()
for (i in units) {
  if (i == UNIT) {
    colors[i] = "red"
  }else{
    colors[i] = "grey70"
  }
}

res.synth2 %>% 
  ggplot(aes(x = time, y = gap.original, color = unit)) +
  geom_line() +
  scale_color_manual(name = NULL, values = colors)

mean(df.quant$sig.original)

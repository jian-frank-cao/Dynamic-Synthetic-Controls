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


## Grid Search Basque ----------------------------------------------------------
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

args.synth = list(predictors = NULL,
                  special.predictors =
                    expression(list(
                      list(dep.var, 1960:1969, c("mean")),
                      list("invest_ratio", 1964:1969, c("mean")),
                      list("popdens", 1969, c("mean")),
                      list("sec.agriculture", 1961:1969, c("mean")),
                      list("sec.energy", 1961:1969, c("mean")),
                      list("sec.industry", 1961:1969, c("mean")),
                      list("sec.construction", 1961:1969, c("mean")),
                      list("sec.services.venta", 1961:1969, c("mean")),
                      list("sec.services.nonventa", 1961:1969, c("mean")),
                      list("school.illit", 1964:1969, c("mean")),
                      list("school.prim", 1964:1969, c("mean")),
                      list("school.med", 1964:1969, c("mean")),
                      list("school.high", 1964:1969, c("mean")),
                      list("school.post.high", 1964:1969, c("mean"))
                    )),
                  time.predictors.prior = 1955:1969,
                  time.optimize.ssr = 1955:1969)

args.TFDTW.synth = list(start.time = 1955, end.time = 1997, treat.time = 1970,
                        args.TFDTW = args.TFDTW, args.synth = args.synth,
                        ## 2nd
                        n.mse = 10,
                        ## other
                        plot.figures = FALSE,
                        plot.path = "./figures/",
                        legend.pos = c(0.3, 0.3))

for (index in 1:length(data.list)) {
  args.TFDTW.synth.all.units = list(target = data.list[[index]]$target,
                                    # data = data,
                                    args.TFDTW.synth = args.TFDTW.synth,
                                    ## 2nd
                                    detailed.output = TRUE,
                                    all.units.parallel = FALSE)
  
  args.TFDTW.synth.all.units[["data"]] = data.list[[index]]$data
  results = SimDesign::quiet(
    grid.search(filter.width.range = filter.width.range,
                k.range = k.range,
                step.pattern.range = step.pattern.range,
                args.TFDTW.synth.all.units = args.TFDTW.synth.all.units,
                grid.search.parallel = grid.search.parallel)
  )
  
  saveRDS(results, paste0("./data/placebo/basque/res_basque_", index, ".Rds"))
  print(index)
}


## Placebo ---------------------------------------------------------------------
folder = "./data/placebo/basque/"
file.list = as.list(list.files(folder))

pre.start = 7
pre.end = 16
post.start = 17
post.end = 26

# df.mse
df.mse = future_map2(
  file.list,
  as.list(1:length(file.list)),
  ~{
    file.name = .x
    data.id = .y
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
              ###
              scales = df.rescale %>% filter(unit == task$dependent)
              value.min = scales$value.min
              multiplier = scales$multiplier
              value.true = task$res.synth.TFDTW$value/multiplier + value.min
              value.sc = task$res.synth.raw$synthetic/multiplier + value.min
              value.dsc = task$res.synth.TFDTW$synthetic/multiplier + value.min
              gap.raw = value.true - value.true
              gap.TFDTW = value.dsc - value.sc
              ###
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

saveRDS(df.mse, "./data/df.mse_basque.Rds")

# t.test for log(MSEdsc/MSEsc)
df.mse = readRDS("./data/df.mse_basque.Rds")
df.mse = df.mse %>% 
  mutate(log.ratio = log(mse.postT.TFDTW/mse.postT.raw),
         data.id = as.character(data.id))

res.aov = aov(log.ratio ~ data.id*unit, df.mse)
summary.aov = summary(res.aov)

BMS = summary.aov[[1]]$`Mean Sq`[1]
JMS = summary.aov[[1]]$`Mean Sq`[2]
EMS = summary.aov[[1]]$`Mean Sq`[3]
n = summary.aov[[1]]$`Df`[1] + 1
k = summary.aov[[1]]$`Df`[2] + 1
res.icc = (BMS - EMS)/(BMS + (k - 1)*EMS + k*(JMS - EMS)/n)
res.vif = 1 + (k - 1)*res.icc
DF = nrow(df.mse)/res.vif

t.value = t.test(df.mse$log.ratio)$statistic
p.value = pt(t.value, df = DF, lower.tail = TRUE)*2


## Plot results ----------------------------------------------------------------
# df.target
results.target = readRDS("./data/placebo/basque/res_basque_1.Rds")
target = "Basque Country (Pais Vasco)"

pre.start = 7
pre.end = 16
post.start = 17
post.end = 26

mse = future_map2(
  results.target,
  as.list(names(results.target)),
  ~{
    result.synth = .x[["results.TFDTW.synth"]][[target]]
    grid.id = .y
    gap.raw = result.synth$gap.raw
    gap.TFDTW = result.synth$gap.TFDTW
    data.frame(grid.id = grid.id,
               mse.preT.raw = mean(gap.raw[pre.start:pre.end]^2, na.rm = T),
               mse.preT.TFDTW = mean(gap.TFDTW[pre.start:pre.end]^2, na.rm = T),
               mse.postT.raw = mean(gap.raw[post.start:post.end]^2, na.rm = T),
               mse.postT.TFDTW = mean(gap.TFDTW[post.start:post.end]^2, na.rm = T))
  }
) %>% do.call("rbind", .)

opt.grid.id = mse %>% 
  top_n(-1, mse.preT.TFDTW) %>% 
  top_n(-1, grid.id) %>% 
  .[["grid.id"]]

df.target = data.frame(
  time = 1955:1997,
  unit = target,
  data.id = 0,
  grid.id = opt.grid.id,
  value = results.target[[opt.grid.id]][[4]][[target]][[3]][["value"]],
  synth.sc = results.target[[opt.grid.id]][[4]][[target]][[3]][["synthetic"]],
  synth.dsc = results.target[[opt.grid.id]][[4]][[target]][[4]][["synthetic"]]
)

df.target = df.target %>% 
  mutate(
    gap.sc = value - synth.sc,
    gap.dsc = value - synth.dsc,
    group = paste0(data.id, "-", grid.id, "-", unit)
  )

saveRDS(df.target, "./data/df.target_basque.Rds")


# df.gap
df.mse = readRDS("./data/df.mse_basque.Rds")
folder = "./data/placebo/basque/"
file.list = as.list(list.files(folder))
results = file.list %>%
  future_map(
    ~{
      file.name = .
      readRDS(paste0(folder, file.name))
    }
  )


df.gap = NULL
for (i in 1:nrow(df.mse)) {
  unit = df.mse$unit[i]
  data.id = df.mse$data.id[i]
  grid.id = df.mse$grid.id[i]
  df.gap[[i]] = data.frame(
    time = 1955:1997,
    unit = unit,
    data.id = data.id,
    grid.id = grid.id,
    value = results[[data.id]][[grid.id]][[4]][[unit]][[3]][["value"]],
    synth.sc = results[[data.id]][[grid.id]][[4]][[unit]][[3]][["synthetic"]],
    synth.dsc = results[[data.id]][[grid.id]][[4]][[unit]][[4]][["synthetic"]]
  )
  print(i)
}

df.gap = df.gap %>%
  do.call("rbind", .) %>%
  mutate(
    gap.sc = value - synth.sc,
    gap.dsc = value - synth.dsc,
    group = paste0(data.id, "-", grid.id, "-", unit)
  )

saveRDS(df.gap, "./data/df.gap_basque.Rds")

avg.mse.sc = mean(df.mse$mse.postT.raw, na.rm = TRUE)
avg.mse.dsc = mean(df.mse$mse.postT.TFDTW, na.rm = TRUE)

# plot
df.target = readRDS("./data/df.target_basque.Rds")
df.gap = readRDS("./data/df.gap_basque.Rds")

df.quantile = df.gap %>%
  group_by(time) %>%
  summarise(quantile.sc.975 = quantile(gap.sc, 0.975, na.rm = T),
            quantile.sc.025 = quantile(gap.sc, 0.025, na.rm = T),
            quantile.dsc.975 = quantile(gap.dsc, 0.975, na.rm = T),
            quantile.dsc.025 = quantile(gap.dsc, 0.025, na.rm = T)) %>% 
  mutate(group = "quantile")

color.sc = "#2ab7ca"
color.dsc = "#fe4a49"
# color.sc = "grey70"
# color.dsc = "grey30"
  
colors = c("TE (SC)" = color.sc,
           "TE (DSC)" = color.dsc)

fills = c("95% Quantile (SC)" = color.sc,
          "95% Quantile (DSC)" = color.dsc)

set.seed(20230812)
group.sample = sample(unique(df.gap$group), 100)
  
fig_basque = df.gap %>%
  filter(group %in% group.sample) %>% 
  ggplot(aes(x = time, group = group)) +
  annotate("rect", xmin = 1970, xmax = 1980,
           ymin = -3, ymax = 3, alpha = .3) +
  geom_line(aes(y = gap.sc), col = color.sc, alpha = 0.4) +
  geom_line(aes(y = gap.dsc), col = color.dsc, alpha = 0.4) +
  geom_ribbon(aes(ymin = quantile.sc.025,
                  ymax = quantile.sc.975,
                  fill="95% Quantile (SC)"),
              data = df.quantile, alpha=0.5) +
  geom_ribbon(aes(ymin = quantile.dsc.025,
                  ymax = quantile.dsc.975,
                  fill="95% Quantile (DSC)"),
              data = df.quantile, alpha=0.5) +
  geom_line(aes(y = gap.sc, color = "TE (SC)"),
            data = df.target, size = 1) +
  geom_line(aes(y = gap.dsc, color = "TE (DSC)"), 
            data = df.target, size = 1) +
  scale_color_manual(name = NULL, values = colors) +
  scale_fill_manual(name = NULL, values = fills) +
  geom_vline(xintercept = 1970, linetype="dashed", col = "grey20") +
  geom_hline(yintercept = 0, linetype="dashed", col = "grey20") +
  annotate("text", x = 1975, y = 0.85,
           label = "t = -7.6121\nP < 0.0001", col = "grey20") +
  annotate("text", x = 1969, y = 0.6, angle = 90,
           label = "Treatment", col = "grey20") +
  # annotate("text", x = 1960, y = -0.5, label = "bar(MSE)[SC]==0.0469", parse = TRUE,
  #          col = "grey20", size = 4, fontface = "bold") +
  # annotate("text", x = 1960, y = -0.8, label = "bar(MSE)[DSC]==0.032", parse = TRUE,
  #          col = "grey20", size = 4, fontface = "bold") +
  coord_cartesian(xlim=c(1950, 1990), ylim = c(-1, 1)) +
  xlab("Year") +
  ylab("TE(y - Synthetic Control)") +
  theme_bw() +
  theme(legend.position = "none",
        legend.box = "horizontal",
        legend.background = element_rect(fill=NA))
  
saveRDS(fig_basque, "./data/placebo_basque.Rds")








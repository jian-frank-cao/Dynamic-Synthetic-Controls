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


## Germany Reunification Data --------------------------------------------------
data = readRDS("./data/repgermany.Rds")
colnames(data)[1:4] = c("id", "unit", "time", "value")
data = data %>% mutate(value_raw = value)

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

# ## Grid Search Germany ---------------------------------------------------------
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
                        n.mse = 10,
                        ## other
                        plot.figures = FALSE,
                        plot.path = "./figures/",
                        legend.pos = c(0.3, 0.8))

args.TFDTW.synth.all.units = list(target = data.list[[index]]$target,
                                  # data = data,
                                  args.TFDTW.synth = args.TFDTW.synth,
                                  ## 2nd
                                  detailed.output = TRUE,
                                  all.units.parallel = FALSE)

for (index in 1:length(data.list)) {
  args.TFDTW.synth.all.units[["data"]] = data.list[[index]]$data
  results = SimDesign::quiet(
    grid.search(filter.width.range = filter.width.range,
                k.range = k.range,
                step.pattern.range = step.pattern.range,
                args.TFDTW.synth.all.units = args.TFDTW.synth.all.units,
                grid.search.parallel = grid.search.parallel)
  )
  
  saveRDS(results, paste0("./data/placebo/germany/res_germany_", index, ".Rds"))
  print(index)
}


## Placebo ---------------------------------------------------------------------
folder = "./data/placebo/germany/"
file.list = as.list(list.files(folder))

pre.start = 22
pre.end = 31
post.start = 32
post.end = 41

# de.mse
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
              # value.min = scales$value.min
              multiplier = scales$multiplier
              gap.raw = task$gap.raw/multiplier
              gap.TFDTW = task$gap.TFDTW/multiplier
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

saveRDS(df.mse, "./data/df.mse_germany2.Rds")

# t.test for log(MSEdsc/MSEsc)
df.mse = readRDS("./data/df.mse_germany2.Rds")
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
results.target = readRDS("./data/res_germany_1204_1.Rds")
target = "West Germany"
scales = df.rescale %>% filter(unit == target)
multiplier = scales$multiplier

pre.start = 22
pre.end = 31
post.start = 32
post.end = 41

mse = future_map2(
  results.target,
  as.list(names(results.target)),
  ~{
    result.synth = .x[["results.TFDTW.synth"]][[target]]
    grid.id = .y
    gap.raw = result.synth$gap.raw/multiplier
    gap.TFDTW = result.synth$gap.TFDTW/multiplier
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
  time = 1960:2003,
  unit = target,
  data.id = 0,
  grid.id = opt.grid.id,
  value = results.target[[opt.grid.id]][[4]][[target]][[3]][["value"]],
  synth.sc = results.target[[opt.grid.id]][[4]][[target]][[3]][["synthetic"]],
  synth.dsc = results.target[[opt.grid.id]][[4]][[target]][[4]][["synthetic"]]
)

df.target = df.target %>% 
  mutate(
    gap.sc = (value - synth.sc)/multiplier,
    gap.dsc = (value - synth.dsc)/multiplier,
    group = paste0(data.id, "-", grid.id, "-", unit)
  )

saveRDS(df.target, "./data/df.target_germany2.Rds")


# df.gap
df.mse = readRDS("./data/df.mse_germany2.Rds")
folder = "./data/placebo/germany/"
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
  data.id = as.numeric(df.mse$data.id[i])
  grid.id = df.mse$grid.id[i]
  scales = df.rescale %>% filter(unit == df.mse$unit[i])
  value.min = scales$value.min
  multiplier = scales$multiplier
  value = results[[data.id]][[grid.id]][[4]][[unit]][[3]][["value"]]
  synth.sc = results[[data.id]][[grid.id]][[4]][[unit]][[3]][["synthetic"]]
  synth.dsc = results[[data.id]][[grid.id]][[4]][[unit]][[4]][["synthetic"]]
  
  df.gap[[i]] = data.frame(
    time = 1960:2003,
    unit = unit,
    data.id = data.id,
    grid.id = grid.id,
    value = value/multiplier + value.min,
    synth.sc = synth.sc/multiplier + value.min,
    synth.dsc = synth.dsc/multiplier + value.min
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

saveRDS(df.gap, "./data/df.gap_germany2.Rds")

avg.log.mse.sc = mean(log(df.mse$mse.postT.raw), na.rm = TRUE)
avg.log.mse.dsc = mean(log(df.mse$mse.postT.TFDTW), na.rm = TRUE)

# plot
set.seed(20230901)
df.target = readRDS("./data/df.target_germany2.Rds")
df.gap = readRDS("./data/df.gap_germany2.Rds")

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

colors = c("Target TE (SC)" = color.sc,
           "Target TE (DSC)" = color.dsc)

fills = c("95% Quantile (SC)" = color.sc,
          "95% Quantile (DSC)" = color.dsc)

set.seed(20230812)
group.sample = sample(unique(df.gap$group), 100)

fig_germany = df.gap %>%
  filter(group %in% group.sample) %>% 
  ggplot(aes(x = time, group = group)) +
  annotate("rect", xmin = 1990, xmax = 2000, 
           ymin = -9000, ymax = 9000, alpha = .3) +
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
  geom_line(aes(y = gap.sc, color = "Target TE (SC)"),
            data = df.target, size = 1) +
  geom_line(aes(y = gap.dsc, color = "Target TE (DSC)"), 
            data = df.target, size = 1) +
  scale_color_manual(name = NULL, values = colors) +
  scale_fill_manual(name = NULL, values = fills) +
  geom_vline(xintercept = 1990, linetype="dashed", col = "grey20") +
  geom_hline(yintercept = 0, linetype="dashed", col = "grey20") +
  annotate("text", x = 1995, y = 6800,
           label = "t = -4.0951\nP < 0.0001", col = "grey20") +
  annotate("text", x = 1989, y = 4800, angle = 90,
           label = "Treatment", col = "grey20") +
  # annotate("text", x = 2007, y = -3000, label = "bar(log(MSE))[SC]==13.90", parse = TRUE,
  #          col = "grey20", size = 4, fontface = "bold") +
  # annotate("text", x = 2007, y = -5000, label = "bar(log(MSE))[DSC]==13.79", parse = TRUE,
  #          col = "grey20", size = 4, fontface = "bold") +
  coord_cartesian(xlim=c(1970, 2010),ylim=c(-8000,8000)) +
  xlab("Year") +
  ylab("Treatment Effect") +
  theme_bw() +
  theme(legend.position = c(0.25, 0.2),
        legend.box = "horizontal",
        legend.background = element_rect(fill=NA))

saveRDS(fig_germany, "./data/placebo_germany.Rds")
  


## Setup -----------------------------------------------------------------------
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
  
  # args.synth[["df"]] = df.synth
  # args.synth[["dependent.id"]] = dependent.id
  # 
  # # Synth
  # if (is.null(res.synth.raw)) {
  #   args.synth[["dep.var"]] = "value_raw"
  #   res.synth.raw = do.call(do.synth, args.synth)
  # }
  # 
  # args.synth[["dep.var"]] = "value_warped"
  # res.synth.TFDTW = do.call(do.synth, args.synth)
  
  args.synth[["dependent.id"]] = dependent.id
  
  # Synth
  if (is.null(res.synth.raw)) {
    args.synth[["df"]] = data
    args.synth[["dep.var"]] = "value_raw"
    res.synth.raw = do.call(do.synth, args.synth)
  }
  
  args.synth[["df"]] = df.synth
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

search.grid = expand.grid(filter.width.range, k.range,
                          names(step.pattern.range)) %>% 
  `colnames<-`(c("filter.width", "k", "step.pattern"))
grid.opt = readRDS("./data/Figure_A3_3_gridOpt.Rds")

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

results = as.list(1:nrow(grid.opt)) %>%
  future_map(
    ~{
      i = .
      set.seed(20220407)
      data.id = as.numeric(grid.opt$data.id[i])
      grid.id = as.numeric(grid.opt$grid.id[i])
      unit = grid.opt$unit[i]
      id = grid.opt$id[i]
      args.TFDTW.synth.target.only = list(target = unit, id = id,
                                          data = data.list[[data.id]][["data"]],
                                          args.TFDTW.synth = args.TFDTW.synth)
      SimDesign::quiet(
        grid.search.opt(filter.width.range = search.grid[grid.id,]$filter.width,
                        k.range = search.grid[grid.id,]$k,
                        step.pattern.range =
                          step.pattern.range[search.grid[grid.id,]$step.pattern],
                        args.TFDTW.synth.target.only = args.TFDTW.synth.target.only,
                        grid.search.parallel = grid.search.parallel)
      )
    }
  )

set.seed(20220407)
data.id = 1
grid.id = 109
unit = "West Germany"
id = 7
args.TFDTW.synth.target.only = list(target = unit, id = id,
                                    data = data.list[[data.id]][["data"]],
                                    args.TFDTW.synth = args.TFDTW.synth)
result_target = SimDesign::quiet(
  grid.search.opt(filter.width.range = search.grid[grid.id,]$filter.width,
                  k.range = search.grid[grid.id,]$k,
                  step.pattern.range =
                    step.pattern.range[search.grid[grid.id,]$step.pattern],
                  args.TFDTW.synth.target.only = args.TFDTW.synth.target.only,
                  grid.search.parallel = grid.search.parallel)
)


## Placebo ---------------------------------------------------------------------
pre.start = 22
pre.end = 31
post.start = 32
post.end = 41

# de.mse
df.mse = lapply(results, "[[", 1) %>%
  future_map(
    ~{
      task = .[["results.TFDTW.synth"]]
      unit = task$dependent
      scales = df.rescale %>% filter(unit == task$dependent)
      multiplier = scales$multiplier
      gap.raw = task$gap.raw/multiplier
      gap.TFDTW = task$gap.TFDTW/multiplier
      data.frame(unit = unit,
                 mse.preT.raw = mean(gap.raw[pre.start:pre.end]^2, na.rm = T),
                 mse.preT.TFDTW = mean(gap.TFDTW[pre.start:pre.end]^2, na.rm = T),
                 mse.postT.raw = mean(gap.raw[post.start:post.end]^2, na.rm = T),
                 mse.postT.TFDTW = mean(gap.TFDTW[post.start:post.end]^2, na.rm = T))
    }
  ) %>% do.call("rbind", .)

df.mse = cbind(df.mse, grid.opt[,c(2,3)])

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
results.target = result_target[[1]]
target = "West Germany"
scales = df.rescale %>% filter(unit == target)
multiplier = scales$multiplier

pre.start = 22
pre.end = 31
post.start = 32
post.end = 41

df.target = data.frame(
  time = 1960:2003,
  unit = target,
  data.id = 1,
  grid.id = 109,
  value = results.target[["results.TFDTW.synth"]][[3]][["value"]],
  synth.sc = results.target[["results.TFDTW.synth"]][[3]][["synthetic"]],
  synth.dsc = results.target[["results.TFDTW.synth"]][[4]][["synthetic"]]
)

df.target = df.target %>% 
  mutate(
    gap.sc = (value - synth.sc)/multiplier,
    gap.dsc = (value - synth.dsc)/multiplier,
    group = paste0(data.id, "-", grid.id, "-", unit)
  )

# df.gap
df.gap = NULL
for (i in 1:nrow(df.mse)) {
  unit = df.mse$unit[i]
  data.id = df.mse$data.id[i]
  grid.id = df.mse$grid.id[i]
  scales = df.rescale %>% filter(unit == df.mse$unit[i])
  value.min = scales$value.min
  multiplier = scales$multiplier
  value = results[[i]][[1]][["results.TFDTW.synth"]][[3]][["value"]]
  synth.sc = results[[i]][[1]][["results.TFDTW.synth"]][[3]][["synthetic"]]
  synth.dsc = results[[i]][[1]][["results.TFDTW.synth"]][[4]][["synthetic"]]
  
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

avg.log.mse.sc = mean(log(df.mse$mse.postT.raw), na.rm = TRUE)
avg.log.mse.dsc = mean(log(df.mse$mse.postT.TFDTW), na.rm = TRUE)

# plot
df.quantile = df.gap %>%
  group_by(time) %>%
  summarise(quantile.sc.975 = quantile(gap.sc, 0.975, na.rm = T),
            quantile.sc.025 = quantile(gap.sc, 0.025, na.rm = T),
            quantile.dsc.975 = quantile(gap.dsc, 0.975, na.rm = T),
            quantile.dsc.025 = quantile(gap.dsc, 0.025, na.rm = T)) %>% 
  mutate(group = "quantile")

color.sc = "#2ab7ca"
color.dsc = "#fe4a49"

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
           label = "t = -5.4611\nP < 0.0001", col = "grey20") +
  annotate("text", x = 1989, y = 4800, angle = 90,
           label = "Treatment", col = "grey20") +
  coord_cartesian(xlim=c(1970, 2010),ylim=c(-8000,8000)) +
  xlab("Year") +
  ylab("Treatment Effect") +
  theme_bw() +
  theme(legend.position = c(0.25, 0.2),
        legend.box = "horizontal",
        legend.background = element_rect(fill=NA))

saveRDS(fig_germany, "./data/Figure_A3_3.Rds")



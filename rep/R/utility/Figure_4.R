## Setup -----------------------------------------------------------------------
library(parallel)
n.cores = detectCores()
library(tidyverse)
library(furrr)
plan(multisession, workers = n.cores)
options(future.rng.onMisuse="ignore")
options(stringsAsFactors = FALSE)
source("./code/utility/misc.R")
source("./code/utility/TFDTW.R")
source("./code/utility/synth.R")
source("./code/utility/implement.R")
source("./code/utility/simulate.R")
source("./code/utility/grid.search.R")
set.seed(20220407)


## Data Simulation -------------------------------------------------------------
n.simulation = 100
length = 100
n = 10
beta = 1
shock = 10

# simulate
data.list = NULL
for (i in 1:n.simulation) {
  data.list[[i]] = sim.data(n = n, length = length,
                            t.treat = 60, shock = shock, ar.x = 0.6,
                            n.SMA = 1, n.diff = 1,
                            speed.upper = 2,
                            speed.lower = 0.5,
                            reweight = TRUE,
                            rescale = TRUE,
                            rescale.multiplier = 20,
                            beta = beta)
}


## Run -------------------------------------------------------------------------
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
grid.opt = readRDS("./data/Figure_4_gridOpt.Rds")

args.TFDTW = list(buffer = 20, match.method = "open.end",
                  dist.quant = 0.95, 
                  window.type = "sakoechiba",
                  ## other
                  norm.method = "t",
                  step.pattern2 = dtw::asymmetricP2,
                  n.burn = 3, n.IQR = 3,
                  ma = 3, ma.na = "original",
                  default.margin = 3,
                  n.q = 1, n.r = 1)

args.synth = list(predictors = NULL,
                  special.predictors = 
                    expression(list(list(dep.var, 50:59, c("mean")),
                                    list(dep.var, 40:49, c("mean")),
                                    list(dep.var, 30:39, c("mean")))),
                  time.predictors.prior = 1:59,
                  time.optimize.ssr = 1:59)

args.TFDTW.synth = list(start.time = 1, end.time = 100, treat.time = 60,
                        args.TFDTW = args.TFDTW, args.synth = args.synth,
                        ## 2nd
                        n.mse = 10, 
                        ## other
                        plot.figures = FALSE,
                        plot.path = "./results/",
                        legend.pos = c(0.3, 0.7))

args.TFDTW.synth.target.only = list(target = "A", id = 1,
                                    args.TFDTW.synth = args.TFDTW.synth)


results = as.list(1:nrow(grid.opt)) %>%
  future_map(
    ~{
      i = .
      set.seed(20220407)
      args.TFDTW.synth.target.only[["data"]] = data.list[[grid.opt$data.id[i]]]
      SimDesign::quiet(
        grid.search.opt(filter.width.range = 
                          search.grid[grid.opt$grid.id[i],]$filter.width,
                        k.range = search.grid[grid.opt$grid.id[i],]$k,
                        step.pattern.range =
                          step.pattern.range[search.grid[grid.opt$grid.id[i],]$step.pattern],
                        args.TFDTW.synth.target.only = args.TFDTW.synth.target.only,
                        grid.search.parallel = grid.search.parallel)
      )
    }
  )


## Test result -----------------------------------------------------------------
length = 100
treat_time = 60
n_mse = 10
treatment = c(rep(0, treat_time),
              seq(0, shock, length.out = round(0.1*length)),
              rep(shock, round(0.9*length - treat_time)))

pre.start = 51
pre.end = 60
post.start = 61
post.end = 70

# de.mse
df.mse = lapply(results, "[[", 1) %>%
  future_map(
    ~{
      result = .
      value = result$res.synth.target.raw$value
      synth.raw = result$res.synth.target.raw$synthetic
      synth.TFDTW = result$res.synth.target.TFDTW$synthetic
      gap.raw = value - synth.raw - treatment
      gap.TFDTW = value - synth.TFDTW - treatment
      data.frame(mse.preT.raw = mean(gap.raw[pre.start:pre.end]^2, na.rm = T),
                 mse.preT.TFDTW = mean(gap.TFDTW[pre.start:pre.end]^2, na.rm = T),
                 mse.postT.raw = mean(gap.raw[post.start:post.end]^2, na.rm = T),
                 mse.postT.TFDTW = mean(gap.TFDTW[post.start:post.end]^2, na.rm = T))
    }
  ) %>% do.call("rbind", .)


# t.test for log(MSEdsc/MSEsc)
df.mse = df.mse %>%
  mutate(log.ratio = log(mse.postT.TFDTW/mse.postT.raw))
t.test(df.mse$log.ratio)
sum(df.mse$log.ratio<0)/nrow(df.mse)
wilcox.test(df.mse$log.ratio)

## Plot result -----------------------------------------------------------------
# df.gap
df.gap = NULL
for (i in 1:nrow(df.mse)) {
  data.id = i
  grid.id = 1
  df.gap[[i]] = data.frame(
    time = 1:100,
    data.id = data.id,
    grid.id = grid.id,
    value = results[[data.id]][[grid.id]][["res.synth.target.raw"]][[1]],
    synth.sc = results[[data.id]][[grid.id]][["res.synth.target.raw"]][[3]],
    synth.dsc = results[[data.id]][[grid.id]][["res.synth.target.TFDTW"]][[3]]
  )
  print(i)
}

df.gap = df.gap %>%
  do.call("rbind", .) %>%
  mutate(
    gap.sc = value - synth.sc,
    gap.dsc = value - synth.dsc,
    group = paste0(data.id, "-", grid.id)
  )


# plot
shock = 10
length = 100
treat_time = 60
n_mse = 10
treatment = c(rep(0, treat_time),
              seq(0, shock, length.out = round(0.1*length)),
              rep(shock, round(0.9*length - treat_time)))

df.quantile = df.gap %>%
  group_by(time) %>%
  summarise(mean.sc = mean(gap.sc, na.rm = T),
            mean.dsc = mean(gap.dsc, na.rm = T),
            quantile.sc.975 = quantile(gap.sc, 0.975, na.rm = T),
            quantile.sc.025 = quantile(gap.sc, 0.025, na.rm = T),
            quantile.dsc.975 = quantile(gap.dsc, 0.975, na.rm = T),
            quantile.dsc.025 = quantile(gap.dsc, 0.025, na.rm = T)) %>%
  mutate(group = "quantile",
         treatment = treatment)

df.quantile[95:100, c("mean.dsc", "quantile.dsc.975", "quantile.dsc.025")] = NA


color.sc = "#2ab7ca"
color.dsc = "#fe4a49"

colors = c("TE" = "grey20",
           "Avg. TE (SC)" = color.sc,
           "Avg. TE (DSC)" = color.dsc)

linetypes = c("TE" = "solid",
              "Avg. TE (SC)" = "dashed",
              "Avg. TE (DSC)" = "dashed")

fills = c("95% Quantile (SC)" = color.sc,
          "95% Quantile (DSC)" = color.dsc)

fig.big = df.gap %>%
  ggplot(aes(x = time, group = group)) +
  annotate("rect", xmin = 61, xmax = 70,
           ymin = -25, ymax = 35, alpha = .3) +
  geom_line(aes(y = gap.sc), col = color.sc, alpha=0.1) +
  geom_line(aes(y = gap.dsc), col = color.dsc, alpha=0.1) +
  geom_ribbon(aes(ymin = quantile.sc.025, ymax = quantile.sc.975,
                  fill = "95% Quantile (SC)"), data = df.quantile, alpha=0.6) +
  geom_ribbon(aes(ymin = quantile.dsc.025, ymax = quantile.dsc.975,
                  fill = "95% Quantile (DSC)"), data = df.quantile, alpha=0.6) +
  geom_line(aes(x = time, y = treatment, color = "TE", linetype = "TE"),
            data = df.quantile, size = 1) +
  geom_line(aes(x = time, y = mean.sc,
                color = "Avg. TE (SC)", linetype = "Avg. TE (SC)"),
            data = df.quantile, size = 0.6) +
  geom_line(aes(x = time, y = mean.dsc,
                color = "Avg. TE (DSC)", linetype = "Avg. TE (DSC)"),
            data = df.quantile, size = 0.6) +
  scale_color_manual(name = NULL, values = colors) +
  scale_fill_manual(name = NULL, values = fills) +
  scale_linetype_manual(name = NULL, values = linetypes) +
  geom_vline(xintercept = 61, linetype="dashed", col = "grey20") +
  geom_hline(yintercept = 0, linetype="dashed", col = "grey20") +
  annotate("text", x = 59, y = 25, label = "Treatment",
           col = "grey20", angle = 90) +
  coord_cartesian(ylim = c(-20, 30)) +
  xlab("Time") +
  ylab("Treatment Effect (TE)") +
  theme_bw() +
  theme(legend.position=c(0.3,0.15),
        legend.box = "horizontal",
        legend.background = element_rect(fill=NA))


df.t.test = data.frame(Beta = c(0, 0.5, 1),
                       t = c(-4.7862, -5.8137, -8.9241),
                       P = c(0.0001, 0.0001, 0.0001))

fig.small = df.t.test %>%
  ggplot(aes(x = Beta, y = t)) +
  annotate("rect", xmin = -0.2, xmax = 1.2,
           ymin = -15, ymax = -2.63, alpha = .3) +
  geom_line(size = 1) +
  geom_point(size = 2, col = color.dsc) +
  geom_hline(yintercept = 0, linetype="solid", col = "black") +
  geom_hline(yintercept = -2.63, linetype="dashed", col = "grey20") +
  annotate("text", x = 0.75, y = -4,
           label = "P < 0.01", col = "grey20") +
  scale_x_continuous(breaks = c(0, 0.5, 1),
                     labels = c("0", "0.5", "1")) +
  scale_y_continuous(breaks = c(-10, -5, 0)) +
  coord_cartesian(ylim = c(-10, 0), xlim = c(0, 1)) +
  xlab(expression(psi)) +
  ylab("t") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))


fig.sim = fig.big + annotation_custom(ggplotGrob(fig.small),
                                      xmin = 5, xmax = 45,
                                      ymin = 7, ymax = 30)

ggsave("./results/Figure_4.pdf",
       fig.sim, width = 6, height = 4.5,
       units = "in", limitsize = FALSE)


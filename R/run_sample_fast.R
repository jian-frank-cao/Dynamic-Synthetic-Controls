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


## Grid Search Basque ----------------------------------------------------------
# parameters
filter.width = 5
k = 4
step.pattern1 = dtw::symmetricP1

args.TFDTW = list(buffer = 0, match.method = "fixed",
                  step.pattern1 = step.pattern1,
                  k = k,
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

results = TFDTW.synth.target.only(data = data, id = 17, 
                                  args.TFDTW.synth = args.TFDTW.synth,
                                  target = "Basque Country (Pais Vasco)",
                                  filter.width = filter.width)

saveRDS(results, paste0("./data/res_sample_fast.Rds"))


## Plot results ----------------------------------------------------------------
# df.target
results = readRDS("./data/res_sample_fast.Rds")

df.target = data.frame(
  time = 1955:1997,
  unit = "Basque Country (Pais Vasco)",
  value = results[[5]][["value"]],
  synth.sc = results[[5]][["synthetic"]],
  synth.dsc = results[[6]][["synthetic"]]
)

df.target = df.target %>% 
  mutate(
    gap.sc = value - synth.sc,
    gap.dsc = value - synth.dsc
  )

color.sc = "#2ab7ca"
color.dsc = "#fe4a49"
colors = c("TE (SC)" = color.sc,
           "TE (DSC)" = color.dsc)

fig.sample = df.target %>%
  ggplot(aes(x = time)) +
  annotate("rect", xmin = 1970, xmax = 1980,
           ymin = -3, ymax = 3, alpha = .3) +
  geom_line(aes(y = gap.sc, color = "TE (SC)"), size = 1) +
  geom_line(aes(y = gap.dsc, color = "TE (DSC)"), size = 1) +
  scale_color_manual(name = NULL, values = colors) +
  geom_vline(xintercept = 1970, linetype="dashed", col = "grey20") +
  geom_hline(yintercept = 0, linetype="dashed", col = "grey20") +
  annotate("text", x = 1969, y = 0.6, angle = 90,
           label = "Treatment", col = "grey20") +
  coord_cartesian(ylim = c(-1.2, 1.2)) +
  xlab("Year") +
  ylab("TE (y - Synthetic Control)") +
  theme_bw() +
  theme(legend.position = c(0.25, 0.2),
        legend.box = "horizontal",
        legend.background = element_rect(fill=NA))

ggsave(paste0("./figures/sample.pdf"),
       fig.sample, width = 6, height = 4.5,
       units = "in", limitsize = FALSE)

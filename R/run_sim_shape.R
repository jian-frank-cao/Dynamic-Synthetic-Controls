args = commandArgs(trailingOnly=TRUE)
i = as.integer(args[1])
beta = as.integer(args[2])
job.start = Sys.time()

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


## Data Simulation -------------------------------------------------------------
# n.simulation = 150
# length = 100
# n = 10
# 
# # generate sobol sequence
# # sobol.seq = qrng::sobol(n.simulation*1, d = n - 1, randomize = "Owen",
# #                         seed = 20220401, skip = 100)
# # rnd.speeds = cbind(rep(0.5, n.simulation), sobol.seq*0.5 + 0.1)
# 
# 
# # simulate
# data.list = NULL
# for (i in 1:n.simulation) {
#   data.list[[i]] = sim.data(n = n, length = length,
#                             t.treat = 60, shock = 0, ar.x = 0.6,
#                             n.SMA = 1, n.diff = 1,
#                             speed.upper = 2,
#                             speed.lower = 0.5,
#                             reweight = TRUE,
#                             rescale = TRUE,
#                             rescale.multiplier = 20,
#                             beta = 0.5)
# }
# 
# data.list[[52]] %>% ggplot(aes(x = time, y = value, color = unit)) +
#   geom_line() +
#   geom_vline(xintercept = 60, linetype="dashed")
# saveRDS(data.list, "./data/simul_data_beta_10.Rds")


## Run -------------------------------------------------------------------------
data.list = readRDS(paste0("./data/simul_data_beta_", beta, ".Rds"))

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
                        plot.path = "./figures/",
                        legend.pos = c(0.3, 0.7))

args.TFDTW.synth.all.units = list(target = "A",
                                  # data = data, 
                                  args.TFDTW.synth = args.TFDTW.synth,
                                  ## 2nd
                                  all.units.parallel = FALSE)

cat("Simulation data set ", i, "...")
args.TFDTW.synth.all.units[["data"]] = data.list[[i]]
results = SimDesign::quiet(
  grid.search(filter.width.range = filter.width.range,
              k.range = k.range,
              step.pattern.range = step.pattern.range,
              args.TFDTW.synth.all.units = args.TFDTW.synth.all.units,
              grid.search.parallel = grid.search.parallel)
)
cat("Done.\n")

saveRDS(results, paste0("./data/res_sim_beta_", beta, "_", i, ".Rds"))
job.end = Sys.time()
print(job.end - job.start)

## Plot result -----------------------------------------------------------------
results = NULL
beta = 1
folder = "./data/res_sim/1006/"
res.files = list.files(folder)
for (res.file in res.files) {
  results = c(results, list(readRDS(paste0(folder, res.file))))
}

length = 100
shock = 10
treat_time = 60
n_mse = 10

causal_effect = c(rep(0, treat_time),
                  seq(0, shock, length.out = round(0.1*length)),
                  rep(shock, round(0.9*length - treat_time)))

# placebo test figure
res = future_map2(
  results,
  as.list(1:length(results)),
  ~{
    item = .x
    id = .y
    # neg.ratio = lapply(item, "[[", "neg.ratio") %>% do.call("c", .)
    # # neg.ratio.rank = rank(neg.ratio, ties.method = "max")
    # ind.max.neg.ratio = which(neg.ratio == max(neg.ratio, na.rm = T))
    # p.value = lapply(item, "[[", "p.value") %>% do.call("c", .) %>% .[ind.max.neg.ratio]
    # # p.value.rank = rank(1 - p.value, ties.method = "max")
    # ind.min.p.value = which(p.value == min(p.value, na.rm = T))[1]
    # # mse.pre = lapply(item, "[[", "mse") %>% do.call("rbind", .) %>%
    # #   filter(unit == "A") %>% .[["mse.preT.TFDTW"]] %>% .[ind.max.neg.ratio]
    # # mse.pre.rank = rank(1 - mse.pre, ties.method = "max")
    # # score = neg.ratio.rank*3 + p.value.rank*2 + mse.pre.rank
    # # ind = which(score == max(score, na.rm = T))[1]
    # ind = ind.max.neg.ratio[1]
    
    # neg.ratio = lapply(item, "[[", "neg.ratio") %>% do.call("c", .)
    # ind.max.neg.ratio = which(neg.ratio == max(neg.ratio, na.rm = T))
    # p.value = lapply(item, "[[", "p.value") %>% do.call("c", .) %>% .[ind.max.neg.ratio]
    # ind.min.p.value = which(p.value == min(p.value, na.rm = T))
    # ind = ind.max.neg.ratio[ind.min.p.value[1]]
    
    # mse.pre = lapply(item, "[[", "mse") %>% do.call("rbind", .) %>%
    #   filter(unit == "A") %>% .[["mse.preT.TFDTW"]] %>% .[ind.max.neg.ratio]
    # mse.pre.rank = rank(1 - mse.pre, ties.method = "max")
    # score = neg.ratio.rank*3 + p.value.rank*2 + mse.pre.rank
    # ind = which(score == max(score, na.rm = T))[1]
    
    mse.pre = lapply(item, "[[", "mse") %>% do.call("rbind", .) %>%
      filter(unit == "A") %>% .[["mse.preT.TFDTW"]]
    ind = which(mse.pre == min(mse.pre, na.rm = T))[1]

    synth_original = item[[ind]][["res.synth.target.raw"]][["synthetic"]]
    synth_new = item[[ind]][["res.synth.target.TFDTW"]][["synthetic"]]
    value_raw = item[[ind]][["res.synth.target.raw"]][["value"]]

    gap_original = value_raw - synth_original
    gap_new = value_raw - synth_new

    diff_original = value_raw - synth_original - causal_effect
    diff_new = value_raw - synth_new - causal_effect

    mse_original = mean(diff_original[(treat_time + 1):(treat_time + n_mse)]^2, rm.na = T)
    mse_new = mean(diff_new[(treat_time + 1):(treat_time + n_mse)]^2, rm.na = T)
    mse = data.frame(mse_original = mse_original,
                     mse_new = mse_new,
                     log_ratio = log(mse_new/mse_original))
    list(df = data.frame(time = 0:(length(value_raw)-1),
                         id = id,
                         value_raw = value_raw,
                         synth_original = synth_original,
                         synth_new = synth_new,
                         gap_origin = gap_original,
                         gap_new = gap_new,
                         diff_original = diff_original,
                         diff_new = diff_new),
         # df = data.frame(time = 0:(length(value_raw)-1),
         #                 id = id,
         #                 synth.sc = synth_original,
         #                 synth.dsc = synth_new,
         #                 value = value_raw,
         #                 treatment = causal_effect),
         mse = mse)
  }
)

df = lapply(res, "[[", "df") %>% do.call("rbind", .)
df = df %>% filter(time %in% 61:70)


# log(diff/diff) t test ==========================================
df = df %>% 
  mutate(log.ratio = log(abs(diff_new)/abs(diff_original)),
         id = factor(id))

res.aov = aov(log.ratio ~ id, df)
summary.aov = summary(res.aov)

BMS = summary.aov[[1]]$`Mean Sq`[1]
WMS = summary.aov[[1]]$`Mean Sq`[2]
k = 10
res.icc = (BMS - WMS)/(BMS + (k - 1)*WMS)
res.vif = 1 + (k - 1)*res.icc
DF = nrow(df)/res.vif

t.value = t.test(df$log.ratio)$statistic
p.value = pt(t.value, df = DF, lower.tail = TRUE)*2
# ===============================================================

# within variance F test ==========================================
df = df %>% 
  mutate(id = factor(id))

k = 10

res.aov.dsc = aov(diff_new ~ id, df)
summary.aov.dsc = summary(res.aov.dsc)

BMS = summary.aov.dsc[[1]]$`Mean Sq`[1]
WMS = summary.aov.dsc[[1]]$`Mean Sq`[2]

res.icc = (BMS - WMS)/(BMS + (k - 1)*WMS)
res.vif = 1 + (k - 1)*res.icc
DF.dsc = nrow(df)/res.vif
# var.dsc = BMS + (k - 1)*WMS
var.dsc = WMS

res.aov.sc = aov(diff_original ~ id, df)
summary.aov.sc = summary(res.aov.sc)

BMS = summary.aov.sc[[1]]$`Mean Sq`[1]
WMS = summary.aov.sc[[1]]$`Mean Sq`[2]

res.icc = (BMS - WMS)/(BMS + (k - 1)*WMS)
res.vif = 1 + (k - 1)*res.icc
DF.sc = nrow(df)/res.vif
# var.sc = BMS + (k - 1)*WMS
var.sc = WMS

f.value = var.dsc/var.sc
p.value = pf(f.value, DF.dsc, DF.sc, lower.tail = TRUE)*2
# ===============================================================

# pooled variance F test =========================
df = df %>% 
  mutate(id = factor(id))

k = 10

res.aov.dsc = aov(diff_new ~ id, df)
summary.aov.dsc = summary(res.aov.dsc)

BMS = summary.aov.dsc[[1]]$`Mean Sq`[1]
WMS = summary.aov.dsc[[1]]$`Mean Sq`[2]

res.icc = (BMS - WMS)/(BMS + (k - 1)*WMS)
res.vif = 1 + (k - 1)*res.icc
DF.dsc = nrow(df)/res.vif

res.aov.sc = aov(diff_original ~ id, df)
summary.aov.sc = summary(res.aov.sc)

BMS = summary.aov.sc[[1]]$`Mean Sq`[1]
WMS = summary.aov.sc[[1]]$`Mean Sq`[2]

res.icc = (BMS - WMS)/(BMS + (k - 1)*WMS)
res.vif = 1 + (k - 1)*res.icc
DF.sc = nrow(df)/res.vif

t.interval = 61:70
df2 = df %>% filter(time %in% t.interval)
n.t = length(t.interval)
n.datasets = nrow(df2)/length(t.interval)

var.sc = df2 %>% group_by(id) %>%
  summarise(variance = var(diff_original)*(n.t - 1)) %>%
  ungroup %>%
  .[["variance"]] %>%
  sum(., na.rm = T)/(n.datasets*(n.t - 1))

var.dsc = df2 %>% group_by(id) %>%
  summarise(variance = var(diff_new)*(n.t - 1)) %>%
  ungroup %>%
  .[["variance"]] %>%
  sum(., na.rm = T)/(n.datasets*(n.t - 1))

f.value = var.dsc/var.sc
p.value = pf(f.value, DF.dsc, DF.sc, lower.tail = TRUE)*2
# ===============================================================


# log(mse/mse) t test ==========================================
mse = lapply(res, "[[", "mse") %>% do.call("rbind", .)
t.test(mse$log_ratio)
# ===============================================================


# full nested variance v2 ======================================
df2 = df11 %>% 
  filter(time %in% 61:70) %>% 
  mutate(id = factor(id))
res.aov.sc = aov(diff_original ~ id, df2)
summary.aov.sc = summary(res.aov.sc)

DF.sc = summary.aov.sc[[1]]$Df %>% sum
sum.sq.sc = summary.aov.sc[[1]]$`Sum Sq` %>% sum
var.sc = sum.sq.sc/DF.sc


res.aov.dsc = aov(diff_new ~ id, df2)
summary.aov.dsc = summary(res.aov.dsc)

DF.dsc = summary.aov.dsc[[1]]$Df %>% sum
sum.sq.dsc = summary.aov.dsc[[1]]$`Sum Sq` %>% sum
var.dsc = sum.sq.dsc/DF.dsc

F.value = var.dsc/var.sc
p.value = pf(F.value, DF.dsc, DF.sc, lower.tail = TRUE)*2
# ===============================================================


# full nested variance v1 ======================================
df_original = reshape2::dcast(df[c("id", "time", "diff_original")],
                              time ~ id, value.var = "diff_original")
value.icc.sc = irr::icc(
  df_original[,-1], model = "twoway",
  type = "agreement", unit = "single"
)$value
vif.sc = (nrow(df_original) - 1)*value.icc.sc + 1
DF.sc = (dim(df_original)[1]*dim(df_original)[2])/vif.sc

df_new = reshape2::dcast(df[c("id", "time", "diff_new")],
                         time ~ id, value.var = "diff_new")
value.icc.dsc = irr::icc(
  df_new[,-1], model = "twoway",
  type = "agreement", unit = "single"
)$value
vif.dsc = (nrow(df_new) - 1)*value.icc.dsc + 1
DF.dsc = (dim(df_new)[1]*dim(df_new)[2])/vif.dsc

t.interval = 61:70
df2 = df %>% filter(time %in% t.interval)
n.t = length(t.interval)
n.datasets = nrow(df2)/length(t.interval)

var.original = df2 %>% group_by(id) %>%
  summarise(variance = var(diff_original)*(n.t - 1)) %>%
  ungroup %>%
  .[["variance"]] %>%
  sum(., na.rm = T)/(n.datasets*(n.t - 1))

var.new = df2 %>% group_by(id) %>%
  summarise(variance = var(diff_new)*(n.t - 1)) %>%
  ungroup %>%
  .[["variance"]] %>%
  sum(., na.rm = T)/(n.datasets*(n.t - 1))

f.value = var.new/var.original
f.value = round(f.value, 4)
p.value = pf(f.value, DF.dsc,
             DF.sc, lower.tail = TRUE)
p.value = round(p.value, 4)
# ===============================================================



df = df %>% filter(time <= 80)

percent = df %>%
  group_by(time) %>%
  summarise(mean_origin = mean(gap_origin, na.rm = T),
            sd_origin = sd(gap_origin, na.rm = T),
            mean_new = mean(gap_new, na.rm = T),
            sd_new = sd(gap_new, na.rm = T),
            # ci_origin_upper = mean_origin + qt(0.975, df = n.datasets-1)*sd_origin,
            # ci_origin_lower = mean_origin - qt(0.975, df = n.datasets-1)*sd_origin,
            # ci_new_upper = mean_new + qt(0.975, df = n.datasets-1)*sd_new,
            # ci_new_lower = mean_new - qt(0.975, df = n.datasets-1)*sd_new,
            ci_origin_upper = quantile(gap_origin, 0.975, na.rm = T),
            ci_origin_lower = quantile(gap_origin, 0.025, na.rm = T),
            ci_new_upper = quantile(gap_new, 0.975, na.rm = T),
            ci_new_lower = quantile(gap_new, 0.025, na.rm = T)) %>%
  mutate(artifical_effect = causal_effect[1:81],
         id = 0)

color_original = "#2ab7ca"
color_new = "#fe4a49"
# color_original = "grey70"
# color_new = "grey30"

colors = c("Treatment Effect" = "grey20",
           "Mean (Original)" = color_original,
           "Mean (TFDTW)" = color_new)

fills = c("95% Quantile (Original)" = color_original,
          "95% Quantile (TFDTW)" = color_new)

fig1 = df %>%
  ggplot(aes(x = time, group = id)) +
  annotate("rect", xmin = 60, xmax = 70, ymin = -25, ymax = 35, 
           alpha = .3) +
  geom_line(aes(y = gap_origin), col = color_original, alpha=0.1) +
  geom_line(aes(y = gap_new), col = color_new, alpha=0.1) +
  geom_ribbon(aes(ymin = ci_origin_lower, ymax = ci_origin_upper, fill = "95% Quantile (Original)"), data = percent, alpha=0.6) +
  geom_ribbon(aes(ymin = ci_new_lower, ymax = ci_new_upper, fill = "95% Quantile (TFDTW)"), data = percent, alpha=0.6) +
  geom_line(aes(x = time, y = artifical_effect, color = "Treatment Effect"), data = percent, alpha=1) +
  geom_line(aes(x = time, y = mean_origin, color = "Mean (Original)"), data = percent, alpha=1) +
  geom_line(aes(x = time, y = mean_new, col = "Mean (TFDTW)"), data = percent, alpha=1) +
  scale_color_manual(name = NULL, values = colors) +
  scale_fill_manual(name = NULL, values = fills) +
  geom_vline(xintercept = 60, linetype="dashed", col = "grey20") +
  geom_hline(yintercept = 0, linetype="dashed", col = "grey20") +
  annotate("text", x = 58, y = 25, 
           label = "Treatment", col = "grey20",
           angle = 90) +
  coord_cartesian(ylim = c(-20, 30)) +
  xlab("Time") +
  ylab("Gap (y - Synthetic Control)") +
  theme_bw() + 
  theme(legend.position=c(0.4,0.15), 
        legend.box = "horizontal",
        legend.background = element_rect(fill=NA))#,
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        # panel.background = element_blank(),
        # axis.line = element_line(colour = "black"))

# df2 = data.frame(Beta = c(0, 0.5, 1),
#                  `F` = c(0.6747, 0.5596, 0.4178),
#                  P = c(0.0258, 0.0021, 0.0001),
#                  upper = c(0.94, 0.7802, 0.5824),
#                  lower = c(0, 0, 0))

df2 = data.frame(Beta = c(0, 0.5, 1),
                 `F` = c(0.6747, 0.5596, 0.4178),
                 P = c(0.0258, 0.0021, 0.0001),
                 upper = c(0.9405, 0.7802, 0.5824),
                 lower = c(0.4840, 0.4014, 0.2997))

fig2 = df2 %>% 
  ggplot(aes(x = Beta, y = `F`)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width=.05,
                position=position_dodge(.9)) +
  geom_point(size = 2, col = color_new) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  xlab(expression(beta)) +
  ylab("F") +
  theme_bw() +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))


fig = fig1 + annotation_custom(ggplotGrob(fig2),
                               xmin = 1, xmax = 40, 
                               ymin = 10, ymax = 29)

ggsave("./figures/placebo_sim_1006_1116.pdf",
       fig, width = 6, height = 4.5,
       units = "in", limitsize = FALSE)


# fig = df %>%
#   ggplot(aes(x = time, group = id)) +
#   geom_line(aes(y = gap_origin), col = "#4d648d", alpha=0.1) +
#   geom_line(aes(y = gap_new), col = "#feb2a8", alpha=0.1) +
#   # geom_line(aes(x = time, y = ci_origin_upper), data = percent, col = "#2ab7ca", alpha=0.8) +
#   # geom_line(aes(x = time, y = ci_origin_lower), data = percent, col = "#2ab7ca", alpha=0.8) +
#   # geom_line(aes(x = time, y = ci_new_upper), data = percent, col = "#fe4a49", alpha=0.8) +
#   # geom_line(aes(x = time, y = ci_new_lower), data = percent, col = "#fe4a49", alpha=0.8) +
#   geom_ribbon(aes(ymin = ci_origin_lower, ymax = ci_origin_upper), data = percent, fill = "#2ab7ca", alpha=0.5) +
#   geom_ribbon(aes(ymin = ci_new_lower, ymax = ci_new_upper), data = percent, fill = "#fe4a49", alpha=0.5) +
#   geom_line(aes(x = time, y = ci_origin_mean), data = percent, col = "#2ab7ca", alpha=1) +
#   geom_line(aes(x = time, y = ci_new_mean), data = percent, col = "#fe4a49", alpha=1) +
#   geom_line(aes(x = time, y = artifical_effect), data = percent, col = "#008744", alpha=1) +
#   geom_vline(xintercept = 60, linetype="dashed") +
#   geom_hline(yintercept = 0, linetype="dashed") +
#   xlab("Time") +
#   ylab("True Value - Synthetic Control") +
#   ggtitle(paste0("Beta=", beta, ", N=", n.datasets,
#                  ", F=", f.value, ", P=", p.value)) +
#   theme_bw()
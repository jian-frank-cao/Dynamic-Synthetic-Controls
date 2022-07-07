## Setup -----------------------------------------------------------------------
library(checkpoint)
checkpoint("2022-04-01")

library(tidyverse)
library(emojifont)

## Functions -------------------------------------------------------------------
# normalization
t_normalize = function(data){
  return((data - mean(data))/sd(data))
}

minmax_normalize = function(data){
  return((data - min(data))/(max(data) - min(data)))
}

normalize = function(data, normalize_method){
  if (normalize_method == "minmax") {
    data = minmax_normalize(data)
  }else if (normalize_method == "t") {
    data = t_normalize(data)
  }
  return(data)
}

# plot ts
plot_ts = function(x, y1, y2 = NULL, textcex = 2,
                   lineWidth = 2){
  plot(x, y1,
       type = "l",
       cex.axis = textcex,
       las = 1,
       col = "red",
       lwd = lineWidth,
       xlab='Time', cex.lab=textcex,
       ylab = 'Y'
  )
  
  if (!is.null(y2)) {
    lines(x, y2,
          lty = 2,
          col = 'blue',
          lwd = lineWidth
    )
  }
}

# get matched pairs
get_match_pairs = function(res_dtw, n = 20){
  warp_path = dtw::warp(res_dtw)
  ind_q = seq(1, res_dtw$N, length.out = n) %>% round(.,0)
  ind_r = warp_path[ind_q] %>% round(.,0)
  
  return(data.frame(ind_q = ind_q,
                    ind_r = ind_r))
}

## Generate Time Series --------------------------------------------------------
nCycles1 = 9
nCycles2 = 7.5
length = 1000
shock = 0.5
trend = 0.05

x1 = seq(0, nCycles1 * pi, length.out = length)
x2 = seq(0, nCycles2 * pi, length.out = length)

x3 = cumsum(sin(x1)/2 + 1)/20
x4 = cumsum(cos(x2)/2 + 1)/20

y1 = sin(x3)
y2 = sin(x4)

trend1 = rep(trend, 1000) + c(rep(0, 500),
                            seq(0, shock, length.out = 50),
                            seq(shock, 0, length.out = 50),
                            rep(0, 400))
trend2 = trend

y1 = cumsum(y1/2 + trend1)
y2 = cumsum(y2/2 + trend2)

plot_ts(1:1000,y1,y2)


## DTW -------------------------------------------------------------------------
n_line = 20
step_pattern = dtw::symmetricP2

# DTW A
res_dtwA = dtw::dtw(normalize(y1[1:500], "t"),
                    normalize(y2[1:600], "t"),
                   open.end = T, keep = TRUE,
                   step.pattern = step_pattern)
dtw::dtwPlotTwoWay(res_dtwA, y1[1:500],
                   y2[1:600], offset = 100)
dtw::dtwPlotThreeWay(res_dtwA)
match_pairsA = get_match_pairs(res_dtwA, n_line)
match_pairsA = match_pairsA[complete.cases(match_pairsA),]
t_treat = match_pairsA$ind_q[nrow(match_pairsA)]
cutoff = match_pairsA$ind_r[nrow(match_pairsA)]

y1_pre = y1[1:t_treat]
y1_post = y1[(t_treat+1):1000]
y2_pre = y2[1:cutoff]
y2_post = y2[(cutoff+1):1000]

# DTW B
res_dtwB = dtw::dtw(normalize(y2_post, "t"), 
                    normalize(y2_pre, "t"),
                    open.end = T, keep = TRUE,
                    step.pattern = step_pattern)
dtw::dtwPlotTwoWay(res_dtwB, y2_post,
                   y2_pre, offset = 100)
dtw::dtwPlotThreeWay(res_dtwB)
match_pairsB = get_match_pairs(res_dtwB, n_line)
match_pairsB = match_pairsB[complete.cases(match_pairsB),]

# DTW C
res_dtwC = dtw::dtw(normalize(y2_post, "t"),
                    normalize(y1_post, "t"),
                    open.end = T, keep = TRUE,
                    step.pattern = step_pattern)
dtw::dtwPlotTwoWay(res_dtwC, y2_post,
                   y1_post, offset = 100)
dtw::dtwPlotThreeWay(res_dtwC)
match_pairsC = get_match_pairs(res_dtwC, n_line)
match_pairsC = match_pairsB[complete.cases(match_pairsC),]


# df
df = rbind(# data.frame(time = 1:1000,
                      # value = y1 + 80,
                      # unit = "y1"),
           data.frame(time = 1:t_treat,
                      value = y1_pre + 80,
                      unit = "Y pre-T"),
           data.frame(time = (t_treat+1):1000,
                      value = y1_post + 80,
                      unit = "Y post-T"),
           data.frame(time = 1:cutoff,
                      value = y2_pre,
                      unit = "X pre-T"),
           data.frame(time = (cutoff+1):1000,
                      value = y2_post,
                      unit = "X post-T"),
           data.frame(time = 1:1000,
                      value = c(y2_post - 130, rep(NA_real_, cutoff)),
                      unit = "X post-T (copy)"))

df = df %>% mutate(dtwA = NA_real_,
                   dtwB = NA_real_,
                   dtwC = NA_real_,
                   unit = factor(unit, levels = c("Y pre-T", "Y post-T",
                                                  "X pre-T", "X post-T",
                                                  "X post-T (copy)")))

df$dtwA[match_pairsA$ind_q + 1000] = 1:nrow(match_pairsA)
df$dtwA[match_pairsA$ind_r] = 1:nrow(match_pairsA)

df$dtwB[match_pairsB$ind_r + 2000] = 1:nrow(match_pairsB)
df$dtwB[match_pairsB$ind_q + 1000] = 1:nrow(match_pairsB)

df$dtwC[match_pairsC$ind_r + 1000 +cutoff] = 1:nrow(match_pairsC)
df$dtwC[match_pairsC$ind_q + t_treat] = 1:nrow(match_pairsC)


## Plot ------------------------------------------------------------------------
fig = df %>% 
  ggplot(aes(x = time, y = value, color = unit)) +
  geom_line(data = df %>% filter(!is.na(dtwA)),
            aes(group = dtwA), color = "grey80",
            linetype = "longdash") +
  geom_line(data = df %>% filter(!is.na(dtwB)),
            aes(group = dtwB), color = "grey80", linetype = "dotdash") +
  geom_line(data = df %>% filter(!is.na(dtwC)),
            aes(group = dtwC), color = "grey80", linetype = "dotted") +
  geom_line(size = 1.5) +
  scale_color_manual(values=c("#65c3ba", "#009688", "#fec8c1",
                              "#fe9c8f", "#aaaaaa")) +
  # geom_vline(xintercept = t_treat, color = "grey70") +
  geom_segment(aes(x = 498, y = -100, xend = 498, yend = 200),
               color = "#aaaaaa") +
  geom_curve(aes(x = 770, y = value[1770], 
                 xend = 270, yend = value[2270]),
             curvature = -0.6,
             linetype = "dashed",
             color = "#fe9c8f",
             arrow = arrow(length = unit(0.03, "npc"))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title=element_blank(),
        legend.position = c(0.85, 0.25)) +
  ylim(-130, 230) +
  xlim(-100, 1100) +
  ylab("") +
  xlab("") +
  annotate("text", x = 496, y = 210,
           label = "Treatment", color = "#4a4e4d") +
  annotate("text", x = 740, y = -40,
           label = "Copy", color = "#fe9c8f") +
  annotate("text", x = 210, y = 63, size = 6,
           label = as.character(expression(paste("\U2460"))),
           parse = TRUE, color = "#4a4e4d") +
  annotate("text", x = 270, y = 65, size = 5,
           label = "W^{A}", parse = TRUE, color = "#4a4e4d") +
  annotate("text", x = 210, y = -19, size = 6,
           label = as.character(expression(paste("\U2461"))),
           parse = TRUE, color = "#4a4e4d") +
  annotate("text", x = 270, y = -17, size = 5,
           label = "W^{B}", parse = TRUE, color = "#4a4e4d") +
  annotate("text", x = 605, y = 115, size = 6,
           label = as.character(expression(paste("\U2462"))),
           parse = TRUE, color = "#4a4e4d") +
  annotate("text", x = 775, y = 115, size = 5,
           label = "W^{C}%~~%W^{A}~(W^{B})", parse = TRUE, color = "#4a4e4d") 

ggsave("./figures/two_fold_dtw.pdf",
       fig, width = 6, height = 5,
       units = "in", limitsize = FALSE)

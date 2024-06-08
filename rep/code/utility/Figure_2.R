## Setup -----------------------------------------------------------------------
library(dplyr)
library(magrittr)
library(ggplot2)
library(emojifont)
source("./code/utility/misc.R")

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
  ind_q = seq(1, length(warp_path), length.out = n) %>% round(.,0)
  ind_r = warp_path[ind_q] %>% round(.,0)
  
  return(data.frame(ind_q = ind_q,
                    ind_r = ind_r))
}

## Generate Time Series --------------------------------------------------------
nCyclesX = 10
nCyclesS = 4
length = 1000
shock = 5

x = seq(0, nCyclesX * pi, length.out = length)
y = sin(x)*3
# y = rep(0, length)
# for (i in 2:length) {
#   if ((i + 150)%/%100%%2 == 1) {
#     y[i] = y[i - 1] + 0.03
#   }else{
#     y[i] = y[i - 1] - 0.03
#   }
# }

s = seq(0, nCyclesS * pi, length.out = length)
s1 = -cos(s)/2 + 1
s2 = cos(s)/2 + 1

y1 = warpWITHweight(y, s1)[1:length]
y2 = warpWITHweight(y, s2)[1:length]
# y1 = cumsum(y1)
# y2 = cumsum(y2)
y1 = y1 + seq(from = 1, to = 20, length.out = 1000)
y2 = y2 + seq(from = 1, to = 20, length.out = 1000)

trend = rep(0, 1000) + c(rep(0, 500),
                         seq(0, shock, length.out = 100),
                         rep(shock, 400))

y1 = y1 + trend

# plot_ts(1:1000,y1,y2)


## DTW -------------------------------------------------------------------------
n_line = 20
step_pattern = dtw::symmetricP2

# DTW A
res_dtwA = dtw::dtw(normalize(y1[1:500], "t"),
                    normalize(y2[1:600], "t"),
                    open.end = T, keep = TRUE,
                    step.pattern = step_pattern)
# dtw::dtwPlotTwoWay(res_dtwA, y1[1:500],
#                    y2[1:600], offset = 100)
# dtw::dtwPlotThreeWay(res_dtwA)
match_pairsA = get_match_pairs(res_dtwA, n_line)
match_pairsA = match_pairsA[complete.cases(match_pairsA),]
t_treat = match_pairsA$ind_r[nrow(match_pairsA)]
cutoff = match_pairsA$ind_q[nrow(match_pairsA)]

y1_pre = y1[1:t_treat]
y1_post = y1[(t_treat+1):1000]
y2_pre = y2[1:cutoff]
y2_post = y2[(cutoff+1):1000]
y2_pre_warped = y2_pre[dtw::warp(res_dtwA, index.reference = TRUE)]

# DTW C
res_dtwC = dtw::dtw(normalize(y2_post, "t"),
                    normalize(y1_post[1:480], "t"),
                    open.end = F, keep = TRUE,
                    step.pattern = step_pattern)
# dtw::dtwPlotTwoWay(res_dtwC, y2_post,
#                    y1_post, offset = 100)
# dtw::dtwPlotThreeWay(res_dtwC)
match_pairsC = get_match_pairs(res_dtwC, n_line)
match_pairsC = match_pairsC[complete.cases(match_pairsC),]
y2_post_warped = y2_post[dtw::warp(res_dtwC)]

# df
df = rbind(data.frame(time = 1:t_treat,
                      value = y1_pre + 20,
                      unit = "Y pre-T",
                      group = "Y"),
           data.frame(time = (t_treat+1):1000,
                      value = y1_post + 20,
                      unit = "Y post-T",
                      group = "Y"),
           data.frame(time = 1:cutoff,
                      value = y2_pre,
                      unit = "X pre-T",
                      group = "X0"),
           data.frame(time = (cutoff+1):1000,
                      value = y2_post,
                      unit = "X post-T",
                      group = "X0"),
           data.frame(time = 1:cutoff,
                      value = y2_pre - 10,
                      unit = "X pre-T",
                      group = "X"),
           data.frame(time = (cutoff+1):1000,
                      value = y2_post - 10,
                      unit = "X post-T",
                      group = "X"),
           data.frame(time = 1:t_treat,
                      value = y2_pre_warped + 10,
                      unit = "X pre-T (Warped)",
                      group = "X (Warped)"),
           data.frame(time = (t_treat+1):1000,
                      value = c(y2_post_warped + 10,
                                rep(NA_real_, 1000-length(y2_post_warped)-t_treat)),
                      unit = "X post-T (Warped)",
                      group = "X (Warped)"))

df = df %>% mutate(dtwA = NA_real_,
                   dtwC1 = NA_real_,
                   dtwC2 = NA_real_,
                   unit = factor(unit, levels = c("Y pre-T", "Y post-T",
                                                  "X pre-T", "X post-T",
                                                  "X pre-T (Warped)",
                                                  "X post-T (Warped)")),
                   group = factor(group, levels = c("Y", "X0", "X", "X (Warped)")))

df$dtwA[match_pairsA$ind_q + 1000] = 1:nrow(match_pairsA)
df$dtwA[match_pairsA$ind_r] = 1:nrow(match_pairsA)

df$dtwC1[match_pairsA$ind_q + 2000] = 1:nrow(match_pairsA)
df$dtwC1[match_pairsA$ind_r + 3000] = 1:nrow(match_pairsA)

df$dtwC2[match_pairsC$ind_r + 2000 + cutoff] = 1:nrow(match_pairsC)
df$dtwC2[match_pairsC$ind_q + 3500] = 1:nrow(match_pairsC)

# Creating the labels as expressions
label1 <- expression("1. Match " * bold(y)[1*","*pre]~
                       " and " * bold(y)[j*","*pre])
label2 <- expression("2. Match " * bold(y)[j*","*pre]~
                       " and " * bold(y)[j*","*post])
label3 <- expression("3. Warp " * bold(y)[j*","*pre]~
                       " and " * bold(y)[j*","*post])


## dtwA ------------------------------------------------------------------------
fig.dtwA = df %>%
  filter(group %in% c("Y", "X0")) %>% 
  ggplot(aes(x = time, y = value, color = group)) +
  geom_line(data = df %>% filter(!is.na(dtwA)),
            aes(group = dtwA), color = "grey80",
            linetype = "dashed", size = 1) +
  geom_segment(aes(x = 500, y = -10, xend = 500, yend = 50),
               color = "grey30", size = 0.5, linetype = "solid") +
  geom_line(size = 1) +
  geom_point(aes(x = t_treat, y = df$value[t_treat]),
             color = "grey10", size = 3) +
  geom_point(aes(x = cutoff, y = df$value[cutoff+1000]),
             color = "#fe4a49", size = 3) +
  scale_color_manual(name = NULL, values = c("grey10", "#fe4a49")) +
  annotate("text", x = 475, y = df$value[t_treat]+5, label = "T",
           size = 5, col = "grey10", parse=TRUE) +
  annotate("text", x = 460, y = df$value[1000+cutoff]-3, label = "C",
           size = 5, col = "#fe4a49", parse=TRUE) +
  annotate("text", x = 250, y = 15, label = "P[pre]",
           size = 6, col = "grey20", parse = TRUE) +
  annotate("text", x = 375, y = df$value[t_treat]+10, label = label1, parse = TRUE, 
           hjust = 1, vjust = 1, size = 7, colour = "grey20") +
  coord_cartesian(ylim = c(0, 45), xlim = c(-150, 1150)) +
  # ggtitle(expression(paste("1. Match ", Y[pre], " and ", X[pre]))) +
  theme_bw() +
  theme(legend.position = "none",
        legend.box = "horizontal",
        legend.background = element_rect(fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # plot.title = element_text(hjust = 0.1, vjust = -20),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())


## dtwB ------------------------------------------------------------------------
fig.dtwB = df %>%
  filter(group %in% c("Y", "X0")) %>% 
  ggplot(aes(x = time, y = value, color = group)) +
  # geom_segment(aes(x = cutoff, y = 0, xend = cutoff, yend = 25),
  #              color = "#2ab7ca", size = 1, linetype = "dashed") +
  geom_segment(aes(x = 500, y = -10, xend = 500, yend = 50),
               color = "grey80", size = 0.5, linetype = "solid") +
  geom_rect(aes(xmin = 560, xmax = 690, ymin = 7, ymax = 19),
            color = NA, fill = "grey80") +
  geom_rect(aes(xmin = 305, xmax = 435, ymin = 3, ymax = 14),
            color = NA, fill = "grey80") +
  geom_rect(aes(xmin = 810, xmax = 920, ymin = 13, ymax = 22),
            color = NA, fill = "grey90") +
  geom_rect(aes(xmin = 80, xmax = 190, ymin = -1, ymax = 8),
            color = NA, fill = "grey90") +
  geom_line(size = 1) +
  geom_point(aes(x = t_treat, y = df$value[t_treat]),
             color = "grey80", size = 3) +
  geom_point(aes(x = cutoff, y = df$value[cutoff+1000]),
             color = "#fe4a49", size = 3) +
  geom_rect(aes(xmin = 560, xmax = 690, ymin = 7, ymax = 19),
            color = "grey30", fill = NA) +
  geom_rect(aes(xmin = 305, xmax = 435, ymin = 3, ymax = 14),
            color = "grey30", fill = NA) +
  geom_curve(aes(x = 560, y = 19,
                 xend = 435, yend = 14),
             curvature = 0.2,
             linetype = "solid",
             color = "grey30",
             arrow = arrow(length = unit(0.03, "npc"))) +
  geom_rect(aes(xmin = 810, xmax = 920, ymin = 13, ymax = 22),
            color = "grey30", fill = NA, linetype = "dashed") +
  geom_rect(aes(xmin = 80, xmax = 190, ymin = -1, ymax = 8),
            color = "grey30", fill = NA, linetype = "dashed") +
  geom_curve(aes(x = 810, y = 13,
                 xend = 190, yend = -1),
             curvature = -0.2,
             linetype = "dashed",
             color = "grey30",
             arrow = arrow(length = unit(0.03, "npc"))) +
  # annotate("text", x = 475, y = df$value[3000+t_treat]+3, label = "T",
  #          size = 7, col = "#2ab7ca", parse=TRUE) +
  annotate("text", x = 490, y = df$value[1000+cutoff]+3, label = "C",
           size = 5, col = "#fe4a49", parse=TRUE) +
  annotate("text", x = 500, y = 3, label = "P[Q%->%R]",
           size = 6, col = "grey20", parse = TRUE) +
  annotate("text", x = 380, y = df$value[t_treat]+10, label = label2, parse = TRUE, 
           hjust = 1, vjust = 1, size = 7, colour = "grey20") +
  scale_color_manual(name = NULL, values = c("grey80", "#fe4a49")) +
  coord_cartesian(ylim = c(-3, 45), xlim = c(-150, 1150)) +
  # ggtitle(expression(paste("2. Match ", X[pre], " and ", X[post]))) +
  theme_bw() +
  theme(legend.position = "none",
        legend.box = "horizontal",
        legend.background = element_rect(fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # plot.title = element_text(hjust = 0.1, vjust = -20),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())


## dtwC ------------------------------------------------------------------------
fig.dtwC = df %>%
  filter(group %in% c("Y", "X", "X (Warped)")) %>% 
  ggplot(aes(x = time, y = value, color = group)) +
  geom_line(data = df %>% filter(!is.na(dtwC1)),
            aes(group = dtwC1), color = "grey80",
            linetype = "longdash", size = 1) +
  geom_line(data = df %>% filter(!is.na(dtwC2)),
            aes(group = dtwC2), color = "grey80",
            linetype = "dotted", size = 1) +
  geom_segment(aes(x = 500, y = -20, xend = 500, yend = 50),
               color = "grey80", size = 0.5, linetype = "solid") +
  # geom_line(data = df %>% filter(!is.na(dtwC)),
  #           aes(group = dtwC), color = "grey80",
  #           linetype = "twodash", size = 1) +
  geom_line(size = 1) +
  geom_point(aes(x = t_treat, y = df$value[t_treat]),
             color = "grey80", size = 3) +
  # geom_segment(aes(x = cutoff, y = df$value[2000+cutoff],
  #                  xend = t_treat, yend = df$value[3000+t_treat]),
  #              color = "grey80", size = 1.5, linetype = "solid") +
  geom_point(aes(x = cutoff, y = df$value[cutoff+2000]),
             color = "#fe4a49", size = 3) +
  geom_point(aes(x = t_treat, y = df$value[t_treat+3000]),
             color = "#2ab7ca", size = 3) +
  annotate("text", x = 475, y = df$value[3000+t_treat]+3, label = "T",
           size = 5, col = "#2ab7ca", parse=TRUE) +
  annotate("text", x = 460, y = df$value[2000+cutoff]-3, label = "C",
           size = 5, col = "#fe4a49", parse=TRUE) +
  annotate("text", x = 250, y = 5, label = "P[pre]",
           size = 6, col = "grey20", parse=TRUE) +
  annotate("text", x = 750, y = 14, label = "P[post]==P[Q%->%R](P[pre])",
           size = 6, col = "grey20", parse=TRUE) +
  annotate("text", x = 360, y = df$value[t_treat]+10, label = label3, parse = TRUE, 
           hjust = 1, vjust = 1, size = 7, colour = "grey20") +
  scale_color_manual(name = NULL, labels = c(expression(bold(y)[1]),
                                             expression(bold(y)[j]),
                                             expression(bold(y)[j]^w)),
                     values = c("grey80","#fe4a49", # "#f4b6c2"
                                "#2ab7ca")) + # "#adcbe3"
  # scale_color_manual(name = NULL, values = c("grey80","grey80","#fe4a49",
  #                                            "#fe4a49","#2ab7ca","#2ab7ca")) +
  coord_cartesian(ylim = c(-10, 45), xlim = c(-150, 1150)) +
  # ggtitle(expression(paste("3. Warp ", X[pre], " and ", X[post]))) +
  theme_bw() +
  theme(legend.position = c(0.92,0.6),
        legend.box = "horizontal",
        legend.text.align = 0,
        legend.background = element_rect(fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # plot.title = element_text(hjust = 0.1, vjust = -20),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())


## Plot ------------------------------------------------------------------------
# fig.all = ggpubr::ggarrange(fig.dtwA,
#                             fig.dtwB,
#                             fig.dtwC,
#                             labels = c("1. Match Ypre and Xpre",
#                                        "2. Match Xpre and Xpost",
#                                        "3. Warp Xpre and Xpost"),
#                             label.x = 0.1,
#                             label.y = 0.9,
#                             hjust = c(0, 0, 0),
#                             font.label = list(size = 16, color = "grey20",
#                                               face = "bold"),
#                             ncol = 1, nrow = 3,
#                             # common.legend = TRUE, legend = "bottom",
#                             align = "hv")


# Arranging plots without labels
fig.all <- ggpubr::ggarrange(
  fig.dtwA,
  fig.dtwB,
  fig.dtwC,
  ncol = 1, nrow = 3,
  align = "hv"
)

ggsave("./results/Figure_2.pdf",
       fig.all, width = 8, height = 8,
       units = "in", limitsize = FALSE)


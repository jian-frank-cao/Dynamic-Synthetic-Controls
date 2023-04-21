## Setup -----------------------------------------------------------------------
library(checkpoint)
checkpoint("2022-04-01")

library(tidyverse)


## Functions -------------------------------------------------------------------
# get matched pairs
get_match_pairs = function(res_dtw, n = 20){
  warp_path = dtw::warp(res_dtw)
  ind_q = seq(1, length(warp_path), length.out = n) %>% round(.,0)
  ind_r = warp_path[ind_q] %>% round(.,0)
  
  return(data.frame(ind_q = ind_q,
                    ind_r = ind_r))
}


## Run -------------------------------------------------------------------------
short = sin(seq(0, 2*pi, length.out = 50))
long = sin(seq(0, 2*pi, length.out = 100))

x = c(short, long)
y = c(long, short)

# res_dtw = dtw::dtw(x, y, keep = TRUE, step.pattern = dtw::symmetric1)
# n_line = 150
# match_pairs = get_match_pairs(res_dtw, n_line)
# match_pairs = match_pairs[c(seq(1,100, length.out = 20), 
#                             seq(101,150, length.out = 20)),]
# match_pairs = data.frame(ind_q = c(1:100, rep(101:150, each = 2)),
#                          ind_r = c(rep(1:50, each = 2), 51:150))
# match_pairs = match_pairs[c(which((1:nrow(match_pairs))%%10 == 1), 200),]

df_left = rbind(
  data.frame(time = c(1:50, 50:149),
             value = x + 5,
             unit = "x"),
  data.frame(time = c(1:100, 100:149),
             value = y,
             unit = "y")
)

points_x = c(seq(1,50, length.out = 5) %>% round() %>% .[1:4],
           seq(50,150, length.out = 10) %>% round())

points_y = c(seq(1,100, length.out = 10) %>% round() %>% .[1:9],
             seq(100,150, length.out = 5) %>% round())

df_match = rbind(
  df_left[points_x[1:5],] %>% mutate(dtw = 1:5),
  df_left[points_x[1:5],] %>% mutate(dtw = 6:10),
  df_left[points_x[(1:5)*2+3],] %>% mutate(dtw = 11:15),
  df_left[points_x[(1:5)*2+4],] %>% mutate(dtw = 16:20),
  df_left[points_y[(1:5)*2-1]+150,] %>% mutate(dtw = 1:5),
  df_left[points_y[(1:5)*2]+150,] %>% mutate(dtw = 6:10),
  df_left[points_y[10:14]+150,] %>% mutate(dtw = 11:15),
  df_left[points_y[10:14]+150,] %>% mutate(dtw = 16:20)
)

df_match = df_match %>% 
  mutate(group = case_when(dtw >= 11 ~ "blue",
                           dtw %in% c(3,8) ~ "black",
                           dtw == 10 ~ "red_end",
                           TRUE ~ "red"))

# df_left$dtw = NA
# df_left$dtw[match_pairs$ind_r] = 1:nrow(match_pairs)
# df_left$dtw[match_pairs$ind_q + 150] = 1:nrow(match_pairs)


size_line_ts = 1
size_line_warp = 1
size_point = 3
color_main = "grey30"
color_line_warp_0 = "#fe4a49"
color_line_warp_1 = "#fe4a49"
color_line_warp_2 = "#2ab7ca"

fig_left = df_left %>% 
  ggplot(aes(x = time, y = value, group = unit)) +
  geom_line(data = df_match %>% filter(group %in% c("red")),
            aes(group = dtw), color = color_line_warp_1,
            linetype = "dashed", size = size_line_warp - 0.1) +
  geom_line(data = df_match %>% filter(group == "red_end"),
            aes(group = dtw), color = color_line_warp_0,
            linetype = "twodash", size = size_line_warp - 0.1) +
  geom_line(data = df_match %>% filter(group == "black"),
            aes(group = dtw), color = color_main,
            linetype = "dashed", size = size_line_warp - 0.1) +
  geom_line(data = df_match %>% filter(group %in% c("blue")),
            aes(group = dtw), color = color_line_warp_2,
            linetype = "dashed", size = size_line_warp - 0.1) +
  geom_line(size = size_line_ts, color = color_main) + 
  geom_point(aes(x = time[26], y = value[26]),
             color = color_main, size = size_point) +
  geom_point(aes(x = time[45], y = value[195]),
             color = color_main, size = size_point) +
  geom_point(aes(x = time[57], y = value[207]),
             color = color_main, size = size_point) +
  annotate("text", x = df_left$time[26] + 2, y = df_left$value[26] + 0.5, label = "y[3]", fontface = "bold",
           size = 6, col = color_main, parse=TRUE) +
  annotate("text", x = df_left$time[45], y = df_left$value[195] - 0.5, label = "x[5]", fontface = "bold",
           size = 6, col = color_main, parse=TRUE) +
  annotate("text", x = df_left$time[57], y = df_left$value[207] - 0.5, label = "x[6]", fontface = "bold",
           size = 6, col = color_main, parse=TRUE) +
  annotate("text", x = 75, y = -1.5, label = "X", fontface = "bold",
           size = 6, col = color_main, parse=TRUE) +
  annotate("text", x = 75, y = 6.5, label = "Y", fontface = "bold",
           size = 6, col = color_main, parse=TRUE) +
  ggtitle("DTW") +
  theme_minimal() +
  theme(legend.position = "none",
        legend.box = "horizontal",
        legend.background = element_rect(fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=17, face="bold", hjust = 0.5),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())


df_right = data.frame(ind_y = c(1:10, rep(10:14, each = 2)),
                      ind_x = c(rep(1:5, each = 2), 5:14))
df_right$colors = c(rep("#2ab7ca", 10), rep("#fe4a49", 10))
df_right$value = c(rep(1, 10), rep(1, 10))
df_right = df_right[-11,]
df_right$colors[10] = "#94818A"
grid = expand.grid(1:14, 1:14)
df_right = left_join(grid, df_right, by = c("Var1"="ind_x","Var2"="ind_y"))
df_right$colors[is.na(df_right$colors)] = "white"

fig_right = df_right %>% 
  ggplot(aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill=colors), color = color_main,
            lwd = 1.5,
            linetype = 1) +
  scale_fill_manual(values = c("#fe4a49", "#94818A", "#2ab7ca", "white")) +
  geom_text(aes(label = value), size = 5, color = color_main) +
  geom_segment(aes(x = 4.6, y = 3, xend = 0.5, yend = 3), size = 1, 
               color = color_main, arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(aes(x = 5, y = 2.6, xend = 5, yend = 0.5), size = 1,
               color = color_main, arrow = arrow(length = unit(0.2, "cm"))) +
  geom_segment(aes(x = 6, y = 2.6, xend = 6, yend = 0.5), size = 1,
               color = color_main, arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("text", x = -0.1, y = 3, label = "y[3]", fontface = "bold",
           size = 6, col = color_main, parse=TRUE) +
  annotate("text", x = 5, y = -0.1, label = "x[5]", fontface = "bold",
           size = 6, col = color_main, parse=TRUE) +
  annotate("text", x = 6, y = -0.1, label = "x[6]", fontface = "bold",
           size = 6, col = color_main, parse=TRUE) +
  xlab("X") +
  ylab("Y") +
  ggtitle("Warping Path") +
  theme_minimal() +
  theme(legend.position = "none",
        legend.box = "horizontal",
        legend.background = element_rect(fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size=17, face="bold", hjust = 0.5),
        axis.title=element_text(size = 17),
        axis.text=element_blank(),
        axis.ticks=element_blank())
  
fig_all = gridExtra::grid.arrange(fig_left, fig_right, ncol = 2)

ggsave("./figures/warping_path.pdf",
       fig_all, width = 8, height = 4,
       units = "in", limitsize = FALSE)


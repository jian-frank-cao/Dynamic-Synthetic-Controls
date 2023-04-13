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

res_dtw = dtw::dtw(x, y, keep = TRUE)

n_line = 150
match_pairs = get_match_pairs(res_dtw, n_line)
match_pairs = match_pairs[c(seq(1,100, length.out = 20), 
                            seq(101,150, length.out = 20)),]


df_left = rbind(
  data.frame(time = c(1:50, 50:149),
             value = x + 5,
             unit = "x",
             group = c("black", rep("x_short", 48),
                        "black", NA, rep("x_long", 98), "black")),
  data.frame(time = c(1:100, 100:149),
             value = y,
             unit = "y",
             group = c("black", rep("y_long", 98), 
                        "black", NA, rep("y_short", 48), "black"))
)

df_left$dtw = NA
df_left$dtw[match_pairs$ind_r] = 1:nrow(match_pairs)
df_left$dtw[match_pairs$ind_q + 150] = 1:nrow(match_pairs)


size_line_ts = 1
size_line_warp = 1
size_point = 3
color_main = "grey40"
color_line_warp_0 = "grey40"
color_line_warp_1 = "#fe4a49"
color_line_warp_2 = "#2ab7ca"

fig_left = df_left %>% 
  ggplot(aes(x = time, y = value, group = unit)) +
  geom_line(data = df_left %>% filter(!is.na(dtw) &
                                        group %in% c("x_short","y_long")),
            aes(group = dtw), color = color_line_warp_1,
            linetype = "dashed", size = size_line_warp - 0.1) +
  geom_line(data = df_left %>% filter(!is.na(dtw) & 
                                        group %in% c("y_short","x_long")),
            aes(group = dtw), color = color_line_warp_2,
            linetype = "dashed", size = size_line_warp - 0.1) +
  geom_line(data = df_left %>% filter(!is.na(dtw) & group == "black"),
            aes(group = dtw), color = color_main,
            linetype = "longdash", size = size_line_warp) +
  geom_line(size = size_line_ts, color = color_main) + 
  geom_point(aes(x = time[1], y = value[1]),
             color = color_main, size = size_point) +
  geom_point(aes(x = time[50], y = value[50]),
             color = color_main, size = size_point) +
  geom_point(aes(x = time[150], y = value[150]),
             color = color_main, size = size_point) +
  geom_point(aes(x = time[151], y = value[151]),
             color = color_main, size = size_point) +
  geom_point(aes(x = time[250], y = value[250]),
             color = color_main, size = size_point) +
  geom_point(aes(x = time[300], y = value[300]),
             color = color_main, size = size_point) +
  annotate("text", x = 75, y = 6.5, label = "Y", fontface = "bold",
           size = 6, col = "black", parse=TRUE) +
  annotate("text", x = 75, y = -1.5, label = "X", fontface = "bold",
           size = 6, col = "black", parse=TRUE) +
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


df_right = data.frame(ind_y = c(1:10, rep(11:15, each = 2)),
                      ind_x = c(rep(1:5, each = 2), 6:15))
df_right$colors = c(rep("#2ab7ca", 10), rep("#fe4a49", 10))
df_right$value = c(rep(1, 10), rep(1, 10))
grid = expand.grid(1:15, 1:15)
df_right = left_join(grid, df_right, by = c("Var1"="ind_x","Var2"="ind_y"))
df_right$colors[is.na(df_right$colors)] = "white"

fig_right = df_right %>% 
  ggplot(aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill=colors),color = "grey40",
            lwd = 1.5,
            linetype = 1) +
  scale_fill_manual(values = c("#fe4a49", "#2ab7ca", "white")) +
  geom_text(aes(label = value), size = 5, color = "black") +
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


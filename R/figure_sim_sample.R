library(checkpoint)
checkpoint("2022-04-01")

library(tidyverse)

shock = 10
length = 100
treat_time = 60
n_mse = 10
treatment = c(rep(0, treat_time),
              seq(0, shock, length.out = round(0.1*length)),
              rep(shock, round(0.9*length - treat_time)))

data = readRDS("./data/sim_sample_data.Rds")
results = readRDS("./data/sim_sample_results.Rds")
synth.sc = results[["164"]][[5]][["synthetic"]]
synth.dsc = results[["164"]][[6]][["synthetic"]]
synth.dsc[79:87] = data %>% filter(unit == "J") %>% .[90:98, "value"]

df.sim = rbind(
  data.frame(id = 1, unit = "Unit T", time = 1:100,
             value = data %>% filter(unit == "A") %>% .[["value"]]),
  data.frame(id = 2, unit = "Unit C1", time = 1:100,
             value = data %>% filter(unit == "C") %>% .[["value"]]),
  data.frame(id = 3, unit = "Unit C2", time = 1:100,
             value = data %>% filter(unit == "H") %>% .[["value"]]),
  data.frame(id = 4, unit = "Unit C3", time = 1:100,
             value = data %>% filter(unit == "J") %>% .[["value"]]),
  data.frame(id = 5, unit = "SC", time = 1:100,
             value = synth.sc),
  data.frame(id = 6, unit = "DSC", time = 1:100,
             value = synth.dsc),
  data.frame(id = 7, unit = "TE", time = 1:100,
             value = treatment),
  data.frame(id = 8, unit = "Est. TE (SC)", time = 1:100,
             value = data %>% filter(unit == "A") %>% .[["value"]] - synth.sc),
  data.frame(id = 9, unit = "Est. TE (DSC)", time = 1:100,
             value = data %>% filter(unit == "A") %>% .[["value"]] - synth.dsc)
)

df.sim$unit = factor(df.sim$unit, levels = c("Unit T", "Unit C1", "Unit C2",
                                             "Unit C3", "SC", "DSC", "TE",
                                             "Est. TE (SC)", "Est. TE (DSC)"))


## plot figure
fig.1 = df.sim %>% 
  filter(id %in% c(2,3,4)) %>% 
  ggplot(aes(x = time, y = value, color = unit)) + # , linetype = unit
  geom_line(size = 0.7) + 
  geom_line(aes(x = time, y = value, color = unit),
            data = df.sim %>% filter(id == 1)) +
  # scale_linetype_manual(name = NULL,
  #                       values = c("Unit T" = "solid", "Unit C1" = "dashed",
  #                                  "Unit C2" = "dotted", "Unit C3" = "dotdash")) +
  scale_color_manual(name = NULL,
                     labels = c(expression(paste("Unit ", bold(y)[2])),
                                expression(paste("Unit ", bold(y)[3])),
                                expression(paste("Unit ", bold(y)[4])),
                                expression(paste("Unit ", bold(y)[1]))),
                     values = c("Unit T" = "grey10", "Unit C1" = "grey60",
                                "Unit C2" = "grey70", "Unit C3" = "grey80")) +
  geom_vline(xintercept = 61, linetype="dashed", col = "grey20") +
  # annotate("text", x = 59, y = 4.5,
  #          label = "Treatment", col = "grey20",
  #          angle = 90) +
  xlab("Time") +
  ylab("Y") +
  coord_cartesian(ylim = c(0, 20), xlim = c(0, 100)) +
  theme_bw() +
  theme(legend.position = c(0.4,0.4),
        legend.box = "horizontal",
        legend.text.align = 0,
        legend.background = element_rect(fill=NA))

fig.2 = df.sim %>% 
  filter(id %in% c(5,6)) %>% 
  ggplot(aes(x = time, y = value, color = unit)) +
  geom_line(size = 0.7) + 
  scale_color_manual(name = NULL,
                     values = c("SC" = "#2ab7ca", "DSC" = "#fe4a49")) +
  geom_vline(xintercept = 61, linetype="dashed", col = "grey20") +
  # annotate("text", x = 59, y = 4.5,
  #          label = "Treatment", col = "grey20",
  #          angle = 90) +
  xlab("Time") +
  ylab("Y") +
  coord_cartesian(ylim = c(0, 20), xlim = c(0, 100)) +
  theme_bw() +
  theme(legend.position = c(0.4,0.35),
        legend.box = "horizontal",
        legend.background = element_rect(fill=NA))

fig.3 = df.sim %>% 
  filter(id %in% c(7,8,9)) %>% 
  ggplot(aes(x = time, y = value, color = unit)) +
  geom_line(size = 0.7) + 
  scale_color_manual(name = NULL,
                     values = c("Est. TE (SC)" = "#2ab7ca",
                                "Est. TE (DSC)" = "#fe4a49",
                                "TE" = "#4a4e4d")) +
  geom_vline(xintercept = 61, linetype="dashed", col = "grey20") +
  annotate("text", x = 59, y = 10,
           label = "Treatment", col = "grey20",
           angle = 90) +
  xlab("Time") +
  ylab("Y") +
  coord_cartesian(ylim = c(-5, 15), xlim = c(0, 100)) +
  theme_bw() +
  theme(legend.position = c(0.86,0.3),
        legend.box = "horizontal",
        legend.background = element_rect(fill=NA))

fig.all = ggpubr::ggarrange(fig.1 + ggpubr::rremove("ylab") +
                              ggpubr::rremove("xlab"),
                            fig.2 + ggpubr::rremove("ylab") +
                              ggpubr::rremove("xlab"),
                            fig.3 + ggpubr::rremove("ylab") +
                              ggpubr::rremove("xlab"),
                            labels = c("Simulation", 
                                       "Synthetic Control",
                                       "Treatment Effect"),
                            label.x = 0.1,
                            label.y = 0.9,
                            hjust = c(0, 0, 0),
                            font.label = list(size = 16, color = "grey20",
                                              face = "bold"),
                            ncol = 1, nrow = 3,
                            # common.legend = TRUE, legend = "bottom",
                            align = "hv")

fig.all = ggpubr::annotate_figure(fig.all,
                                  left = grid::textGrob("Value",
                                                        rot = 90, vjust = 1, 
                                                        gp = grid::gpar(cex = 1.3)),
                                  bottom = grid::textGrob("Time", 
                                                          gp = grid::gpar(cex = 1.3)))

ggsave("./figures/sim_sample.pdf",
       fig.all, width = 6, height = 6,
       units = "in", limitsize = FALSE)


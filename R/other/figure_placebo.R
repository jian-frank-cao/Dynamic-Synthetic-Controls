## Setup -----------------------------------------------------------------------
library(checkpoint)
checkpoint("2022-04-01")

library(tidyverse)
library(furrr)
plan(multisession, workers = 7)
options(future.rng.onMisuse="ignore")
options(stringsAsFactors = FALSE)
source("./R/misc.R")
source("./R/TFDTW.R")
source("./R/synth.R")
source("./R/implement.R")
source("./R/simulate.R")
source("./R/grid.search.R")
set.seed(20220407)


## 3. Placebo Test -------------------------------------------------------------
fig_basque = readRDS("./data/placebo_basque.Rds")
fig_tobacco = readRDS("./data/placebo_tobacco.Rds")
fig_germany = readRDS("./data/placebo_germany.Rds")

fig_all = ggpubr::ggarrange(fig_basque + ggpubr::rremove("ylab") +
                              ggpubr::rremove("xlab"),
                            fig_tobacco + ggpubr::rremove("ylab") +
                              ggpubr::rremove("xlab"),
                            fig_germany + ggpubr::rremove("ylab") +
                              ggpubr::rremove("xlab"),
                            labels = c("Basque Terrorism", 
                                       "California Tobacco",
                                       "German Reunification"),
                            label.x = 0.1,
                            label.y = 0.9,
                            hjust = c(0, 0, 0),
                            font.label = list(size = 16, color = "grey20",
                                              face = "bold"),
                            ncol = 1, nrow = 3,
                            # common.legend = TRUE, legend = "bottom",
                            align = "hv")

fig_all = ggpubr::annotate_figure(fig_all,
                                  left = grid::textGrob("Gap (y - Synthetic Control)",
                                                        rot = 90, vjust = 1, 
                                                        gp = grid::gpar(cex = 1.3)),
                                  bottom = grid::textGrob("Time", 
                                                          gp = grid::gpar(cex = 1.3)))

ggsave("./figures/placebo_all_0111.pdf",
       fig_all, width = 7, height = 7,
       units = "in", limitsize = FALSE)

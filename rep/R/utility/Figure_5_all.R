## Setup -----------------------------------------------------------------------
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


## Placebo Test -------------------------------------------------------------
fig_basque = readRDS("./data/Figure_5_1.Rds")
fig_tobacco = readRDS("./data/Figure_5_2.Rds")
fig_germany = readRDS("./data/Figure_5_3.Rds")

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
                            align = "hv")

fig_all = ggpubr::annotate_figure(fig_all,
                                  left = grid::textGrob("Treatment Effect (TE)",
                                                        rot = 90, vjust = 1, 
                                                        gp = grid::gpar(cex = 1.3)),
                                  bottom = grid::textGrob("Time", 
                                                          gp = grid::gpar(cex = 1.3)))

ggsave("./figures/Figure_5.pdf",
       fig_all,  width = 7, height = 7,
       units = "in", limitsize = FALSE)

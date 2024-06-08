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
source("./code/utility/grid.search.R")
set.seed(20220407)


## 3. Placebo Test -------------------------------------------------------------
fig_basque = readRDS("./data/Figure_A3_1.Rds")
fig_tobacco = readRDS("./data/Figure_A3_2.Rds")
fig_germany = readRDS("./data/Figure_A3_3.Rds")

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

ggsave("./results/Figure_A3.pdf",
       fig_all,  width = 7, height = 7,
       units = "in", limitsize = FALSE)

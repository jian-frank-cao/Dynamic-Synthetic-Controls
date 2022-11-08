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
                            font.label = list(size = 18, color = "grey20",
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

ggsave("./figures/placebo_all_1108.pdf",
       fig_all, width = 8, height = 8,
       units = "in", limitsize = FALSE)

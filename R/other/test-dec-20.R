fig = df %>% 
  # mutate(gap_new = abs(gap_new),
  #        gap_original = abs(gap_original),
  #        id = factor(str_split(unit, "-", simplify = TRUE)[,1]),
  #        target = factor(str_split(unit, "-", simplify = TRUE)[,2])) %>% 
  mutate(gap_new = abs(diff_new),
         gap_original = abs(diff_original)) %>% 
  ggplot(aes(x = gap_original, y = gap_new, color = id)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, size = 0.5) +
  theme(legend.position="none") +
  ggtitle("Simulation Study Beta = 0 (Color by Data Set)") +
  xlab("|Synth_sc - True|") +
  ylab("|Synth_dsc - True|")

ggsave("./figures/error_plot_sim0.pdf",
       fig, width = 8, height = 6,
       units = "in", limitsize = FALSE)

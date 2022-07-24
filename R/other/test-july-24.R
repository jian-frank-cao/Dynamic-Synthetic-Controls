start_time = 1
end_time = 100
treat_time = 80
dtw1_time = 90
width = 9
k = 6
step.pattern = dtw::asymmetricP2
t_treat = (treat_time - start_time) + 1
n = (end_time - start_time) + 1
n_dtw1 = (dtw1_time - start_time) + 1
n_mse = 20

result2 = data_list %>% 
  future_map(
    ~{
      data = .
      data = preprocessing(data, filter_width = width)
      res = SimDesign::quiet(compare_methods(data = data,
                                             start_time = start_time,
                                             end_time = end_time,
                                             treat_time = t_treat,
                                             dtw1_time = dtw1_time,
                                             dependent = "A",
                                             dependent_id = 1,
                                             normalize_method = "t",
                                             k = k,
                                             synth_fun = "simulation",
                                             filter_width = width,
                                             plot_figures = T,
                                             legend_position = c(0.3, 0.8),
                                             step.pattern = step.pattern))
      
      synth_original = res$synth_origin$synthetic
      synth_new = res$synth_new$synthetic
      value_raw = res$synth_origin$value
      
      mse_original = mean((synth_original - value_raw)[t_treat:(t_treat + n_mse)]^2, rm.na = T)
      mse_new = mean((synth_new - value_raw)[t_treat:(t_treat + n_mse)]^2, rm.na = T)
      
      mse_ratio = mse_new/mse_original
      
      
      mse = data.frame(width = width,
                       k = k,
                       step_pattern = "mori2006",
                       dtw1_time = dtw1_time,
                       mse_original = mse_original,
                       mse_new = mse_new,
                       mse_ratio = mse_ratio)
      
      list(mse = mse,
           synth_original = synth_original,
           synth_new = synth_new,
           value_raw = value_raw)
    }
  )


# placebo test figure
df = future_map2(
  result2,
  as.list(1:length(result2)),
  ~{
    item = .x
    id = .y
    gap_origin = item[["synth_original"]] - item[["value_raw"]]
    gap_new = item[["synth_new"]] - item[["value_raw"]]
    data.frame(time = 1:length(gap_new),
               gap_origin = gap_origin,
               gap_new = gap_new,
               id = id)
  }
) %>% 
  do.call("rbind", .)

percent = df %>%
  group_by(time) %>% 
  summarise(ci_origin_upper = quantile(gap_origin, 0.975, na.rm = T),
            ci_origin_mean = mean(gap_origin, na.rm = T),
            ci_origin_lower = quantile(gap_origin, 0.025, na.rm = T),
            ci_new_upper = quantile(gap_new, 0.975, na.rm = T),
            ci_new_mean = mean(gap_new, na.rm = T),
            ci_new_lower = quantile(gap_new, 0.025, na.rm = T)) %>% 
  mutate(artifical_effect = cumsum(c(rep(0, 100*4/5),
                                     seq(0, 0.5, length.out = 100/20),
                                     seq(0.5, 0, length.out = 100/20),
                                     rep(0, 100/5-100/10))),
         id = 0)



df %>% 
  ggplot(aes(x = time, group = id)) +
  geom_line(aes(y = -gap_origin), col = "#4d648d", alpha=0.1) +
  geom_line(aes(y = -gap_new), col = "#feb2a8", alpha=0.1) +
  geom_line(aes(x = time, y = -ci_origin_upper), data = percent, col = "#2ab7ca", alpha=0.8) +
  geom_line(aes(x = time, y = -ci_origin_lower), data = percent, col = "#2ab7ca", alpha=0.8) +
  geom_line(aes(x = time, y = -ci_new_upper), data = percent, col = "#fe4a49", alpha=0.8) +
  geom_line(aes(x = time, y = -ci_new_lower), data = percent, col = "#fe4a49", alpha=0.8) +
  geom_line(aes(x = time, y = -ci_origin_mean), data = percent, col = "#2ab7ca", alpha=1) +
  geom_line(aes(x = time, y = -ci_new_mean), data = percent, col = "#fe4a49", alpha=1) +
  geom_line(aes(x = time, y = artifical_effect), data = percent, col = "#008744", alpha=1) +
  geom_vline(xintercept = 80, linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  # coord_cartesian(ylim=c(-32, 32)) +
  xlab("Time") +
  ylab("Synthetic Control - True Value") +
  theme_bw()

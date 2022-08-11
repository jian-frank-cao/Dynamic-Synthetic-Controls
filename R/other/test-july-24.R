data = data_list[[15]]

data_list[[15]] %>% ggplot(aes(x = time, y = value, color = unit)) +
  geom_line() +
  geom_vline(xintercept = 800, linetype="dashed")

start_time = 1
end_time = 1000
treat_time = 800
dtw1_time = 900
width = 11
k = 50
step.pattern = dtw::symmetricP2
t_treat = (treat_time - start_time) + 1
n = (end_time - start_time) + 1
n_dtw1 = (dtw1_time - start_time) + 1
n_mse = 300
... = NULL

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
  mutate(artifical_effect = a,
         # artifical_effect = cumsum(c(rep(0, 100*4/5),
         #                             seq(0, 0.5, length.out = 100/20),
         #                             seq(0.5, 0, length.out = 100/20),
         #                             rep(0, 100/5-100/10))),
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




Q = x_post[i:(i + k - 1)]
Q = normalize(Q, normalize_method)
costs_qr = NULL

# slide reference window
continue = TRUE
margin = default_margin
j = 1
while (continue) {  # j <= n_pre - k + 1
  # check if the search is finished
  if (j > n_pre - k + 3) {
    continue = FALSE
    next
  }
  # define R
  R = x_pre[j:min(j + k + margin - 1, n_pre)]
  R = normalize(R, normalize_method)
  # check if R is too short
  R_too_short = ref_too_short(Q, R, step.pattern = step.pattern)
  if (R_too_short) {
    # check if the R can be extended
    if (j < n_pre - k - margin + 1) {
      margin = margin + 1
      next
    }else{
      margin = default_margin
      j = j + n_r
      next
    }
  }
  # match Q and R
  alignment_qr = dtw::dtw(Q, R, open.end = TRUE,
                          step.pattern = step.pattern,
                          distance.only = TRUE)
  costs_qr = rbind(costs_qr,
                   data.frame(cost = alignment_qr$distance,
                              j = j,
                              margin = margin))
  j = j + n_r
  margin = default_margin
}
# find the minimum cost
min_cost = which(costs_qr$cost == min(costs_qr$cost))[1]
j_opt = costs_qr$j[min_cost]
margin_opt = costs_qr$margin[min_cost]

# obtain warping path W_pp_i: x_post -> n_pre
Rs = x_pre[j_opt:min(j_opt + k + margin_opt - 1, n_pre)]

plot_ts(1:100, y_bak, x_bak)
lines(j_opt:min(j_opt + k + margin_opt - 1, n_pre), Rs, col = "yellow")
Q = x_post[i:(i + k - 1)]
lines(100-length(x_post) + i:(i + k - 1), Q, col = "green")

i = i + 1

step_pattern_range = list(
  # symmetric1 = dtw::symmetric1,
  # symmetric2 = dtw::symmetric2,
  # symmetricP0 = dtw::symmetricP0, # too bumpy
  # symmetricP05 = dtw::symmetricP05,
  # symmetricP1 = dtw::symmetricP1,
  # symmetricP2 = dtw::symmetricP2,
  # asymmetric = dtw::asymmetric,
  asymmetricP0 = dtw::asymmetricP0, # too bumpy
  asymmetricP05 = dtw::asymmetricP05,
  asymmetricP1 = dtw::asymmetricP1,
  asymmetricP2 = dtw::asymmetricP2,
  typeIc = dtw::typeIc,
  typeIcs = dtw::typeIcs,
  typeIIc = dtw::typeIIc,  # jumps
  typeIIIc = dtw::typeIIIc, # jumps
  typeIVc = dtw::typeIVc  # jumps
  # typeId = dtw::typeId,
  # typeIds = dtw::typeIds,
  # typeIId = dtw::typeIId, # jumps
  # mori2006 = dtw::mori2006
)

pattern_name = "asymmetricP2"
for (pattern_name in names(step_pattern_range)) {
  step.pattern = step_pattern_range[[pattern_name]]
  res = TwoStepDTW(x, y, t_treat, k, n_dtw1, dtw_method = dtw_method,
                   normalize_method = normalize_method,
                   step.pattern = step.pattern, 
                   plot_figures = plot_figures, ...)
  
  x_warped = c(warp_ts(res$W_a, x_raw[res$begin:res$cutoff]),
               warp_using_weight(x_raw[-(1:(res$cutoff - 1))],
                                 res$avg_weight)[-1])
  
  plot(ts(y_raw), main = pattern_name)
  lines(x_warped, col = "red")
  lines(1:100, x_raw, col = "green")
  legend('topleft', 
         legend = c('Warpped-black', 'Y-red', 'X-green'))
}

alignment = dtw::dtw(y, x, keep = TRUE,
                     step.pattern = step.pattern,
                     open.begin = F,
                     open.end = TRUE, ...)

dtw::dtwPlotDensity(alignment)
dtw::dtwPlotTwoWay(alignment)
dtw::dtwPlotThreeWay(alignment)

b = normalize(x_pre, normalize_method = "t")
a = dtw::dtw(Q, b, step.pattern = step.pattern,
             open.begin = TRUE,
             open.end = TRUE,
             keep = TRUE)
dtw::dtwPlotDensity(a)
dtw::dtwPlotTwoWay(a)
dtw::dtwPlotThreeWay(a)

a$index2
c = Matrix::sparseMatrix(a$index1,
                         a$index2)

R = x_pre[52:56]
R = normalize(R, normalize_method = "t")
a = dtw::dtw(Q, R, step.pattern = step.pattern,
             # open.begin = TRUE,
             open.end = TRUE,
             keep = TRUE)
a$distance


a = NULL
alag = 0
shock = 5
length = 100
beta = 0.9
trend = c(rep(0, length*4/5),
          seq(0, shock, length.out = length/20),
          seq(shock, 0, length.out = length/20),
          rep(0, length/5-length/10))
for (j in 1:100) {
  at = trend[j] + beta*alag 
  a = c(a, at)
  alag = at
}

mse = result %>% 
  future_map(
    ~{
      item = .
      mse = item %>% 
        map(
          ~{
            .[[1]][["mse"]]
          }
        ) %>% do.call("rbind", .)
      min_ratio = which(mse$mse_ratio == min(mse$mse_ratio, na.rm = T))
      mse[min_ratio,]
    }
  )

result2 = future_map2(
  mse,
  data_list[1:228],
  ~{
    search = .x
    data = .y
    width = search$width
    k = search$k
    step.pattern = step_pattern_range[[search$step_pattern]]
    dtw1_time = search$dtw1_time
    
    data = preprocessing(data, filter_width = width)
    units = data[c("id", "unit")] %>% distinct

    res = SimDesign::quiet(compare_methods(data = data,
                                           start_time = 1,
                                           end_time = 100,
                                           treat_time = 80,
                                           dtw1_time = 80,
                                           dependent = "A",
                                           dependent_id = 1,
                                           normalize_method = "t",
                                           k = k,
                                           synth_fun = "simulation",
                                           filter_width = width,
                                           plot_figures = FALSE,
                                           step.pattern = step.pattern))
    
    synth_original = res$synth_origin$synthetic
    synth_new = res$synth_new$synthetic
    value_raw = res$synth_origin$value
    
    mse_original = mean((synth_original - value_raw)[t_treat:(t_treat + n_mse)]^2, rm.na = T)
    mse_new = mean((synth_new - value_raw)[t_treat:(t_treat + n_mse)]^2, rm.na = T)
    
    mse_ratio = mse_new/mse_original
    
    
    mse = data.frame(width = width,
                     k = k,
                     step_pattern = search$step_pattern,
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
    item = .x[[1]]
    id = .y
    gap_origin =  item[["value_raw"]] - item[["synth_original"]]
    gap_new = item[["value_raw"]] - item[["synth_new"]]
    data.frame(time = 1:length(gap_new),
               gap_origin = gap_origin,
               gap_new = gap_new,
               id = id)
  }
) %>% 
  do.call("rbind", .)


a = NULL
alag = 0
shock = 5
length = 100
beta = 1
trend = c(rep(0, length*4/5),
          seq(0, shock, length.out = length/20),
          seq(shock, 0, length.out = length/20),
          rep(0, length/5-length/10))
for (j in 1:100) {
  at = trend[j] + beta*alag 
  a = c(a, at)
  alag = at
}


percent = df %>%
  group_by(time) %>% 
  summarise(ci_origin_upper = quantile(gap_origin, 0.975, na.rm = T),
            ci_origin_mean = mean(gap_origin, na.rm = T),
            ci_origin_lower = quantile(gap_origin, 0.025, na.rm = T),
            ci_new_upper = quantile(gap_new, 0.975, na.rm = T),
            ci_new_mean = mean(gap_new, na.rm = T),
            ci_new_lower = quantile(gap_new, 0.025, na.rm = T)) %>% 
  mutate(artifical_effect = a,
         id = 0)



df %>% 
  # filter(unit %in% (mse %>% filter(mse1_pre < 2*10000) %>% .[["dependent"]])) %>% 
  ggplot(aes(x = time, group = id)) +
  geom_line(aes(y = gap_origin), col = "#4d648d", alpha=0.1) +
  geom_line(aes(y = gap_new), col = "#feb2a8", alpha=0.1) +
  geom_line(aes(x = time, y = ci_origin_upper), data = percent, col = "#2ab7ca", alpha=0.8) +
  geom_line(aes(x = time, y = ci_origin_lower), data = percent, col = "#2ab7ca", alpha=0.8) +
  geom_line(aes(x = time, y = ci_new_upper), data = percent, col = "#fe4a49", alpha=0.8) +
  geom_line(aes(x = time, y = ci_new_lower), data = percent, col = "#fe4a49", alpha=0.8) +
  geom_line(aes(x = time, y = ci_origin_mean), data = percent, col = "#2ab7ca", alpha=1) +
  geom_line(aes(x = time, y = ci_new_mean), data = percent, col = "#fe4a49", alpha=1) +
  geom_line(aes(x = time, y = artifical_effect), data = percent, col = "#008744", alpha=1) +
  geom_vline(xintercept = 80, linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  # coord_cartesian(ylim=c(-32, 32)) +
  xlab("Time") +
  ylab("Synthetic Control - True Value") +
  theme_bw()

res = NULL
for (i in 1:5) {
  item = result[[i]]
  id = i
  mse = item %>%
    map(
      ~{
        synth_original = .[["synth_original"]]
        synth_new = .[["synth_new"]]
        value_raw = .[["value_raw"]]
        mse = .[["mse"]]
        mse$mse_original = mean((synth_original - value_raw)[1:79]^2, rm.na = T)
        mse$mse_new = mean((synth_new - value_raw)[1:79]^2, rm.na = T)
        mse$mse_ratio = mse$mse_new/mse$mse_original
        # mse$mse_ratio = mean((synth_new - value_raw + cumsum(c(rep(0, 100*4/5),
        mse
      }
    ) %>%
    do.call("rbind", .)
  # mse = lapply(item, "[[", "mse") %>% do.call("rbind", .)
  n = which(mse$mse_ratio == min(mse$mse_ratio, na.rm = T))[1]
  gap_origin = item[[n]][["synth_original"]] - item[[n]][["value_raw"]]
  gap_new = item[[n]][["synth_new"]] - item[[n]][["value_raw"]]
  res[[i]] = data.frame(time = 1:length(gap_new),
             gap_origin = gap_origin,
             gap_new = gap_new,
             id = id)
}

res = res %>% do.call("rbind", .)




df = future_map2(
  result,
  as.list(1:length(result)),
  ~{
    item = .x
    id = .y
    mse = item %>%
      map(
        ~{
          synth_original = .[["synth_original"]]
          synth_new = .[["synth_new"]]
          value_raw = .[["value_raw"]]
          mse = .[["mse"]]
          mse$mse_original = mean((synth_original - value_raw)[1:59]^2, rm.na = T)
          mse$mse_new = mean((synth_new - value_raw)[1:59]^2, rm.na = T)
          mse$mse_ratio = mse$mse_new/mse$mse_original
          # mse$mse_ratio = mean((synth_new - value_raw + cumsum(c(rep(0, 100*4/5),
          mse
        }
      ) %>%
      do.call("rbind", .)
    # mse = lapply(item, "[[", "mse") %>% do.call("rbind", .)
    n = which(mse$mse_ratio == min(mse$mse_ratio, na.rm = T))[1]
    gap_origin = item[[n]][["synth_original"]] - item[[n]][["value_raw"]]
    gap_new = item[[n]][["synth_new"]] - item[[n]][["value_raw"]]
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
  mutate(artifical_effect = cumsum(c(rep(0, 100*3/5),
                                     seq(0, 1, length.out = 100/20),
                                     seq(1, 0, length.out = 100/20),
                                     rep(0, 100*2/5-100/10))),
         id = 0)



df %>% 
  # filter(unit %in% (mse %>% filter(mse1_pre < 2*10000) %>% .[["dependent"]])) %>% 
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
  geom_vline(xintercept = 60, linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  # coord_cartesian(ylim=c(-32, 32)) +
  xlab("Time") +
  ylab("Synthetic Control - True Value") +
  theme_bw()


simulate_data = function(n = 3,
                         length = 100,
                         rnd_nCycles = seq(0.1, 0.9, length.out = n),
                         rnd_shift = seq(0.9, 0.1, length.out = n),
                         rnd_trend = seq(0.2, 0.8, length.out = n),
                         nCycles_min = 6,
                         nCycles_max = 12,
                         trend_min = 0.01,
                         trend_max = 0.1,
                         beta = 1,
                         t_treat = 60,
                         length_divide = 50,
                         shock = 0.5){
  
  # prepare random numbers
  nCycles = rnd_nCycles * (nCycles_max - nCycles_min) + nCycles_min
  shifts = rnd_shift * 2 * pi
  trends = rnd_trend * (trend_max - trend_min) + trend_min
  
  # simulate
  data = NULL
  rnd_ind = round(runif(1, length - t_treat, t_treat - 1), 0)
  for (i in 1:n) {
    x = seq(0, nCycles[i] * pi, length.out = length)
    x = cumsum(sin(x + shifts[i])/2 + 1)/(length/length_divide)
    y = sin(x)
    if (i == 1) {
      trend = rep(trends[i], length) + c(rep(0, length*3/5),
                                         seq(0, shock, length.out = length/20),
                                         seq(shock, 0, length.out = length/20),
                                         rep(0, length*2/5-length/10))
    }else{
      trend = rep(trends[i], length)
    }
    
    y1 = NULL
    y1lag = 1
    for (j in 1:t_treat) {
      y1t = trend[j] + beta*y1lag + y[j]
      y1 <- c(y1, y1t)
      y1lag = y1t
    }
    
    y2 = NULL
    y2lag = y1lag
    for (j in (t_treat+1):length) {
      y2t = trend[j] + beta*y2lag + y[j - rnd_ind]
      y2 <- c(y2, y2t)
      y2lag = y2t
    }
    
    y = c(y1, y2)
    
    data = rbind(data,
                 data.frame(id = i,
                            unit = LETTERS[i],
                            time = 1:length,
                            value = y,
                            value_raw = y))
  }
  return(data)
}


n_simulation = 1000
n = 5
length = 100
t_treat = 60
shock = 1

# generate sobol sequence
sobol_seq = qrng::sobol(n_simulation*3, d = n, randomize = "Owen",
                        seed = 20220401, skip = 100)
rnd_nCycles = sobol_seq[1:n_simulation,]
rnd_shift = sobol_seq[(n_simulation + 1):(2*n_simulation),]
rnd_trend = sobol_seq[(2*n_simulation + 1):(3*n_simulation),]

# generate data
data_list = NULL
for (i in 1:n_simulation) {
  data_list[[i]] = simulate_data(n = n,
                                 length = length,
                                 rnd_nCycles = rnd_nCycles[i,],
                                 rnd_shift = rnd_shift[i,],
                                 rnd_trend = rnd_trend[i,],                         
                                 t_treat = t_treat,
                                 shock = shock)
}


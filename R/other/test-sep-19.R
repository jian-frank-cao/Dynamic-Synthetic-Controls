results = readRDS("./data/res_sim_0916.Rds")
data_list = readRDS("./data/simul_data_list_0915.Rds")

res = results[[1]][[35]]

data_list[[2]] %>% ggplot(aes(x = time, y = value, color = unit)) +
  geom_line() +
  geom_vline(xintercept = 80, linetype="dashed")

log(result[["results.TFDTW.synth"]][["A"]][["mse"]][["mse.postT.TFDTW"]]/result[["results.TFDTW.synth"]][["A"]][["mse"]][["mse.postT.raw"]])

mse = future_map2(
  results[[2]],
  search.grid.list,
  ~{
    item = .x
    task = .y
    task$neg.ratio = item$neg.ratio
    task$p.value = item$p.value
    task
  }
) %>% do.call("rbind", .)

n = 10
length = 100
rnd_nCycles = seq(0.1, 0.9, length.out = n)
rnd_shift = seq(0.9, 0.1, length.out = n)
rnd_lag = seq(0.1, 0.9, length.out = n)
nCycles_min = 5
nCycles_max = 15
noise_mean = 0
noise_sd = 0.01
n_lag_min = 5
n_lag_max = 15
extra_x = 20
beta = 0.9
ar_x = 0.9
t_treat = 80
shock = 50

data %>%  ggplot(aes(x = time, y = value, color = unit)) +
  geom_line() +
  geom_vline(xintercept = 80, linetype="dashed")

y = NULL
ylag = 1
for (j in 1:length) {
  yt = trend[j] + beta*ylag + x[j]
  y <- c(y, yt)
  ylag = yt
}

i = 3
phi = sin(seq(0, nCycles[i] * pi, length.out = length) + shifts[i])
phi = phi/2+0.5
# trend
if (i == 1) {
  trend = rep(0, length) + c(rep(0, length*4/5),
                             seq(0, shock, length.out = length/20),
                             seq(shock, 0, length.out = length/20),
                             rep(0, length/5-length/10))
}else{
  trend = rep(0, length)
}
y = NULL
ylag = 1
for (j in 1:length) {
  yt = trend[j] + beta*ylag + phi[j]*x[j + n_lags[i]] +
    (1 - phi[j])*x[j] + 
    rnorm(n = 1, mean = noise_mean, sd = noise_sd)
  y <- c(y, yt)
  ylag = yt
}

par(mfrow=c(2,2)) 
plot(ts(Y), ylim = c(-300, 300))
lines(y, col = "red")
plot(ts(phi))
plot(ts(x))

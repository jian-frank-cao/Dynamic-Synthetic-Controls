library(tidyverse)
set.seed(1)

n = 100
burn_in = 20
n_lag = 20
noise_mean = 0
noise_sd = 0.1
beta1 = 0.9

x1 = seq(0, 5 * pi, length.out = n + burn_in)
x2 = seq(0, 7 * pi, length.out = n + burn_in)
phi1 = 0.5+sin(x1)/2
phi2 = 0.5+cos(x2)/2

# initialize y and z
y <- z <- yt <- zt <- NULL
# starting values of y and z
ylag <- zlag <- 1
# Create a sequence of n values for y and z
# x = arima.sim(list(order = c(1,0,0), ar = ar_x), n = n + burn_in)
x = cumsum(sin(x1+0.5*pi)+0.3)

# time series
for (i in (1 + burn_in):(n + burn_in)) {
  yt = beta1*ylag + phi1[i]*x[i] +
    (1 - phi1[i])*x[i-n_lag] + 
    rnorm(n = 1, mean = noise_mean, sd = noise_sd)
  zt = beta1*zlag + phi2[i]*x[i] + 
    (1 - phi2[i])*x[i-n_lag] + 
    rnorm(n = 1, mean = noise_mean, sd = noise_sd)
  y <- c(y, yt)
  z <- c(z, zt)
  ylag = yt
  zlag = zt
}

y = y + 100
y[51:100] = y[51:100] + (1:50)*5
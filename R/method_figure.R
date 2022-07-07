library(checkpoint)
checkpoint("2022-04-01")


library(tidyverse)

plot_ts = function(x, y1, y2 = NULL, textcex = 2,
                   lineWidth = 2){
  plot(x, y1,
       type = "l",
       cex.axis = textcex,
       las = 1,
       col = "red",
       lwd = lineWidth,
       xlab='Time', cex.lab=textcex,
       ylab = 'Y'
  )
  
  if (!is.null(y2)) {
    lines(x, y2,
          lty = 2,
          col = 'blue',
          lwd = lineWidth
    )
  }
}


nCycles1 = 9
nCycles2 = 7.5
length = 1000
shock = 0.5
trend = 0.05

x1 = seq(0, nCycles1 * pi, length.out = length)
x2 = seq(0, nCycles2 * pi, length.out = length)

x3 = cumsum(sin(x1)/2 + 1)/20
x4 = cumsum(cos(x2)/2 + 1)/20

y1 = sin(x3)
y2 = sin(x4)

trend1 = rep(trend, 1000) + c(rep(0, 500),
                            seq(0, shock, length.out = 50),
                            seq(shock, 0, length.out = 50),
                            rep(0, 400))
trend2 = trend

y1 = cumsum(y1/2 + trend1)
y2 = cumsum(y2/2 + trend2)

plot_ts(1:1000,y1,y2)


res_dtw = dtw::dtw(y1[1:500], y2[1:500],
                   open.end = TRUE, keep = TRUE,
                   step.pattern = dtw::symmetricP2)
dtw::dtwPlotTwoWay(res_dtw)


# Parameters
noise = 0.05
nCycles = 3.5
speedUp = 1.2
speedDown = 0.8
par(bty = 'n')
trend = 0.3
length = 2000

x <- seq(0, nCycles * pi, length.out = length)
x0 = seq(0, 10*nCycles * pi, length.out = length)
speed1 = sin(x0)
#normalize speed to be between 0.8 and 1.2 (otherwise get swings that are too extreme)
speed1 =  0.5 + (speed1/2 )
# plot(speed1, main = 'sample speed 1')

speed2 = sin(-x)
#normalize speed 
speed2 =  0.5 + (speed2/2 )
# plot(speed2, main = 'sample speed 2')

x1 = cumsum(speed1)/10
x2 = cumsum(speed2)/20

y1 = sin(x1)
y2 = sin(x2)

y1 <- sin(x) #+ rnorm(length, sd = noise)
y2 <- sin(speed1 * x) #+ rnorm(length, sd = noise)
y3 <- sin(speed2 * x) #+ rnorm(length, sd = noise)
dat = data.frame(unit = c(rep('A', length),
                          rep('B', length),
                          rep('C',length)),
                 time = rep(1:length, 3),
                 Y = c(y1,y2,y3) )
write.csv(dat, './data/simulData.csv', row.names = FALSE)


# Plot curves
plot(x, y2,
     type = "l",
     cex.axis = textcex,
     las = 1,
     lwd = lineWidth,
     xlab='Time', cex.lab=textcex,
     ylab = 'Y'
)

lines(x, y2,
      lty = 2,
      #col = 'darkgrey',
      lwd = lineWidth
)

x <- seq(0, nCycles * pi, length.out = length)
lines(x, y3,
      lty = 3,
      #col = 'lightgrey',
      lwd = lineWidth
)

legend('topleft', 
       legend = c('Unit T (constant speed)', 'Unit U1 (speed 1)', 'Unit U2 (Speed 2)'),
       cex=textcex,
       lty=1:3)













# n = 100
# burn_in = 20
# n_lag = 20
# noise_mean = 0
# noise_sd = 0.1
# beta1 = 0.9
# 
# x1 = seq(0, 5 * pi, length.out = n + burn_in)
# x2 = seq(0, 7 * pi, length.out = n + burn_in)
# phi1 = 0.5+sin(x1 + 1.5*pi)/2
# phi2 = 0.5+cos(x2)/2
# 
# 
# plot( phi1,
#       type = 'l',
#       main = 'Speed profiles, i.e., \\phi_{i,t}',
#       ylab = 'phi',
#       xlab = 't'
# )
# lines(phi2, col = 2)
# legend('bottom',
#        c('phi1', 'phi2'),
#        col = c(1, 2),
#        lty = 1)
# 
# 
# # initialize y and z
# y <- z <- yt <- zt <- NULL
# # starting values of y and z
# ylag <- zlag <- 1
# # Create a sequence of n values for y and z
# # x = arima.sim(list(order = c(1,0,0), ar = ar_x), n = n + burn_in)
# x = cumsum(sin(x1+0.5*pi)+0.3)
# 
# # time series
# for (i in (1 + burn_in):(n + burn_in)) {
#   yt = beta1*ylag + phi1[i]*x[i] +
#     (1 - phi1[i])*x[i-n_lag] + 
#     rnorm(n = 1, mean = noise_mean, sd = noise_sd)
#   zt = beta1*zlag + phi2[i]*x[i] + 
#     (1 - phi2[i])*x[i-n_lag] + 
#     rnorm(n = 1, mean = noise_mean, sd = noise_sd)
#   y <- c(y, yt)
#   z <- c(z, zt)
#   ylag = yt
#   zlag = zt
# }
# 
# y = y + 100
# y[51:100] = y[51:100] + (1:50)*5
# 
# # plot the sequences
# mini = min(c(y, z))
# maxi = max(c(y, z))
# plot(y,
#      type = 'l',
#      ylim = c(mini, maxi),
#      main = 'Time series')
# lines(z, col = 2)
# legend('bottomright',
#        c('unit 1', 'unit 2'),
#        col = c(1, 2),
#        lty = 1)


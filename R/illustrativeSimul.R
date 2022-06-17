# Parameters
noise = 0.05
nCycles = 3.5
speedUp = 1.2
speedDown = 0.8
par(bty = 'n')
lineWidth = 2
textcex=2
trend = 0.3
length = 150
set.seed(20220407)

x <- seq(0, nCycles * pi, length.out = length)
speed1 = sin(x)
#normalize speed to be between 0.8 and 1.2 (otherwise get swings that are too extreme)
speed1 =  1 + (speed1/5 )
plot(speed1, main = 'sample speed 1')

speed2 = sin(-x)
#normalize speed 
speed2 =  1 + (speed2/5 )
plot(speed2, main = 'sample speed 2')

y1 <- sin(x) + rnorm(length, sd = noise)
y2 <- sin(speed1 * x) + rnorm(length, sd = noise)
y3 <- sin(speed2 * x) + rnorm(length, sd = noise)
dat = data.frame(unit = c(rep('A', length),
                          rep('B', length),
                          rep('C',length)),
                 time = rep(1:length, 3),
                 Y = c(y1,y2,y3) )
write.csv(dat, './data/simulData.csv', row.names = FALSE)


# Plot curves
plot(x, y1,
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


## Ver. 2 ----------------------------------------------------------------------
set.seed(3)
nCycles = 2.5
x <- seq(0, nCycles * pi, length.out = 100)
phi1 = sin(x)

#normalize speed to be between 0.9 and 1.1
# (otherwise we get swings that are too extreme) 
phi1= 1+(phi1/10)

phi2 = sin(-x) 
#normalize speed
phi2= 1+(phi2/10)
# Plot the speed profiles
plot( phi1,
      type = 'l',
      main = 'Speed profiles, i.e., \\phi_{i,t}',
      ylab = 'phi',
      xlab = 't'
)
lines(phi2, col = 2)
legend('bottom',
       c('phi1', 'phi2'),
       col = c(1, 2),
       lty = 1)


# initialize y and z
y <- z <- u <- NULL
# starting values of y and z
ylag <- zlag <- ulag <- ulag2 <- 1
# Create a sequence of 100 values for y and z
noise = 100:1/200 - 0.25
for (i in 1:100) {
  trend = 0
  # noise = rnorm(n = 1, mean = 0)
  yt = noise[i] + 0.9 * ylag * phi1[i]
  zt = noise[i] + 0.9 * zlag * phi2[i]
  y <- c(y, yt)
  z <- c(z, zt)
  u <- c(u, ut)
  ylag = yt
  zlag = zt
}
# plot the sequences
mini = min(c(y, z))
maxi = max(c(y, z))
plot(y,
     type = 'l',
     ylim = c(mini, maxi),
     main = 'Time series')
lines(z, col = 2)
legend('bottomright',
       c('unit 1', 'unit 2'),
       col = c(1, 2),
       lty = 1)

## Ver. 3 ----------------------------------------------------------------------
set.seed(3)
n = 100
nCycles = 2.5
x1 <- seq(0, nCycles * pi, length.out = n)
x2 <- seq(0, nCycles * pi, length.out = n)
x3 <- seq(0, nCycles * pi, length.out = n)


phi1 =  0.5 + (sin(x1) / 2)
phi2 = 0.5 + (sin(x2 + 0.66*pi) / 2)
phi3 = 0.5 + (sin(x3 + 0.33*pi) /2)

phi1 =  1/(1+exp(cumsum(rnorm(100))))
phi2 = 1/(1+exp(cumsum(rnorm(100))))
phi3 = 1/(1+exp(cumsum(rnorm(100))))



# Plot the speed profiles
plot( phi1,
      type = 'l',
      main = 'Speed profiles, i.e., \\phi_{i,t}',
      ylab = 'phi',
      xlab = 't',
      col = 3
)
lines(phi2, col = 2)
lines(phi3, col = 1)
legend('bottom',
       c('phi1', 'phi2', "phi3"),
       col = c(3, 2, 1),
       lty = 1)

# initialize y and z
y <- z <- yt <- zt <- u <- ut <- NULL
# starting values of y and z
ylag <- zlag <- ulag <- 1
# Create a sequence of 100 values for y and z
set.seed(9)
# x =  cumsum(rep(1, 100))
x = cumsum(rnorm(n = 100, mean = 0))
# x =  cumsum(sin(x1+0.5*pi)/2+0.3)

for (i in 21:100) {
  yt = 0.9*ylag + phi1[i]*x[i] + (1-phi1[i])*x[i-10] + rnorm(n = 1, mean = 0, sd = 0.1)
  zt = 0.9*zlag + phi3[i]*x[i] + (1-phi3[i])*x[i-10] + rnorm(n = 1, mean = 0, sd = 0.1)
  ut = 0.9*ulag + phi2[i]*x[i] + (1-phi2[i])*x[i-10] + rnorm(n = 1, mean = 0, sd = 0.1)
  # ut =  0.9* ulag + x[i] + rnorm(n = 1, mean = 0, sd = 0.1)
  y <- c(y, yt)
  z <- c(z, zt)
  u <- c(u, ut)
  ylag = yt
  zlag = zt
  ulag = ut
}

# plot the sequences
mini = min(c(y, z, u))
maxi = max(c(y, z, u))
#pdf('sample.pdf')
plot(y, lwd=2, col = 3,
     type = 'l',
     ylim = c(mini, maxi),
     main = 'Time series')
lines(z, col = 2, lwd=2)
lines(u, col = 1, lwd=2)
legend('bottomright',
       c('unit 1', 'unit 2', 'unit 3'),
       col = c(3, 2, 1),
       lty = 1,
       lwd=2, cex=2)

length = length(y)
dat = data.frame(unit = c(rep('A', length),
                          rep('B', length),
                          rep('C',length)),
                 time = rep(1:length, 3),
                 Y = c(u,y,z) )
write.csv(dat, './data/simulData_v3.csv', row.names = FALSE)



